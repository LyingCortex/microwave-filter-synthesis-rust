mod backend;

use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyPlan};
use crate::matrix::CouplingMatrix;

/// One sampled point of the synthesized S-parameter response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseSample {
    /// Physical frequency associated with the sample.
    pub frequency_hz: f64,
    /// Matching normalized frequency used internally by the solver.
    pub normalized_omega: f64,
    /// Approximate group delay derived from the solved response matrix.
    pub group_delay: f64,
    /// Real part of the reflection coefficient.
    pub s11_re: f64,
    /// Imaginary part of the reflection coefficient.
    pub s11_im: f64,
    /// Real part of the forward transmission coefficient.
    pub s21_re: f64,
    /// Imaginary part of the forward transmission coefficient.
    pub s21_im: f64,
}

/// Collection of sampled S-parameter data.
#[derive(Debug, Clone, PartialEq)]
pub struct SParameterResponse {
    /// Response samples in the same order as the evaluation grid.
    pub samples: Vec<ResponseSample>,
}

/// Front-end API for response evaluation.
#[derive(Debug, Default, Clone, Copy)]
pub struct ResponseSolver;

/// Source/load terminations used by the response backend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseSettings {
    /// Source resistance in normalized units.
    pub source_resistance: f64,
    /// Load resistance in normalized units.
    pub load_resistance: f64,
}

impl Default for ResponseSettings {
    fn default() -> Self {
        Self {
            source_resistance: 1.0,
            load_resistance: 1.0,
        }
    }
}

impl ResponseSolver {
    /// Evaluates a matrix directly on a grid that is already normalized.
    pub fn evaluate(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
    ) -> Result<SParameterResponse> {
        self.evaluate_with_settings(matrix, grid, ResponseSettings::default())
    }

    /// Evaluates a matrix on a physical grid after mapping it through a frequency plan.
    pub fn evaluate_with_plan(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        plan: &impl FrequencyPlan,
    ) -> Result<SParameterResponse> {
        self.evaluate_with_plan_and_settings(matrix, grid, plan, ResponseSettings::default())
    }

    /// Evaluates a normalized grid with explicitly supplied source/load settings.
    pub fn evaluate_with_settings(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        settings: ResponseSettings,
    ) -> Result<SParameterResponse> {
        backend::evaluate_normalized_response(matrix, grid.as_slice(), settings)
    }

    /// Evaluates a physical grid with explicit frequency mapping and terminations.
    pub fn evaluate_with_plan_and_settings(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        plan: &impl FrequencyPlan,
        settings: ResponseSettings,
    ) -> Result<SParameterResponse> {
        let normalized_omegas = plan
            .map_grid_hz_to_normalized(grid)?
            .into_iter()
            .map(|sample| sample.omega)
            .collect::<Vec<_>>();
        // Keep both physical and normalized axes so callers can inspect either one.
        backend::evaluate_response(matrix, grid.as_slice(), &normalized_omegas, settings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::{BandPassPlan, LowPassPlan};
    use crate::matrix::CouplingMatrixBuilder;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn response_with_plan_varies_across_frequency() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.8)?
            .set(1, 1, -0.3)?
            .set_symmetric(1, 2, 0.5)?
            .set(2, 2, 0.3)?
            .set_symmetric(2, 3, 0.9)?
            .build()?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 5)?;
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;

        let response = ResponseSolver::default().evaluate_with_plan(&matrix, &grid, &plan)?;
        assert_eq!(response.samples.len(), 5);
        assert_ne!(response.samples[0].s21_re, response.samples[2].s21_re);
        assert_eq!(response.samples[2].normalized_omega, 0.0);
        assert!(response.samples[2].group_delay.is_finite());
        Ok(())
    }

    #[test]
    fn response_default_evaluate_uses_normalized_grid_values() -> Result<()> {
        let matrix = CouplingMatrix::identity(2)?;
        let grid = FrequencyGrid::linspace(1.0, 2.0, 3)?;
        let response = ResponseSolver::default().evaluate(&matrix, &grid)?;

        assert_eq!(response.samples[0].normalized_omega, 1.0);
        assert!(response.samples[0].group_delay.is_finite());
        Ok(())
    }

    #[test]
    fn response_with_lowpass_plan_uses_normalized_mapping() -> Result<()> {
        let matrix = CouplingMatrix::identity(1)?;
        let grid = FrequencyGrid::linspace(1.0e9, 2.0e9, 2)?;
        let plan = LowPassPlan::new(1.0e9)?;

        let response = ResponseSolver::default().evaluate_with_plan(&matrix, &grid, &plan)?;
        assert_eq!(response.samples[0].normalized_omega, 1.0);
        assert_eq!(response.samples[1].normalized_omega, 2.0);
        Ok(())
    }

    #[test]
    fn lossless_response_preserves_power_for_symmetric_case() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(1)?
            .set_symmetric(0, 1, 1.0)?
            .set(1, 1, 0.0)?
            .set_symmetric(1, 2, 1.0)?
            .build()?;
        let grid = FrequencyGrid::linspace(-1.0, 1.0, 5)?;
        let response = ResponseSolver::default().evaluate(&matrix, &grid)?;

        for sample in response.samples {
            let s11_mag_sq = sample.s11_re.powi(2) + sample.s11_im.powi(2);
            let s21_mag_sq = sample.s21_re.powi(2) + sample.s21_im.powi(2);
            approx_eq(s11_mag_sq + s21_mag_sq, 1.0, 1e-9);
        }
        Ok(())
    }
}
