mod backend;
pub mod pole_expansion;

use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyMapping};
use crate::matrix::CouplingMatrix;

/// One sampled point of the synthesized S-parameter response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseSample {
    /// Physical frequency associated with the sample.
    pub frequency_hz: f64,
    /// Matching normalized frequency used internally by the solver.
    pub normalized_omega: f64,
    /// Approximate group delay derived from the solved response matrix.
    /// `None` when the solver cannot compute it (e.g. the pole-expansion
    /// fast path) — never a silent 0.0.
    pub group_delay: Option<f64>,
    /// Real part of the reflection coefficient.
    pub s11_re: f64,
    /// Imaginary part of the reflection coefficient.
    pub s11_im: f64,
    /// Real part of the forward transmission coefficient.
    pub s21_re: f64,
    /// Imaginary part of the forward transmission coefficient.
    pub s21_im: f64,
}

impl ResponseSample {
    /// S11 magnitude (linear).
    pub fn s11_mag(&self) -> f64 {
        (self.s11_re.powi(2) + self.s11_im.powi(2)).sqrt()
    }

    /// S21 magnitude (linear).
    pub fn s21_mag(&self) -> f64 {
        (self.s21_re.powi(2) + self.s21_im.powi(2)).sqrt()
    }

    /// S11 magnitude in dB.
    pub fn s11_db(&self) -> f64 {
        20.0 * self.s11_mag().log10()
    }

    /// S21 magnitude in dB.
    pub fn s21_db(&self) -> f64 {
        20.0 * self.s21_mag().log10()
    }

    /// S11 phase in degrees.
    pub fn s11_phase_deg(&self) -> f64 {
        self.s11_im.atan2(self.s11_re).to_degrees()
    }

    /// S21 phase in degrees.
    pub fn s21_phase_deg(&self) -> f64 {
        self.s21_im.atan2(self.s21_re).to_degrees()
    }
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

/// Source/load terminations and loss parameters used by the response backend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseSettings {
    /// Source resistance in normalized units.
    pub source_resistance: f64,
    /// Load resistance in normalized units.
    pub load_resistance: f64,
    /// Unloaded Q factor for all resonators (0 or infinity = lossless).
    pub unloaded_q: f64,
}

impl Default for ResponseSettings {
    fn default() -> Self {
        Self {
            source_resistance: 1.0,
            load_resistance: 1.0,
            unloaded_q: f64::INFINITY,
        }
    }
}

impl ResponseSolver {
    /// Evaluates a matrix on a physical grid after mapping it through a frequency mapping.
    pub fn evaluate(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        mapping: &dyn FrequencyMapping,
    ) -> Result<SParameterResponse> {
        self.evaluate_with_settings(matrix, grid, mapping, ResponseSettings::default())
    }

    /// Evaluates a matrix on a physical grid with explicit frequency mapping and terminations.
    pub fn evaluate_with_settings(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        mapping: &dyn FrequencyMapping,
        settings: ResponseSettings,
    ) -> Result<SParameterResponse> {
        let normalized_omegas = mapping
            .map_grid_hz_to_normalized(grid)?
            .into_iter()
            .map(|sample| sample.omega)
            .collect::<Vec<_>>();
        // Keep both physical and normalized axes so callers can inspect either one.
        backend::evaluate_response(matrix, grid.as_slice(), &normalized_omegas, settings)
    }

    /// Evaluates a matrix directly on a grid that is already normalized.
    pub fn evaluate_normalized(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
    ) -> Result<SParameterResponse> {
        self.evaluate_normalized_with_settings(matrix, grid, ResponseSettings::default())
    }

    /// Evaluates a normalized grid with explicitly supplied source/load settings.
    pub fn evaluate_normalized_with_settings(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        settings: ResponseSettings,
    ) -> Result<SParameterResponse> {
        backend::evaluate_normalized_response(matrix, grid.as_slice(), settings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::{BandPassMapping, LowPassMapping};
    use crate::matrix::CouplingMatrixBuilder;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn response_with_mapping_varies_across_frequency() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.8)?
            .set(1, 1, -0.3)?
            .set_symmetric(1, 2, 0.5)?
            .set(2, 2, 0.3)?
            .set_symmetric(2, 3, 0.9)?
            .build()?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 5)?;
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let response = ResponseSolver {}.evaluate(&matrix, &grid, &mapping)?;
        assert_eq!(response.samples.len(), 5);
        assert_ne!(response.samples[0].s21_re, response.samples[2].s21_re);
        assert_eq!(response.samples[2].normalized_omega, 0.0);
        assert!(response.samples[2].group_delay.unwrap().is_finite());
        Ok(())
    }

    #[test]
    fn response_normalized_evaluate_uses_normalized_grid_values() -> Result<()> {
        // A loaded 2-resonator chain (non-degenerate) so the LU backend's
        // group-delay expression has a non-zero denominator.
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.5)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.5)?
            .build()?;
        let grid = FrequencyGrid::linspace(1.0, 2.0, 3)?;
        let response = ResponseSolver {}.evaluate_normalized(&matrix, &grid)?;

        assert_eq!(response.samples[0].normalized_omega, 1.0);
        assert!(response.samples[0].group_delay.unwrap().is_finite());
        Ok(())
    }

    /// The LU backend must mark group delay `None` (not a fake 0.0) when its
    /// expression is numerically degenerate — e.g. an identity matrix whose
    /// transmission column is zero.
    #[test]
    fn group_delay_is_none_when_lu_expression_is_degenerate() -> Result<()> {
        let matrix = CouplingMatrix::identity(2)?;
        let grid = FrequencyGrid::linspace(1.0, 2.0, 3)?;
        let response = ResponseSolver.evaluate_normalized(&matrix, &grid)?;
        assert!(
            response.samples[0].group_delay.is_none(),
            "degenerate LU group delay must be None, not 0.0"
        );
        Ok(())
    }

    /// T4.3: group delay must be explicit per solver path, never a silent 0.0.
    /// The LU backend computes it (`Some`); the pole-expansion fast path does
    /// not and marks it `None`.
    #[test]
    fn group_delay_is_explicit_per_solver_path() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.5)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.5)?
            .build()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 5)?;

        let lu = ResponseSolver.evaluate_normalized(&matrix, &grid)?;
        assert!(
            lu.samples[0].group_delay.is_some(),
            "LU backend must provide a group delay"
        );

        let pole = crate::response::pole_expansion::PoleExpansionData::from_matrix(&matrix)?;
        let pe = pole.evaluate_normalized(&grid)?;
        assert!(
            pe.samples.iter().all(|s| s.group_delay.is_none()),
            "pole-expansion fast path must mark group_delay None, not a silent 0.0"
        );
        Ok(())
    }

    /// P1-3 (value level): the LU group delay must equal `-d/dω arg(S21(ω))`,
    /// not just be finite/some. Derivation: S21 = -2√(R_S R_L)·x[L] where
    /// x = A(ω)⁻¹e₀; by x' = -A⁻¹A'x and A' = I_r (resonator block only),
    /// (S21'/S21) = -Σ_{k=1}^{L-1}(A⁻¹A')_{L,k}x[k]/x[L] = -numerator/denominator,
    /// so group delay = (numerator/denominator).im = -d/dω arg(S21). A central
    /// difference of S21 over the (uniform) grid is an independent reference;
    /// the historical `0..side` fold (which added the two port-node terms)
    /// deviates by O(1) from this (measured 2.679 vs 0.905 at ω=0 on an
    /// order-2 matrix), so the value comparison pins the fix.
    #[test]
    fn lu_group_delay_matches_s21_phase_derivative() -> Result<()> {
        use num_complex::Complex64;
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.5)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.5)?
            .build()?;
        // Fine uniform grid (h = 0.01): the central difference's O(h²)
        // truncation error (amplified past the passband edge where S21 varies
        // rapidly) must stay well under the assert tolerance. At h = 0.1 the
        // worst (edge) gap measured 1.3e-2; at h = 0.01 it drops ~100× to
        // ~1e-4, while the O(1) port-node bug the fix removes is ~1.8.
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 401)?;
        let resp = ResponseSolver.evaluate_normalized(&matrix, &grid)?;
        let step = 4.0 / 400.0; // linspace(-2, 2, 401)
        let s21c = |i: usize| Complex64::new(resp.samples[i].s21_re, resp.samples[i].s21_im);

        let mut checked = 0usize;
        for i in 1..resp.samples.len() - 1 {
            let s21 = s21c(i);
            // Skip points where S21 is near a null: the phase derivative is
            // ill-conditioned there and carries no information about the fix.
            if s21.norm_sqr() < 1e-6 {
                continue;
            }
            // S21' ≈ (S21[i+1] - S21[i-1]) / 2h (central difference).
            let ds21 = (s21c(i + 1) - s21c(i - 1)) / (2.0 * step);
            let numerical = -(ds21 / s21).im;
            let actual = resp.samples[i].group_delay.expect("LU group delay present");
            checked += 1;
            // Loose enough for O(h²)≈1e-4, tight enough to reject the O(1)
            // port-node error the P1-3 fix removes.
            assert!(
                (actual - numerical).abs() < 1e-2 + 1e-3 * numerical.abs(),
                "group delay at ω={}: LU {actual} vs S21 phase derivative {numerical}",
                resp.samples[i].normalized_omega
            );
        }
        assert!(checked > 5, "expected to check > 5 interior points, checked {checked}");
        Ok(())
    }
    #[test]
    fn response_with_lowpass_mapping_uses_normalized_mapping() -> Result<()> {
        let matrix = CouplingMatrix::identity(1)?;
        let grid = FrequencyGrid::linspace(1.0e9, 2.0e9, 2)?;
        let mapping = LowPassMapping::new(1.0e9)?;

        let response = ResponseSolver {}.evaluate(&matrix, &grid, &mapping)?;
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
        let response = ResponseSolver {}.evaluate_normalized(&matrix, &grid)?;

        for sample in response.samples {
            let s11_mag_sq = sample.s11_re.powi(2) + sample.s11_im.powi(2);
            let s21_mag_sq = sample.s21_re.powi(2) + sample.s21_im.powi(2);
            approx_eq(s11_mag_sq + s21_mag_sq, 1.0, 1e-9);
        }
        Ok(())
    }
}
