//! Microwave filter synthesis core library.
//!
//! This crate defines the base domain model for a future Rust rewrite of the
//! Python prototype. The current implementation focuses on a stable module
//! layout, typed specifications, and a minimal end-to-end synthesis pipeline.

pub mod approx;
pub mod error;
pub mod freq;
pub mod matrix;
pub mod response;
pub mod spec;
pub mod synthesis;

pub use approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet, PrototypePoint};
pub use error::{MfsError, Result};
pub use freq::{BandPassPlan, FrequencyGrid, FrequencyPlan, LowPassPlan, NormalizedSample};
pub use matrix::{CouplingMatrix, CouplingMatrixBuilder, CouplingMatrixSynthesizer, MatrixShape};
pub use response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use spec::{
    ApproximationFamily, FilterClass, FilterSpec, FilterType, PerformanceSpec, TransmissionZero,
    TransmissionZeroDomain,
};
pub use synthesis::{ChebyshevSynthesis, EvaluationOutcome, SynthesisOutcome};

/// High-level helper that mirrors the intended synthesis workflow.
pub fn synthesize_chebyshev(
    spec: &FilterSpec,
    plan: &impl FrequencyPlan,
) -> Result<(PolynomialSet, CouplingMatrix)> {
    let outcome = ChebyshevSynthesis::default().synthesize(spec, plan)?;
    Ok((outcome.polynomials, outcome.matrix))
}

/// High-level helper that synthesizes a prototype and evaluates its response.
pub fn synthesize_and_evaluate_chebyshev(
    spec: &FilterSpec,
    plan: &impl FrequencyPlan,
    grid: &FrequencyGrid,
) -> Result<(PolynomialSet, CouplingMatrix, SParameterResponse)> {
    let outcome = ChebyshevSynthesis::default().synthesize_and_evaluate(spec, plan, grid)?;
    Ok((outcome.polynomials, outcome.matrix, outcome.response))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_synthesis_pipeline_compiles_and_runs() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::finite(-1.25)]);
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;

        let (polynomials, matrix) = synthesize_chebyshev(&spec, &plan)?;
        assert_eq!(polynomials.order, 4);
        assert_eq!(polynomials.transmission_zeros_normalized.len(), 1);
        assert_eq!(matrix.order(), 4);
        assert!(matrix.at(0, 1).unwrap_or_default() > 0.0);
        Ok(())
    }

    #[test]
    fn response_solver_returns_matching_grid_length() -> Result<()> {
        let matrix = CouplingMatrix::identity(3)?;
        let grid = FrequencyGrid::linspace(6.0e9, 7.0e9, 11)?;
        let response = ResponseSolver::default().evaluate(&matrix, &grid)?;

        assert_eq!(response.samples.len(), 11);
        assert_eq!(response.samples[0].frequency_hz, 6.0e9);
        Ok(())
    }

    #[test]
    fn high_level_synthesis_and_evaluation_pipeline_runs() -> Result<()> {
        let spec = FilterSpec::chebyshev(3, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 7)?;

        let (polynomials, matrix, response) =
            synthesize_and_evaluate_chebyshev(&spec, &plan, &grid)?;
        assert_eq!(polynomials.order, 3);
        assert_eq!(matrix.order(), 3);
        assert_eq!(response.samples.len(), 7);
        assert_ne!(response.samples[0].s21_re, response.samples[3].s21_re);
        Ok(())
    }
}
