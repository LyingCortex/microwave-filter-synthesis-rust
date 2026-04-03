//! Core library for microwave filter synthesis experiments in Rust.
//!
//! The crate is organized around a small end-to-end workflow:
//! build a [`FilterParameter`], map physical frequencies through a [`FrequencyMapping`],
//! synthesize prototype polynomials, convert them into a coupling matrix, and
//! optionally evaluate the resulting S-parameter response on a frequency grid.

pub mod approx;
pub mod error;
pub mod freq;
pub mod matrix;
pub mod response;
pub mod spec;
pub mod synthesis;

pub use approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet, PrototypePoint};
pub use error::{MfsError, Result};
pub use freq::{BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping, NormalizedSample};
pub use matrix::{
    AdmittancePolynomials, BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder,
    CouplingMatrixSynthesizer, MatrixShape, MatrixSynthesisMethod, MatrixSynthesisOutcome,
    MatrixTopology, ResidueExpansion, ResiduePole,
};
pub use response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use spec::{
    ApproximationFamily, FilterClass, FilterParameter, ReturnLossSpec, TransmissionZero, TransmissionZeroDomain,
};
pub use synthesis::{ChebyshevSynthesis, EvaluationOutcome, SynthesisOutcome};

/// Convenience helper that runs the current Chebyshev synthesis flow.
pub fn synthesize_chebyshev(
    spec: &FilterParameter,
    mapping: &impl FrequencyMapping,
) -> Result<(PolynomialSet, CouplingMatrix)> {
    let outcome = ChebyshevSynthesis::default().synthesize(spec, mapping)?;
    Ok((outcome.polynomials, outcome.matrix))
}

/// Convenience helper that synthesizes a design and evaluates it on a physical grid.
pub fn synthesize_and_evaluate_chebyshev_with_mapping(
    spec: &FilterParameter,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<(PolynomialSet, CouplingMatrix, SParameterResponse)> {
    let outcome =
        ChebyshevSynthesis::default().synthesize_and_evaluate_with_mapping(spec, mapping, grid)?;
    Ok((outcome.polynomials, outcome.matrix, outcome.response))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_synthesis_pipeline_compiles_and_runs() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::finite(-1.25)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let (polynomials, matrix) = synthesize_chebyshev(&spec, &mapping)?;
        assert_eq!(polynomials.order, 4);
        assert_eq!(polynomials.transmission_zeros_normalized.len(), 1);
        assert!(polynomials.generalized.is_some());
        assert_eq!(matrix.order(), 4);
        assert!(matrix.at(0, 1).unwrap_or_default() > 0.0);
        Ok(())
    }

    #[test]
    fn response_solver_returns_matching_grid_length() -> Result<()> {
        let matrix = CouplingMatrix::identity(3)?;
        let grid = FrequencyGrid::linspace(6.0e9, 7.0e9, 11)?;
        let response = ResponseSolver::default().evaluate_normalized(&matrix, &grid)?;

        assert_eq!(response.samples.len(), 11);
        assert_eq!(response.samples[0].frequency_hz, 6.0e9);
        Ok(())
    }

    #[test]
    fn high_level_synthesis_and_evaluation_pipeline_runs() -> Result<()> {
        let spec = FilterParameter::chebyshev(3, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 7)?;

        let (polynomials, matrix, response) =
            synthesize_and_evaluate_chebyshev_with_mapping(&spec, &mapping, &grid)?;
        assert_eq!(polynomials.order, 3);
        assert_eq!(matrix.order(), 3);
        assert_eq!(response.samples.len(), 7);
        assert_ne!(response.samples[0].s21_re, response.samples[3].s21_re);
        Ok(())
    }
}
