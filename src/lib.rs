//! Core library for microwave filter synthesis experiments in Rust.
//!
//! The crate is organized around a small end-to-end workflow:
//! build a [`FilterSpec`], map physical frequencies through a [`FrequencyMapping`],
//! synthesize prototype polynomials, convert them into a coupling matrix, and
//! optionally evaluate the resulting S-parameter response on a frequency grid.

pub mod approx;
pub mod error;
pub mod fixtures;
pub mod freq;
pub mod matrix;
pub mod prelude;
pub mod response;
pub mod spec;
pub mod synthesis;
pub mod transform;
pub mod verify;

pub use approx::{
    ApproximationEngine, ChebyshevApproximation, GeneralizedChebyshevApproximation, PolynomialSet,
    PrototypePoint,
};
pub use error::{MfsError, Result};
pub use freq::{BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping, NormalizedSample};
pub use matrix::{
    BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder, MatrixShape, MatrixTopology,
};
pub use response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use spec::{
    ApproximationFamily, FilterClass, FilterSpec, FilterSpecBuilder, ReturnLossSpec,
    TransmissionZero, TransmissionZeroDomain,
};
pub use synthesis::{
    AdmittancePolynomials, ApproximationStageKind, CanonicalMatrixSynthesis, ChebyshevSynthesis,
    EvaluationOutcome, MatrixSynthesisEngine, MatrixSynthesisMethod, MatrixSynthesisOutcome,
    ResidueExpansion, ResiduePole, SectionSynthesis, SynthesisOutcome,
    VerifiedSectionSynthesis,
};
pub use transform::{
    extract_quadruplet_section, extract_quadruplet_section_with_report,
    extract_quadruplet_section_with_response_check, extract_triplet_section,
    extract_triplet_section_with_report, extract_triplet_section_with_response_check,
    extract_trisection_section, extract_trisection_section_with_report,
    extract_trisection_section_with_response_check, to_wheel, transform_matrix,
    transform_matrix_with_response_check, SectionTransformEngine, SectionTransformOutcome,
    TopologyKind, TransformEngine, TransformOutcome, TransformReport,
};
pub use verify::{
    compare_responses, matches_arrow_pattern, matches_folded_pattern, matches_topology_pattern,
    verify_quadruplet_extraction, verify_triplet_extraction, verify_trisection_extraction,
    MatrixPatternTolerance, ResponseCheckReport, ResponseComparison, ResponseTolerance,
    SectionVerificationReport,
};

/// Synthesizes a canonical coupling matrix from precomputed approximation polynomials.
pub fn synthesize_canonical_matrix(polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
    CanonicalMatrixSynthesis::default().synthesize(polynomials)
}

/// Synthesizes a canonical matrix and reports which internal synthesis path was used.
pub fn synthesize_canonical_matrix_with_details(
    polynomials: &PolynomialSet,
) -> Result<MatrixSynthesisOutcome> {
    CanonicalMatrixSynthesis::default().synthesize_with_details(polynomials)
}

/// Synthesizes a canonical matrix and converts it into the requested topology.
pub fn synthesize_matrix_with_topology(
    polynomials: &PolynomialSet,
    topology: TopologyKind,
) -> Result<CouplingMatrix> {
    CanonicalMatrixSynthesis::default().synthesize_with_topology(polynomials, topology)
}

/// Convenience helper that runs the current Chebyshev synthesis flow.
pub fn synthesize_chebyshev(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
) -> Result<(PolynomialSet, CouplingMatrix)> {
    let outcome = ChebyshevSynthesis::default().synthesize(spec, mapping)?;
    Ok((outcome.polynomials, outcome.matrix))
}

/// Convenience helper that exposes approximation and matrix-path metadata.
pub fn synthesize_chebyshev_with_details(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
) -> Result<SynthesisOutcome> {
    ChebyshevSynthesis::default().synthesize_with_details(spec, mapping)
}

/// Convenience helper that synthesizes a design and evaluates it on a physical grid.
pub fn synthesize_and_evaluate_chebyshev_with_mapping(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<(PolynomialSet, CouplingMatrix, SParameterResponse)> {
    let outcome =
        ChebyshevSynthesis::default().synthesize_and_evaluate_with_mapping(spec, mapping, grid)?;
    Ok((outcome.polynomials, outcome.matrix, outcome.response))
}

/// Convenience helper that evaluates a design and preserves approximation and matrix-path metadata.
pub fn synthesize_and_evaluate_chebyshev_with_details(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<EvaluationOutcome> {
    ChebyshevSynthesis::default().synthesize_and_evaluate_with_details(spec, mapping, grid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::load_filter_database_end_to_end_fixture;

    #[test]
    fn basic_synthesis_pipeline_compiles_and_runs() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
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
        let spec = FilterSpec::chebyshev(3, 20.0)?
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

    #[test]
    fn canonical_helper_functions_follow_new_high_level_naming() -> Result<()> {
        let spec = FilterSpec::builder()
            .order(4)
            .return_loss_db(20.0)
            .transmission_zeros(vec![
                TransmissionZero::normalized(-2.0),
                TransmissionZero::normalized(1.5),
            ])
            .build()?;
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let outcome = synthesize_canonical_matrix_with_details(&polynomials)?;
        let arrow = synthesize_matrix_with_topology(&polynomials, TopologyKind::Arrow)?;

        assert_eq!(outcome.matrix.order(), 4);
        assert_eq!(arrow.order(), 4);
        assert!(arrow.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn high_level_chebyshev_helpers_can_report_generalized_main_flow_details() -> Result<()> {
        let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;

        let outcome = synthesize_chebyshev_with_details(&fixture.spec, &fixture.mapping)?;
        assert_eq!(outcome.approximation_kind, ApproximationStageKind::GeneralizedChebyshev);
        assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::ResidueExpansion);
        Ok(())
    }
}
