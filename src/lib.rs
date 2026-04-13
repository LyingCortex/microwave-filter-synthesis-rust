//! Core library for microwave filter synthesis experiments in Rust.
//!
//! The crate is organized around a small end-to-end workflow:
//! build a [`FilterSpec`] with normalized transmission zeros, synthesize
//! prototype polynomials, convert them into a coupling matrix, and optionally
//! evaluate the resulting S-parameter response on a physical frequency grid.
//!
//! Frequency convention:
//! - [`FilterSpec`] stores transmission zeros in normalized prototype coordinates
//! - physical-frequency zeros should be converted up front with
//!   [`normalize_transmission_zeros_hz`]
//! - [`FrequencyMapping`] is used for physical-frequency evaluation and reporting,
//!   not for normalizing zeros inside a spec
//! - approximation reads and validates the normalized zeros already present in
//!   the spec; it does not normalize physical-frequency zeros itself

pub mod approx;
pub mod error;
pub mod fixtures;
pub mod freq;
pub mod matrix;
pub mod output;
pub mod prelude;
pub mod response;
pub mod spec;
pub mod synthesis;
pub mod transform;
pub mod verify;

pub use approx::{
    ApproximationEngine, GeneralizedChebyshevApproximation, PolynomialSet, PrototypePoint,
};
pub use error::{MfsError, Result};
pub use freq::{
    BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping, NormalizedSample,
};
pub use matrix::{
    BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder, MatrixShape, MatrixTopology,
};
pub use response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use spec::{
    FilterSpec, FilterSpecBuilder, OutOfBandAttenuationSpec, OutOfBandAttenuationWindow,
    TransmissionZero,
};
pub use synthesis::{
    AdmittancePolynomials, CanonicalMatrixSynthesis, EvaluationOutcome, MatrixSynthesisEngine,
    MatrixSynthesisMethod, MatrixSynthesisOutcome, ResidueExpansion, ResiduePole,
    SectionSynthesis, SynthesisOutcome, VerifiedSectionSynthesis,
    synthesize_and_evaluate_with_mapping,
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

/// Default unloaded Q used by `filter_spec` when none is specified.
pub const DEFAULT_UNLOADED_Q: f64 = 2000.0;

/// Builds a normalized prototype filter specification from order, return loss,
/// and transmission zeros.
///
/// Usage:
/// - pass normalized prototype zeros such as `[f64, ...]`
/// - no zeros: pass `None`
/// - if your zeros start in physical Hz, pre-normalize them with
///   `normalize_transmission_zeros_hz(...)` before building the spec
/// - unloaded Q defaults to `2000.0`
///
/// ```rust
/// use mfs::prelude::*;
///
/// let spec = filter_spec(4, 20.0, [-1.5, 2.0], None)?;
/// # Ok::<(), mfs::MfsError>(())
/// ```
///
/// ```rust
/// use mfs::prelude::*;
///
/// let spec = filter_spec::<Vec<f64>>(4, 20.0, None, None)?;
/// # Ok::<(), mfs::MfsError>(())
/// ```
///
/// ```rust
/// use mfs::prelude::*;
///
/// let mapping = bandpass(6.75e9, 300.0e6)?;
/// let zeros = normalize_transmission_zeros_hz([6.5e9, 7.0e9], &mapping)?;
/// let spec = filter_spec(4, 20.0, zeros, None)?;
/// # Ok::<(), mfs::MfsError>(())
/// ```
pub fn filter_spec<T>(
    order: usize,
    return_loss_db: f64,
    zeros: impl Into<Option<T>>,
    unloaded_q: impl Into<Option<f64>>,
) -> Result<FilterSpec>
where
    T: IntoIterator<Item = f64>,
{
    let transmission_zeros = zeros
        .into()
        .into_iter()
        .flat_map(|iter| iter.into_iter());
    let mut spec = FilterSpec::new(order, return_loss_db)?
        .with_normalized_transmission_zeros(transmission_zeros);
    spec.unloaded_q = Some(unloaded_q.into().unwrap_or(DEFAULT_UNLOADED_Q));
    Ok(spec)
}

/// Builds a normalized low-pass mapping helper.
pub fn lowpass(cutoff: f64) -> Result<LowPassMapping> {
    LowPassMapping::new(cutoff)
}

/// Normalizes physical-frequency transmission zeros into prototype coordinates.
pub fn normalize_transmission_zeros_hz(
    zeros_hz: impl IntoIterator<Item = f64>,
    mapping: &impl FrequencyMapping,
) -> Result<Vec<f64>> {
    zeros_hz
        .into_iter()
        .map(|zero_hz| mapping.map_hz_to_normalized(zero_hz).map(|sample| sample.omega))
        .collect()
}

/// Builds a band-pass mapping helper.
pub fn bandpass(center_hz: f64, bandwidth_hz: f64) -> Result<BandPassMapping> {
    BandPassMapping::new(center_hz, bandwidth_hz)
}

/// Synthesizes a normalized generalized-Chebyshev prototype and canonical matrix.
///
/// `spec.transmission_zeros` must already be expressed in normalized prototype
/// coordinates. If you start from physical Hz zeros, convert them first with
/// [`normalize_transmission_zeros_hz`].
pub fn generalized_chebyshev(spec: &FilterSpec) -> Result<SynthesisOutcome> {
    synthesis::synthesize_generalized_chebyshev(spec)
}

/// Exposes the intermediate generalized-Chebyshev prototype polynomials for debugging or inspection.
pub fn generalized_chebyshev_polynomials(spec: &FilterSpec) -> Result<PolynomialSet> {
    GeneralizedChebyshevApproximation.synthesize(spec)
}

/// Synthesizes and evaluates a design on a physical-frequency grid.
///
/// The supplied `spec` still expects normalized transmission zeros. `mapping`
/// is used here for the physical-frequency evaluation grid, not to normalize
/// zeros inside the spec.
pub fn generalized_chebyshev_with_response(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<EvaluationOutcome> {
    synthesize_and_evaluate_generalized_chebyshev_with_details(spec, mapping, grid)
}

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

/// Convenience helper that runs the normalized generalized-Chebyshev synthesis flow.
pub fn synthesize_generalized_chebyshev(spec: &FilterSpec) -> Result<(PolynomialSet, CouplingMatrix)> {
    let outcome = synthesis::synthesize_generalized_chebyshev(spec)?;
    Ok((outcome.polynomials, outcome.matrix))
}

/// Convenience helper that synthesizes a design and evaluates it on a physical grid.
pub fn synthesize_and_evaluate_generalized_chebyshev(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<(PolynomialSet, CouplingMatrix, SParameterResponse)> {
    let outcome = synthesize_and_evaluate_with_mapping(spec, mapping, grid)?;
    Ok((outcome.synthesis.polynomials, outcome.synthesis.matrix, outcome.response))
}

/// Convenience helper that evaluates a design and preserves approximation and matrix-path metadata.
pub fn synthesize_and_evaluate_generalized_chebyshev_with_details(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<EvaluationOutcome> {
    synthesize_and_evaluate_with_mapping(spec, mapping, grid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::load_filter_database_end_to_end_fixture;

    #[test]
    fn basic_synthesis_pipeline_compiles_and_runs() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::finite(-1.25)]);

        let (polynomials, matrix) = synthesize_generalized_chebyshev(&spec)?;
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
        let spec = FilterSpec::new(3, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::normalized(2.0)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 7)?;

        let (polynomials, matrix, response) =
            synthesize_and_evaluate_generalized_chebyshev(&spec, &mapping, &grid)?;
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
            .normalized_transmission_zeros(vec![-2.0, 1.5])
            .build()?;
        let polynomials = GeneralizedChebyshevApproximation.synthesize(&spec)?;

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

        let outcome = generalized_chebyshev(&fixture.spec)?;
        assert_eq!(outcome.approximation_kind(), "GeneralizedChebyshev");
        assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::ResidueExpansion);
        Ok(())
    }

    #[test]
    fn shortcut_helpers_support_default_normalized_prototype_flow() -> Result<()> {
        let spec = filter_spec(4, 20.0, [2.0, -1.5], None)?;
        let outcome = generalized_chebyshev(&spec)?;

        assert_eq!(spec.order, 4);
        assert_eq!(spec.transmission_zeros.len(), 2);
        assert_eq!(outcome.polynomials.order, 4);
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn shortcut_helpers_support_physical_zero_specs() -> Result<()> {
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let zeros = normalize_transmission_zeros_hz([7.0e9], &mapping)?;
        let spec = filter_spec(3, 20.0, zeros, None)?;
        let polynomials = generalized_chebyshev(&spec)?.polynomials;

        assert_eq!(spec.transmission_zeros.len(), 1);
        assert!(spec.transmission_zeros[0].value > 1.0);
        assert_eq!(polynomials.order, 3);
        Ok(())
    }

    #[test]
    fn shortcut_helpers_use_default_unloaded_q() -> Result<()> {
        let spec = filter_spec(4, 20.0, [-1.5], None)?;

        assert_eq!(spec.unloaded_q, Some(DEFAULT_UNLOADED_Q));
        Ok(())
    }
}
