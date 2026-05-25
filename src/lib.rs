//! Microwave filter synthesis library.
//!
//! # Quick Start
//!
//! ```rust
//! use mfs::prelude::*;
//!
//! let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
//!     .zeros_hz([6.4e9, 6.5e9, 7.0e9])
//!     .synthesize()?;
//!
//! let folded = design.to_folded()?;
//! let response = design.response(6.0e9, 7.5e9, 201)?;
//! # Ok::<(), MfsError>(())
//! ```
//!
//! See [`design::FilterDesign`] for the full API.

// ─── Public modules ──────────────────────────────────────────────────────────

/// High-level filter design API (start here).
pub mod design;
/// Error types.
pub mod error;
/// Coupling matrix types and operations.
pub mod matrix;
/// Prelude for convenient imports.
pub mod prelude;
/// S-parameter response evaluation.
pub mod response;
/// Touchstone (.s2p) file export and import.
pub mod touchstone;

// ─── Advanced modules (for power users) ──────────────────────────────────────

/// Polynomial approximation internals.
pub mod approx;
/// Frequency mapping and grid helpers.
pub mod freq;
/// Filter specification types.
pub mod spec;
/// Matrix synthesis engine.
pub mod synthesis;
/// Topology transform facades.
pub mod transform;
/// Response verification helpers.
pub mod verify;
/// Output formatting.
pub mod output;
/// Pipeline orchestration (JSON/CLI).
pub mod pipeline;
/// Literature fixtures for testing.
pub mod fixtures;

/// Python bindings (requires `python` feature).
#[cfg(feature = "python")]
pub mod python;

// ─── Convenience re-exports at crate root ────────────────────────────────────

pub use design::FilterDesign;
pub use error::{MfsError, Result};
pub use matrix::{CouplingMatrix, MatrixTopology};
pub use response::{ResponseSample, SParameterResponse};

// ─── Legacy free functions (kept for backward compatibility) ─────────────────

use crate::synthesis::{
    synthesize_and_evaluate_with_mapping, EvaluationOutcome, MatrixSynthesisEngine,
    MatrixSynthesisOutcome, SynthesisOutcome,
};

pub use freq::{BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping};
pub use matrix::CouplingMatrixBuilder;
pub use response::ResponseSolver;
pub use spec::{FilterSpec, FilterSpecBuilder, TransmissionZero};
pub use synthesis::MatrixSynthesisMethod;
pub use transform::TopologyKind;
pub use verify::ResponseTolerance;
pub use approx::PolynomialSet;

/// Default unloaded Q used by `filter_spec` when none is specified.
pub const DEFAULT_UNLOADED_Q: f64 = 2000.0;

/// Legacy: builds a filter spec. Prefer `FilterDesign::bandpass()` or `FilterDesign::prototype()`.
pub fn filter_spec<T>(
    order: usize,
    return_loss_db: f64,
    zeros: impl Into<Option<T>>,
    unloaded_q: impl Into<Option<f64>>,
) -> Result<FilterSpec>
where
    T: IntoIterator<Item = f64>,
{
    let transmission_zeros = zeros.into().into_iter().flat_map(|iter| iter.into_iter());
    let spec = FilterSpec::new(order, return_loss_db)?
        .with_normalized_transmission_zeros(transmission_zeros)
        .with_unloaded_q(unloaded_q.into().unwrap_or(DEFAULT_UNLOADED_Q));
    Ok(spec)
}

/// Legacy: builds a low-pass mapping.
pub fn lowpass(cutoff: f64) -> Result<LowPassMapping> {
    LowPassMapping::new(cutoff)
}

/// Legacy: normalizes Hz zeros to prototype coordinates.
pub fn normalize_transmission_zeros_hz(
    zeros_hz: impl IntoIterator<Item = f64>,
    mapping: &impl FrequencyMapping,
) -> Result<Vec<f64>> {
    zeros_hz
        .into_iter()
        .map(|hz| mapping.map_hz_to_normalized(hz).map(|s| s.omega))
        .collect()
}

/// Legacy: builds a band-pass mapping.
pub fn bandpass(center_hz: f64, bandwidth_hz: f64) -> Result<BandPassMapping> {
    BandPassMapping::new(center_hz, bandwidth_hz)
}

/// Legacy: synthesizes prototype and matrix. Prefer `FilterDesign`.
pub fn generalized_chebyshev(spec: &FilterSpec) -> Result<SynthesisOutcome> {
    synthesis::synthesize_generalized_chebyshev(spec)
}

/// Legacy: exposes intermediate polynomials.
pub fn generalized_chebyshev_polynomials(spec: &FilterSpec) -> Result<PolynomialSet> {
    approx::generalized_chebyshev_polynomials(spec)
}

/// Legacy: synthesizes and evaluates on a physical grid.
pub fn generalized_chebyshev_with_response(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<EvaluationOutcome> {
    synthesize_and_evaluate_with_mapping(spec, mapping, grid)
}

/// Legacy: synthesizes canonical matrix from polynomials.
pub fn synthesize_canonical_matrix(polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
    MatrixSynthesisEngine.synthesize(polynomials)
}

/// Legacy: synthesizes with details.
pub fn synthesize_canonical_matrix_with_details(
    polynomials: &PolynomialSet,
) -> Result<MatrixSynthesisOutcome> {
    MatrixSynthesisEngine.synthesize_with_details(polynomials)
}

/// Legacy: synthesizes with topology.
pub fn synthesize_matrix_with_topology(
    polynomials: &PolynomialSet,
    topology: TopologyKind,
) -> Result<CouplingMatrix> {
    MatrixSynthesisEngine.synthesize_with_topology(polynomials, topology)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::load_filter_database_end_to_end_fixture;
    use crate::transform::{transform_matrix_with_response_check, extract_triplet_section_with_response_check};

    #[test]
    fn new_api_basic_flow() -> Result<()> {
        let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
        assert_eq!(design.order(), 4);
        let r = design.response(6.5e9, 7.0e9, 11)?;
        assert_eq!(r.samples.len(), 11);
        Ok(())
    }

    #[test]
    fn new_api_with_zeros() -> Result<()> {
        let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
            .zeros_hz([6.4e9, 7.0e9])
            .synthesize()?;
        assert_eq!(design.order(), 6);
        let folded = design.to_folded()?;
        assert_eq!(folded.order(), 6);
        Ok(())
    }

    #[test]
    fn legacy_api_still_works() -> Result<()> {
        let spec = filter_spec(4, 20.0, [-1.5, 2.0], None)?;
        let outcome = generalized_chebyshev(&spec)?;
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn legacy_response_solver() -> Result<()> {
        let matrix = CouplingMatrix::identity(3)?;
        let grid = FrequencyGrid::linspace(6.0e9, 7.0e9, 11)?;
        let response = ResponseSolver::default().evaluate_normalized(&matrix, &grid)?;
        assert_eq!(response.samples.len(), 11);
        Ok(())
    }

    #[test]
    fn legacy_evaluation_pipeline() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::normalized(2.0)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 7)?;
        let outcome = generalized_chebyshev_with_response(&spec, &mapping, &grid)?;
        assert_eq!(outcome.response.samples.len(), 7);
        Ok(())
    }

    #[test]
    fn legacy_topology_and_details() -> Result<()> {
        let spec = FilterSpec::builder()
            .order(4)
            .return_loss_db(20.0)
            .normalized_transmission_zeros(vec![-2.0, 1.5])
            .build()?;
        let polynomials = approx::generalized_chebyshev_polynomials(&spec)?;
        let outcome = synthesize_canonical_matrix_with_details(&polynomials)?;
        let arrow = synthesize_matrix_with_topology(&polynomials, TopologyKind::Arrow)?;
        assert_eq!(outcome.matrix.order(), 4);
        assert_eq!(arrow.order(), 4);
        Ok(())
    }

    #[test]
    fn legacy_fixture_flow() -> Result<()> {
        let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
        let outcome = generalized_chebyshev(&fixture.spec)?;
        assert_eq!(outcome.approximation, crate::synthesis::ApproximationKind::GeneralizedChebyshev);
        assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::ResidueExpansion);
        Ok(())
    }
}
