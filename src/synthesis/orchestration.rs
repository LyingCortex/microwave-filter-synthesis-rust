use crate::approx::{generalized_chebyshev_polynomials, PolynomialSet};
use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyMapping};
use crate::matrix::CouplingMatrix;
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterSpec;
use crate::synthesis::{MatrixSynthesisEngine, MatrixSynthesisMethod};

/// Identifies the approximation algorithm used during polynomial generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproximationKind {
    /// Generalized Chebyshev approximation with prescribed transmission zeros.
    GeneralizedChebyshev,
}

impl std::fmt::Display for ApproximationKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ApproximationKind::GeneralizedChebyshev => write!(f, "GeneralizedChebyshev"),
        }
    }
}

/// Shared synthesis data produced before optional response evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisOutcome {
    /// Approximation-stage polynomial data.
    pub polynomials: PolynomialSet,
    /// Coupling matrix synthesized from those polynomials.
    pub matrix: CouplingMatrix,
    /// Matrix-construction path used by the synthesis engine.
    pub matrix_method: MatrixSynthesisMethod,
    /// The approximation algorithm that produced the polynomials.
    pub approximation: ApproximationKind,
}

/// Result of the full synthesize-then-evaluate flow.
#[derive(Debug, Clone, PartialEq)]
pub struct EvaluationOutcome {
    /// Shared synthesis-stage output.
    pub synthesis: SynthesisOutcome,
    /// Sampled S-parameter response over the requested grid.
    pub response: SParameterResponse,
}

/// Canonical `FilterSpec` construction shared by every synthesis entry point
/// (high-level `FilterDesign`, the JSON pipeline, and the CLI). Centralizes the
/// order/return-loss/normalized-zeros/unloaded-Q assembly so the convention
/// lives in one place instead of being copy-pasted at each call site.
pub fn build_filter_spec(
    order: usize,
    return_loss_db: f64,
    normalized_transmission_zeros: &[f64],
    unloaded_q: Option<f64>,
) -> Result<FilterSpec> {
    let mut spec = FilterSpec::new(order, return_loss_db)?;
    if !normalized_transmission_zeros.is_empty() {
        spec =
            spec.with_normalized_transmission_zeros(normalized_transmission_zeros.iter().copied());
    }
    if let Some(q) = unloaded_q {
        spec = spec.with_unloaded_q(q);
    }
    Ok(spec)
}

/// Synthesizes using the default Chebyshev approximation flow.
pub fn synthesize_generalized_chebyshev(spec: &FilterSpec) -> Result<SynthesisOutcome> {
    let polynomials = generalized_chebyshev_polynomials(spec)?;
    let matrix_outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
    Ok(SynthesisOutcome {
        polynomials,
        matrix: matrix_outcome.matrix,
        matrix_method: matrix_outcome.method,
        approximation: ApproximationKind::GeneralizedChebyshev,
    })
}

/// Synthesizes a design and immediately evaluates its response over a physical grid.
pub fn synthesize_and_evaluate_with_mapping(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
    grid: &FrequencyGrid,
) -> Result<EvaluationOutcome> {
    let synthesis = synthesize_generalized_chebyshev(spec)?;
    let response = ResponseSolver.evaluate(&synthesis.matrix, grid, mapping)?;
    Ok(EvaluationOutcome {
        synthesis,
        response,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::load_filter_database_end_to_end_fixture;
    use crate::freq::BandPassMapping;
    use crate::synthesis::MatrixSynthesisEngine;

    #[test]
    fn orchestration_object_runs_full_flow() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_transmission_zeros(vec![crate::spec::TransmissionZero::normalized(2.0)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 9)?;

        let outcome = synthesize_and_evaluate_with_mapping(&spec, &mapping, &grid)?;
        assert_eq!(outcome.synthesis.polynomials.order, 4);
        assert_eq!(outcome.synthesis.matrix.order(), 4);
        assert_eq!(outcome.response.samples.len(), 9);
        assert_eq!(
            outcome.synthesis.approximation,
            ApproximationKind::GeneralizedChebyshev
        );
        Ok(())
    }

    #[test]
    fn canonical_matrix_synthesis_wraps_matrix_synthesizer() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn orchestration_reports_generalized_path_when_helper_data_is_used() -> Result<()> {
        let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;

        let outcome = synthesize_generalized_chebyshev(&fixture.spec)?;
        assert_eq!(
            outcome.approximation,
            ApproximationKind::GeneralizedChebyshev
        );
        assert_eq!(
            outcome.matrix_method,
            MatrixSynthesisMethod::ResidueExpansion
        );
        Ok(())
    }

    #[test]
    fn orchestration_supports_all_pole_case() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;

        let outcome = synthesize_generalized_chebyshev(&spec)?;
        assert_eq!(
            outcome.approximation,
            ApproximationKind::GeneralizedChebyshev
        );
        assert!(outcome.polynomials.generalized.is_some());
        Ok(())
    }
}
