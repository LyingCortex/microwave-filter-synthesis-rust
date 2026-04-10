use crate::approx::{
    ApproximationEngine, ChebyshevApproximation, GeneralizedChebyshevApproximation, PolynomialSet,
};
use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyMapping};
use crate::matrix::CouplingMatrix;
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::{ApproximationFamily, FilterSpec};
use crate::synthesis::{MatrixSynthesisEngine, MatrixSynthesisMethod};

/// Indicates whether approximation output included generalized Chebyshev helper data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproximationStageKind {
    /// Baseline Chebyshev polynomial bundle without generalized helper attachments.
    ClassicalChebyshev,
    /// Polynomial bundle augmented with generalized Chebyshev helper data.
    GeneralizedChebyshev,
}

/// Result of the synthesis stage before response evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisOutcome {
    /// Approximation-stage polynomial data.
    pub polynomials: PolynomialSet,
    /// Coupling matrix synthesized from those polynomials.
    pub matrix: CouplingMatrix,
    /// Approximation-stage path reflected in the polynomial bundle.
    pub approximation_kind: ApproximationStageKind,
    /// Matrix-construction path used by the synthesis engine.
    pub matrix_method: MatrixSynthesisMethod,
}

/// Result of the full synthesize-then-evaluate flow.
#[derive(Debug, Clone, PartialEq)]
pub struct EvaluationOutcome {
    /// Approximation-stage polynomial data.
    pub polynomials: PolynomialSet,
    /// Coupling matrix synthesized from those polynomials.
    pub matrix: CouplingMatrix,
    /// Approximation-stage path reflected in the polynomial bundle.
    pub approximation_kind: ApproximationStageKind,
    /// Matrix-construction path used by the synthesis engine.
    pub matrix_method: MatrixSynthesisMethod,
    /// Sampled S-parameter response over the requested grid.
    pub response: SParameterResponse,
}

/// High-level orchestrator that ties together approximation, matrix synthesis, and response solving.
#[derive(Debug, Default, Clone, Copy)]
pub struct ChebyshevSynthesis {
    approximation: ChebyshevApproximation,
    generalized_approximation: GeneralizedChebyshevApproximation,
    matrix_synthesizer: MatrixSynthesisEngine,
    response_solver: ResponseSolver,
}

impl ChebyshevSynthesis {
    /// Synthesizes prototype polynomials and a coupling matrix, exposing the chosen internal path.
    pub fn synthesize_with_details(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
    ) -> Result<SynthesisOutcome> {
        let polynomials = match spec.approximation_family {
            ApproximationFamily::GeneralizedChebyshev => {
                self.generalized_approximation.synthesize(spec, mapping)?
            }
            ApproximationFamily::Chebyshev => self.approximation.synthesize(spec, mapping)?,
        };
        let approximation_kind = if polynomials.generalized.is_some() {
            ApproximationStageKind::GeneralizedChebyshev
        } else {
            ApproximationStageKind::ClassicalChebyshev
        };
        let matrix_outcome = self.matrix_synthesizer.synthesize_with_details(&polynomials)?;
        Ok(SynthesisOutcome {
            polynomials,
            matrix: matrix_outcome.matrix,
            approximation_kind,
            matrix_method: matrix_outcome.method,
        })
    }

    /// Synthesizes prototype polynomials and a coupling matrix from the input spec.
    pub fn synthesize(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
    ) -> Result<SynthesisOutcome> {
        self.synthesize_with_details(spec, mapping)
    }

    /// Synthesizes a design and immediately evaluates its response over a physical grid.
    pub fn synthesize_and_evaluate_with_mapping(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
        grid: &FrequencyGrid,
    ) -> Result<EvaluationOutcome> {
        let synthesis = self.synthesize(spec, mapping)?;
        let response = self.response_solver.evaluate(&synthesis.matrix, grid, mapping)?;

        Ok(EvaluationOutcome {
            polynomials: synthesis.polynomials,
            matrix: synthesis.matrix,
            approximation_kind: synthesis.approximation_kind,
            matrix_method: synthesis.matrix_method,
            response,
        })
    }

    /// Synthesizes and evaluates a design while preserving approximation and matrix-path metadata.
    pub fn synthesize_and_evaluate_with_details(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
        grid: &FrequencyGrid,
    ) -> Result<EvaluationOutcome> {
        self.synthesize_and_evaluate_with_mapping(spec, mapping, grid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::load_filter_database_end_to_end_fixture;
    use crate::freq::BandPassMapping;
    use crate::spec::TransmissionZero;
    use crate::synthesis::CanonicalMatrixSynthesis;

    #[test]
    fn orchestration_object_runs_full_flow() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 9)?;

        let outcome = ChebyshevSynthesis::default()
            .synthesize_and_evaluate_with_mapping(&spec, &mapping, &grid)?;
        assert_eq!(outcome.polynomials.order, 4);
        assert_eq!(outcome.matrix.order(), 4);
        assert_eq!(outcome.response.samples.len(), 9);
        assert_eq!(outcome.approximation_kind, ApproximationStageKind::ClassicalChebyshev);
        Ok(())
    }

    #[test]
    fn canonical_matrix_synthesis_wraps_matrix_synthesizer() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::normalized(-2.0)]);
        let mapping = crate::freq::LowPassMapping::new(1.0)?;
        let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let outcome = CanonicalMatrixSynthesis::default().synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn orchestration_reports_generalized_path_when_helper_data_is_used() -> Result<()> {
        let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;

        let outcome = ChebyshevSynthesis::default()
            .synthesize_with_details(&fixture.spec, &fixture.mapping)?;
        assert_eq!(outcome.approximation_kind, ApproximationStageKind::GeneralizedChebyshev);
        assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::ResidueExpansion);
        Ok(())
    }

    #[test]
    fn orchestration_reports_fallback_path_for_all_pole_case() -> Result<()> {
        let spec = FilterSpec::chebyshev(3, 20.0)?;
        let mapping = crate::freq::LowPassMapping::new(1.0)?;

        let outcome = ChebyshevSynthesis::default().synthesize_with_details(&spec, &mapping)?;
        assert_eq!(outcome.approximation_kind, ApproximationStageKind::ClassicalChebyshev);
        assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::PlaceholderFallback);
        Ok(())
    }

    #[test]
    fn orchestration_allows_explicit_generalized_family_without_finite_zeros() -> Result<()> {
        let spec = FilterSpec::generalized_chebyshev(3, 20.0)?;
        let mapping = crate::freq::LowPassMapping::new(1.0)?;

        let outcome = ChebyshevSynthesis::default().synthesize_with_details(&spec, &mapping)?;
        assert_eq!(outcome.approximation_kind, ApproximationStageKind::GeneralizedChebyshev);
        assert!(outcome.polynomials.generalized.is_some());
        Ok(())
    }
}
