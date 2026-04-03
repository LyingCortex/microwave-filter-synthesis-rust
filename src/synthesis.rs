use crate::approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet};
use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyMapping};
use crate::matrix::{CouplingMatrix, CouplingMatrixSynthesizer};
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterParameter;

/// Result of the synthesis stage before response evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisOutcome {
    /// Approximation-stage polynomial data.
    pub polynomials: PolynomialSet,
    /// Coupling matrix synthesized from those polynomials.
    pub matrix: CouplingMatrix,
}

/// Result of the full synthesize-then-evaluate flow.
#[derive(Debug, Clone, PartialEq)]
pub struct EvaluationOutcome {
    /// Approximation-stage polynomial data.
    pub polynomials: PolynomialSet,
    /// Coupling matrix synthesized from those polynomials.
    pub matrix: CouplingMatrix,
    /// Sampled S-parameter response over the requested grid.
    pub response: SParameterResponse,
}

/// High-level orchestrator that ties together approximation, matrix synthesis, and response solving.
#[derive(Debug, Default, Clone, Copy)]
pub struct ChebyshevSynthesis {
    approximation: ChebyshevApproximation,
    matrix_synthesizer: CouplingMatrixSynthesizer,
    response_solver: ResponseSolver,
}

impl ChebyshevSynthesis {
    /// Synthesizes prototype polynomials and a coupling matrix from the input spec.
    pub fn synthesize(
        &self,
        spec: &FilterParameter,
        mapping: &impl FrequencyMapping,
    ) -> Result<SynthesisOutcome> {
        let polynomials = self.approximation.synthesize(spec, mapping)?;
        let matrix = self.matrix_synthesizer.synthesize(&polynomials)?;
        Ok(SynthesisOutcome {
            polynomials,
            matrix,
        })
    }

    /// Synthesizes a design and immediately evaluates its response over a physical grid.
    pub fn synthesize_and_evaluate_with_mapping(
        &self,
        spec: &FilterParameter,
        mapping: &impl FrequencyMapping,
        grid: &FrequencyGrid,
    ) -> Result<EvaluationOutcome> {
        let synthesis = self.synthesize(spec, mapping)?;
        let response = self.response_solver.evaluate(&synthesis.matrix, grid, mapping)?;

        Ok(EvaluationOutcome {
            polynomials: synthesis.polynomials,
            matrix: synthesis.matrix,
            response,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::BandPassMapping;
    use crate::spec::TransmissionZero;

    #[test]
    fn orchestration_object_runs_full_flow() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 9)?;

        let outcome = ChebyshevSynthesis::default()
            .synthesize_and_evaluate_with_mapping(&spec, &mapping, &grid)?;
        assert_eq!(outcome.polynomials.order, 4);
        assert_eq!(outcome.matrix.order(), 4);
        assert_eq!(outcome.response.samples.len(), 9);
        Ok(())
    }
}
