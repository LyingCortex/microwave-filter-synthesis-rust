use crate::approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet};
use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyPlan};
use crate::matrix::{CouplingMatrix, CouplingMatrixSynthesizer};
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterSpec;

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
        spec: &FilterSpec,
        plan: &impl FrequencyPlan,
    ) -> Result<SynthesisOutcome> {
        let polynomials = self.approximation.synthesize(spec, plan)?;
        let matrix = self.matrix_synthesizer.synthesize(&polynomials)?;
        Ok(SynthesisOutcome {
            polynomials,
            matrix,
        })
    }

    /// Synthesizes a design and immediately evaluates its response over a grid.
    pub fn synthesize_and_evaluate(
        &self,
        spec: &FilterSpec,
        plan: &impl FrequencyPlan,
        grid: &FrequencyGrid,
    ) -> Result<EvaluationOutcome> {
        let synthesis = self.synthesize(spec, plan)?;
        let response = self
            .response_solver
            .evaluate_with_plan(&synthesis.matrix, grid, plan)?;

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
    use crate::freq::BandPassPlan;
    use crate::spec::TransmissionZero;

    #[test]
    fn orchestration_object_runs_full_flow() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 9)?;

        let outcome = ChebyshevSynthesis::default().synthesize_and_evaluate(&spec, &plan, &grid)?;
        assert_eq!(outcome.polynomials.order, 4);
        assert_eq!(outcome.matrix.order(), 4);
        assert_eq!(outcome.response.samples.len(), 9);
        Ok(())
    }
}
