use crate::approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet};
use crate::error::Result;
use crate::freq::{FrequencyGrid, FrequencyPlan};
use crate::matrix::{CouplingMatrix, CouplingMatrixSynthesizer};
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterSpec;

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisOutcome {
    pub polynomials: PolynomialSet,
    pub matrix: CouplingMatrix,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EvaluationOutcome {
    pub polynomials: PolynomialSet,
    pub matrix: CouplingMatrix,
    pub response: SParameterResponse,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct ChebyshevSynthesis {
    approximation: ChebyshevApproximation,
    matrix_synthesizer: CouplingMatrixSynthesizer,
    response_solver: ResponseSolver,
}

impl ChebyshevSynthesis {
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
