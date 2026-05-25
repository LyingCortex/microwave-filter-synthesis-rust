//! Stage trait definition and concrete pipeline stage implementations.
//!
//! Each stage is a standalone callable unit with explicit input/output types.
//! Stages delegate to existing library implementations rather than reimplementing
//! core algorithms.

use crate::approx::{generalized_chebyshev_polynomials, PolynomialSet};
use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::CouplingMatrix;
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterSpec;
use crate::synthesis::MatrixSynthesisEngine;
use crate::transform::{transform_matrix, TopologyKind, TransformOutcome};

/// A single pipeline stage with explicit input and output types.
///
/// Each stage accepts its declared input by reference and produces its declared
/// output type wrapped in a `Result`. Stages are composable, testable, and
/// replaceable.
pub trait Stage {
    /// The input type consumed by this stage.
    type Input;
    /// The output type produced by this stage.
    type Output;

    /// Executes the stage logic on the given input.
    fn execute(&self, input: &Self::Input) -> Result<Self::Output>;

    /// Returns a human-readable name identifying this stage.
    fn name(&self) -> &'static str;
}

// ---------------------------------------------------------------------------
// ApproximationStage
// ---------------------------------------------------------------------------

/// Generates prototype polynomials from a filter specification.
///
/// Delegates to [`generalized_chebyshev_polynomials`] which produces a
/// [`PolynomialSet`] containing the E, F, P polynomials and generalized
/// Chebyshev helper data.
#[derive(Debug, Default, Clone, Copy)]
pub struct ApproximationStage;

impl Stage for ApproximationStage {
    type Input = FilterSpec;
    type Output = PolynomialSet;

    fn execute(&self, input: &Self::Input) -> Result<Self::Output> {
        generalized_chebyshev_polynomials(input)
    }

    fn name(&self) -> &'static str {
        "approximation"
    }
}

// ---------------------------------------------------------------------------
// MatrixSynthesisStage
// ---------------------------------------------------------------------------

/// Synthesizes a canonical coupling matrix from prototype polynomials.
///
/// Delegates to [`MatrixSynthesisEngine::synthesize`] which performs residue
/// expansion to build the transversal coupling matrix.
#[derive(Debug, Default, Clone, Copy)]
pub struct MatrixSynthesisStage;

impl Stage for MatrixSynthesisStage {
    type Input = PolynomialSet;
    type Output = CouplingMatrix;

    fn execute(&self, input: &Self::Input) -> Result<Self::Output> {
        MatrixSynthesisEngine.synthesize(input)
    }

    fn name(&self) -> &'static str {
        "matrix_synthesis"
    }
}

// ---------------------------------------------------------------------------
// TopologyTransformStage
// ---------------------------------------------------------------------------

/// Input for the topology transform stage: a coupling matrix paired with the
/// desired output topology.
#[derive(Debug, Clone)]
pub struct TopologyTransformInput {
    /// The coupling matrix to transform.
    pub matrix: CouplingMatrix,
    /// The target topology for the transformation.
    pub topology: TopologyKind,
}

/// Transforms a coupling matrix into a requested topology.
///
/// Delegates to [`transform_matrix`] which applies similarity rotations to
/// achieve the target coupling pattern (folded, arrow, wheel, or identity for
/// transversal).
#[derive(Debug, Default, Clone, Copy)]
pub struct TopologyTransformStage;

impl Stage for TopologyTransformStage {
    type Input = TopologyTransformInput;
    type Output = TransformOutcome;

    fn execute(&self, input: &Self::Input) -> Result<Self::Output> {
        transform_matrix(&input.matrix, input.topology)
    }

    fn name(&self) -> &'static str {
        "topology_transform"
    }
}

// ---------------------------------------------------------------------------
// ResponseEvaluationStage
// ---------------------------------------------------------------------------

/// Input for the response evaluation stage: a coupling matrix paired with a
/// frequency grid for normalized evaluation.
#[derive(Debug, Clone)]
pub struct ResponseEvaluationInput {
    /// The coupling matrix to evaluate.
    pub matrix: CouplingMatrix,
    /// The frequency grid (normalized) on which to evaluate the response.
    pub grid: FrequencyGrid,
}

/// Evaluates the S-parameter response of a coupling matrix on a frequency grid.
///
/// Delegates to [`ResponseSolver::evaluate_normalized`] which computes S11 and
/// S21 at each grid point using the normalized prototype frequency axis.
#[derive(Debug, Default, Clone, Copy)]
pub struct ResponseEvaluationStage;

impl Stage for ResponseEvaluationStage {
    type Input = ResponseEvaluationInput;
    type Output = SParameterResponse;

    fn execute(&self, input: &Self::Input) -> Result<Self::Output> {
        ResponseSolver::default().evaluate_normalized(&input.matrix, &input.grid)
    }

    fn name(&self) -> &'static str {
        "response_evaluation"
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::MfsError;
    use crate::spec::FilterSpec;

    #[test]
    fn approximation_stage_produces_polynomials_from_spec() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);

        let stage = ApproximationStage;
        assert_eq!(stage.name(), "approximation");

        let polynomials = stage.execute(&spec)?;
        assert_eq!(polynomials.order, 4);
        assert!(polynomials.generalized.is_some());
        Ok(())
    }

    #[test]
    fn matrix_synthesis_stage_produces_matrix_from_polynomials() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let stage = MatrixSynthesisStage;
        assert_eq!(stage.name(), "matrix_synthesis");

        let matrix = stage.execute(&polynomials)?;
        assert_eq!(matrix.order(), 4);
        assert_eq!(matrix.side(), 6);
        // Source coupling must be non-zero
        assert!(matrix.at(0, 1).unwrap_or_default().abs() > 1e-12);
        Ok(())
    }

    #[test]
    fn matrix_synthesis_stage_returns_error_without_generalized_data() -> Result<()> {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![-1.5, 1.5],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, 0.5, -2.25],
        )?;

        let stage = MatrixSynthesisStage;
        let error = stage.execute(&polynomials).unwrap_err();
        assert!(matches!(error, MfsError::PreconditionViolation(_)));
        Ok(())
    }

    #[test]
    fn topology_transform_stage_converts_to_folded() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;

        let stage = TopologyTransformStage;
        assert_eq!(stage.name(), "topology_transform");

        let input = TopologyTransformInput {
            matrix,
            topology: TopologyKind::Folded,
        };
        let outcome = stage.execute(&input)?;
        assert_eq!(outcome.topology, TopologyKind::Folded);
        assert!(outcome.report.pattern_verified);
        Ok(())
    }

    #[test]
    fn topology_transform_stage_converts_to_arrow() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;

        let stage = TopologyTransformStage;
        let input = TopologyTransformInput {
            matrix,
            topology: TopologyKind::Arrow,
        };
        let outcome = stage.execute(&input)?;
        assert_eq!(outcome.topology, TopologyKind::Arrow);
        assert!(outcome.report.pattern_verified);
        Ok(())
    }

    #[test]
    fn response_evaluation_stage_produces_samples() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 21)?;

        let stage = ResponseEvaluationStage;
        assert_eq!(stage.name(), "response_evaluation");

        let input = ResponseEvaluationInput { matrix, grid };
        let response = stage.execute(&input)?;
        assert_eq!(response.samples.len(), 21);
        // Response should vary across frequency
        assert_ne!(response.samples[0].s21_re, response.samples[10].s21_re);
        Ok(())
    }

    #[test]
    fn response_evaluation_stage_preserves_power_conservation() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        let grid = FrequencyGrid::linspace(-1.5, 1.5, 11)?;

        let stage = ResponseEvaluationStage;
        let input = ResponseEvaluationInput { matrix, grid };
        let response = stage.execute(&input)?;

        for sample in &response.samples {
            let s11_mag_sq = sample.s11_re.powi(2) + sample.s11_im.powi(2);
            let s21_mag_sq = sample.s21_re.powi(2) + sample.s21_im.powi(2);
            let power_sum = s11_mag_sq + s21_mag_sq;
            assert!(
                (power_sum - 1.0).abs() < 1e-9,
                "power conservation violated: |S11|² + |S21|² = {power_sum}"
            );
        }
        Ok(())
    }
}
