//! Stage trait definition and concrete pipeline stage implementations.
//!
//! Each stage is a standalone callable unit with explicit input/output types.
//! Stages delegate to existing library implementations rather than reimplementing
//! core algorithms.

use super::context::MappingConfig;
use crate::approx::{generalized_chebyshev_polynomials, PolynomialSet};
use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::CouplingMatrix;
use crate::response::{ResponseSettings, ResponseSolver, SParameterResponse};
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
/// frequency grid for evaluation. When `mapping` is present the grid is
/// treated as physical frequencies and evaluated in the physical domain;
/// otherwise it is treated as a normalized-prototype grid.
#[derive(Debug, Clone)]
pub struct ResponseEvaluationInput {
    /// The coupling matrix to evaluate.
    pub matrix: CouplingMatrix,
    /// The frequency grid on which to evaluate the response.
    pub grid: FrequencyGrid,
    /// Optional physical-frequency mapping; when set, `grid` is physical Hz.
    pub mapping: Option<MappingConfig>,
    /// Physical unloaded Q (Hz) for a lossy response, or `None`/`∞` for
    /// lossless (review P2-8). Only meaningful with a band-pass mapping: a
    /// physical Q is incommensurate with a normalized-prototype grid, so the
    /// loss term is applied only in the physical domain.
    pub unloaded_q: Option<f64>,
}

/// Evaluates the S-parameter response of a coupling matrix on a frequency grid.
///
/// Delegates to [`ResponseSolver::evaluate_with_settings`] (physical grid,
/// mapped through the given [`MappingConfig`], lossy when a positive finite
/// `unloaded_q` is supplied) or [`ResponseSolver::evaluate_normalized`]
/// (already-normalized grid, always lossless — a physical Q is
/// incommensurate with a normalized-prototype grid).
#[derive(Debug, Default, Clone, Copy)]
pub struct ResponseEvaluationStage;

impl Stage for ResponseEvaluationStage {
    type Input = ResponseEvaluationInput;
    type Output = SParameterResponse;

    fn execute(&self, input: &Self::Input) -> Result<Self::Output> {
        match &input.mapping {
            Some(mapping) => {
                let settings = lossy_settings(input.unloaded_q, mapping)?;
                let mapping = mapping.build_mapping()?;
                ResponseSolver.evaluate_with_settings(&input.matrix, &input.grid, mapping.as_ref(), settings)
            }
            None => {
                validate_response_q(input.unloaded_q)?;
                ResponseSolver.evaluate_normalized(&input.matrix, &input.grid)
            }
        }
    }

    fn name(&self) -> &'static str {
        "response_evaluation"
    }
}

/// Validates a user-supplied physical unloaded Q for a lossy response.
///
/// The response backend silently drops a non-finite or non-positive Q (it
/// only applies the dissipation term when `is_finite() && > 0.0`), which
/// would return a *lossless* response for a "lossy" request with no warning
/// (the same class the P2-4 fix closed on `FilterDesign::response_lossy*`).
/// `f64::INFINITY` is the documented lossless sentinel and is allowed through;
/// everything else must be a finite positive value (review round-2 F5/F6).
fn validate_response_q(unloaded_q: Option<f64>) -> Result<()> {
    if let Some(q) = unloaded_q {
        if q == f64::INFINITY {
            return Ok(());
        }
        if !q.is_finite() || q <= 0.0 {
            return Err(crate::error::MfsError::InvalidFrequency(format!(
                "unloaded_q must be a positive finite value (or f64::INFINITY for \
                 lossless), got {q}"
            )));
        }
    }
    Ok(())
}

/// Builds the [`ResponseSettings`] for a physical unloaded Q.
///
/// A physical Q is incommensurate with a normalized-prototype grid, so the
/// loss term is applied only for a band-pass mapping, converting physical Q to
/// the normalized dissipation exactly as `FilterDesign::response_lossy` does
/// (`normalized_q = Q_u * bandwidth / center`). `None` and the lossless
/// `INFINITY` sentinel stay lossless; any other non-finite or non-positive Q
/// is rejected (it would otherwise silently degrade to a lossless response,
/// review round-2 F5). A low-pass mapping with a finite Q is likewise
/// rejected: the value is incommensurate with a low-pass grid and applying it
/// silently would drop it.
fn lossy_settings(
    unloaded_q: Option<f64>,
    mapping: &MappingConfig,
) -> Result<ResponseSettings> {
    validate_response_q(unloaded_q)?;
    let mut settings = ResponseSettings::default();
    if let Some(q) = unloaded_q {
        if q.is_finite() && q > 0.0 {
            if mapping.kind != "bandpass" {
                return Err(crate::error::MfsError::InvalidFrequency(
                    "a finite unloaded_q is only meaningful with a bandpass \
                     mapping (the normalized-prototype grid is incommensurate \
                     with a physical Q)".to_string(),
                ));
            }
            let center_hz = mapping.center_hz.ok_or_else(|| {
                crate::error::MfsError::InvalidFrequency(
                    "bandpass mapping requires center_hz for a lossy response".to_string(),
                )
            })?;
            let bandwidth_hz = mapping.bandwidth_hz.ok_or_else(|| {
                crate::error::MfsError::InvalidFrequency(
                    "bandpass mapping requires bandwidth_hz for a lossy response".to_string(),
                )
            })?;
            settings.unloaded_q = q * bandwidth_hz / center_hz;
        }
    }
    Ok(settings)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::PolynomialSetParams;
    use crate::error::MfsError;
    use crate::spec::FilterSpec;

    #[test]
    fn approximation_stage_produces_polynomials_from_spec() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);

        let stage = ApproximationStage;
        assert_eq!(stage.name(), "approximation");

        let polynomials = stage.execute(&spec)?;
        assert_eq!(polynomials.order, 4);
        assert!(polynomials.generalized.is_some());
        Ok(())
    }

    #[test]
    fn matrix_synthesis_stage_produces_matrix_from_polynomials() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
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
        let polynomials = PolynomialSet::new(&PolynomialSetParams {
            order: 3,
            ripple_factor: 0.1,
            eps: 0.1,
            eps_r: 1.0,
            transmission_zeros_normalized: vec![-1.5, 1.5],
            e: vec![1.0, 0.2, 0.3, 0.4],
            f: vec![0.8, 0.6, 0.4, 0.2],
            p: vec![1.0, 0.5, -2.25],
        })?;

        let stage = MatrixSynthesisStage;
        let error = stage.execute(&polynomials).unwrap_err();
        assert!(matches!(error, MfsError::PreconditionViolation(_)));
        Ok(())
    }

    #[test]
    fn topology_transform_stage_converts_to_folded() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 21)?;

        let stage = ResponseEvaluationStage;
        assert_eq!(stage.name(), "response_evaluation");

        let input = ResponseEvaluationInput {
            matrix,
            grid,
            mapping: None,
            unloaded_q: None,
        };
        let response = stage.execute(&input)?;
        assert_eq!(response.samples.len(), 21);
        // Response should vary across frequency
        assert_ne!(response.samples[0].s21_re, response.samples[10].s21_re);
        Ok(())
    }

    #[test]
    fn response_evaluation_stage_preserves_power_conservation() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        let grid = FrequencyGrid::linspace(-1.5, 1.5, 11)?;

        let stage = ResponseEvaluationStage;
        let input = ResponseEvaluationInput {
            matrix,
            grid,
            mapping: None,
            unloaded_q: None,
        };
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

    /// Regression for review P2-8: the pipeline previously stored `unloaded_q`
    /// in the spec but never consumed it — the response stage always ran
    /// lossless. A finite physical Q on a band-pass mapping must now add
    /// dissipation, so the peak in-band |S21| drops below the lossless value.
    #[test]
    fn response_evaluation_stage_applies_lossy_unloaded_q() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;

        let mapping = MappingConfig {
            kind: "bandpass".to_string(),
            center_hz: Some(6.75e9),
            bandwidth_hz: Some(300.0e6),
            cutoff_hz: None,
        };
        let grid = FrequencyGrid::linspace(6.6e9, 6.9e9, 51)?;

        let stage = ResponseEvaluationStage;
        let peak_s21 = |q: Option<f64>| {
            let input = ResponseEvaluationInput {
                matrix: matrix.clone(),
                grid: grid.clone(),
                mapping: Some(mapping.clone()),
                unloaded_q: q,
            };
            stage
                .execute(&input)
                .unwrap()
                .samples
                .iter()
                .map(|s| s.s21_mag())
                .fold(0.0_f64, f64::max)
        };

        let lossless = peak_s21(None);
        let lossy = peak_s21(Some(3000.0));
        assert!(
            lossless > 0.999,
            "lossless in-band |S21| should be near unity, got {lossless}"
        );
        assert!(
            lossy < lossless - 1e-4,
            "finite unloaded_q must add insertion loss (lossy {lossy} >= lossless {lossless})"
        );
        assert!(
            lossy.is_finite() && lossy > 0.0,
            "lossy in-band |S21| must stay finite and positive, got {lossy}"
        );
        Ok(())
    }

    #[test]
    fn lossy_settings_rejects_invalid_q() -> Result<()> {
        // A "lossy" request with a non-positive or non-finite Q (other than
        // the lossless INFINITY sentinel) must be rejected, not silently
        // degraded to a lossless response (review round-2 F5).
        let bandpass = MappingConfig {
            kind: "bandpass".to_string(),
            center_hz: Some(6.75e9),
            bandwidth_hz: Some(300.0e6),
            cutoff_hz: None,
        };
        for q in [0.0, -1.0, f64::NAN, f64::NEG_INFINITY] {
            let err = lossy_settings(Some(q), &bandpass).unwrap_err();
            assert!(
                err.to_string().contains("unloaded_q"),
                "q={q} should be rejected, got {err:?}"
            );
        }
        // The lossless sentinel and a valid finite Q stay accepted.
        assert!(lossy_settings(Some(f64::INFINITY), &bandpass).is_ok());
        assert!(lossy_settings(Some(3000.0), &bandpass).is_ok());
        assert!(lossy_settings(None, &bandpass).is_ok());
        Ok(())
    }

    #[test]
    fn lossy_settings_rejects_finite_q_with_lowpass() -> Result<()> {
        // A finite physical Q is incommensurate with a low-pass grid; applying
        // it would be silently dropped, so reject it (review round-2 F5).
        let lowpass = MappingConfig {
            kind: "lowpass".to_string(),
            center_hz: None,
            bandwidth_hz: None,
            cutoff_hz: Some(1.0e9),
        };
        let err = lossy_settings(Some(3000.0), &lowpass).unwrap_err();
        assert!(err.to_string().contains("bandpass"));
        Ok(())
    }

    #[test]
    fn response_stage_rejects_invalid_q_on_normalized_path() -> Result<()> {
        // The normalized (no-mapping) path evaluates lossless, but an invalid
        // Q must still be rejected rather than silently ignored (review
        // round-2 F6).
        let spec = FilterSpec::new(2, 20.0)?;
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 9)?;
        let input = ResponseEvaluationInput {
            matrix,
            grid,
            mapping: None,
            unloaded_q: Some(0.0),
        };
        let err = ResponseEvaluationStage.execute(&input).unwrap_err();
        assert!(err.to_string().contains("unloaded_q"));
        Ok(())
    }
}
