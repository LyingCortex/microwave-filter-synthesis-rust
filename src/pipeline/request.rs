//! JSON-based synthesis request handling and full pipeline execution.
//!
//! This module provides:
//! - [`SynthesisRequest`]: a JSON-serializable input format for the synthesis pipeline
//! - [`run_full_pipeline`]: executes the complete pipeline from a request
//! - [`run_from_json`]: parses JSON, executes, and serializes the result
//! - [`ValidationError`]: structured field-level validation errors

use std::time::Instant;

use serde::{Deserialize, Serialize};

use crate::error::{MfsError, Result};
use crate::freq::FrequencyGrid;
use crate::spec::FilterSpec;
use crate::transform::TopologyKind;

use super::context::{GridConfig, MappingConfig, PipelineOptions, SynthesisContext};
use super::stages::{
    ApproximationStage, MatrixSynthesisStage, ResponseEvaluationInput, ResponseEvaluationStage,
    Stage, TopologyTransformInput, TopologyTransformStage,
};

/// A structured field-level validation error.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ValidationError {
    /// The field path that failed validation (e.g. "order", "grid.points").
    pub field: String,
    /// A human-readable description of the validation failure.
    pub message: String,
}

/// A collection of validation errors returned when input is malformed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ValidationErrors {
    /// Individual field-level errors.
    pub errors: Vec<ValidationError>,
}

impl std::fmt::Display for ValidationErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "validation failed: ")?;
        for (i, err) in self.errors.iter().enumerate() {
            if i > 0 {
                write!(f, "; ")?;
            }
            write!(f, "{}: {}", err.field, err.message)?;
        }
        Ok(())
    }
}

/// JSON input format for the synthesis pipeline.
///
/// All fields except `order` and `return_loss_db` are optional with sensible
/// defaults. A minimal request needs only the filter order and return loss:
///
/// ```json
/// { "order": 4, "return_loss_db": 20.0 }
/// ```
///
/// A full request can specify transmission zeros, topology, frequency mapping,
/// and evaluation grid:
///
/// ```json
/// {
///   "order": 4,
///   "return_loss_db": 20.0,
///   "transmission_zeros": [-2.0, 1.5],
///   "unloaded_q": 2000.0,
///   "topology": "Folded",
///   "mapping": { "kind": "bandpass", "center_hz": 6.75e9, "bandwidth_hz": 3e8 },
///   "grid": { "start": 6.0e9, "stop": 7.5e9, "points": 201 }
/// }
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SynthesisRequest {
    /// Filter order (number of resonators). Must be > 0.
    pub order: usize,
    /// Minimum passband return loss in dB. Must be > 0 and finite.
    pub return_loss_db: f64,
    /// Normalized transmission zeros (prototype coordinates). Optional.
    #[serde(default)]
    pub transmission_zeros: Vec<f64>,
    /// Unloaded Q factor. Optional.
    pub unloaded_q: Option<f64>,
    /// Requested output topology for the matrix transform stage.
    pub topology: Option<TopologyKind>,
    /// Frequency mapping configuration for physical-frequency evaluation.
    pub mapping: Option<MappingConfig>,
    /// Frequency grid configuration for response evaluation.
    pub grid: Option<GridConfig>,
}

impl SynthesisRequest {
    /// Validates the request and returns structured field-level errors if invalid.
    pub fn validate(&self) -> std::result::Result<(), ValidationErrors> {
        let mut errors = Vec::new();

        if self.order == 0 {
            errors.push(ValidationError {
                field: "order".to_string(),
                message: "must be greater than 0".to_string(),
            });
        }

        if self.return_loss_db <= 0.0 || !self.return_loss_db.is_finite() {
            errors.push(ValidationError {
                field: "return_loss_db".to_string(),
                message: "must be greater than 0 and finite".to_string(),
            });
        }

        if let Some(ref grid) = self.grid {
            if grid.points < 2 {
                errors.push(ValidationError {
                    field: "grid.points".to_string(),
                    message: "must be at least 2".to_string(),
                });
            }
            if !grid.start.is_finite() || !grid.stop.is_finite() || grid.stop <= grid.start {
                errors.push(ValidationError {
                    field: "grid.start/grid.stop".to_string(),
                    message: "start must be less than stop, both must be finite".to_string(),
                });
            }
        }

        if errors.is_empty() {
            Ok(())
        } else {
            Err(ValidationErrors { errors })
        }
    }
}

/// Executes the full synthesis pipeline from a validated request.
///
/// The pipeline runs stages in sequence:
/// 1. Approximation (always)
/// 2. Matrix synthesis (always)
/// 3. Topology transform (if `topology` is specified)
/// 4. Response evaluation (if `grid` is specified)
///
/// Each stage's execution time is tracked in the context metadata.
pub fn run_full_pipeline(request: SynthesisRequest) -> Result<SynthesisContext> {
    // Validate the request, converting validation errors to MfsError
    request.validate().map_err(|ve| {
        MfsError::PreconditionViolation(ve.to_string())
    })?;

    // Build FilterSpec from request
    let mut spec = FilterSpec::new(request.order, request.return_loss_db)?;
    if !request.transmission_zeros.is_empty() {
        spec = spec.with_normalized_transmission_zeros(request.transmission_zeros.iter().copied());
    }
    if let Some(q) = request.unloaded_q {
        spec = spec.with_unloaded_q(q);
    }

    // Build PipelineOptions from request
    let options = PipelineOptions {
        topology: request.topology,
        grid: request.grid.clone(),
        mapping: request.mapping.clone(),
        response_tolerance: None,
    };

    // Create context
    let mut ctx = SynthesisContext::with_options(spec, options);

    // Stage 1: Approximation
    {
        let start = Instant::now();
        let stage = ApproximationStage;
        let polynomials = stage.execute(&ctx.spec)?;
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        ctx.polynomials = Some(std::sync::Arc::new(polynomials));
        ctx.metadata.stages_executed.push(stage.name().to_string());
        ctx.metadata.stage_timings_ms.push(elapsed_ms);
    }

    // Stage 2: Matrix synthesis
    {
        let start = Instant::now();
        let stage = MatrixSynthesisStage;
        let polynomials = ctx.polynomials().ok_or_else(|| {
            MfsError::PreconditionViolation(
                "approximation stage must complete before matrix synthesis".to_string(),
            )
        })?;
        let matrix = stage.execute(polynomials)?;
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        ctx.matrix = Some(std::sync::Arc::new(matrix));
        ctx.metadata.stages_executed.push(stage.name().to_string());
        ctx.metadata.stage_timings_ms.push(elapsed_ms);
    }

    // Stage 3: Topology transform (optional)
    if let Some(topology) = request.topology {
        let start = Instant::now();
        let stage = TopologyTransformStage;
        let matrix = ctx.matrix().ok_or_else(|| {
            MfsError::PreconditionViolation(
                "matrix synthesis must complete before topology transform".to_string(),
            )
        })?;
        let input = TopologyTransformInput {
            matrix: matrix.clone(),
            topology,
        };
        let outcome = stage.execute(&input)?;
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        ctx.transform = Some(outcome);
        ctx.metadata.stages_executed.push(stage.name().to_string());
        ctx.metadata.stage_timings_ms.push(elapsed_ms);
    }

    // Stage 4: Response evaluation (optional, requires grid)
    if let Some(ref grid_config) = request.grid {
        let start = Instant::now();
        let stage = ResponseEvaluationStage;
        // Use the transformed matrix if available, otherwise the synthesis matrix
        let eval_matrix = if let Some(ref transform) = ctx.transform {
            transform.matrix.clone()
        } else {
            ctx.matrix()
                .ok_or_else(|| {
                    MfsError::PreconditionViolation(
                        "matrix must be available for response evaluation".to_string(),
                    )
                })?
                .clone()
        };
        let grid = FrequencyGrid::linspace(grid_config.start, grid_config.stop, grid_config.points)?;
        let input = ResponseEvaluationInput {
            matrix: eval_matrix,
            grid,
        };
        let response = stage.execute(&input)?;
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        ctx.response = Some(response);
        ctx.metadata.stages_executed.push(stage.name().to_string());
        ctx.metadata.stage_timings_ms.push(elapsed_ms);
    }

    Ok(ctx)
}

/// Parses a JSON string into a [`SynthesisRequest`], executes the full pipeline,
/// and serializes the resulting context back to JSON.
///
/// Returns the JSON string representation of the completed [`SynthesisContext`].
pub fn run_from_json(json: &str) -> Result<String> {
    let request: SynthesisRequest = serde_json::from_str(json).map_err(|e| {
        MfsError::PreconditionViolation(format!("failed to parse JSON request: {e}"))
    })?;

    let ctx = run_full_pipeline(request)?;

    serde_json::to_string(&ctx).map_err(|e| {
        MfsError::NumericalFailure(format!("failed to serialize pipeline result: {e}"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimal_request_runs_approximation_and_matrix() -> Result<()> {
        let request = SynthesisRequest {
            order: 4,
            return_loss_db: 20.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: None,
        };

        let ctx = run_full_pipeline(request)?;
        assert!(ctx.polynomials().is_some());
        assert!(ctx.matrix().is_some());
        assert!(ctx.transform().is_none());
        assert!(ctx.response().is_none());
        assert_eq!(ctx.metadata.stages_executed.len(), 2);
        assert_eq!(ctx.metadata.stages_executed[0], "approximation");
        assert_eq!(ctx.metadata.stages_executed[1], "matrix_synthesis");
        Ok(())
    }

    #[test]
    fn full_request_runs_all_stages() -> Result<()> {
        let request = SynthesisRequest {
            order: 4,
            return_loss_db: 20.0,
            transmission_zeros: vec![-2.0, 1.5],
            unloaded_q: Some(2000.0),
            topology: Some(TopologyKind::Folded),
            mapping: Some(MappingConfig {
                kind: "bandpass".to_string(),
                center_hz: Some(6.75e9),
                bandwidth_hz: Some(300.0e6),
                cutoff_hz: None,
            }),
            grid: Some(GridConfig {
                start: -3.0,
                stop: 3.0,
                points: 201,
            }),
        };

        let ctx = run_full_pipeline(request)?;
        assert!(ctx.polynomials().is_some());
        assert!(ctx.matrix().is_some());
        assert!(ctx.transform().is_some());
        assert!(ctx.response().is_some());
        assert_eq!(ctx.metadata.stages_executed.len(), 4);
        assert_eq!(ctx.metadata.stages_executed[0], "approximation");
        assert_eq!(ctx.metadata.stages_executed[1], "matrix_synthesis");
        assert_eq!(ctx.metadata.stages_executed[2], "topology_transform");
        assert_eq!(ctx.metadata.stages_executed[3], "response_evaluation");
        // Verify timing data is populated
        assert_eq!(ctx.metadata.stage_timings_ms.len(), 4);
        for timing in &ctx.metadata.stage_timings_ms {
            assert!(*timing >= 0.0);
        }
        Ok(())
    }

    #[test]
    fn validation_rejects_zero_order() {
        let request = SynthesisRequest {
            order: 0,
            return_loss_db: 20.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: None,
        };

        let err = request.validate().unwrap_err();
        assert_eq!(err.errors.len(), 1);
        assert_eq!(err.errors[0].field, "order");
    }

    #[test]
    fn validation_rejects_negative_return_loss() {
        let request = SynthesisRequest {
            order: 4,
            return_loss_db: -5.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: None,
        };

        let err = request.validate().unwrap_err();
        assert_eq!(err.errors.len(), 1);
        assert_eq!(err.errors[0].field, "return_loss_db");
    }

    #[test]
    fn validation_rejects_infinite_return_loss() {
        let request = SynthesisRequest {
            order: 4,
            return_loss_db: f64::INFINITY,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: None,
        };

        let err = request.validate().unwrap_err();
        assert_eq!(err.errors.len(), 1);
        assert_eq!(err.errors[0].field, "return_loss_db");
    }

    #[test]
    fn validation_rejects_grid_with_less_than_2_points() {
        let request = SynthesisRequest {
            order: 4,
            return_loss_db: 20.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: Some(GridConfig {
                start: -2.0,
                stop: 2.0,
                points: 1,
            }),
        };

        let err = request.validate().unwrap_err();
        assert!(err.errors.iter().any(|e| e.field == "grid.points"));
    }

    #[test]
    fn validation_collects_multiple_errors() {
        let request = SynthesisRequest {
            order: 0,
            return_loss_db: -1.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: Some(GridConfig {
                start: -2.0,
                stop: 2.0,
                points: 0,
            }),
        };

        let err = request.validate().unwrap_err();
        assert!(err.errors.len() >= 3);
    }

    #[test]
    fn run_full_pipeline_rejects_invalid_request() {
        let request = SynthesisRequest {
            order: 0,
            return_loss_db: 20.0,
            transmission_zeros: vec![],
            unloaded_q: None,
            topology: None,
            mapping: None,
            grid: None,
        };

        let err = run_full_pipeline(request).unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
    }

    #[test]
    fn run_from_json_minimal_request() -> Result<()> {
        let json = r#"{"order": 4, "return_loss_db": 20.0, "transmission_zeros": [-1.5, 1.5]}"#;
        let result = run_from_json(json)?;
        // Result should be valid JSON
        let _: serde_json::Value = serde_json::from_str(&result).expect("output should be valid JSON");
        Ok(())
    }

    #[test]
    fn run_from_json_full_request() -> Result<()> {
        let json = r#"{
            "order": 4,
            "return_loss_db": 20.0,
            "transmission_zeros": [-2.0, 1.5],
            "unloaded_q": 2000.0,
            "topology": "Folded",
            "grid": { "start": -3.0, "stop": 3.0, "points": 51 }
        }"#;
        let result = run_from_json(json)?;
        let value: serde_json::Value = serde_json::from_str(&result).expect("output should be valid JSON");
        // Check metadata is present
        assert!(value.get("metadata").is_some());
        Ok(())
    }

    #[test]
    fn run_from_json_rejects_malformed_json() {
        let json = r#"{ not valid json }"#;
        let err = run_from_json(json).unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
    }

    #[test]
    fn synthesis_request_deserializes_from_minimal_json() {
        let json = r#"{"order": 4, "return_loss_db": 20.0}"#;
        let request: SynthesisRequest = serde_json::from_str(json).expect("should parse");
        assert_eq!(request.order, 4);
        assert_eq!(request.return_loss_db, 20.0);
        assert!(request.transmission_zeros.is_empty());
        assert!(request.topology.is_none());
        assert!(request.grid.is_none());
        assert!(request.mapping.is_none());
        assert!(request.unloaded_q.is_none());
    }

    #[test]
    fn synthesis_request_deserializes_from_full_json() {
        let json = r#"{
            "order": 4,
            "return_loss_db": 20.0,
            "transmission_zeros": [-2.0, 1.5],
            "unloaded_q": 2000.0,
            "topology": "Folded",
            "mapping": {
                "kind": "bandpass",
                "center_hz": 6.75e9,
                "bandwidth_hz": 3e8
            },
            "grid": {
                "start": 6.0e9,
                "stop": 7.5e9,
                "points": 201
            }
        }"#;
        let request: SynthesisRequest = serde_json::from_str(json).expect("should parse");
        assert_eq!(request.order, 4);
        assert_eq!(request.return_loss_db, 20.0);
        assert_eq!(request.transmission_zeros, vec![-2.0, 1.5]);
        assert_eq!(request.unloaded_q, Some(2000.0));
        assert_eq!(request.topology, Some(TopologyKind::Folded));
        assert_eq!(request.mapping.as_ref().unwrap().kind, "bandpass");
        assert_eq!(request.grid.as_ref().unwrap().points, 201);
    }
}
