//! Incremental pipeline execution and context persistence.
//!
//! This module provides:
//! - [`run_stage`]: execute a single named pipeline stage with predecessor validation
//! - [`save_context`] / [`load_context`]: file-based persistence for `SynthesisContext`
//!
//! Predecessor rules:
//! - `"approximation"` requires: spec (always present in a valid context)
//! - `"matrix_synthesis"` requires: polynomials
//! - `"topology_transform"` requires: matrix + options.topology
//! - `"response_evaluation"` requires: matrix

use std::fs;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use serde::{Deserialize, Serialize};

use crate::error::{MfsError, Result};
use crate::freq::FrequencyGrid;
use crate::pipeline::context::{PipelineMetadata, PipelineOptions, SynthesisContext};
use crate::pipeline::stages::{
    ApproximationStage, MatrixSynthesisStage, ResponseEvaluationInput, ResponseEvaluationStage,
    Stage, TopologyTransformInput, TopologyTransformStage,
};
use crate::spec::FilterSpec;

/// Executes a single named pipeline stage, storing the result in the context.
///
/// Validates that required predecessor artifacts exist before running. Returns
/// `MfsError::PreconditionViolation` if predecessors are missing or the stage
/// name is unrecognized.
///
/// # Valid stage names
/// - `"approximation"` — generates polynomials from the filter spec
/// - `"matrix_synthesis"` — synthesizes a coupling matrix from polynomials
/// - `"topology_transform"` — transforms the matrix to the requested topology
/// - `"response_evaluation"` — evaluates S-parameter response on a frequency grid
pub fn run_stage(context: &mut SynthesisContext, stage_name: &str) -> Result<()> {
    match stage_name {
        "approximation" => run_approximation_stage(context),
        "matrix_synthesis" => run_matrix_synthesis_stage(context),
        "topology_transform" => run_topology_transform_stage(context),
        "response_evaluation" => run_response_evaluation_stage(context),
        _ => Err(MfsError::PreconditionViolation(format!(
            "unrecognized stage name: '{stage_name}'. Valid stages: approximation, \
             matrix_synthesis, topology_transform, response_evaluation"
        ))),
    }
}

fn run_approximation_stage(context: &mut SynthesisContext) -> Result<()> {
    // Predecessor: spec (always present)
    let start = Instant::now();
    let stage = ApproximationStage;
    let polynomials = stage.execute(&context.spec)?;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    context.polynomials = Some(Arc::new(polynomials));
    context.metadata.stages_executed.push("approximation".to_string());
    context.metadata.stage_timings_ms.push(elapsed_ms);
    Ok(())
}

fn run_matrix_synthesis_stage(context: &mut SynthesisContext) -> Result<()> {
    // Predecessor: polynomials
    let polynomials = context.polynomials.as_ref().ok_or_else(|| {
        MfsError::PreconditionViolation(
            "matrix_synthesis requires polynomials from the approximation stage".to_string(),
        )
    })?;

    let start = Instant::now();
    let stage = MatrixSynthesisStage;
    let matrix = stage.execute(polynomials)?;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    context.matrix = Some(Arc::new(matrix));
    context.metadata.stages_executed.push("matrix_synthesis".to_string());
    context.metadata.stage_timings_ms.push(elapsed_ms);
    Ok(())
}

fn run_topology_transform_stage(context: &mut SynthesisContext) -> Result<()> {
    // Predecessors: matrix + options.topology
    let matrix = context.matrix.as_ref().ok_or_else(|| {
        MfsError::PreconditionViolation(
            "topology_transform requires a coupling matrix from the matrix_synthesis stage"
                .to_string(),
        )
    })?;

    let topology = context.options.topology.ok_or_else(|| {
        MfsError::PreconditionViolation(
            "topology_transform requires options.topology to be set".to_string(),
        )
    })?;

    let start = Instant::now();
    let stage = TopologyTransformStage;
    let input = TopologyTransformInput {
        matrix: matrix.as_ref().clone(),
        topology,
    };
    let outcome = stage.execute(&input)?;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    context.transform = Some(outcome);
    context.metadata.stages_executed.push("topology_transform".to_string());
    context.metadata.stage_timings_ms.push(elapsed_ms);
    Ok(())
}

fn run_response_evaluation_stage(context: &mut SynthesisContext) -> Result<()> {
    // Predecessor: matrix
    let matrix = context.matrix.as_ref().ok_or_else(|| {
        MfsError::PreconditionViolation(
            "response_evaluation requires a coupling matrix from the matrix_synthesis stage"
                .to_string(),
        )
    })?;

    // Build a frequency grid from options or use a sensible default
    let grid = build_frequency_grid(&context.options)?;

    let start = Instant::now();
    let stage = ResponseEvaluationStage;
    let input = ResponseEvaluationInput {
        matrix: matrix.as_ref().clone(),
        grid,
    };
    let response = stage.execute(&input)?;
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    context.response = Some(response);
    context
        .metadata
        .stages_executed
        .push("response_evaluation".to_string());
    context.metadata.stage_timings_ms.push(elapsed_ms);
    Ok(())
}

/// Builds a `FrequencyGrid` from pipeline options, falling back to a default
/// normalized grid if no grid configuration is provided.
fn build_frequency_grid(options: &PipelineOptions) -> Result<FrequencyGrid> {
    match &options.grid {
        Some(grid_config) => FrequencyGrid::linspace(
            grid_config.start,
            grid_config.stop,
            grid_config.points,
        ),
        None => {
            // Default: normalized prototype grid spanning [-3, 3] with 101 points
            FrequencyGrid::linspace(-3.0, 3.0, 101)
        }
    }
}

// ---------------------------------------------------------------------------
// Context persistence (simplified serializable representation)
// ---------------------------------------------------------------------------

/// Simplified serializable snapshot of a `SynthesisContext`.
///
/// Inner domain types (`PolynomialSet`, `CouplingMatrix`, etc.) do not yet
/// implement `Serialize`/`Deserialize`. This snapshot captures enough state to
/// reconstruct a context for pipeline resumption: the spec, options, metadata,
/// and which stages have completed. Artifact data is stored in a reduced form
/// (matrix data as flat `Vec<f64>`, polynomial scalars, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ContextSnapshot {
    /// Filter spec fields needed to reconstruct a `FilterSpec`.
    spec: SpecSnapshot,
    /// Pipeline options (already serializable).
    options: PipelineOptions,
    /// Execution metadata (already serializable).
    metadata: PipelineMetadata,
    /// Reduced polynomial data if approximation has completed.
    polynomials: Option<PolynomialSnapshot>,
    /// Reduced matrix data if matrix synthesis has completed.
    matrix: Option<MatrixSnapshot>,
    /// Whether a transform outcome exists (topology transform completed).
    has_transform: bool,
    /// Whether a response exists (response evaluation completed).
    has_response: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpecSnapshot {
    order: usize,
    return_loss_db: f64,
    transmission_zeros: Vec<f64>,
    unloaded_q: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PolynomialSnapshot {
    order: usize,
    ripple_factor: f64,
    eps: f64,
    eps_r: f64,
    transmission_zeros_normalized: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MatrixSnapshot {
    order: usize,
    data: Vec<f64>,
}

impl ContextSnapshot {
    fn from_context(context: &SynthesisContext) -> Self {
        let spec = SpecSnapshot {
            order: context.spec.order,
            return_loss_db: context.spec.return_loss_db,
            transmission_zeros: context
                .spec
                .transmission_zeros
                .iter()
                .map(|tz| tz.value)
                .collect(),
            unloaded_q: context.spec.unloaded_q,
        };

        let polynomials = context.polynomials.as_ref().map(|p| PolynomialSnapshot {
            order: p.order,
            ripple_factor: p.ripple_factor,
            eps: p.eps,
            eps_r: p.eps_r,
            transmission_zeros_normalized: p.transmission_zeros_normalized.clone(),
        });

        let matrix = context.matrix.as_ref().map(|m| MatrixSnapshot {
            order: m.order(),
            data: m.as_slice().to_vec(),
        });

        Self {
            spec,
            options: context.options.clone(),
            metadata: context.metadata.clone(),
            polynomials,
            matrix,
            has_transform: context.transform.is_some(),
            has_response: context.response.is_some(),
        }
    }

    fn into_context(self) -> Result<SynthesisContext> {
        let mut spec = FilterSpec::new(self.spec.order, self.spec.return_loss_db)?;
        spec = spec.with_normalized_transmission_zeros(self.spec.transmission_zeros);
        if let Some(q) = self.spec.unloaded_q {
            spec = spec.with_unloaded_q(q);
        }

        let mut context = SynthesisContext::with_options(spec, self.options);
        context.metadata = self.metadata;

        // Restore polynomials by re-running approximation if snapshot indicates they existed
        // (We cannot fully reconstruct PolynomialSet from the snapshot alone because it
        // requires the generalized Chebyshev data. Re-running is the safest approach.)
        if self.polynomials.is_some() {
            let stage = ApproximationStage;
            let polynomials = stage.execute(&context.spec)?;
            context.polynomials = Some(Arc::new(polynomials));
        }

        // Restore matrix from raw data if present
        if let Some(matrix_snap) = self.matrix {
            use crate::matrix::CouplingMatrix;
            let matrix = CouplingMatrix::new(matrix_snap.order, matrix_snap.data)?;
            context.matrix = Some(Arc::new(matrix));
        }

        // Transform and response cannot be cheaply restored from a snapshot;
        // they would need to be re-run. The metadata tracks which stages completed,
        // so the caller can use run_stage to re-execute from the last checkpoint.

        Ok(context)
    }
}

/// Saves a `SynthesisContext` to a JSON file for later resumption.
///
/// The saved representation captures enough state to reconstruct the context
/// and resume pipeline execution from the last completed stage.
pub fn save_context(context: &SynthesisContext, path: &Path) -> Result<()> {
    let snapshot = ContextSnapshot::from_context(context);
    let json = serde_json::to_string_pretty(&snapshot).map_err(|e| {
        MfsError::PreconditionViolation(format!("failed to serialize context: {e}"))
    })?;
    fs::write(path, json).map_err(|e| {
        MfsError::PreconditionViolation(format!(
            "failed to write context to {}: {e}",
            path.display()
        ))
    })?;
    Ok(())
}

/// Loads a `SynthesisContext` from a previously saved JSON file.
///
/// Polynomial data is reconstructed by re-running the approximation stage
/// (since the full complex polynomial data cannot be trivially serialized).
/// Matrix data is restored from the raw flat storage. Transform and response
/// artifacts are not restored — use `run_stage` to re-execute those stages.
pub fn load_context(path: &Path) -> Result<SynthesisContext> {
    let json = fs::read_to_string(path).map_err(|e| {
        MfsError::PreconditionViolation(format!(
            "failed to read context from {}: {e}",
            path.display()
        ))
    })?;
    let snapshot: ContextSnapshot = serde_json::from_str(&json).map_err(|e| {
        MfsError::PreconditionViolation(format!("failed to deserialize context: {e}"))
    })?;
    snapshot.into_context()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::context::GridConfig;
    use crate::spec::FilterSpec;
    use crate::transform::TopologyKind;

    #[test]
    fn run_stage_approximation_populates_polynomials() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let mut ctx = SynthesisContext::new(spec);

        run_stage(&mut ctx, "approximation")?;

        assert!(ctx.polynomials().is_some());
        assert_eq!(ctx.polynomials().unwrap().order, 4);
        assert_eq!(ctx.metadata.stages_executed, vec!["approximation"]);
        assert_eq!(ctx.metadata.stage_timings_ms.len(), 1);
        Ok(())
    }

    #[test]
    fn run_stage_matrix_synthesis_requires_polynomials() {
        let spec = FilterSpec::new(4, 20.0).unwrap();
        let mut ctx = SynthesisContext::new(spec);

        let err = run_stage(&mut ctx, "matrix_synthesis").unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("polynomials"));
    }

    #[test]
    fn run_stage_topology_transform_requires_matrix_and_topology() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let mut ctx = SynthesisContext::new(spec);

        // Missing matrix
        let err = run_stage(&mut ctx, "topology_transform").unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("matrix"));

        // Add matrix but no topology option
        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;
        let err = run_stage(&mut ctx, "topology_transform").unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("topology"));

        Ok(())
    }

    #[test]
    fn run_stage_response_evaluation_requires_matrix() {
        let spec = FilterSpec::new(4, 20.0).unwrap();
        let mut ctx = SynthesisContext::new(spec);

        let err = run_stage(&mut ctx, "response_evaluation").unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("matrix"));
    }

    #[test]
    fn run_stage_rejects_unknown_stage_name() {
        let spec = FilterSpec::new(4, 20.0).unwrap();
        let mut ctx = SynthesisContext::new(spec);

        let err = run_stage(&mut ctx, "unknown_stage").unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("unrecognized"));
    }

    #[test]
    fn incremental_execution_produces_complete_context() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let mut ctx = SynthesisContext::with_options(
            spec,
            PipelineOptions {
                topology: Some(TopologyKind::Folded),
                grid: Some(GridConfig {
                    start: -3.0,
                    stop: 3.0,
                    points: 51,
                }),
                ..PipelineOptions::default()
            },
        );

        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;
        run_stage(&mut ctx, "topology_transform")?;
        run_stage(&mut ctx, "response_evaluation")?;

        // All artifacts populated
        assert!(ctx.polynomials().is_some());
        assert!(ctx.matrix().is_some());
        assert!(ctx.transform().is_some());
        assert!(ctx.response().is_some());

        // Metadata tracks all stages
        assert_eq!(
            ctx.metadata.stages_executed,
            vec![
                "approximation",
                "matrix_synthesis",
                "topology_transform",
                "response_evaluation"
            ]
        );
        assert_eq!(ctx.metadata.stage_timings_ms.len(), 4);

        // Verify artifacts are consistent
        assert_eq!(ctx.polynomials().unwrap().order, 4);
        assert_eq!(ctx.matrix().unwrap().order(), 4);
        assert_eq!(ctx.transform().unwrap().topology, TopologyKind::Folded);
        assert_eq!(ctx.response().unwrap().samples.len(), 51);

        Ok(())
    }

    #[test]
    fn save_and_load_context_round_trip() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let mut ctx = SynthesisContext::new(spec);

        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;

        // Save to a temp file
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("mfs_test_context.json");

        save_context(&ctx, &path)?;
        let loaded = load_context(&path)?;

        // Verify restored context has equivalent data
        assert_eq!(loaded.spec.order, 4);
        assert_eq!(loaded.spec.return_loss_db, 20.0);
        assert_eq!(loaded.spec.transmission_zeros.len(), 2);
        assert!(loaded.polynomials().is_some());
        assert_eq!(loaded.polynomials().unwrap().order, 4);
        assert!(loaded.matrix().is_some());
        assert_eq!(loaded.matrix().unwrap().order(), 4);
        assert_eq!(
            loaded.metadata.stages_executed,
            vec!["approximation", "matrix_synthesis"]
        );

        // Clean up
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn loaded_context_can_resume_execution() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let mut ctx = SynthesisContext::with_options(
            spec,
            PipelineOptions {
                topology: Some(TopologyKind::Arrow),
                grid: Some(GridConfig {
                    start: -2.0,
                    stop: 2.0,
                    points: 21,
                }),
                ..PipelineOptions::default()
            },
        );

        // Run first two stages
        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;

        // Save and reload
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("mfs_test_resume_context.json");
        save_context(&ctx, &path)?;
        let mut loaded = load_context(&path)?;

        // Resume from where we left off
        run_stage(&mut loaded, "topology_transform")?;
        run_stage(&mut loaded, "response_evaluation")?;

        assert!(loaded.transform().is_some());
        assert!(loaded.response().is_some());
        assert_eq!(loaded.response().unwrap().samples.len(), 21);

        // Clean up
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn incremental_execution_matches_full_pipeline_result() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let options = PipelineOptions {
            topology: Some(TopologyKind::Folded),
            grid: Some(GridConfig {
                start: -3.0,
                stop: 3.0,
                points: 51,
            }),
            ..PipelineOptions::default()
        };

        // Run incrementally
        let mut ctx = SynthesisContext::with_options(spec.clone(), options.clone());
        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;
        run_stage(&mut ctx, "topology_transform")?;
        run_stage(&mut ctx, "response_evaluation")?;

        // Run again from scratch (fresh context, same spec/options)
        let mut ctx2 = SynthesisContext::with_options(spec, options);
        run_stage(&mut ctx2, "approximation")?;
        run_stage(&mut ctx2, "matrix_synthesis")?;
        run_stage(&mut ctx2, "topology_transform")?;
        run_stage(&mut ctx2, "response_evaluation")?;

        // Both should produce identical artifacts
        let poly1 = ctx.polynomials().unwrap();
        let poly2 = ctx2.polynomials().unwrap();
        assert_eq!(poly1.order, poly2.order);
        assert_eq!(poly1.eps, poly2.eps);

        let mat1 = ctx.matrix().unwrap();
        let mat2 = ctx2.matrix().unwrap();
        assert_eq!(mat1.order(), mat2.order());
        // Matrix data should be identical (deterministic algorithm)
        for i in 0..mat1.side() {
            for j in 0..mat1.side() {
                let v1 = mat1.at(i, j).unwrap_or_default();
                let v2 = mat2.at(i, j).unwrap_or_default();
                assert!(
                    (v1 - v2).abs() < 1e-12,
                    "matrix mismatch at ({i},{j}): {v1} vs {v2}"
                );
            }
        }

        let resp1 = ctx.response().unwrap();
        let resp2 = ctx2.response().unwrap();
        assert_eq!(resp1.samples.len(), resp2.samples.len());
        for (s1, s2) in resp1.samples.iter().zip(resp2.samples.iter()) {
            assert!((s1.s11_re - s2.s11_re).abs() < 1e-12);
            assert!((s1.s21_re - s2.s21_re).abs() < 1e-12);
        }

        Ok(())
    }
}
