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
use crate::synthesis::build_filter_spec;

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
    context
        .metadata
        .stages_executed
        .push("approximation".to_string());
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
    context
        .metadata
        .stages_executed
        .push("matrix_synthesis".to_string());
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
    context
        .metadata
        .stages_executed
        .push("topology_transform".to_string());
    context.metadata.stage_timings_ms.push(elapsed_ms);
    Ok(())
}

fn run_response_evaluation_stage(context: &mut SynthesisContext) -> Result<()> {
    // Use the transformed (folded) matrix when a transform exists, otherwise
    // the synthesis matrix — exactly mirroring `run_full_pipeline`, which picks
    // the same `eval_matrix`. Otherwise a resumed context would evaluate the
    // untransformed (transversal) matrix and diverge from a fresh full run at
    // machine precision (review: matrix-source consistency between the two
    // response-evaluation drivers).
    let eval_matrix = match context.transform.as_ref() {
        Some(transform) => transform.matrix.clone(),
        None => context
            .matrix()
            .ok_or_else(|| {
                MfsError::PreconditionViolation(
                    "response_evaluation requires a coupling matrix from the matrix_synthesis stage"
                        .to_string(),
                )
            })?
            .clone(),
    };

    // Build a frequency grid from options or use a sensible default
    let grid = build_frequency_grid(&context.options)?;

    let start = Instant::now();
    let stage = ResponseEvaluationStage;
    let input = ResponseEvaluationInput {
        matrix: eval_matrix,
        grid,
        mapping: context.options.mapping.clone(),
        unloaded_q: context.spec.unloaded_q,
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
        Some(grid_config) => {
            FrequencyGrid::linspace(grid_config.start, grid_config.stop, grid_config.points)
        }
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
}

/// A finite normalized value, or `None` for a non-finite one (a zero at
/// infinity is `f64::INFINITY`). `serde_json` cannot round-trip a non-finite
/// `f64` (it serializes to `null`, which fails to deserialize back to `f64`),
/// so the snapshot stores them as `Option` — the same convention the output
/// serializer (`context.rs` `SpecSummary`) uses (review P2-10 / round-2 F1).
type FiniteZero = Option<f64>;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpecSnapshot {
    order: usize,
    return_loss_db: f64,
    transmission_zeros: Vec<FiniteZero>,
    unloaded_q: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PolynomialSnapshot {
    order: usize,
    ripple_factor: f64,
    eps: f64,
    eps_r: f64,
    transmission_zeros_normalized: Vec<FiniteZero>,
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
            // A zero at infinity is `f64::INFINITY`; store it as `None` so the
            // JSON round-trips (serde_json cannot serialize non-finite f64).
            transmission_zeros: context
                .spec
                .transmission_zeros
                .iter()
                .map(|tz| tz.value.is_finite().then_some(tz.value))
                .collect(),
            // A non-finite Q (e.g. the lossless `INFINITY` sentinel) would
            // break the round trip the same way; drop it to `None` (lossless).
            unloaded_q: context
                .spec
                .unloaded_q
                .filter(|q| q.is_finite()),
        };

        let polynomials = context.polynomials.as_ref().map(|p| PolynomialSnapshot {
            order: p.order,
            ripple_factor: p.ripple_factor,
            eps: p.eps,
            eps_r: p.eps_r,
            transmission_zeros_normalized: p
                .transmission_zeros_normalized
                .iter()
                .copied()
                .map(|z| z.is_finite().then_some(z))
                .collect(),
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
        }
    }

    fn into_context(self) -> Result<SynthesisContext> {
        // A `null` zero was a zero at infinity (`f64::INFINITY`); a finite
        // value is used as-is. Reconstruct the original `Vec<f64>`.
        let transmission_zeros: Vec<f64> = self
            .spec
            .transmission_zeros
            .iter()
            .map(|z| z.unwrap_or(f64::INFINITY))
            .collect();
        let spec = build_filter_spec(
            self.spec.order,
            self.spec.return_loss_db,
            &transmission_zeros,
            self.spec.unloaded_q,
        )?;

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

        // A saved context can be hand-crafted or partially written: reject any
        // snapshot whose metadata disagrees with its artifacts instead of
        // letting run_remaining_stages skip a stage on bad metadata.
        //
        // Transform and response artifacts are never restored from a snapshot,
        // so the metadata must not claim they completed either: otherwise
        // run_remaining_stages would skip them and a resumed fully-completed
        // context would come back with no transform/response while its metadata
        // says those stages ran (review P2-9). Keep only the stages whose
        // artifacts were actually restored, dropping the parallel per-stage
        // timings so stage_timings_ms stays aligned with stages_executed.
        let mut stages = Vec::new();
        let mut timings = Vec::new();
        // A hand-crafted or truncated snapshot can claim more stages than it
        // has per-stage timings; reject it up front instead of indexing out of
        // bounds in the loop below (review round-2 F2).
        if context.metadata.stage_timings_ms.len()
            < context.metadata.stages_executed.len()
        {
            return Err(MfsError::PreconditionViolation(format!(
                "corrupt context snapshot: stage_timings_ms has {} entries \
                 but stages_executed has {}",
                context.metadata.stage_timings_ms.len(),
                context.metadata.stages_executed.len()
            )));
        }
        for (i, name) in context.metadata.stages_executed.iter().enumerate() {
            let has_artifact = match name.as_str() {
                "approximation" => context.polynomials.is_some(),
                "matrix_synthesis" => context.matrix.is_some(),
                _ => true,
            };
            if !has_artifact {
                return Err(MfsError::PreconditionViolation(format!(
                    "corrupt context snapshot: metadata claims stage '{name}' ran, but its \
                     artifact is missing"
                )));
            }
            if matches!(name.as_str(), "approximation" | "matrix_synthesis") {
                stages.push(name.clone());
                timings.push(context.metadata.stage_timings_ms[i]);
            }
        }
        context.metadata.stages_executed = stages;
        context.metadata.stage_timings_ms = timings;

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

/// Runs every stage that has not yet executed, in pipeline order.
///
/// This is the resume driver used after [`load_context`]: it skips stages
/// whose names appear in [`PipelineMetadata::stages_executed`] and runs the
/// rest. `topology_transform` is optional and is skipped (rather than
/// failed) when `options.topology` is not set; every other stage —
/// including `response_evaluation` — is always attempted once its
/// predecessors are present.
pub fn run_remaining_stages(ctx: &mut SynthesisContext) -> Result<()> {
    const ALL_STAGES: [&str; 4] = [
        "approximation",
        "matrix_synthesis",
        "topology_transform",
        "response_evaluation",
    ];

    for stage_name in ALL_STAGES {
        if ctx
            .metadata
            .stages_executed
            .iter()
            .any(|name| name == stage_name)
        {
            continue;
        }
        // Skip the optional transform when no topology was requested.
        if stage_name == "topology_transform" && ctx.options.topology.is_none() {
            continue;
        }
        run_stage(ctx, stage_name)?;
    }

    Ok(())
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
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
    fn save_and_load_context_round_trip_with_infinite_zero() -> Result<()> {
        // A zero at infinity (`f64::INFINITY`) and the lossless `INFINITY` Q
        // are legal spec values; serde_json cannot round-trip non-finite f64,
        // so the snapshot must store them as `null`/`None` and restore them
        // faithfully (review round-2 F1).
        let spec = FilterSpec::new(4, 20.0)?
            .with_transmission_zeros(vec![
                crate::spec::TransmissionZero::normalized(2.0),
                crate::spec::TransmissionZero::infinite(),
            ])
            .with_unloaded_q(f64::INFINITY);
        let mut ctx = SynthesisContext::new(spec);
        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;

        let path = std::env::temp_dir().join("mfs_test_context_inf_zero.json");
        save_context(&ctx, &path)?;
        let loaded = load_context(&path)?;

        assert_eq!(loaded.spec.transmission_zeros.len(), 2);
        assert_eq!(loaded.spec.transmission_zeros[0].value, 2.0);
        assert_eq!(loaded.spec.transmission_zeros[1].value, f64::INFINITY);
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn context_serialize_emits_null_for_infinite_normalized_zero() -> Result<()> {
        // The public `Serialize` impl for `SynthesisContext` (what CLI
        // `output_context` and any JSON dump use) must emit a zero at
        // infinity in the polynomial summary as `null` — matching the
        // schema (`transmission_zeros_normalized` items are `["number",
        // "null"]`) and the spec summary, not a bare `inf` (review round-2
        // F4 follow-up: the spec summary was fixed in P2-10 but the
        // polynomial summary was missed).
        let spec = FilterSpec::new(4, 20.0)?
            .with_transmission_zeros(vec![
                crate::spec::TransmissionZero::normalized(-1.5),
                crate::spec::TransmissionZero::infinite(),
            ]);
        let mut ctx = SynthesisContext::new(spec);
        run_stage(&mut ctx, "approximation")?;

        let json = serde_json::to_string(&ctx).expect("serialization must succeed");
        let value: serde_json::Value =
            serde_json::from_str(&json).expect("output must be valid JSON");
        let zeros = &value["polynomials"]["transmission_zeros_normalized"];
        let zeros = zeros
            .as_array()
            .expect("polynomials.transmission_zeros_normalized must be an array");
        assert_eq!(zeros.len(), 2, "expected 2 normalized zeros, got {zeros:?}");
        assert!(
            zeros[1].is_null(),
            "a zero at infinity must serialize as null, got {:?}",
            zeros[1]
        );
        assert!(
            zeros[0].as_f64() == Some(-1.5),
            "finite zero must be preserved, got {:?}",
            zeros[0]
        );
        Ok(())
    }

    #[test]
    fn into_context_rejects_mismatched_stage_timings() -> Result<()> {
        // A hand-crafted snapshot that claims more stages than it has
        // per-stage timings must be rejected, not panic on an index out of
        // bounds (review round-2 F2).
        let spec = FilterSpec::new(4, 20.0)?;
        let mut ctx = SynthesisContext::new(spec);
        run_stage(&mut ctx, "approximation")?;
        run_stage(&mut ctx, "matrix_synthesis")?;
        // 3 claimed stages but only 2 timings recorded.
        ctx.metadata.stages_executed.push("response_evaluation".to_string());
        let snapshot = ContextSnapshot::from_context(&ctx);
        let err = snapshot.into_context().unwrap_err();
        assert!(matches!(err, MfsError::PreconditionViolation(_)));
        assert!(err.to_string().contains("stage_timings_ms"));
        Ok(())
    }

    #[test]
    fn loaded_context_can_resume_execution() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
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
    fn resume_of_completed_context_reproduces_full_pipeline() -> Result<()> {
        // Regression for review P2-9: a fully-completed context, saved and
        // resumed via run_remaining_stages (the `--resume` driver), must come
        // back with a transform and a response that match the original run —
        // not a context whose metadata claims those stages ran but whose
        // artifacts are missing.
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let options = PipelineOptions {
            topology: Some(TopologyKind::Folded),
            grid: Some(GridConfig {
                start: -2.0,
                stop: 2.0,
                points: 31,
            }),
            ..PipelineOptions::default()
        };

        // Full pipeline reference run.
        let mut full = SynthesisContext::with_options(spec.clone(), options.clone());
        run_remaining_stages(&mut full)?;
        assert!(full.transform().is_some());
        assert!(full.response().is_some());

        // Persist and reload.
        let path = std::env::temp_dir().join("mfs_test_resume_full_context.json");
        save_context(&full, &path)?;
        let mut resumed = load_context(&path)?;

        // The loaded context must NOT claim transform/response completed,
        // since their artifacts were not restored.
        assert!(!resumed
            .metadata
            .stages_executed
            .iter()
            .any(|s| s == "topology_transform"));
        assert!(!resumed
            .metadata
            .stages_executed
            .iter()
            .any(|s| s == "response_evaluation"));

        // Resuming runs the still-missing stages on top of the restored
        // approximation + matrix artifacts.
        run_remaining_stages(&mut resumed)?;

        assert!(resumed.transform().is_some());
        assert!(resumed.response().is_some());
        assert_eq!(
            resumed.transform().unwrap().topology,
            full.transform().unwrap().topology
        );

        let r_full = full.response().unwrap();
        let r_res = resumed.response().unwrap();
        assert_eq!(r_full.samples.len(), r_res.samples.len());
        for (a, b) in r_full.samples.iter().zip(r_res.samples.iter()) {
            assert!((a.s11_re - b.s11_re).abs() < 1e-9, "S11 re mismatch");
            assert!((a.s11_im - b.s11_im).abs() < 1e-9, "S11 im mismatch");
            assert!((a.s21_re - b.s21_re).abs() < 1e-9, "S21 re mismatch");
            assert!((a.s21_im - b.s21_im).abs() < 1e-9, "S21 im mismatch");
        }

        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn incremental_execution_matches_full_pipeline_result() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
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
