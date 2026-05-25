# Implementation Plan: Architecture Cleanup

## Overview

本实现计划将 MFS 库的架构清理分为 4 个阶段、15 个任务。Phase 1-2 是内部质量改进（无新功能），Phase 3-4 是新增 Pipeline 架构和 CLI 适配层。

## Tasks

### Phase 1: Foundation (错误类型 + 矩阵拆分)

- [x] 1. Refactor `MfsError` enum: replace `Unsupported` with `NumericalFailure`, `NotImplemented`, `PreconditionViolation`
  - Add `#[derive(Serialize, Deserialize)]` to `MfsError`
  - Update all existing `MfsError::Unsupported(...)` call sites to use the appropriate new variant
  - Preserve `Display` formatting for unchanged variants
  - Update tests that match on `Unsupported`

- [x] 2. Split `coupling_matrix.rs` into submodules [depends on: 1]
  - Create `src/matrix/core.rs` with `CouplingMatrix` struct, constructors, accessors, `MatrixTopology`, `MatrixShape`, `BandPassScaledCouplingMatrix`
  - Create `src/matrix/rotations.rs` with `rotate_matrix`, `rotate_matrix_with_indices`, `rotation_matrix_basic`, `safe_angle`, `diagonal_rotation_angle`, `RotationAxis`
  - Create `src/matrix/sections.rs` with `extract_triplet`, `extract_quadruplet`, `extract_trisection` and their validators
  - Create `src/matrix/scaling.rs` with `denormalize_bandpass`, `normalize_bandpass`, `denormalize_bandpass_with_external_q`, `normalize_bandpass_with_external_q` and helper functions
  - Update `src/matrix/mod.rs` to re-export all public types unchanged
  - Verify all existing tests pass without modification

- [x] 3. Replace `unwrap_or_default()` with debug-panicking accessor [depends on: 2]
  - Add `pub(crate) fn get(&self, row: usize, col: usize) -> f64` that uses `debug_assert!` + unchecked indexing
  - Replace all internal `self.at(x, y).unwrap_or_default()` calls in rotations and sections with `self.get(x, y)`
  - Keep public `at()` returning `Option<f64>` unchanged
  - Run full test suite in debug mode to verify no panics on valid inputs

- [x] 4. Replace hand-written matrix multiply/transpose with nalgebra [depends on: 3]
  - Rewrite `multiply()` to use `DMatrix<f64>` multiplication
  - Rewrite `transpose()` to use nalgebra transpose
  - Add a comparison test that verifies numerical equivalence (< 1e-12 deviation) against the old implementation for several known matrices
  - Remove the old triple-loop implementation

### Phase 2: Cleanup (Placeholder 移除 + Transform 统一)

- [x] 5. Remove placeholder synthesis path [depends on: 1]
  - Delete `src/synthesis/placeholder.rs`
  - Remove `pub(crate) use placeholder::synthesize_placeholder_matrix` from `synthesis/mod.rs`
  - Remove `MatrixSynthesisMethod::PlaceholderFallback` variant
  - Update `MatrixSynthesisEngine::synthesize_with_details` to return `PreconditionViolation` when generalized data is missing
  - Update any tests that relied on placeholder fallback to either provide valid generalized data or assert the expected error

- [x] 6. Unify `SectionTransformOutcome` and `VerifiedSectionSynthesis` into single type [depends on: 2]
  - Keep `SectionTransformOutcome` in `transform/sections.rs` as the canonical type
  - Remove `VerifiedSectionSynthesis` from `synthesis/sections.rs`
  - Update `SectionSynthesis` methods to return `SectionTransformOutcome` directly
  - Update all call sites and tests

- [x] 7. Extract shared response-check logic into single helper [depends on: 3, 6]
  - Identify the duplicated `attach_response_check` pattern in `transform/mod.rs` and `transform/sections.rs`
  - Create a single `pub(crate) fn attach_response_invariance_check(...)` in a shared location
  - Replace both duplicates with calls to the shared helper
  - Verify response-check tests still pass

- [x] 8. Deprecate Wheel topology and clean up SynthesisOutcome [depends on: 1]
  - Add `#[deprecated]` attribute to `MatrixTopology::Wheel` with explanatory message
  - Add doc comment on `Wheel` variant explaining current behavior
  - Remove `SynthesisOutcome::used_generalized_approximation()` method
  - Remove `SynthesisOutcome::approximation_kind()` method
  - Add `pub approximation: ApproximationKind` enum field to `SynthesisOutcome`
  - Update all call sites that used the removed methods

- [x] 9. Demote `CouplingMatrix::transform_topology` to `pub(crate)` [depends on: 1]
  - Change `pub fn transform_topology` to `pub(crate) fn transform_topology`
  - Ensure all external callers use `transform::transform_matrix(...)` instead
  - Update any examples or tests that called `transform_topology` directly

### Phase 3: Pipeline Architecture (上下文 + 阶段)

- [x] 10. Create `src/pipeline/` module with `SynthesisContext` and `PipelineMetadata` [depends on: 4, 7]
  - Create `src/pipeline/mod.rs` with module declarations and re-exports
  - Create `src/pipeline/context.rs` with `SynthesisContext`, `PipelineOptions`, `PipelineMetadata` structs
  - All structs derive `Debug, Clone, Serialize, Deserialize`
  - Add typed accessor methods (`polynomials()`, `matrix()`, `response()`) returning `Option<&T>`
  - Add `pub mod pipeline` to `src/lib.rs`

- [x] 11. Define `Stage` trait and implement four concrete stages [depends on: 10]
  - Create `src/pipeline/stages.rs` with `Stage` trait definition
  - Implement `ApproximationStage` (FilterSpec → PolynomialSet)
  - Implement `MatrixSynthesisStage` (PolynomialSet → CouplingMatrix)
  - Implement `TopologyTransformStage` (CouplingMatrix + TopologyKind → TransformOutcome)
  - Implement `ResponseEvaluationStage` (CouplingMatrix + FrequencyGrid → SParameterResponse)
  - Add unit tests for each stage in isolation

- [x] 12. Implement `SynthesisRequest` JSON input and `run_full_pipeline` [depends on: 11]
  - Define `SynthesisRequest` struct with all optional fields and defaults
  - Implement `run_full_pipeline(request: SynthesisRequest) -> Result<SynthesisContext>`
  - Implement `run_from_json(json: &str) -> Result<String>` that parses, executes, and serializes
  - Add validation that returns structured field-level errors for malformed input
  - Add tests for minimal request (order + return_loss only) and full request

- [x] 13. Implement `describe_schema()` and JSON output format [depends on: 11]
  - Create `src/pipeline/schema.rs`
  - Implement `describe_schema()` returning JSON Schema for `SynthesisRequest` and output format
  - Ensure output JSON includes: all artifacts, metadata (version, stages_executed, stage_timings_ms, warnings)
  - Add test that validates a known-good request against the generated schema

- [x] 14. Implement incremental execution and context resume [depends on: 11]
  - Add `run_stage(context: &mut SynthesisContext, stage_name: &str) -> Result<()>`
  - Validate predecessor artifacts exist before running each stage
  - Return `PreconditionViolation` if predecessors are missing
  - Add `save_context` / `load_context` helpers for file-based persistence
  - Add test for incremental execution producing same result as full pipeline

### Phase 4: CLI Adapter

- [x] 15. Create CLI binary with JSON stdin/file support [depends on: 12, 13, 14]
  - Add `src/bin/mfs_cli.rs` with argument parsing (clap or similar)
  - Accept `--input <file>` or read from stdin when no file specified
  - Output full pipeline result as JSON to stdout
  - Add `--format table` flag for human-readable terminal output

- [x] 16. Implement `--stage` and `--resume` flags [depends on: 15]
  - `--stage <name>` limits execution to a single stage, outputs only that stage's artifact
  - `--resume <context_file>` loads a saved `SynthesisContext` and continues from last completed stage
  - Structured JSON error output to stderr on failure

- [x] 17. Add integration tests for CLI [depends on: 16]
  - Test JSON file input → JSON stdout output
  - Test stdin input → JSON stdout output
  - Test `--stage` flag produces partial output
  - Test invalid input produces structured error on stderr

## Notes

- Phase 1 和 Phase 2 可以部分并行执行（1 是所有后续任务的前置）
- Phase 3 依赖 Phase 1-2 完成，因为 Pipeline 需要使用重构后的类型
- Phase 4 依赖 Phase 3 完成，因为 CLI 是 Pipeline 的薄包装
- 每个任务完成后应运行 `cargo test` 确保无回归
- `MfsError::Unsupported` 的移除是 breaking change，建议在同一个 PR 中完成所有 Phase 1 任务
