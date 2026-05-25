//! Pipeline orchestration: context, stages, and schema generation.
//!
//! This module defines the synthesis pipeline architecture:
//! - `SynthesisContext`: accumulates stage artifacts across the workflow
//! - `PipelineOptions`: configures pipeline behavior
//! - `PipelineMetadata`: tracks execution progress and diagnostics
//! - `Stage`: trait for composable, testable pipeline stages
//! - Concrete stages: `ApproximationStage`, `MatrixSynthesisStage`,
//!   `TopologyTransformStage`, `ResponseEvaluationStage`
//! - `describe_schema`: JSON Schema generation for structured I/O
//! - `run_stage`: incremental single-stage execution with predecessor validation
//! - `save_context` / `load_context`: file-based context persistence

pub mod context;
pub mod execution;
pub mod request;
pub mod schema;
pub mod stages;

pub use context::{
    GridConfig, MappingConfig, PipelineMetadata, PipelineOptions, ResponseToleranceConfig,
    SynthesisContext,
};
pub use execution::{load_context, run_stage, save_context};
pub use request::{
    run_from_json, run_full_pipeline, SynthesisRequest, ValidationError, ValidationErrors,
};
pub use schema::describe_schema;
pub use stages::{
    ApproximationStage, MatrixSynthesisStage, ResponseEvaluationInput, ResponseEvaluationStage,
    Stage, TopologyTransformInput, TopologyTransformStage,
};
