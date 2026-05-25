# Requirements Document

## Introduction

本文档定义了 MFS（微波滤波器综合）Rust 库的架构清理与重构需求。核心目标是：
1. 优化高层数据流组织，使 pipeline 各阶段的数据传递更高效、更清晰
2. 使 API 适合 AI agent 调用（结构化输入/输出、可序列化、自描述）
3. 使 API 适合 CLI 工具集成（单次调用完成完整流程、JSON 输入/输出）
4. 消除代码冗余，改善模块边界
5. 统一错误处理语义

重构不改变核心算法行为，侧重架构实现质量。

## Glossary

- **Pipeline**: 完整的滤波器综合流水线：FilterSpec → Approximation → CouplingMatrix → Transform → Response
- **Stage_Artifact**: 流水线每个阶段产出的不可变数据对象（PolynomialSet、CouplingMatrix、SParameterResponse 等）
- **Pipeline_Context**: 贯穿整个综合流程的上下文对象，携带所有阶段产出和配置
- **Synthesis_Engine**: 矩阵综合引擎，负责从多项式生成耦合矩阵
- **Transform_Module**: `src/transform/` 目录，拓扑变换和截面提取的公共门面
- **MfsError**: 库级错误枚举类型
- **Coupling_Matrix_Module**: 耦合矩阵数据结构及其操作
- **Section_Extraction**: triplet/quadruplet/trisection 截面提取操作
- **Rotation_Helpers**: 矩阵旋转辅助函数
- **Placeholder_Path**: 残余展开不可用时的链式矩阵回退构建器
- **AI_Agent**: 通过结构化 API 调用库功能的自动化程序
- **CLI_Adapter**: 命令行接口适配层，接受 JSON 输入并产出 JSON 输出

## Requirements

### Requirement 1: 统一 Pipeline 上下文与数据流

**User Story:** As a library integrator, I want a single pipeline context that accumulates stage artifacts, so that I can inspect any intermediate result without re-running earlier stages.

#### Acceptance Criteria

1. THE Pipeline SHALL define a `SynthesisContext` struct that holds all stage artifacts: spec, frequency mapping, polynomials, matrix, transform result, and response
2. WHEN a pipeline stage completes, THE Pipeline SHALL store its output in the context without consuming or invalidating earlier artifacts
3. THE Pipeline SHALL support both full-pipeline execution (spec → response in one call) and incremental stage-by-stage execution using the same context type
4. WHEN incremental execution is used, THE Pipeline SHALL validate that required predecessor artifacts exist before running a stage
5. THE `SynthesisContext` SHALL derive `Serialize` and `Deserialize` so that intermediate state can be saved to disk and resumed later
6. THE `SynthesisContext` SHALL provide typed accessor methods for each stage artifact that return `Option<&T>` rather than requiring the caller to know internal field names

### Requirement 2: 流水线阶段的显式类型化

**User Story:** As a library maintainer, I want each pipeline stage to be a standalone callable unit with explicit input/output types, so that stages can be composed, tested, and called independently.

#### Acceptance Criteria

1. THE Pipeline SHALL define a `Stage` trait with associated `Input` and `Output` types and a single `execute` method
2. WHEN a stage is executed, THE Stage SHALL accept only its declared input type and return only its declared output type
3. THE Pipeline SHALL implement the following stages as distinct types: `ApproximationStage`, `MatrixSynthesisStage`, `TopologyTransformStage`, `ResponseEvaluationStage`
4. WHEN stages are composed into a full pipeline, THE Pipeline SHALL use the `Stage` trait to chain them without type erasure
5. THE Pipeline SHALL allow users to replace any stage implementation with a custom one that satisfies the same trait bounds

### Requirement 3: 结构化输入/输出适配 AI 调用

**User Story:** As an AI agent developer, I want all pipeline inputs and outputs to be JSON-serializable with self-describing schemas, so that an AI agent can construct valid requests and parse responses without hardcoded knowledge.

#### Acceptance Criteria

1. THE Pipeline SHALL accept a single JSON object as input that fully describes a synthesis request (spec, mapping, grid, topology, options)
2. THE Pipeline SHALL produce a single JSON object as output that contains all requested artifacts (polynomials, matrix entries, response samples, verification reports)
3. WHEN the input JSON is malformed or incomplete, THE Pipeline SHALL return a structured error response with field-level validation messages
4. THE Pipeline SHALL provide a `describe_schema()` function that returns the JSON Schema for valid input and output objects
5. THE Pipeline input format SHALL support optional fields with documented defaults, so that minimal requests (order + return_loss only) produce valid results
6. THE Pipeline output SHALL include metadata fields: pipeline version, stages executed, execution time per stage, and warnings

### Requirement 4: CLI 适配层设计

**User Story:** As a CLI user, I want to run the full synthesis pipeline from a single command with JSON input/output, so that I can integrate MFS into scripts and automation workflows.

#### Acceptance Criteria

1. THE CLI_Adapter SHALL accept a JSON file path or stdin JSON as the synthesis request
2. THE CLI_Adapter SHALL output the full pipeline result as JSON to stdout
3. WHEN the `--format` flag is set to `table`, THE CLI_Adapter SHALL output a human-readable summary instead of JSON
4. THE CLI_Adapter SHALL support a `--stage` flag that limits execution to a specific pipeline stage and outputs only that stage's artifact
5. THE CLI_Adapter SHALL support a `--resume` flag that loads a previously saved `SynthesisContext` and continues from the last completed stage
6. IF an error occurs, THEN THE CLI_Adapter SHALL output a JSON error object to stderr with error type, message, and context

### Requirement 5: 消除数据传递中的冗余克隆

**User Story:** As a performance-conscious user, I want the pipeline to avoid unnecessary data copies between stages, so that high-order matrix operations remain efficient.

#### Acceptance Criteria

1. THE Pipeline SHALL pass stage artifacts by reference between stages rather than cloning them at each boundary
2. WHEN a stage needs to mutate data (such as topology transform), THE Pipeline SHALL clone only the specific artifact being mutated, not the entire context
3. THE `CouplingMatrix` SHALL implement a zero-copy view for read-only operations (response evaluation, verification) that does not require cloning the underlying data
4. WHEN the `SynthesisContext` is serialized, THE Pipeline SHALL serialize lazily (on demand) rather than eagerly converting all artifacts to serializable form at each stage
5. THE Pipeline SHALL use `Arc<T>` for shared immutable artifacts when multiple downstream consumers need the same data simultaneously

### Requirement 6: 错误类型语义细分

**User Story:** As a library consumer, I want distinct error variants for different failure modes, so that I can programmatically distinguish numerical failures from unimplemented features and invalid inputs.

#### Acceptance Criteria

1. THE MfsError enum SHALL replace the single `Unsupported` variant with at least three distinct variants: `NumericalFailure`, `NotImplemented`, and `PreconditionViolation`
2. WHEN a root solver fails to converge, THE MfsError SHALL use the `NumericalFailure` variant with a descriptive message
3. WHEN a feature is recognized but not yet implemented (such as real Wheel topology), THE MfsError SHALL use the `NotImplemented` variant
4. WHEN a precondition check fails (such as wrong input topology for trisection extraction), THE MfsError SHALL use the `PreconditionViolation` variant
5. THE MfsError enum SHALL preserve backward-compatible `Display` formatting for all existing variants that remain unchanged
6. THE MfsError SHALL derive `Serialize` so that error details can be included in JSON output for AI and CLI consumers

### Requirement 7: 耦合矩阵模块拆分

**User Story:** As a library maintainer, I want the coupling matrix module split into focused submodules, so that each concern is independently navigable and testable.

#### Acceptance Criteria

1. THE Coupling_Matrix_Module SHALL separate core data structure and accessor methods into a dedicated `core` submodule
2. THE Coupling_Matrix_Module SHALL separate Section_Extraction operations (extract_triplet, extract_quadruplet, extract_trisection) into a dedicated `sections` submodule
3. THE Coupling_Matrix_Module SHALL separate Rotation_Helpers into a dedicated `rotations` submodule
4. THE Coupling_Matrix_Module SHALL separate denormalization/scaling logic into a dedicated `scaling` submodule
5. WHEN the module is split, THE Coupling_Matrix_Module SHALL preserve all existing public and `pub(crate)` API signatures

### Requirement 8: 移除 Placeholder 回退路径

**User Story:** As a library maintainer, I want the placeholder matrix fallback removed, so that synthesis failures are reported explicitly rather than producing incorrect matrices silently.

#### Acceptance Criteria

1. THE Synthesis_Engine SHALL remove the `synthesize_placeholder_matrix` function and its module
2. WHEN the residue expansion path fails, THE Synthesis_Engine SHALL return an explicit `NumericalFailure` error
3. WHEN generalized data is missing, THE Synthesis_Engine SHALL return an explicit `PreconditionViolation` error
4. THE Synthesis_Engine SHALL remove the `MatrixSynthesisMethod::PlaceholderFallback` variant
5. IF removing the placeholder path causes existing tests to fail, THEN THE Synthesis_Engine SHALL update those tests to provide valid generalized data or assert the expected error

### Requirement 9: 拓扑变换职责边界与样板消除

**User Story:** As a library maintainer, I want a clear single-owner for transform logic with minimal boilerplate, so that adding new topologies or section types is straightforward.

#### Acceptance Criteria

1. THE Transform_Module SHALL be the sole public entry point for all topology conversion and section extraction operations
2. THE Coupling_Matrix_Module SHALL expose transform operations only through `pub(crate)` methods
3. THE Transform_Module SHALL extract the shared response-check attachment logic into a single reusable helper (currently duplicated between `mod.rs` and `sections.rs`)
4. THE Transform_Module SHALL unify `VerifiedSectionSynthesis` and `SectionTransformOutcome` into a single canonical outcome type
5. THE Transform_Module SHALL consolidate the `_matrix`, plain, and `_with_response_check` function variants into a builder or options pattern to reduce public API surface

### Requirement 10: 矩阵运算使用 nalgebra 替代手写循环

**User Story:** As a library maintainer, I want matrix multiplication to use the existing nalgebra dependency, so that the hand-rolled triple loop is replaced by an optimized and tested implementation.

#### Acceptance Criteria

1. THE Coupling_Matrix_Module SHALL replace the hand-written `multiply` method with nalgebra `DMatrix<f64>` multiplication
2. THE Coupling_Matrix_Module SHALL replace the hand-written `transpose` method with nalgebra transpose
3. WHEN nalgebra is used, THE Coupling_Matrix_Module SHALL maintain numerical equivalence within floating-point tolerance (max deviation < 1e-12)
4. THE Coupling_Matrix_Module SHALL keep the internal `DMatrix` conversion private and expose only the existing `Vec<f64>` storage to public consumers

### Requirement 11: 消除 unwrap_or_default 滥用

**User Story:** As a library maintainer, I want matrix access to fail explicitly on out-of-bounds indices in internal code, so that bugs are caught early rather than silently producing zero values.

#### Acceptance Criteria

1. THE Coupling_Matrix_Module SHALL provide an internal accessor that panics on out-of-bounds access in debug builds
2. WHEN internal rotation and transform code accesses matrix entries, THE Rotation_Helpers SHALL use the panicking accessor instead of `unwrap_or_default()`
3. THE Coupling_Matrix_Module SHALL retain the existing `at()` method returning `Option<f64>` for public API consumers
4. WHEN compiled in release mode, THE internal accessor SHALL use unchecked indexing for performance parity

### Requirement 12: Wheel 拓扑诚实标记与 SynthesisOutcome 清理

**User Story:** As a library consumer, I want honest API markers and no misleading constant-returning methods, so that I can trust the library's self-reported state.

#### Acceptance Criteria

1. THE Transform_Module SHALL mark `TopologyKind::Wheel` with `#[deprecated]` and a message directing users to Arrow until a real implementation exists
2. THE SynthesisOutcome SHALL remove the `used_generalized_approximation()` method that unconditionally returns `true`
3. THE SynthesisOutcome SHALL remove the `approximation_kind()` method that unconditionally returns a string literal
4. WHEN callers need to identify the approximation path, THE SynthesisOutcome SHALL expose the information through a typed enum field

