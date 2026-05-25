# Design Document

## Overview

本设计文档描述 MFS 库架构清理的实现方案。核心设计理念：

1. **Pipeline-as-Data**: 将综合流水线建模为可序列化的上下文对象，每个阶段是纯函数变换
2. **Stage Composition**: 每个阶段是独立的、可替换的、可测试的单元
3. **Structured I/O**: 所有输入输出均为 JSON 可序列化，适合 AI agent 和 CLI 集成
4. **Zero-Cost Abstractions**: 引用传递 + 按需克隆，避免不必要的数据复制

## Architecture

### Pipeline 数据流架构

```
┌─────────────┐     ┌──────────────────┐     ┌─────────────────┐     ┌──────────────────┐     ┌────────────────┐
│ FilterSpec  │────▶│ ApproximationStage│────▶│MatrixSynthesis  │────▶│TopologyTransform │────▶│ResponseEval    │
│ + Mapping   │     │                  │     │Stage            │     │Stage             │     │Stage           │
│ + Options   │     │ Output:          │     │ Output:         │     │ Output:          │     │ Output:        │
│             │     │ PolynomialSet    │     │ CouplingMatrix  │     │ TransformOutcome │     │ SParameter     │
│             │     │                  │     │                 │     │                  │     │ Response       │
└─────────────┘     └──────────────────┘     └─────────────────┘     └──────────────────┘     └────────────────┘
                                    │                    │                      │                       │
                                    ▼                    ▼                      ▼                       ▼
                              ┌──────────────────────────────────────────────────────────────────────────────┐
                              │                        SynthesisContext                                      │
                              │  spec: FilterSpec                                                            │
                              │  mapping: Option<MappingConfig>                                              │
                              │  polynomials: Option<Arc<PolynomialSet>>                                     │
                              │  matrix: Option<Arc<CouplingMatrix>>                                         │
                              │  transform: Option<TransformOutcome>                                         │
                              │  response: Option<SParameterResponse>                                        │
                              │  metadata: PipelineMetadata                                                  │
                              └──────────────────────────────────────────────────────────────────────────────┘
```

### 模块重组后的结构

```
src/
  pipeline/           # 新增：流水线编排与上下文
    mod.rs            # Pipeline trait, SynthesisContext, stage composition
    context.rs        # SynthesisContext 实现
    stages.rs         # Stage trait 定义和具体阶段实现
    schema.rs         # JSON Schema 生成
  approx/             # 不变：近似多项式生成
  error.rs            # 重构：细分错误类型 + Serialize
  fixtures/           # 不变
  freq.rs             # 不变
  lib.rs              # 精简：保留 facade 函数，新增 pipeline 入口
  matrix/             # 重构：拆分子模块
    mod.rs            # 重新导出
    core.rs           # 数据结构 + 访问器
    rotations.rs      # 旋转辅助函数
    sections.rs       # 截面提取
    scaling.rs        # 归一化/反归一化
    builder.rs        # 不变
  output/             # 不变
  prelude.rs          # 更新：加入 pipeline 相关导出
  response/           # 不变
  spec/               # 不变
  synthesis/          # 重构：移除 placeholder
    mod.rs
    engine.rs
    orchestration.rs  # 简化：委托给 pipeline
    residues.rs
    sections.rs
  transform/          # 重构：统一类型，消除样板
    mod.rs
    arrow.rs
    folded.rs
    sections.rs       # 合并 outcome 类型
  verify/             # 不变
```

## Components and Interfaces

### Component 1: Pipeline 模块 (`src/pipeline/`)

**职责**: 定义流水线上下文、阶段 trait、JSON schema 生成

```rust
// src/pipeline/mod.rs
pub mod context;
pub mod stages;
pub mod schema;

pub use context::SynthesisContext;
pub use stages::{Stage, ApproximationStage, MatrixSynthesisStage, 
                 TopologyTransformStage, ResponseEvaluationStage};
pub use schema::describe_schema;

/// 一次性执行完整流水线
pub fn run_full_pipeline(request: SynthesisRequest) -> Result<SynthesisContext> { ... }

/// 从 JSON 字符串执行流水线
pub fn run_from_json(json: &str) -> Result<String> { ... }
```

**Interface: Stage Trait**

```rust
// src/pipeline/stages.rs
pub trait Stage {
    type Input;
    type Output;
    fn execute(&self, input: &Self::Input) -> Result<Self::Output>;
    fn name(&self) -> &'static str;
}
```

**Interface: SynthesisContext**

```rust
// src/pipeline/context.rs
use std::sync::Arc;
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SynthesisContext {
    pub spec: FilterSpec,
    pub options: PipelineOptions,
    pub polynomials: Option<Arc<PolynomialSet>>,
    pub matrix: Option<Arc<CouplingMatrix>>,
    pub transform: Option<TransformOutcome>,
    pub response: Option<SParameterResponse>,
    pub metadata: PipelineMetadata,
}

impl SynthesisContext {
    pub fn polynomials(&self) -> Option<&PolynomialSet> { ... }
    pub fn matrix(&self) -> Option<&CouplingMatrix> { ... }
    pub fn response(&self) -> Option<&SParameterResponse> { ... }
}
```

### Component 2: 错误类型重构 (`src/error.rs`)

**Interface: MfsError**

```rust
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum MfsError {
    InvalidOrder { order: usize },
    InvalidReturnLoss { return_loss_db: f64 },
    InvalidFrequency(String),
    InvalidGridSize { points: usize },
    InvalidTransmissionZero(String),
    DimensionMismatch { expected: usize, actual: usize },
    NumericalFailure(String),
    NotImplemented(String),
    PreconditionViolation(String),
}
```

### Component 3: 耦合矩阵模块拆分 (`src/matrix/`)

**Interface: Module Re-exports**

```rust
// src/matrix/mod.rs
mod builder;
mod core;
mod rotations;
mod scaling;
mod sections;

pub use builder::CouplingMatrixBuilder;
pub use core::{CouplingMatrix, MatrixShape, MatrixTopology, BandPassScaledCouplingMatrix};
```

**Interface: Internal Accessor**

```rust
// src/matrix/core.rs
impl CouplingMatrix {
    /// Debug-panicking accessor for internal use. Panics on OOB in debug, unchecked in release.
    #[inline]
    pub(crate) fn get(&self, row: usize, col: usize) -> f64 {
        debug_assert!(row < self.side() && col < self.side(),
            "matrix access out of bounds: ({row}, {col}) for side {}",  self.side());
        unsafe { *self.data.get_unchecked(row * self.side() + col) }
    }
}
```

**Interface: nalgebra-based Rotation**

```rust
// src/matrix/rotations.rs
pub(crate) fn similarity_transform(
    matrix: &CouplingMatrix, 
    rotation: &CouplingMatrix
) -> CouplingMatrix { ... }
```

### Component 4: Transform 模块统一 (`src/transform/`)

**Interface: Unified TransformOutcome**

```rust
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TransformOutcome {
    pub matrix: CouplingMatrix,
    pub topology: TopologyKind,
    pub report: TransformReport,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TransformReport {
    pub source_topology: TopologyKind,
    pub requested_topology: TopologyKind,
    pub result_topology: TopologyKind,
    pub pattern_verified: bool,
    pub section_verification: Option<SectionVerificationReport>,
    pub response: ResponseCheckReport,
    pub notes: Vec<String>,
}
```

### Component 5: Synthesis 引擎清理

**Interface: Simplified Engine**

```rust
impl MatrixSynthesisEngine {
    pub fn synthesize_with_details(&self, polynomials: &PolynomialSet) -> Result<MatrixSynthesisOutcome> {
        let _generalized = polynomials.generalized.as_ref()
            .ok_or_else(|| MfsError::PreconditionViolation(
                "residue expansion requires generalized Chebyshev data".to_string()
            ))?;
        
        let (y11, y12, y22) = synthesize_residue_expansions(polynomials)?;
        let matrix = build_transversal_from_residues(polynomials, &y11, &y12, &y22)?;
        
        Ok(MatrixSynthesisOutcome {
            matrix,
            method: MatrixSynthesisMethod::ResidueExpansion,
        })
    }
}
```

### Component 6: CLI Adapter

**Interface: CLI Arguments**

```
mfs_cli [OPTIONS] [INPUT_FILE]

Options:
  --input <FILE>     JSON input file (reads stdin if omitted)
  --format <FORMAT>  Output format: json (default) | table
  --stage <STAGE>    Execute only this stage: approximation | matrix | transform | response
  --resume <FILE>    Resume from saved context file
  --output <FILE>    Write output to file instead of stdout
```

## Data Models

### SynthesisRequest (JSON Input)

```json
{
  "order": 4,
  "return_loss_db": 20.0,
  "transmission_zeros": [-2.0, 1.5],
  "unloaded_q": 2000.0,
  "topology": "folded",
  "mapping": {
    "kind": "bandpass",
    "center_hz": 6.75e9,
    "bandwidth_hz": 300e6
  },
  "grid": {
    "start": 6.0e9,
    "stop": 7.5e9,
    "points": 201
  }
}
```

### SynthesisResponse (JSON Output)

```json
{
  "metadata": {
    "version": "0.1.0",
    "stages_executed": ["approximation", "matrix_synthesis", "topology_transform", "response_evaluation"],
    "stage_timings_ms": [1.2, 0.8, 0.3, 5.1],
    "warnings": []
  },
  "spec": { "order": 4, "return_loss_db": 20.0, "transmission_zeros": [-2.0, 1.5] },
  "polynomials": { "order": 4, "epsilon": 0.123, "transmission_zeros_normalized": [-2.0, 1.5] },
  "matrix": { "order": 4, "topology": "folded", "data": [...] },
  "transform": { "topology": "folded", "pattern_verified": true, "response_invariant": true },
  "response": { "samples": [{ "frequency_hz": 6.0e9, "s11_db": -5.2, "s21_db": -0.3 }, ...] }
}
```

### PipelineOptions

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineOptions {
    pub topology: Option<TopologyKind>,
    pub grid: Option<GridConfig>,
    pub mapping: Option<MappingConfig>,
    pub response_tolerance: Option<ResponseTolerance>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MappingConfig {
    pub kind: String,
    pub center_hz: Option<f64>,
    pub bandwidth_hz: Option<f64>,
    pub cutoff_hz: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridConfig {
    pub start: f64,
    pub stop: f64,
    pub points: usize,
}
```

## Error Handling

### Error Categories

| Error Variant | When Used | Example |
|---|---|---|
| `InvalidOrder` | Order is 0 | `filter_spec(0, 20.0, ...)` |
| `InvalidReturnLoss` | RL ≤ 0 or non-finite | `filter_spec(4, -5.0, ...)` |
| `InvalidFrequency` | Frequency validation fails | Negative bandwidth |
| `InvalidGridSize` | Grid < 2 points | `FrequencyGrid::linspace(a, b, 1)` |
| `InvalidTransmissionZero` | Zero placement invalid | Triplet center out of range |
| `DimensionMismatch` | Size mismatch | Matrix data length ≠ (N+2)² |
| `NumericalFailure` | Algorithm diverges | Root solver non-convergence |
| `NotImplemented` | Feature recognized but absent | Real Wheel topology |
| `PreconditionViolation` | Required state missing | Trisection on non-Arrow matrix |

### Error Propagation Strategy

- All public functions return `Result<T, MfsError>`
- Pipeline stages propagate errors upward without wrapping
- JSON output includes full error details when serialized
- CLI outputs error JSON to stderr with exit code 1

## Testing Strategy

### Unit Tests
- Each new submodule (`matrix/core.rs`, `matrix/rotations.rs`, etc.) carries its own `#[cfg(test)]` module
- Stage implementations tested in isolation with mock inputs

### Integration Tests
- Full pipeline round-trip: JSON input → execute → JSON output → validate
- Incremental vs full execution equivalence
- Context serialization round-trip

### Regression Tests
- All existing tests must pass after refactoring (no behavior change)
- Numerical equivalence tests for nalgebra replacement

### Property-Based Tests
- Context serialization round-trip
- nalgebra vs hand-written multiplication equivalence
- Error serialization round-trip

## Correctness Properties

### Property 1: Pipeline Context Round-Trip Serialization
- **What**: For any valid `SynthesisContext`, serializing to JSON and deserializing back produces an equivalent context
- **Type**: Round-trip property
- **Validates: Requirements 1.5**

### Property 2: Full Pipeline Equivalence
- **What**: Running the full pipeline in one call produces the same artifacts as running each stage incrementally
- **Type**: Metamorphic property
- **Validates: Requirements 1.3, 2.1**

### Property 3: nalgebra Multiplication Equivalence
- **What**: For any two valid coupling matrices, the nalgebra-based multiply produces results within 1e-12 of the hand-written loop
- **Type**: Model-based testing (optimized vs reference implementation)
- **Validates: Requirements 10.3**

### Property 4: Error Serialization Round-Trip
- **What**: For any `MfsError` variant, serializing to JSON and deserializing back produces an equivalent error
- **Type**: Round-trip property
- **Validates: Requirements 6.6**

### Property 5: Stage Predecessor Validation
- **What**: Running any stage without its required predecessor artifact in the context always returns a `PreconditionViolation` error
- **Type**: Error condition property
- **Validates: Requirements 1.4, 6.4**

### Property 6: JSON Schema Validation
- **What**: For any valid `SynthesisRequest`, the output of `describe_schema()` validates that request as conforming
- **Type**: Metamorphic property
- **Validates: Requirements 3.4, 3.5**

## Implementation Phases

### Phase 1: Foundation (错误类型 + 矩阵拆分)
1. 重构 `MfsError` 枚举，添加 Serialize derive
2. 拆分 `coupling_matrix.rs` 为子模块
3. 替换 `unwrap_or_default()` 为 debug-panicking 访问器
4. 替换手写矩阵乘法为 nalgebra

### Phase 2: Cleanup (Placeholder 移除 + Transform 统一)
5. 移除 `placeholder.rs` 和 `PlaceholderFallback` 变体
6. 统一 `SectionTransformOutcome` 和 `VerifiedSectionSynthesis`
7. 合并 response-check 样板代码
8. 标记 Wheel 为 deprecated，清理 SynthesisOutcome

### Phase 3: Pipeline Architecture (上下文 + 阶段)
9. 创建 `pipeline/` 模块，定义 `Stage` trait 和 `SynthesisContext`
10. 实现四个具体 Stage 类型
11. 实现 `SynthesisRequest` JSON 输入和 `run_from_json`
12. 实现 `describe_schema()` 和 JSON 输出

### Phase 4: CLI Adapter
13. 添加 CLI binary，支持 JSON stdin/file
14. 实现 `--stage`、`--resume`、`--format` 标志
15. 结构化错误输出到 stderr

## Dependencies and Risks

- **nalgebra 已是依赖**: 无新增依赖风险
- **serde 已是依赖**: Serialize/Deserialize 支持已就绪
- **公共 API 变更**: `Unsupported` 移除是 breaking change，需要 semver bump
- **Placeholder 移除**: 可能影响某些边缘测试用例，需要逐一验证
- **性能**: nalgebra 矩阵乘法对小矩阵（N<10）可能有额外开销，需要 benchmark 验证
- **CLI 依赖**: 需要新增 `clap` 或类似 argument parsing crate

