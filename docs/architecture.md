# MFS Rust Library Architecture

## Overview

`mfs` is a microwave filter synthesis library in Rust. It implements the complete
Cameron generalized Chebyshev coupling matrix synthesis pipeline:

```
FilterSpec → Approximation → Coupling Matrix → Topology Transform → S-Parameter Response
```

The library provides both a programmatic Rust API and a JSON-driven CLI for
automated filter design workflows. All numerical algorithms are fully implemented
and tested — the library is production-ready for generalized Chebyshev filter
synthesis with finite transmission zeros.

## Design Principles

- **Immutable value objects** — Specifications, polynomial sets, and coupling matrices
  are immutable once constructed. Each pipeline stage produces new data rather than
  mutating shared state.
- **Explicit types over implicit state** — Domain concepts (poles, residues,
  topologies, frequency mappings) are encoded in the type system.
- **Structured errors** — All failures are typed `MfsError` variants with descriptive
  messages. No panics in the synthesis path.
- **Layered architecture** — Mathematical synthesis, frequency mapping, and response
  evaluation are independent modules with typed boundaries.
- **Property-based testing** — Core algorithms are validated with 10 formal correctness
  properties using `proptest`, covering the full input space.

## Module Structure

```
src/
├── lib.rs              Crate root: re-exports and legacy compatibility
├── design.rs           High-level FilterDesign API (primary entry point)
├── prelude.rs          Re-exports for common workflows
├── error.rs            Structured error enum (MfsError) and Result alias
├── freq.rs             Frequency mappings (LowPass, HighPass, BandPass, BandStop, MultiBand) and grid generation
├── duplexer.rs         2-channel (Tx/Rx) duplexer synthesis + N-port analysis
├── multiplexer.rs      Generalized N-channel multiplexer synthesis + analysis
├── optimize.rs         Coupling-matrix tuning/extraction (Levenberg-Marquardt)
├── spectral.rs         Spectral factorization (cepstrum + companion eigenvalues)
├── touchstone.rs       Touchstone (.s2p) export and import
├── python.rs           Python bindings (optional, behind `python` feature)
├── approx/             Polynomial approximation
│   ├── mod.rs              Module facade
│   ├── polynomial.rs       PolynomialSet and GeneralizedChebyshevData types
│   ├── complex_poly.rs     ComplexPolynomial with arithmetic, roots, evaluation
│   ├── generalized_chebyshev.rs        Main approximation entry point
│   ├── generalized_chebyshev_helpers.rs  E(s), F(s), P(s) construction
│   └── generalized_ops.rs             Recursive Chebyshev operations
├── spec/               Filter specification types
│   ├── mod.rs              Module facade
│   ├── types.rs            FilterSpec, TransmissionZero
│   └── builder.rs          FilterSpecBuilder with validation
├── matrix/             Coupling matrix representation and operations
│   ├── mod.rs              Module facade and re-exports
│   ├── core.rs             CouplingMatrix struct, topology, shape, accessors
│   ├── builder.rs          CouplingMatrixBuilder with dimension validation
│   ├── rotations.rs        Givens rotations, diagonal rotation angles
│   ├── sections.rs         Triplet, quadruplet, trisection extraction
│   └── scaling.rs          Bandpass denormalization, external Q extraction
├── synthesis/          Coupling matrix synthesis engine
│   ├── mod.rs              Module facade and public exports
│   ├── engine.rs           MatrixSynthesisEngine orchestration
│   ├── residues/           Residue synthesis (split by engine seam)
│   │   ├── mod.rs              Facade re-exporting the three seams
│   │   ├── admittance.rs       Y-parameter polynomials from approximation output
│   │   ├── expansion.rs        Residue expansion, classification
│   │   └── transversal.rs      Transversal coupling-matrix assembly
│   ├── orchestration.rs    High-level synthesis + evaluation helpers
│   ├── recovery.rs         Shared duplexer/multiplexer channel recovery
│   └── sections.rs         SectionSynthesis (triplet/quadruplet/trisection)
├── transform/          Topology transformations
│   ├── mod.rs              Transform facade, response invariance checking
│   ├── folded.rs           Transversal → Folded (Givens rotation sequence)
│   ├── arrow.rs            Transversal → Arrow (first-row reduction)
│   ├── wheel.rs            Wheel topology (deprecated, Arrow alias)
│   └── sections.rs         Section extraction with response verification
├── response/           S-parameter response evaluation
│   ├── mod.rs              ResponseSolver facade
│   ├── backend.rs          LU-based complex matrix solver
│   └── pole_expansion.rs   Fast pole-expansion solver (O(N) per point)
├── verify/             Verification and comparison utilities
│   └── mod.rs              Response comparison, topology pattern matching
├── pipeline/           JSON-driven pipeline orchestration
│   ├── mod.rs              Module facade
│   ├── context.rs          SynthesisContext, PipelineOptions, PipelineMetadata
│   ├── stages.rs           Stage trait + 4 concrete stages
│   ├── request.rs          SynthesisRequest, context_from_request, run_full_pipeline
│   ├── execution.rs        Incremental execution, context persistence
│   └── schema.rs           JSON Schema generation
├── output/             Report rendering
│   ├── mod.rs              Module facade
│   ├── format.rs           Number formatting utilities
│   ├── terminal.rs         Terminal table rendering
│   ├── markdown.rs         Markdown report generation
│   └── report.rs           Report data structures
├── fixtures/           Literature-backed test fixtures
│   ├── mod.rs              Fixture loading
│   └── database.rs         Filter database JSON parsing
└── bin/
    └── mfs_cli.rs      CLI binary (--input, --stage, --resume, --format)
```

## Data Flow

### User-Facing API (FilterDesign)

```mermaid
flowchart LR
    A["FilterDesign::bandpass(6, 23, 6.75e9, 300e6)<br/>.zeros_hz([...])<br/>.synthesize()"] --> B[FilterDesign]
    B --> C[".to_folded()"]
    B --> D[".response(start, stop, points)"]
    B --> E[".scale()"]
    C --> F[CouplingMatrix]
    D --> G[SParameterResponse]
    E --> H[CouplingMatrix<br/>in Hz]
```

### Internal Pipeline

```mermaid
flowchart LR
    A[FilterSpec] --> B[generalized_chebyshev_polynomials]
    B --> C[PolynomialSet]
    C --> D[MatrixSynthesisEngine::synthesize]
    D --> E[CouplingMatrix<br/>Transversal]
    E --> F[transform_matrix]
    F --> G[CouplingMatrix<br/>Folded/Arrow]
    G --> H[ResponseSolver::evaluate]
    H --> I[SParameterResponse]
```

### Response Evaluation (Two Paths)

```mermaid
flowchart TD
    A[CouplingMatrix] --> B{Pole Expansion<br/>from matrix}
    B -->|Success + verified| C[Fast O(N) response]
    B -->|Fallback| D[LU solver O(N³)]
    C --> E[SParameterResponse]
    D --> E
```

### Synthesis Engine Detail

```mermaid
flowchart TD
    P[PolynomialSet] --> A[synthesize_admittance_polynomials]
    A --> B[AdmittancePolynomials<br/>y11, y12, y22, denominator]
    B --> C[residue_expansion<br/>Durand-Kerner root finding]
    C --> D[ResidueExpansion × 3]
    D --> E{classify_residues}
    E -->|Real| F[Real-valued extraction]
    E -->|ComplexPair| G[combine_conjugate_pair]
    F --> H[build_transversal_from_residues]
    G --> H
    H --> I[CouplingMatrix]
```

## Module Responsibilities

### `spec` — Filter Specification

Captures engineering intent with validation at construction time.

**Key types:**
- `FilterSpec` — order, return loss, transmission zeros, unloaded Q
- `TransmissionZero` — finite (normalized value) or infinite
- `FilterSpecBuilder` — fluent builder with field-level validation

**Invariants enforced:**
- Order ≥ 1
- Return loss > 0 dB
- Transmission zeros are finite or infinite; passband zeros (`|zero| < 1`)
  are accepted by `spec` but rejected with a typed error at the
  approximation layer (`generalized_ops::safe_sqrt_term`), which currently
  supports stopband-only zeros

### `freq` — Frequency Mapping

Separates physical frequency handling from normalized prototype math.

**Key types:**
- `FrequencyMapping` trait — `map_hz_to_normalized`, `map_normalized_to_hz`
- `LowPassMapping` — single cutoff frequency
- `HighPassMapping` — single cutoff frequency (`ω = -fc/f`); `f = 0` clamps to `-SINGULAR_CLAMP`
- `BandPassMapping` — center frequency + bandwidth
- `BandStopMapping` — stopband center + width (`ω = -1/Ω_BPF(f)`); `f = f0` clamps to `-SINGULAR_CLAMP`
- `MultiBandMapping` — disjoint, center-sorted passbands; piecewise forward
  map (nearest band by prototype-transform distance), band-local inverse
- `FrequencyGrid` — uniform or custom frequency point sets

### `approx` — Polynomial Approximation

Implements the generalized Chebyshev approximation algorithm.

**Key types:**
- `PolynomialSet` — E(s), F(s), P(s) polynomials with metadata (eps, eps_r, order)
- `ComplexPolynomial` — complex-coefficient polynomial with arithmetic, evaluation, root finding
- `AdaptiveRootSolver` — canonical production chain: **faer companion-matrix
  eigenvalues → Durand-Kerner → Aberth-Ehrlich** (fall-through on error). The
  faer stage is gated on a backward-residual check (≤ 1e-3) and a
  non-finiteness check, so an ill-conditioned companion matrix yields a typed
  error and fall-through rather than silently shipped garbage eigenvalues.
- `FaerRootSolver` / `DurandKernerRootSolver` / `AberthRootSolver` — the three
  solvers in the chain
- `CompanionMatrixRootSolver` — nalgebra Schur companion solver; kept as a
  real-coefficient test oracle only (its complex-coefficient branch can return
  conjugate-flipped non-roots — see the dedicated test in
  `approx::complex_poly`)
- `GeneralizedChebyshevData` — helper polynomials (E_s, F_s, P_s) for synthesis

**Algorithm:**
Recursive Chebyshev function generation with transmission zero injection,
producing the characteristic polynomials needed for coupling matrix synthesis.
Denominator roots are found via `AdaptiveRootSolver` (faer → Durand-Kerner →
Aberth), not Durand-Kerner alone.

### `matrix` — Coupling Matrix

Dense (N+2)×(N+2) real-valued coupling matrix with source and load nodes.

**Key types:**
- `CouplingMatrix` — immutable matrix with topology metadata
- `CouplingMatrixBuilder` — symmetric entry setting with bounds checking
- `MatrixTopology` — Transversal, Folded, Arrow, Wheel (deprecated)
- `BandPassScaledCouplingMatrix` — physical-frequency scaled matrix with external Q

**Operations:**
- Givens rotations (arbitrary axis pairs)
- Bandpass denormalization and normalization
- External Q extraction
- Section extraction (triplet, quadruplet, trisection)
- nalgebra-backed matrix multiply and transpose

### `synthesis` — Matrix Synthesis Engine

Orchestrates the conversion from polynomials to coupling matrix.

**Key types:**
- `MatrixSynthesisEngine` — main synthesis entry point
- `ResidueExpansion` — partial-fraction decomposition result
- `ResidueClassification` — Real or ComplexPair classification
- `AdmittancePolynomials` — y11, y12, y22 numerators + shared denominator
- `SectionSynthesis` — triplet/quadruplet/trisection synthesis workflows

**Algorithm:**
1. Construct admittance polynomials from E, F, P (with leading-coefficient normalization)
2. Find denominator roots via Durand-Kerner
3. Compute residues at each pole for y11, y12, y22
4. Classify residues as real or complex-conjugate pairs
5. Build transversal coupling matrix from classified residues
6. Optionally transform to target topology

### `transform` — Topology Transformations

Converts transversal matrices to physically realizable topologies.

**Key types:**
- `TopologyKind` — target topology enum
- `TransformOutcome` — transformed matrix + verification report
- `SectionTransformOutcome` — section extraction result with verification

**Supported topologies:**
- **Folded** — sequential coupling with cross-couplings (Givens rotation sequence)
- **Arrow** — first-row reduction form
- **Wheel** — deprecated, currently Arrow alias

**Features:**
- Response invariance verification (S-parameter comparison before/after)
- Topology pattern matching (structural sparsity validation)
- Section extraction with response checking

### `response` — S-Parameter Evaluation

Evaluates electrical response from coupling matrix over frequency grid.

**Key types:**
- `ResponseSolver` — configurable solver (lossless/lossy)
- `SParameterResponse` — collection of frequency-indexed S-parameter samples
- `ResponseSample` — S11, S21 (real + imaginary) at one frequency point, plus
  `group_delay: Option<f64>` (explicitly `None` when the solver path cannot
  compute it, e.g. the pole-expansion fast path — never a silent `0.0`)

**Algorithm:**
Complex matrix inversion method: `[S] = [I] - 2j[R]^(1/2) ([jΩ] - [M])^(-1) [R]^(1/2)`.
A fast O(N) pole-expansion path is used when the matrix is transversal,
cross-checked against the LU solver with automatic fallback.

### `verify` — Verification Utilities

Provides comparison and pattern-matching tools for validation.

**Key types:**
- `ResponseTolerance` — configurable S-parameter comparison thresholds
- `ResponseCheckReport` — pass/fail with deviation metrics
- `MatrixPatternTolerance` — structural zero detection threshold

**Capabilities:**
- Response comparison (max deviation, RMS error)
- Topology pattern verification (folded, arrow sparsity patterns)
- Matrix symmetry checking

### `pipeline` — JSON-Driven Orchestration

Composable stage-based pipeline with JSON I/O and incremental execution.

**Key types:**
- `Stage` trait — `fn execute(&self, context: &mut SynthesisContext) -> Result<()>`
- `SynthesisContext` — accumulates artifacts across stages
- `SynthesisRequest` — JSON-parseable input specification
- `PipelineMetadata` — execution timing, warnings, version info

**Concrete stages:**
1. `ApproximationStage` — FilterSpec → PolynomialSet
2. `MatrixSynthesisStage` — PolynomialSet → CouplingMatrix
3. `TopologyTransformStage` — CouplingMatrix → TransformOutcome
4. `ResponseEvaluationStage` — CouplingMatrix + Grid (+ optional Mapping)
   → SParameterResponse

**Frequency mapping:** when a request carries `mapping` (`bandpass`/`lowpass`),
the response stage evaluates a *physical-Hz* grid through the mapping
(`MappingConfig::build_mapping`); otherwise the grid is treated as already
normalized. Unknown or incomplete mappings are rejected with
`MfsError::InvalidFrequency` (previously the field was accepted but ignored).

**Features:**
- `context_from_request(request)` — canonical request→context prologue
- `run_full_pipeline(request)` — complete execution from JSON
- `run_stage(context, stage_name)` — incremental single-stage execution
- `save_context` / `load_context` — file-based persistence for resume
- `describe_schema()` — JSON Schema for input/output formats

### `output` — Report Rendering

Terminal and markdown report generation for synthesis results.

### `fixtures` — Test Fixtures

Literature-backed reference filter configurations loaded from JSON database.

### `duplexer` — 2-Channel Duplexer

Synthesizes and analyzes a Tx/Rx duplexer sharing a common antenna port.
Combines two band-pass `FilterDesign`s and computes the 3-port S-parameters
(including an ideal-junction model). Provides `ChannelSpec`, `JunctionType`,
`DuplexerSynthesis` (iterative Macchiarella-style perturbation), and
`Duplexer::response`.

### `multiplexer` — N-Channel Multiplexer

Generalizes the duplexer to an arbitrary number of channels over a common
junction. `MultiplexerSynthesis` runs the iterative mutual-loading loop;
`Multiplexer` provides (N+1)-port analysis (`response`, `channel`,
`channel_folded`, `channel_matrices`). Built on `FilterDesign`.

### `optimize` — Coupling-Matrix Optimization

Levenberg-Marquardt nonlinear least-squares for two workflows: **tuning**
(match a target S-parameter response) and **extraction** (recover a matrix from
measured data). `OptimizeConfig` / `OptimizeResult`.

### `spectral` — Spectral Factorization

Polynomial decomposition `R(s) = D(s)·D*(s)` with `D` all-LHP. Two methods:
cepstrum (FFT-based, `rustfft`) and balanced companion-matrix eigenvalues.
Used by the duplexer/multiplexer full-iteration path.

### `touchstone` — Touchstone I/O & Network Operations

`.s2p` text export/import (v1.0/v2.0, RI/MA/DB, Hz/kHz/MHz/GHz,
configurable reference impedance) via `design.save_touchstone` /
`design.to_touchstone` / `read_touchstone`. Backed by the `touchstone`
crate 0.14.5 (import via `Network::new`).

Network-level operations (0.5.0), wrapping the crate's `NetworkBuilder` /
`Network::cascade` behind mfs's typed-error policy:
- `build_network(response, z0)` / `FilterDesign::build_network()` — build a
  live symmetric-reciprocal 2-port `touchstone::Network` (resample /
  Y-Z-ABCD / save).
- `touchstone::cascade(a, b)` / `FilterDesign::cascade(other)` — cascade two
  2-ports (ABCD). Pre-validates rank / z0 / frequency unit / identical
  frequency grid and returns `MfsError` instead of inheriting upstream
  panics or silent index-wise grid pairing.

## CLI Interface

```
mfs_cli [OPTIONS]

Options:
  --input <FILE>     Read JSON request from file (default: stdin)
  --format <FORMAT>  Output format: json (default) or table
  --stage <NAME>     Execute only the named stage
  --resume <FILE>    Resume from a saved context file
```

## Error Model

```rust
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

All errors are `Serialize + Deserialize` for structured JSON error reporting.

## Numerical Thresholds

| Threshold | Value | Purpose |
|-----------|-------|---------|
| Root solver convergence | 1e-12 | Durand-Kerner iteration termination |
| Conjugate pole matching | 1e-8 | Match poles as conjugate pairs |
| Complex residue detection | 1e-6 | Distinguish real from complex residues |
| Pole real threshold | 1e-10 | Classify poles as essentially real |
| Matrix entry realness | 1e-10 | Final validation of real-valued output |
| Structural zero | 1e-6 | Topology pattern verification |
| Leading coefficient cancellation | 1e-20 | Denominator normalization trigger |

## Testing Strategy

### Test Layers

1. **Unit tests** (202 in `src/`, as of 2026-09-28) — per-function validation with known inputs
2. **Property-based tests** (10 properties, 100 iterations each in `tests/synthesis_properties.rs`)
3. **Regression tests** (36 in `tests/synthesis_regression.rs`) — previously-failing configurations
4. **Integration tests** (15 files in `tests/`) — end-to-end pipeline validation
5. **Literature fixtures** — reference values from Cameron textbook and published papers

### Correctness Properties

| # | Property | Validates |
|---|----------|-----------|
| 1 | Complex conjugate residues correctly paired | Residue classification |
| 2 | Residue data preservation (round-trip) | Pipeline fidelity |
| 3 | Transversal matrix entries are real-valued | Complex-to-real conversion |
| 4 | Matrix structural invariants (dimensions, non-zero couplings) | Builder correctness |
| 5 | Synthesis succeeds for all valid configurations | No NumericalFailure |
| 6 | Topology transformation succeeds | Folded transform robustness |
| 7 | Root solver accuracy | Durand-Kerner precision |
| 8 | Backward compatibility for real-residue configs | No regressions |
| 9 | Admittance polynomial degree constraints | Polynomial construction |
| 10 | Admittance polynomial coefficient parity | Symmetry preservation |

## Dependencies

```toml
[dependencies]
nalgebra = "0.33"        # Matrix operations (LU solve, Schur, eigenvalues)
num-complex = "0.4"      # Complex arithmetic for polynomials and response
faer = "0.20"            # Complex companion-matrix eigenvalues (root finding)
rustfft = "6.2"          # FFT for cepstrum spectral factorization
serde = { version = "1", features = ["derive"] }  # Serialization
serde_json = "1"         # JSON parsing for pipeline requests
touchstone = { version = "0.14.5", default-features = false }  # .s2p parsing + NetworkBuilder
clap = { version = "4", features = ["derive"], optional = true }   # CLI (feature "cli")
pyo3 = { version = "0.23", features = ["extension-module"], optional = true }  # Python (feature "python")

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }  # Benchmarking
proptest = "1.6"         # Property-based testing
```

## Public API Levels

### High-Level (common workflows — `FilterDesign`)

```rust
use mfs::prelude::*;

// Synthesize a 4th-order band-pass filter with 2 transmission zeros
let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
    .zeros_hz([6.5e9, 7.0e9])
    .synthesize()?;

let folded = design.to_folded()?;
let response = design.response(6.5e9, 7.0e9, 201)?;

// Network operations (touchstone 0.14+)
let network = design.build_network()?;          // live touchstone::Network
let cascaded = design.cascade(&other_design)?;  // 2-port cascade
```

> The crate-root legacy free functions (`filter_spec`, `generalized_chebyshev`,
> `bandpass`, `lowpass`, `normalize_transmission_zeros_hz`,
> `synthesize_canonical_matrix*`, `synthesize_matrix_with_topology`) are
> `#[deprecated]` in 0.5.0 — use `FilterDesign` or the module APIs
> (`mfs::approx`, `mfs::synthesis`, `mfs::freq`, `mfs::transform`) instead.

### Low-Level (advanced control)

```rust
use mfs::approx::generalized_chebyshev_polynomials;
use mfs::synthesis::{MatrixSynthesisEngine, synthesize_residue_expansions, classify_residues};

let polynomials = generalized_chebyshev_polynomials(&spec)?;
let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
let classifications = classify_residues(&y11, &y12, &y22, 1e-6)?;
let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
```

### Pipeline (JSON-driven)

```rust
use mfs::pipeline::{run_full_pipeline, SynthesisRequest};

let request = SynthesisRequest { order: 4, return_loss_db: 20.0, .. };
let context = run_full_pipeline(request)?;
```

## Supported Filter Configurations

| Configuration | Status |
|---------------|--------|
| All-pole (order 2–8, fixture-anchored) | ✅ Fully supported |
| Finite TZs, \|zero\| ≥ 1.1 (order 2–8, fixture-anchored) | ✅ Fully supported |
| High order (10–30), all-pole + finite TZs | ✅ Numerically stable — see `examples/high_order_stability.rs` (no literature fixtures above order 8) |
| Symmetric TZ pairs | ✅ Fully supported |
| Asymmetric TZ placement | ✅ Fully supported |
| Complex-residue configurations | ✅ Fixed (conjugate pairing) |
| Passband TZs (\|zero\| < 1) | ❌ Not yet supported |
| Repeated poles | ❌ Not yet supported |

## Known Limitations

1. **Passband transmission zeros** (`|zero| < 1`) are not supported — the `safe_sqrt_term`
   computation restricts zeros to the stopband.
2. **Repeated poles** in the residue expansion are not handled — Laurent expansion
   is not implemented.
3. **Wheel topology** is deprecated and currently an alias of Arrow.
4. **N-channel mutual loading fidelity** — the perturbation now propagates
   into the output coupling matrices (P1 fixed, 0.5.0), but the simplified
   frequency-offset perturbation has limited passband fidelity for adjacent
   channels (it does not converge to machine tolerance for well-separated
   bands, and the loaded channel loses some passband peak). The full
   Macchiarella spectral path behind `DuplexerSynthesis::full_iteration(true)`
   is **converged and gate-anchored** at Type-I orders 2–8 and at Type-II
   (resonating node, incl. the §VI degree-40 GSM design) (2026-10-03, plan
   D1–D4); it remains **opt-in** (default is the simplified path) and is not
   yet generalized to N > 2 channels. See `docs/multiplexer-code-review.md`,
   `docs/internal/architecture-review-2026-09-25.md`, and
   `docs/continuation-plan.md` section A.

## Future Roadmap

Worked plan with sequencing and acceptance criteria:
`docs/continuation-plan.md`.

| Phase | Feature | Priority |
|-------|---------|----------|
| Next | Passband transmission zeros | High |
| Next | Repeated-pole residue expansion | High |
| Next | Generalize the full spectral path to N-channel (multiplexer) mutual loading | High |
| Done | Full Macchiarella spectral path converged for the duplexer (Type-I orders 2–8 + Type-II, gate-anchored, 2026-10-03) | — |
| Later | Butterworth/Elliptic approximation families | Low |
| Later | Multi-band synthesis (true multi-passband, scope-1) — design doc first | Medium |

> Touchstone export, Python bindings, and criterion benchmarks are already
> implemented (see `src/touchstone.rs`, `src/python.rs`, `benches/`).
