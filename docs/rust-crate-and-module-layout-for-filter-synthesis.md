# Suggested Rust Crate and Module Layout for a Filter Synthesis Library

## Purpose

This document proposes a Rust crate and module layout for an open-source
microwave filter synthesis library.

It complements the higher-level architecture document by answering a more
practical question:

- how should the codebase be organized on disk so that approximation,
  coupling-matrix synthesis, topology transforms, analysis, and verification
  can evolve without becoming tangled

This proposal assumes a library-first project. CLI tools, examples, and future
applications should sit on top of a stable core crate.

## Current Repository Note

The current repository is still a single-crate layout rather than the full
workspace sketched below. Within that single crate, however, several of the
intended boundaries have already landed:

- `src/spec/`
- `src/approx/`
- `src/matrix/`
- `src/response/`
- `src/synthesis/`
- `src/transform/`
- `src/verify/`
- `src/fixtures/`

In particular, `approx` is no longer a single flat file. It already contains:

- `chebyshev.rs`
- `complex_poly.rs`
- `generalized_ops.rs`
- `generalized_chebyshev.rs`
- `polynomial.rs`

So this document should be read as:

- a target workspace layout for future growth
- plus a rough explanation of where the current crate has already converged

## Design Principles

The layout is guided by a few simple rules:

- keep the mathematical core independent from front-end tooling
- separate domain types from algorithms
- isolate numerical infrastructure from topology-specific logic
- make verification and test fixtures easy to extend
- keep public APIs small even if internal modules are numerous

The codebase should favor clarity of responsibility over minimizing file
count.

## Recommended Repository Layout

```text
filter-synthesis/
├─ Cargo.toml
├─ README.md
├─ LICENSE
├─ docs/
├─ examples/
├─ benches/
├─ tests/
├─ crates/
│  ├─ filter-synthesis-core/
│  ├─ filter-synthesis-cli/
│  ├─ filter-synthesis-export/
│  └─ filter-synthesis-visualize/        # optional later
└─ .github/
```

For a small first release, the project could begin with a single crate and
split later. But the long-term structure should anticipate multiple crates.

## Recommended Workspace Strategy

Use a Cargo workspace from the beginning.

Suggested roles:

- `filter-synthesis-core`
  The main library crate containing the domain model and algorithms.
- `filter-synthesis-cli`
  A thin command-line interface for running synthesis pipelines and exporting
  results.
- `filter-synthesis-export`
  Optional helpers for report generation, serialization, and integration
  formats.
- `filter-synthesis-visualize`
  Optional plotting or graph-export helpers if visualization becomes a distinct
  concern.

This keeps the core library reusable in notebooks, services, and other tools
without inheriting CLI or visualization dependencies.

## Core Crate Layout

The most important crate is `filter-synthesis-core`.

Suggested structure:

```text
crates/filter-synthesis-core/
├─ Cargo.toml
└─ src/
   ├─ lib.rs
   ├─ prelude.rs
   ├─ errors/
   │  ├─ mod.rs
   │  ├─ synthesis.rs
   │  ├─ transform.rs
   │  ├─ analysis.rs
   │  └─ validation.rs
   ├─ spec/
   │  ├─ mod.rs
   │  ├─ types.rs
   │  └─ builder.rs
   ├─ approx/
   │  ├─ mod.rs
   │  ├─ chebyshev.rs
   │  ├─ generalized_chebyshev.rs
   │  ├─ generalized_ops.rs
   │  ├─ complex_poly.rs
   │  └─ polynomial.rs
   ├─ math/
   │  ├─ mod.rs
   │  ├─ polynomial.rs
   │  ├─ roots.rs
   │  ├─ residues.rs
   │  ├─ tolerance.rs
   │  ├─ real.rs
   │  └─ complex.rs
   ├─ matrix/
   │  ├─ mod.rs
   │  ├─ coupling_matrix.rs
   │  ├─ indexing.rs
   │  ├─ pattern.rs
   │  ├─ invariants.rs
   │  └─ rotation.rs
   ├─ synthesis/
   │  ├─ mod.rs
   │  ├─ orchestration.rs
   │  ├─ engine.rs
   │  ├─ canonical.rs
   │  ├─ placeholder.rs
   │  ├─ residues.rs
   │  └─ sections.rs
   ├─ transform/
   │  ├─ mod.rs
   │  ├─ folded.rs
   │  ├─ arrow.rs
   │  ├─ wheel.rs
   │  └─ sections.rs
   ├─ response/
   │  ├─ mod.rs
   │  └─ backend.rs
   ├─ verify/
   │  ├─ mod.rs
   ├─ freq.rs
   ├─ error.rs
   └─ fixtures/
      ├─ mod.rs
      └─ ...
```

## `lib.rs` Strategy

`lib.rs` should be intentionally small.

It should:

- expose the stable public modules
- re-export the most important public types
- avoid leaking internal helper structure

Example:

```rust
pub mod analysis;
pub mod approx;
pub mod errors;
pub mod export;
pub mod matrix;
pub mod spec;
pub mod synthesis;
pub mod topology;
pub mod transform;
pub mod verify;

pub mod prelude;

pub use approx::Approximation;
pub use matrix::CouplingMatrix;
pub use spec::FilterSpec;
```

This makes the crate approachable even if the internal organization is rich.

## Module Responsibilities

### `errors/`

Use a dedicated error hierarchy rather than generic strings.

Examples:

- invalid specification
- unsupported approximation configuration
- singular or ill-conditioned synthesis step
- invalid topology precondition
- response mismatch after transform

If the project grows, warning types should also be modeled explicitly.

### `spec/`

This module should define the external input vocabulary.

Keep it free from synthesis-specific details wherever possible. A user should
be able to define a valid filter specification without needing to know whether
the library will synthesize a transversal or folded matrix first.

### `approx/`

This module should own the transformation from electrical specification to
mathematical filtering functions.

Important rule:

- approximation code should not directly manipulate topology transforms

That boundary keeps approximation algorithms replaceable.

### `math/`

This should be the shared numerical foundation layer.

It should contain reusable utilities that are domain-aware enough to be useful
but not so domain-specific that they know about trisections or folded forms.

If this module becomes too large, it can later split into a separate internal
crate.

### `matrix/`

This module should define the `N+2` coupling-matrix object and the low-level
operations that preserve its invariants.

Important contents:

- indexing conventions for source, resonators, and load
- symmetry enforcement
- threshold-aware pattern inspection
- rotation primitives

This is the place to be conservative and heavily tested.

### `synthesis/`

This module should build canonical matrices from approximation results.

The first stable implementation should target:

- canonical transversal form

This module should not contain topology-specific shortcuts that skip canonical
construction unless they are explicitly modeled as alternate algorithms.

### `topology/`

This module should define topology concepts independently of specific
transform implementations.

Examples:

- topology kind enum
- structural constraints
- shortest-path rules
- transmission-zero capacity helpers

This separation is useful because many checks belong to topology semantics,
not to matrix rotation code.

### `transform/`

This module should implement matrix reconfiguration.

Recommended split:

- `engine.rs`
  Generic similarity transform machinery
- `annihilation.rs`
  Helpers for solving rotation angles that eliminate a target entry
- topology-specific files
  Structured algorithms such as folded, arrow, trisection, and box conversion

This avoids stuffing every procedure into one giant transform file.

### `analysis/` or `response/`

This module should compute the electrical behavior implied by a matrix.

Treat this as foundational infrastructure, not as a post-processing extra.
Without it, transform verification will be weak.

In the current crate this role is still implemented under `response/`, with a
public facade plus `backend.rs`. A later rename to `analysis/` would be an API
clarity choice, not an architectural necessity.

### `verify/`

This module should centralize all correctness checks.

Avoid scattering invariance tests and precondition logic across many modules.
Centralizing them makes contributor expectations much clearer.

### `export/`

If export logic stays modest, it can live in the core crate initially.
If it grows significantly, it can move to `filter-synthesis-export`.

Exports should focus on:

- machine-readable data
- human-readable reports
- topology visualization support

### `fixtures/`

This module should hold reusable canonical examples for tests and examples.

Examples:

- known low-order canonical matrices
- arrow-form transition cases
- trisection extraction reference cases
- standard tolerance profiles for numerical tests

Keeping fixtures centralized prevents every test file from inventing its own
half-documented examples.

## Prelude Design

A small `prelude.rs` can improve ergonomics for downstream users.

It should re-export:

- `FilterSpec`
- `Approximation`
- `CouplingMatrix`
- `TopologyKind`
- common result types

It should not re-export every internal helper or the crate will become harder
to evolve safely.

## Suggested Internal Traits

Traits should be introduced carefully and only where they create stable
extension points.

Good candidates:

- `ApproximationMethod`
- `MatrixAnalyzer`
- `TopologyTransform`
- `ExportFormat`
- `Verifier`

Bad candidates for early abstraction:

- generic matrix algebra traits for all numeric backends
- highly abstract solver traits before multiple implementations exist

Start concrete, then abstract when repetition proves the need.

## Testing Layout

The test layout should mirror the architecture.

Recommended structure:

```text
tests/
├─ spec_validation.rs
├─ approximation_chebyshev.rs
├─ canonical_transversal.rs
├─ response_analysis.rs
├─ transform_folded.rs
├─ transform_arrow.rs
├─ transform_trisection.rs
├─ topology_capacity.rs
├─ invariance_regression.rs
└─ export_smoke.rs
```

### Unit Tests

Place close to implementation for:

- indexing rules
- polynomial helpers
- rotation math
- error construction

### Integration Tests

Place in `tests/` for:

- end-to-end synthesis
- topology conversion workflows
- invariance regression
- serialization and export behavior

### Golden or Reference Tests

Use carefully for:

- canonical matrices from known literature cases
- expected transform sparsity patterns
- markdown and JSON export output

Golden tests are especially helpful when the library begins to serve external
users who depend on stable output formats.

## Examples Layout

The `examples/` directory should be educational, not just a dumping ground.

Suggested examples:

- `synthesize_transversal.rs`
- `convert_to_folded.rs`
- `convert_to_arrow.rs`
- `extract_trisections.rs`
- `analyze_response.rs`
- `export_markdown_report.rs`

Each example should demonstrate one user-visible workflow with minimal noise.

## Benchmarks Layout

The `benches/` directory should focus on meaningful performance questions.

Examples:

- polynomial evaluation cost
- matrix rotation sequences
- response analysis scaling with order
- repeated verification in optimization loops

Benchmarks should be added after correctness and API shape stabilize.

## CLI Crate Layout

If a CLI crate is added, keep it thin.

Suggested structure:

```text
crates/filter-synthesis-cli/
├─ Cargo.toml
└─ src/
   ├─ main.rs
   ├─ args.rs
   ├─ commands/
   │  ├─ mod.rs
   │  ├─ synthesize.rs
   │  ├─ transform.rs
   │  ├─ analyze.rs
   │  └─ export.rs
   └─ io/
      ├─ mod.rs
      ├─ read_spec.rs
      └─ write_output.rs
```

The CLI should orchestrate the core crate, not duplicate business logic.

## Dependency Policy

The dependency policy should be conservative.

Prefer:

- small, well-maintained crates
- explicit numerical dependencies
- serialization crates only where needed

Avoid:

- introducing heavyweight symbolic-math dependencies too early
- binding core logic to plotting frameworks
- making the core crate depend on CLI or UI concerns

This matters for open-source adoption and long-term maintenance.

## Visibility Policy

Use Rust visibility intentionally.

Recommended approach:

- keep helper modules `pub(crate)` unless they are true extension points
- expose only stable domain types publicly
- avoid letting tests force internal APIs into the public surface

This reduces future breaking changes.

## Evolution Strategy

The on-disk structure should support gradual growth.

Likely future refactors:

- split `math/` into a dedicated internal crate
- move `export/` into its own crate if formatting features expand
- add a `physical/` module or crate for realization-specific mapping
- add optional `serde` and plotting features behind Cargo features

The initial layout should make these changes easy rather than requiring a full
reorganization later.

## Recommended MVP Layout

For an initial release, it is reasonable to start smaller:

```text
src/
├─ lib.rs
├─ errors.rs
├─ spec.rs
├─ approx.rs
├─ math.rs
├─ matrix.rs
├─ synthesis.rs
├─ transform.rs
├─ analysis.rs
├─ verify.rs
└─ export.rs
```

Then split into submodules as soon as:

- files become too long
- independent contributors begin working in parallel
- different algorithms need separate maintenance paths

The point is not to over-engineer on day one, but to have a clear destination.

## Summary

The recommended Rust organization is:

- a Cargo workspace
- a library-first core crate
- small companion crates for CLI and optional export or visualization features
- clear module boundaries between specification, approximation, matrix model,
  synthesis, transform, analysis, verification, and export

The most important implementation rule is:

- do not let topology-specific algorithms, numerical helpers, and user-facing
  I/O collapse into the same module layer

If the project follows that rule, the codebase can grow from a research-grade
prototype into a maintainable open-source library without constant structural
rewrites.
