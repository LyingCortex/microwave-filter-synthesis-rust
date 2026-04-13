# MFS

`mfs` is an open-source Rust library for microwave filter synthesis.

Repository:
<https://github.com/LyingCortex/microwave-filter-synthesis-rust>

It is building a typed, testable synthesis core around this pipeline:

`Physical zeros (optional) -> normalize_transmission_zeros_hz(...) -> FilterSpec -> Approximation -> Canonical Matrix -> Transform -> Response`

Frequency convention:

- `FilterSpec` stores transmission zeros in normalized prototype coordinates
- if your zeros start in physical Hz, convert them first with
  `normalize_transmission_zeros_hz(...)`
- the current generalized helper path expects real normalized zeros with
  `|zero| >= 1`
- frequency mappings are used for physical-grid evaluation and reporting, not
  for implicit zero normalization inside synthesis
- approximation reads and validates the normalized zeros already stored in the
  spec; it does not normalize physical-frequency zeros itself

The project is intended to become a reusable core for:

- native Rust workflows
- future Python bindings
- batch or CLI tooling
- benchmark and verification pipelines

## Why MFS

The goal is not to wrap a script-oriented workflow in Rust syntax. The goal is
to build a maintainable synthesis core with:

- explicit domain types
- validated stage boundaries
- reusable coupling-matrix artifacts
- testable numerical kernels
- room for future Python and CLI adapters

## Status

The crate is usable as an early-stage open-source research and development
core, but it is not a finished production synthesis library yet.

What exists today:

- typed filter specifications and frequency plans
- builder-style filter-spec construction
- explicit normalized transmission-zero handling
- validated polynomial and coupling-matrix artifacts
- generalized Chebyshev helper routines and staged orchestration reporting
- a separated approximation-internal structure with reusable complex-polynomial
  primitives and generalized-domain helper operations
- a dedicated synthesis subsystem with canonical, residue, and section workflows
- topology-aware coupling-matrix metadata and precondition checks
- topology transform facades for folded, arrow, and wheel workflows
- transform and section-transform reports with optional response-invariance checks
- reusable response-invariance, section, and topology-shape verification helpers
- small literature-backed fixture coverage for generalized and section workflows
- lossless coupling-matrix response solving with group-delay extraction
- examples and test coverage for the current implemented stages

What is still in progress:

- deeper generalized Chebyshev approximation fidelity beyond the current helper-backed path
- coupling-matrix recovery from richer prototype artifacts
- broader advanced-topology coverage beyond folded and arrow
- more benchmark and regression fixtures

## Flow Coverage

Relative to the standard generalized-Chebyshev workflow

`E/F/P generation -> Y-parameter synthesis -> residue expansion -> transversal matrix -> topology transforms`

the current codebase is best described as partially complete:

- `P(s)` generation from normalized transmission zeros is implemented
- generalized `F(s)` generation through Cameron-style recurrence is implemented
- generalized `E(s)` recovery through the current helper-domain root workflow is implemented
- generalized `E/F/P` generation is the only approximation path exposed by the crate
- Y-parameter synthesis, residue expansion, and transversal-matrix recovery are implemented for the generalized helper path
- the matrix stage can still fall back to a placeholder chain-style builder when the residue path is unavailable
- topology conversion to folded, arrow, and wheel forms is implemented as an engineering backend, but the public API does not expose the full similarity-rotation sequence explicitly

The broad pipeline now runs through the generalized-helper route only.

## Current Capabilities

- Build normalized prototype specs with the short facade
  `filter_spec(order, return_loss_db, zeros, unloaded_q)`
- Build specs from physical-frequency zeros by normalizing first with
  `normalize_transmission_zeros_hz(...)`
- Synthesize the default normalized generalized-Chebyshev prototype and canonical matrix with `generalized_chebyshev(&spec)`
- Inspect intermediate generalized-Chebyshev prototype polynomials only when needed with `generalized_chebyshev_polynomials(&spec)`
- Add physical mappings only when needed through `lowpass(...)`,
  `bandpass(...)`, and `generalized_chebyshev_with_response(...)`
- Normalize physical frequencies and transmission zeros into prototype space
- Build validated coupling-matrix artifacts and evaluate lossless response
- Default unloaded Q is `2000.0` when `filter_spec` receives `None`.
- Access generalized Chebyshev helper routines derived from the Python core
- Run end-to-end orchestration through pure functions such as `generalized_chebyshev(...)`
  and `generalized_chebyshev_with_response(...)`
- Observe orchestration-stage details such as `approximation_kind()` and
  `MatrixSynthesisMethod`
- Create typed filter specifications through `FilterSpec::builder()` and
  related low-level APIs
- Use higher-level facades such as `MatrixSynthesisEngine`,
  `CanonicalMatrixSynthesis`, `SectionSynthesis`, and `TransformEngine`
- Use crate-level helpers such as `synthesize_canonical_matrix(...)` and
  `synthesize_matrix_with_topology(...)`
- Use `generalized_chebyshev(&spec)` for the detail-preserving synthesis outcome and
  `synthesize_and_evaluate_generalized_chebyshev(...)` when you want the tuple-style
  synthesize-plus-evaluate facade
- Inspect generalized helper stage details such as `a_stage` and `e_stage`
  when the strict generalized path is used
- Attach topology metadata to matrices and reject invalid advanced-transform inputs
- Request transform and section-extraction reports, including optional response
  comparisons on normalized sweeps
- Reuse small literature-backed fixtures from `mfs::fixtures` in tests and
  examples

## Quick Start

Add the crate to your Rust project once it is published:

```toml
[dependencies]
mfs = "0.1"
```

Current local usage example:

```rust
use mfs::prelude::*;

fn main() -> mfs::Result<()> {
    // `filter_spec(...)` expects normalized prototype zeros.
    let spec = filter_spec(6, 23.0, [-2.0, -1.2, 1.5], None)?;
    let mapping = bandpass(6.75e9, 300.0e6)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 201)?;

    let outcome = generalized_chebyshev(&spec)?;
    println!("prototype order: {}", outcome.polynomials.order);
    println!("matrix order: {}", outcome.matrix.order());

    let evaluated = generalized_chebyshev_with_response(&spec, &mapping, &grid)?;
    let synthesis = &evaluated.synthesis;
    println!("samples: {}", evaluated.response.samples.len());
    println!("approximation stage: {}", synthesis.approximation_kind());
    println!("matrix method: {:?}", synthesis.matrix_method);

    let normalized_grid = FrequencyGrid::linspace(-2.0, 2.0, 81)?;
    let transform = transform_matrix_with_response_check(
        &synthesis.matrix,
        TopologyKind::Folded,
        &normalized_grid,
        ResponseTolerance::default(),
    )?;
    println!("folded report passed: {}", transform.report.passes());
    Ok(())
}
```

The shortest facade keeps specification, approximation, and mapping concerns
separate:

- `filter_spec(...)` builds a normalized prototype spec
- Normalize physical zeros with `normalize_transmission_zeros_hz(...)` before `filter_spec(...)`
- `generalized_chebyshev(&spec)` synthesizes the default normalized generalized-Chebyshev prototype and matrix
- `generalized_chebyshev_polynomials(&spec)` exposes the internal polynomial bundle when you want debug visibility
- `bandpass(...)` or `lowpass(...)` become relevant only when you want a
  physical mapping or response evaluation
- `generalized_chebyshev_with_response(...)` evaluates the synthesized matrix on a physical
  grid without changing the normalized-zero contract of `FilterSpec`

Two common zero-input styles are both supported:

```rust
let normalized = filter_spec(4, 20.0, [-1.5, 2.0], None)?;
let mapping = bandpass(6.75e9, 300.0e6)?;
let physical_zeros = normalize_transmission_zeros_hz([6.72e9, 6.84e9], &mapping)?;
let physical = filter_spec(4, 20.0, physical_zeros, None)?;
```

The lower-level `FilterSpec::builder()` entry point still exists when you want
more explicit control over how normalized transmission zeros are attached:

```rust
let mapping = bandpass(6.75e9, 300.0e6)?;
let zeros = normalize_transmission_zeros_hz([6.72e9, 6.84e9], &mapping)?;
let spec = FilterSpec::builder()
    .order(4)
    .return_loss_db(20.0)
    .normalized_transmission_zeros(zeros)
    .build()?;
```

When the generalized helper path is active, the returned helper data now keeps
not only the final `F/P/A/E` artifacts but also detailed `A`-stage and
`E`-stage intermediates. This makes it easier to debug or fixture-anchor
generalized polynomial construction without re-running lower-level steps by
hand.

Transform and section-transform reports now share the same electrical-check
shape:

```rust
let transform = transform_matrix_with_response_check(
    &matrix,
    TopologyKind::Arrow,
    &grid,
    ResponseTolerance::default(),
)?;

if let Some(comparison) = &transform.report.response.comparison {
    println!("max |S21| deviation: {}", comparison.max_s21_magnitude_deviation);
}

let section = extract_triplet_section_with_response_check(
    &matrix,
    -1.3,
    2,
    &grid,
    ResponseTolerance::default(),
)?;

assert!(section.passes());
```

High-level section synthesis now supports the same style directly:

```rust
let sections = SectionSynthesis::default();
let verified = sections.synthesize_trisection_with_response_check(
    &polynomials,
    -1.25,
    (2, 4),
    &grid,
    ResponseTolerance::default(),
)?;

assert!(verified.passes());
```

Example command:

```powershell
cargo run --example chebyshev_bandpass
cargo run --example quickstart_prelude
cargo run --example quickstart_report
cargo run --example literature_fixtures
```

Quickstart notes:

- `quickstart_prelude` prints a compact formatted terminal summary.
- `quickstart_report` generates `docs/quickstart_report_output.md`.

See also:

- [docs/development-guide.md](docs/development-guide.md)
- [docs/refactor-roadmap.md](docs/refactor-roadmap.md)
- [docs/algorithm-roadmap.md](docs/algorithm-roadmap.md)
- [docs/literature-fixture-catalog.md](docs/literature-fixture-catalog.md)
- [examples/chebyshev_bandpass.rs](examples/chebyshev_bandpass.rs)
- [examples/quickstart_prelude.rs](examples/quickstart_prelude.rs)
- [examples/quickstart_report.rs](examples/quickstart_report.rs)
- [examples/literature_fixtures.rs](examples/literature_fixtures.rs)

## Project Layout

```text
src/
  approx/       Approximation engines, complex-polynomial primitives, and generalized helpers
  fixtures/     Small literature-backed and literature-shaped reference cases
  freq.rs       Frequency plans and normalization helpers
  matrix/       Coupling-matrix data structures, topology metadata, and low-level operations
  prelude.rs    Ergonomic re-exports for common workflows
  response/     Response API and internal solver backend
  spec/         Input-domain filter specification types and builder
  synthesis/    Canonical, residue, section, and orchestration workflows
  transform/    Topology and section transform facades with reports
  verify/       Response, section, and topology verification helpers
```

## Open Source Readiness

The repository already includes:

- dual MIT / Apache-2.0 licensing
- CI workflows for linting, cross-platform checks, docs build, packaging, and
  dependency audit
- contribution, support, and code-of-conduct documents
- issue and pull-request templates
- changelog scaffold
- release checklist and first-release notes draft

The repository metadata now points at the public GitHub repository. You can add
badges later if you want, but the library is already structured to work as a
public open-source project.

## Development

Common commands:

```powershell
cargo check --all-targets
cargo test
```

The repository is intentionally architecture-driven. If you want to understand
the direction first, read the development guide before making larger changes.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT license

at your option.
