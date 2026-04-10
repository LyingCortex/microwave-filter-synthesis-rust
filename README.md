# MFS

`mfs` is an open-source Rust library for microwave filter synthesis.

Repository:
<https://github.com/LyingCortex/microwave-filter-synthesis-rust>

It is building a typed, testable synthesis core around this pipeline:

`FilterSpec -> Normalization -> Approximation -> Canonical Matrix -> Transform -> Response`

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
- transmission-zero normalization
- validated polynomial and coupling-matrix artifacts
- explicit Chebyshev and generalized-Chebyshev specification families
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

## Current Capabilities

- Create typed filter specifications with explicit approximation family,
  filter class, performance spec, and transmission zeros
- Normalize physical frequencies and transmission zeros into prototype space
- Build validated coupling-matrix artifacts and evaluate lossless response
- Access generalized Chebyshev helper routines derived from the Python core
- Run end-to-end orchestration through `ChebyshevSynthesis`
- Observe orchestration-stage details such as `ApproximationStageKind` and
  `MatrixSynthesisMethod`
- Use higher-level facades such as `MatrixSynthesisEngine`,
  `CanonicalMatrixSynthesis`, `SectionSynthesis`, and `TransformEngine`
- Use crate-level helpers such as `synthesize_canonical_matrix(...)` and
  `synthesize_matrix_with_topology(...)`
- Use detail-preserving helpers such as `synthesize_chebyshev_with_details(...)`
  and `synthesize_and_evaluate_chebyshev_with_details(...)`
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
    let spec = FilterSpec::generalized_chebyshev(6, 23.0)?
        .with_transmission_zeros(vec![
        TransmissionZero::normalized(-2.0),
        TransmissionZero::normalized(-1.2),
        TransmissionZero::normalized(1.5),
    ]);

    let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 201)?;

    let outcome = synthesize_and_evaluate_chebyshev_with_details(
        &spec,
        &mapping,
        &grid,
    )?;
    println!("matrix order: {}", outcome.matrix.order());
    println!("samples: {}", outcome.response.samples.len());
    println!("approximation stage: {:?}", outcome.approximation_kind);
    println!("matrix method: {:?}", outcome.matrix_method);

    let normalized_grid = FrequencyGrid::linspace(-2.0, 2.0, 81)?;
    let transform = transform_matrix_with_response_check(
        &outcome.matrix,
        TopologyKind::Folded,
        &normalized_grid,
        ResponseTolerance::default(),
    )?;
    println!("folded report passed: {}", transform.report.passes());
    Ok(())
}
```

`FilterSpec::chebyshev(...)` keeps the default path lightweight and only
attaches generalized helper data when the current helper can support the
normalized transmission-zero pattern. `FilterSpec::generalized_chebyshev(...)`
requests the strict generalized path instead, with no silent fallback.

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
cargo run --example literature_fixtures
```

See also:

- [docs/development-guide.md](docs/development-guide.md)
- [docs/refactor-roadmap.md](docs/refactor-roadmap.md)
- [docs/algorithm-roadmap.md](docs/algorithm-roadmap.md)
- [docs/literature-fixture-catalog.md](docs/literature-fixture-catalog.md)
- [examples/chebyshev_bandpass.rs](examples/chebyshev_bandpass.rs)
- [examples/quickstart_prelude.rs](examples/quickstart_prelude.rs)
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
