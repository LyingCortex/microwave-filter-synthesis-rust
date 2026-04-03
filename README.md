# MFS

`mfs` is an open-source Rust library for microwave filter synthesis.

Repository:
<https://github.com/LyingCortex/microwave-filter-synthesis-rust>

It is building a typed, testable synthesis core around this pipeline:

`Spec -> Normalization -> Approximation -> Matrix -> Response`

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
- transmission-zero normalization
- validated polynomial and coupling-matrix artifacts
- generalized Chebyshev helper routines
- lossless coupling-matrix response solving with group-delay extraction
- examples and test coverage for the current implemented stages

What is still in progress:

- full generalized Chebyshev approximation integration into the main output
- coupling-matrix recovery from richer prototype artifacts
- topology transforms
- more benchmark and regression fixtures

## Current Capabilities

- Create typed filter specifications with explicit approximation family,
  filter class, performance spec, and transmission zeros
- Normalize physical frequencies and transmission zeros into prototype space
- Build validated coupling-matrix artifacts and evaluate lossless response
- Access generalized Chebyshev helper routines derived from the Python core
- Run end-to-end orchestration through `ChebyshevSynthesis`

## Quick Start

Add the crate to your Rust project once it is published:

```toml
[dependencies]
mfs = "0.1"
```

Current local usage example:

```rust
use mfs::{
    BandPassPlan, ChebyshevSynthesis, FilterSpec, FrequencyGrid, TransmissionZero,
};

fn main() -> mfs::Result<()> {
    let spec = FilterSpec::chebyshev(6, 23.0)?.with_transmission_zeros(vec![
        TransmissionZero::normalized(-2.0),
        TransmissionZero::normalized(-1.2),
        TransmissionZero::normalized(1.5),
    ]);

    let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 201)?;

    let outcome = ChebyshevSynthesis::default().synthesize_and_evaluate(&spec, &plan, &grid)?;
    println!("matrix order: {}", outcome.matrix.order());
    println!("samples: {}", outcome.response.samples.len());
    Ok(())
}
```

Example command:

```powershell
cargo run --example chebyshev_bandpass
```

See also:

- [docs/development-guide.md](docs/development-guide.md)
- [docs/refactor-roadmap.md](docs/refactor-roadmap.md)
- [docs/algorithm-roadmap.md](docs/algorithm-roadmap.md)
- [examples/chebyshev_bandpass.rs](examples/chebyshev_bandpass.rs)

## Project Layout

```text
src/
  approx/    Approximation and generalized Chebyshev helpers
  freq.rs    Frequency plans and normalization helpers
  matrix/    Coupling-matrix artifacts, builders, and synthesis
  response/  Response API and internal solver backend
  spec.rs    Input-domain filter specification types
  synthesis.rs
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
