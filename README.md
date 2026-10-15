# MFS — Microwave Filter Synthesis

[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE-MIT)

A Rust library for microwave coupled-resonator filter synthesis. Implements the
generalized Chebyshev (Cameron) method for coupling matrix generation with
topology transformation and S-parameter evaluation.

## Quick Start

```toml
[dependencies]
mfs = "0.6"
```

```rust
use mfs::prelude::*;

fn main() -> Result<()> {
    // Design a 6th-order band-pass filter
    let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
        .zeros_hz([6.4e9, 6.5e9, 7.0e9])
        .synthesize()?;

    // Get coupling matrices
    let folded = design.to_folded()?;
    println!("Source coupling: {:.4}", folded.source_coupling());
    println!("Chain couplings: {:?}", folded.chain_couplings());

    // Evaluate S-parameters
    let response = design.response(6.0e9, 7.5e9, 501)?;
    for s in &response.samples {
        println!("{:.3} GHz | S21={:.1} dB | S11={:.1} dB",
            s.frequency_hz / 1e9, s.s21_db(), s.s11_db());
    }

    Ok(())
}
```

## Features

- **Generalized Chebyshev synthesis** — Cameron-style polynomial generation
  with prescribed transmission zeros
- **Coupling matrix topologies** — Transversal, Folded, Arrow
- **Fast S-parameter evaluation** — Pole-expansion method (O(N) per frequency
  point) with automatic LU fallback
- **High-order support** — verified numerically stable through order 28
  (all-pole) and order 24 (finite transmission zeros), RL 20–22 dB, via
  adaptive root-finding (faer companion-matrix eigenvalues → Durand-Kerner →
  Aberth-Ehrlich fallback). Order 30+ may be rejected with a typed
  `NumericalFailure` (the adaptive pole-axis tolerance) — see
  `examples/high_order_stability.rs` for the verified boundary and a
  reproducible rejection.
- **Band-pass scaling** — Convert normalized matrices to physical units
- **Python bindings** (optional) — `pip install .` via maturin

## API Overview

### Band-pass filter (most common)

```rust
let design = FilterDesign::bandpass(order, rl_db, center_hz, bandwidth_hz)
    .zeros_hz([...])       // transmission zeros in Hz
    .unloaded_q(3000.0)    // optional Q factor
    .synthesize()?;

// Or specify by passband edges — center = √(f1·f2), bandwidth = f2 − f1:
// FilterDesign::bandpass_edges(order, rl_db, f_start_hz, f_stop_hz)

let response = design.response(start_hz, stop_hz, points)?;
let folded = design.to_folded()?;
let scaled = design.scale()?;
```

### Normalized prototype

```rust
let design = FilterDesign::prototype(4, 20.0)
    .zeros([-1.5, 1.5])
    .synthesize()?;

let response = design.response_normalized(-3.0, 3.0, 201)?;
```

### Reading S-parameter data

```rust
for s in &response.samples {
    s.frequency_hz      // frequency (Hz)
    s.s21_db()          // |S21| in dB
    s.s11_db()          // |S11| in dB
    s.s21_phase_deg()   // ∠S21 in degrees
    s.group_delay       // Option<f64> normalized-prototype group delay;
                        // computed by both solver paths, None only at a transmission zero
}
```

### Touchstone (.s2p) export

```rust
// One-line save (default: GHz, RI, 50Ω)
design.save_touchstone("filter.s2p")?;

// With configuration
use mfs::touchstone::{FreqUnit, DataFormat};
design.to_touchstone()
    .freq_unit(FreqUnit::MHz)
    .format(DataFormat::DB)
    .save("filter.s2p")?;

// Live touchstone::Network (NetworkBuilder, touchstone 0.14+):
// resample / cascade / Y-Z-ABCD conversion, then save or inspect
let network = design.build_network()?;

// Cascade two designs (port 2 -> port 1) into one 2-port network
let combined = design_a.cascade(&design_b)?;
```

## Python Bindings

Requires Rust ≥ 1.86 (mfs MSRV) and a Python whose ABI matches the build.

```bash
pip install maturin
maturin develop --features python
```

```python
import mfs

design = mfs.bandpass(order=6, rl=23.0, center=6.75e9, bw=300e6,
                      zeros=[6.4e9, 6.5e9, 7.0e9])

freq, s21, s11 = design.response(6.0e9, 7.5e9, 501)
m = design.folded()   # coupling matrix as 2D list

# Touchstone export
design.save_touchstone("filter.s2p")
design.save_touchstone("filter.s2p", freq_unit="MHz", format="DB")
```

## Command-Line Interface

The `mfs_cli` binary (behind the `cli` feature) runs the synthesis pipeline
from a JSON request, or via quick subcommands:

```bash
# Quick design from the command line (band-pass, 6th order, 23 dB, export .s2p)
cargo run --features cli --bin mfs_cli -- design \
    -n 6 -r 23.0 --center 6.75e9 --bandwidth 300e6 \
    --zeros 6.4e9,6.5e9,7.0e9 --topology folded --s2p filter.s2p

# JSON-driven full pipeline (request from a file or stdin)
cargo run --features cli --bin mfs_cli -- -f table request.json
cargo run --features cli --bin mfs_cli -- -s approximation request.json   # one stage only
cargo run --features cli --bin mfs_cli -- --resume context.json           # resume a saved context

# Other subcommands: multiplexer | optimize | extract   (see `mfs_cli --help`)
```

## Examples

```bash
cargo run --example chebyshev_bandpass
cargo run --example high_order_stability
cargo run --example literature_fixtures
cargo run --example multi_band_response
cargo run --example duplexer_macchiarella
```

## Documentation

- [API Reference](docs/api-reference.md) — Complete method reference (Rust)
- [Python API Usage](docs/python-api-usage.md) — `pip install .` / maturin usage
- [Algorithm Reference](docs/algorithms.md) — what each numerical stage actually computes (formulas, conventions, invariants, `file:line` anchors)
- [Architecture](docs/architecture.md) — module layout, data flow, supported configurations, known limitations
- [Numerical Optimization Report](docs/numerical-optimization-report.md) — performance optimizations (O(N) Givens, selective solve, root finders)
- [Offline Algorithm Book](docs/book/mfs-algorithms.html) — 34-chapter + 6-appendix single-file HTML book (math + Rust implementation, fully offline); [verification record](docs/book/verification.md)
- [CHANGELOG](CHANGELOG.md) — per-version change history

## Project Structure

```
src/
  design.rs       High-level FilterDesign API (start here)
  prelude.rs      Convenience re-exports (`use mfs::prelude::*`)
  duplexer.rs     2-channel (Tx/Rx) duplexer synthesis + analysis
  multiplexer.rs  N-channel multiplexer synthesis + analysis
  optimize.rs     Coupling-matrix tuning/extraction (Levenberg-Marquardt)
  spectral.rs     Spectral factorization (cepstrum + companion eigenvalues)
  touchstone.rs   Touchstone (.s2p) export/import + network ops
  error.rs        `MfsError` typed error enum
  approx/         Polynomial approximation (Cameron recurrence, root solvers)
  spec/           Filter specification types
  freq.rs         Frequency mapping helpers
  matrix/         Coupling matrix operations and topology transforms
  synthesis/      Matrix synthesis from polynomials
  transform/      Topology transforms (folded, arrow, sections)
  response/       S-parameter evaluation (LU + O(N) pole expansion)
  verify/         Response verification helpers
  pipeline/       JSON-driven pipeline orchestration (CLI backend)
  output/         Terminal / markdown report rendering
  fixtures/       Literature-backed test fixtures
  bin/mfs_cli/    CLI (design / multiplexer / optimize / extract subcommands)
```

## License

Licensed under either of [Apache License 2.0](LICENSE-APACHE) or
[MIT License](LICENSE-MIT) at your option.
