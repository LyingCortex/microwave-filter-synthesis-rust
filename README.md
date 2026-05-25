# MFS — Microwave Filter Synthesis

[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE-MIT)

A Rust library for microwave coupled-resonator filter synthesis. Implements the
generalized Chebyshev (Cameron) method for coupling matrix generation with
topology transformation and S-parameter evaluation.

## Quick Start

```toml
[dependencies]
mfs = "0.1"
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
- **High-order support** — Numerically stable up to order 30+ with adaptive
  root-finding (Durand-Kerner → Aberth → Companion Matrix fallback)
- **Band-pass scaling** — Convert normalized matrices to physical units
- **Python bindings** (optional) — `pip install .` via maturin

## API Overview

### Band-pass filter (most common)

```rust
let design = FilterDesign::bandpass(order, rl_db, center_hz, bandwidth_hz)
    .zeros_hz([...])       // transmission zeros in Hz
    .unloaded_q(3000.0)    // optional Q factor
    .synthesize()?;

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
    s.group_delay       // group delay
}
```

## Python Bindings

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
```

## Examples

```bash
cargo run --example chebyshev_bandpass
cargo run --example high_order_stability
cargo run --example literature_fixtures
```

## Documentation

- [API Reference](docs/api-reference.md) — Complete method reference
- [Numerical Optimization Report](docs/numerical-optimization-report.md) — Algorithm details

## Project Structure

```
src/
  design.rs     High-level FilterDesign API (start here)
  approx/       Polynomial approximation (Cameron recurrence, root solvers)
  matrix/       Coupling matrix operations and topology transforms
  response/     S-parameter evaluation (LU + pole expansion)
  synthesis/    Matrix synthesis from polynomials
  spec/         Filter specification types
  freq.rs       Frequency mapping helpers
```

## License

Licensed under either of [Apache License 2.0](LICENSE-APACHE) or
[MIT License](LICENSE-MIT) at your option.
