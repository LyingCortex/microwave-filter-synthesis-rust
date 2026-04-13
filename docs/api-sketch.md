# MFS API Sketch

## Purpose

This document sketches the intended public API for the Rust crate and its
future Python-facing wrapper. It is not a stability guarantee yet, but it sets
direction for implementation.

## Rust-Native API

### Core Types

```rust
use mfs::{
    BandPassMapping,
    GeneralizedChebyshevApproximation,
    FrequencyGrid,
    ResponseSolver,
    filter_spec,
    generalized_chebyshev_with_response,
};
```

### Typical Workflow

```rust
let spec = filter_spec(6, 23.0, [-2.0, -1.2, 1.5], None)?;

let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
let approximation = GeneralizedChebyshevApproximation::default();
let polynomials = approximation.synthesize(&spec)?;
let matrix = mfs::synthesize_canonical_matrix(&polynomials)?;
let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 2001)?;
let response = ResponseSolver::default().evaluate(&matrix, &grid, &mapping)?;
```

### Orchestrated Workflow

```rust
let outcome = generalized_chebyshev_with_response(&spec, &mapping, &grid)?;

println!("{}", outcome.response.samples.len());
println!("{}", outcome.response.samples[0].group_delay);
```

## Planned Public Modules

### `mfs::spec`

Expected public items:

- `FilterSpec`
- `TransmissionZero`

### `mfs::freq`

Expected public items:

- `FrequencyMapping`
- `LowPassMapping`
- `BandPassMapping`
- `FrequencyGrid`

### `mfs::approx`

Expected public items:

- `ApproximationEngine`
- `GeneralizedChebyshevApproximation`
- `PolynomialSet`
- generalized Chebyshev helper artifacts for advanced workflows

### `mfs::matrix`

Expected public items:

- `CouplingMatrix`
- `CouplingMatrixBuilder`
- `MatrixSynthesisOutcome`
- `MatrixSynthesisMethod`
- topology transforms and structured extractions

### `mfs::response`

Expected public items:

- `ResponseSolver`
- `ResponseSettings`
- `SParameterResponse`
- `ResponseSample`

### `mfs::synthesis`

Expected public items:

- `synthesize_generalized_chebyshev(...)`
- `synthesize_and_evaluate_with_mapping(...)`
- `SynthesisOutcome`
- `EvaluationOutcome`

## High-Level Helpers

The crate should also expose compact helpers for common workflows:

```rust
let (polynomials, matrix) = mfs::synthesize_generalized_chebyshev(&spec)?;
```

Possible future helpers:

- `synthesize_canonical_matrix(...)`
- `evaluate_bandpass_response(...)`
- `folded_topology(...)`

## Proposed Builder Style

For more advanced configuration, builder-style APIs may be useful:

```rust
let spec = filter_spec(8, 22.0, [-1.8, 1.8], None)?;
```

This style keeps validation localized and avoids long argument lists.

## Error Handling

All fallible public APIs should return:

```rust
mfs::Result<T>
```

Error categories should remain typed:

- invalid specification
- invalid frequency mapping
- dimensional mismatch
- unsupported operation
- numerical failure

## Future Python Binding Shape

The Python wrapper should aim for an ergonomic API close to the original
prototype while delegating all real logic to Rust.

Illustrative Python-facing API:

```python
import mfs

spec = mfs.filter_spec(
    order=6,
    return_loss_db=23.0,
    zeros=[-2.0, -1.2, 1.5],
)

mapping = mfs.BandPassMapping(center_hz=6.75e9, bandwidth_hz=300e6)
polys, matrix = mfs.synthesize_generalized_chebyshev(spec)
response = mfs.ResponseSolver().evaluate(matrix, grid, mapping)
```

## API Stability Notes

The following are likely to remain stable:

- explicit `FilterSpec`
- explicit frequency-mapping types
- explicit response solver stage
- typed error returns

The following may still change:

- internal polynomial representation
- how generalized Chebyshev helper artifacts are surfaced in the main
  approximation output
- coupling-matrix backend details
- naming of high-level convenience helpers
- topology transform placement
