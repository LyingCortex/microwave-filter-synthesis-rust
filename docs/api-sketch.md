# MFS API Sketch

## Purpose

This document sketches the intended public API for the Rust crate and its
future Python-facing wrapper. It is not a stability guarantee yet, but it sets
direction for implementation.

## Rust-Native API

### Core Types

```rust
use mfs::{
    ApproximationFamily,
    BandPassMapping,
    ChebyshevApproximation,
    ChebyshevSynthesis,
    CouplingMatrixSynthesizer,
    FilterClass,
    FilterParameter,
    FrequencyGrid,
    ReturnLossSpec,
    ResponseSettings,
    ResponseSolver,
    TransmissionZero,
};
```

### Typical Workflow

```rust
let spec = FilterParameter::chebyshev(6, 23.0)?
    .with_transmission_zeros(vec![
        TransmissionZero::finite(-2.0),
        TransmissionZero::finite(-1.2),
        TransmissionZero::finite(1.5),
    ]);

let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
let approximation = ChebyshevApproximation::default();
let polynomials = approximation.synthesize(&spec, &mapping)?;

let matrix = CouplingMatrixSynthesizer::default().synthesize(&polynomials)?;
let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 2001)?;
let response = ResponseSolver::default().evaluate_with_settings(
    &matrix,
    &grid,
    &mapping,
    ResponseSettings {
        source_resistance: 1.0,
        load_resistance: 1.0,
    },
)?;
```

### Explicit Semantic Construction

```rust
let spec = FilterParameter::new(
    6,
    FilterClass::BandPass,
    ApproximationFamily::Chebyshev,
    ReturnLossSpec::new(23.0)?,
)?
.with_transmission_zeros(vec![
    TransmissionZero::finite(-2.0),
    TransmissionZero::finite(1.5),
]);
```

### Orchestrated Workflow

```rust
let outcome = ChebyshevSynthesis::default()
    .synthesize_and_evaluate_with_mapping(&spec, &mapping, &grid)?;

println!("{}", outcome.response.samples.len());
println!("{}", outcome.response.samples[0].group_delay);
```

## Planned Public Modules

### `mfs::spec`

Expected public items:

- `FilterParameter`
- `FilterClass`
- `ApproximationFamily`
- `ReturnLossSpec`
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
- `ChebyshevApproximation`
- `PolynomialSet`
- optional generalized Chebyshev helper artifacts for advanced workflows

### `mfs::matrix`

Expected public items:

- `CouplingMatrix`
- `CouplingMatrixBuilder`
- future topology transforms

### `mfs::response`

Expected public items:

- `ResponseSolver`
- `ResponseSettings`
- `SParameterResponse`
- `ResponseSample`

### `mfs::synthesis`

Expected public items:

- `ChebyshevSynthesis`
- `SynthesisOutcome`
- `EvaluationOutcome`

## High-Level Helpers

The crate should also expose compact helpers for common workflows:

```rust
let (polynomials, matrix) = mfs::synthesize_chebyshev(&spec, &mapping)?;
```

Possible future helpers:

- `synthesize_chebyshev_matrix(...)`
- `evaluate_bandpass_response(...)`
- `folded_topology(...)`

## Proposed Builder Style

For more advanced configuration, builder-style APIs may be useful:

```rust
let spec = FilterParameter::chebyshev(8, 22.0)?
    .with_filter_class(mfs::FilterClass::BandPass)
    .with_transmission_zeros(vec![
        TransmissionZero::finite(-1.8),
        TransmissionZero::finite(1.8),
    ]);
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

spec = mfs.FilterParameter.chebyshev(
    order=6,
    return_loss_db=23.0,
    transmission_zeros=[-2.0, -1.2, 1.5],
)

mapping = mfs.BandPassMapping(center_hz=6.75e9, bandwidth_hz=300e6)
polys, matrix = mfs.synthesize_chebyshev(spec, mapping)
response = mfs.ResponseSolver().evaluate(matrix, grid, mapping)
```

## API Stability Notes

The following are likely to remain stable:

- explicit `FilterParameter`
- explicit filter-class / approximation-family split
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
