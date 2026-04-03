# MFS Migration Plan

## Goal

Move the Python prototype into a Rust-first library without losing the core
filter-synthesis workflow or making the API harder to use for engineering and
research work.

## Principles

- Port by domain stage, not by file-by-file translation.
- Keep the Rust core authoritative for new logic.
- Use the Python prototype as a behavioral reference, not as a structural
  template.
- Prefer small validated types over a large mutable object graph.
- Keep the crate compiling and testable at every step.

## Source-to-Target Mapping

- `base.py`
  Move orchestration ideas into crate-level helpers and typed stage outputs.
- `physics/frequency_plan.py`
  Port into `src/freq.rs`.
- `filters/chebyshev.py`
  Split across `src/spec.rs`, `src/approx/`, and top-level helpers.
- `response/sparameters.py`
  Port into `src/response/`.

## Migration Phases

### Phase 1: Architectural Stabilization

Status:
Started

Tasks:

- create crate skeleton
- define module boundaries
- define error model
- define filter specification types
- add basic smoke tests

Exit criteria:

- crate compiles
- public API shape is coherent
- docs explain intended layering

### Phase 2: Frequency Mapping Port

Tasks:

- fully port low-pass and band-pass mapping formulas
- define reverse mapping interfaces where needed
- add validation for invalid physical frequencies
- add unit tests against hand-calculated reference values

Exit criteria:

- `LowPassMapping` and `BandPassMapping` behavior matches Python formulas
- frequency-grid mapping is numerically verified

### Phase 3: Approximation Port

Tasks:

- port generalized Chebyshev synthesis
- define polynomial utilities
- normalize transmission zeros consistently
- add reference test vectors

Exit criteria:

- `PolynomialSet` contents are meaningful, not placeholders
- Rust outputs match Python prototype for selected benchmark cases

### Phase 4: Coupling Matrix Synthesis

Tasks:

- define synthesis entry points
- port coupling-matrix construction
- decide dense matrix backend strategy
- add topology transform placeholders or real implementations

Exit criteria:

- a realistic coupling matrix can be synthesized from prototype data
- dimensional checks and invariants are enforced

### Phase 5: Response Solver

Tasks:

- select numerical backend if needed
- implement complex solve for response evaluation
- add `S11` and `S21`
- add regression tests against prototype output

Exit criteria:

- response solver returns physically meaningful results
- benchmark cases match expected curves within tolerance

### Phase 6: Ecosystem Integration

Tasks:

- add Python bindings with `pyo3`
- add examples
- add serialization or export helpers
- consider CLI tooling for batch synthesis

Exit criteria:

- Python users can call the Rust core without reimplementing logic
- docs cover both Rust-native and Python-facing workflows

## Validation Strategy

At each phase, validation should happen in this order:

1. type and invariant tests
2. formula-level unit tests
3. regression tests against frozen reference values
4. comparison against Python prototype output where practical

## Risks

- translating equations too literally from Python may preserve design mistakes
- choosing a matrix backend too early may distort the public API
- lack of benchmark data may slow algorithm verification
- topology transforms may need a dedicated module if they grow beyond simple
  matrix operations

## Recommended Immediate Next Step

Implement and verify the frequency mapping layer first. It is the cleanest
piece to port, and it provides a solid base for both approximation and response
work.
