# MFS Algorithm Roadmap

## Purpose

This document tracks the order in which the mathematical and numerical pieces
of the library should be implemented.

The current crate already has the right architectural buckets. The remaining
work is to fill those buckets with correct algorithms in a sequence that keeps
verification manageable.

## Recommended Order

1. Frequency mapping
2. Transmission-zero normalization
3. Generalized Chebyshev polynomial synthesis
4. Coupling-matrix synthesis
5. Topology transforms
6. Response solving
7. Cross-validation and optimization

## Stage 1: Frequency Mapping

Source inspiration:

- Python `physics/frequency_plan.py`

Target module:

- `src/freq.rs`

Tasks:

- finalize low-pass mapping
- finalize band-pass mapping
- add reverse mapping where needed
- define numeric tolerances for tests

Why first:

- small surface area
- low dependency on other algorithm modules
- essential for both synthesis and evaluation

## Stage 2: Transmission-Zero Normalization

Target modules:

- `src/spec.rs`
- `src/approx/`
- `src/freq.rs`

Tasks:

- decide whether transmission zeros are stored in physical or normalized domain
- add normalization helpers tied to `FrequencyMapping`
- validate finite and special-case zero handling

Key design choice:

The Rust core should avoid ambiguous zero semantics. If a value is physical,
its type or API entry point should make that explicit.

## Stage 3: Generalized Chebyshev Approximation

Target module:

- `src/approx/`

Tasks:

- replace placeholder coefficients with real synthesis logic
- define `E`, `F`, and `P` polynomial representation
- add polynomial validation helpers
- freeze several benchmark cases

Current status:

- transmission-zero padding helpers exist
- Cameron recursion has been ported
- `P` polynomial and epsilon calculation have been ported
- complex-polynomial root solving and `A/E` helper construction now exist
- the remaining work is connecting those helper outputs to the main
  approximation artifact instead of using placeholder `e/f` vectors

Verification ideas:

- known low-order benchmark filters
- symmetry and degree checks
- numerical comparison against Python outputs

## Stage 4: Coupling-Matrix Synthesis

Target modules:

- `src/matrix/`
- possibly a future `src/synthesis.rs`

Tasks:

- decide whether synthesis logic lives in `matrix` or a dedicated module
- implement matrix construction from polynomial data
- model source/load nodes explicitly
- validate matrix dimensions and symmetry where applicable

Open decision:

If coupling synthesis grows substantially, create a dedicated `synthesis`
module rather than overloading `src/matrix/`.

## Stage 5: Topology Transforms

Likely target:

- `src/matrix/` at first
- later maybe `src/topology.rs`

Tasks:

- folded transform
- arrow form transform
- wheel topology transform
- test invariants across equivalent network responses

Why later:

- transforms are easier to trust once base matrix synthesis exists
- response-level equivalence testing becomes possible at this point

## Stage 6: Response Solving

Target module:

- `src/response/`

Tasks:

- choose complex arithmetic approach
- solve the coupling-matrix response equation
- return `S11` and `S21`
- support dense grid evaluation

Current status:

- lossless response solving is now implemented through explicit complex matrix
  inversion
- group delay is also extracted from the inverted response operator
- the main remaining work is numeric hardening, benchmarking, and eventual
  backend strategy decisions

Possible backend options:

- manual dense implementation
- `nalgebra`
- another complex linear algebra crate if justified later

Key caution:

Do not expose backend-specific types in the public API unless there is a very
strong reason.

## Stage 7: Cross-Validation

Tasks:

- compare Rust outputs to Python prototype outputs
- create frozen regression fixtures
- add tolerance-based tests for amplitude and phase
- benchmark high-order cases

Success criteria:

- reproducible outputs
- acceptable numerical stability
- clear mismatch investigation workflow

## Reference Data Plan

The project should eventually maintain a small benchmark set:

- low-order generalized Chebyshev examples
- band-pass mapping sanity cases
- coupling-matrix reference fixtures
- response curves for representative filters

These fixtures should be versioned and used in regression tests.

## Deferred Topics

These are important, but not first-wave priorities:

- optimization-based synthesis extensions
- sparse matrix acceleration
- automatic topology selection
- graphical visualization helpers
- CAD export/import tooling

## Recommended Next Implementation Task

Integrate the generalized Chebyshev helper outputs into the main approximation
artifact model so the crate can stop relying on placeholder `e/f` coefficient
generation.
