# MFS Rust Development Guide

## Purpose

This document is the standalone development guide for the Rust version of MFS.
It is intentionally written from the perspective of "what architecture best
fits a long-lived Rust library" rather than "how to mirror the Python codebase".

The goal is to make `mfs` a reliable synthesis core that can later support:

- native Rust engineering workflows
- Python bindings
- CLI or batch tools
- benchmark and verification pipelines

The core design principle is simple:

`domain model first, numerical kernels second, integration surfaces last`

If we keep those boundaries clean, algorithms can evolve without rewriting the
public API every time.

## Current Implementation Status

The current crate has already moved in the intended direction:

- transmission-zero normalization now lives in the frequency layer
- matrix synthesis is separated from `CouplingMatrix` storage
- `approx` has been split into directory modules with a polynomial submodule
- `response` has been split into a public facade plus a private backend helper
- the response backend now performs real complex matrix inversion for lossless
  S-parameter evaluation and group-delay extraction
- a generalized Chebyshev helper layer now exists for padded transmission
  zeros, Cameron recursion, `P/F/A/E`-related polynomial steps, and epsilon
  calculation
- the high-level approximation output can carry optional generalized
  polynomial data when the zero set fits the currently supported helper domain

That means the remaining work is no longer about basic architectural cleanup.
The next major gains come from deeper mathematical integration and a more
complete coupling-matrix recovery path.

## Product View

The library should own one narrow responsibility:

`turn validated filter design intent into well-defined synthesis artifacts`

That means the crate is not primarily:

- a plotting toolkit
- a notebook helper package
- a GUI backend
- a direct clone of script-style engineering workflows

Instead, it should expose stable engineering objects and deterministic stages:

1. specification input
2. domain normalization
3. approximation synthesis
4. coupling-matrix synthesis
5. response evaluation
6. export or adapter layers

## Architectural Direction

The most reasonable long-term shape is a layered architecture with explicit
domain ownership.

### Layer 1: Domain

This layer defines the engineering concepts and invariants.

Examples:

- filter order
- filter class
- return loss
- passband definition
- transmission zeros
- topology intent
- frequency units and normalized variables

Rules:

- no solver logic
- no matrix-backend details
- constructors enforce invariants
- types should be cheap to inspect and hard to misuse

### Layer 2: Algorithm Kernels

This layer contains the mathematical transformations.

Examples:

- frequency mapping
- generalized Chebyshev polynomial generation
- polynomial normalization
- coupling-matrix extraction
- topology transforms
- response solving

Rules:

- accept typed inputs from the domain layer
- return typed outputs, not side effects
- avoid hidden mutable state
- keep numerical backend choices behind local abstractions

### Layer 3: Workflow / Orchestration

This layer composes kernels into ergonomic use cases.

Examples:

- `ChebyshevSynthesis`
- convenience end-to-end helpers
- future recipe-style builders

Rules:

- minimal logic
- no duplicate math
- mostly wiring, validation, and error propagation

### Layer 4: Adapters

This layer is optional and should stay outside the core math model when
possible.

Examples:

- `pyo3` bindings
- serialization helpers
- Touchstone export
- CLI entry points

Rules:

- adapters depend on the core, never the reverse
- adapter-friendly DTOs are fine, but core types should not be polluted by
  transport concerns

## Recommended Module Layout

The current single-file-per-area layout is fine for bootstrapping, but the
best long-term structure is a directory layout with deeper ownership.

Recommended target shape:

```text
src/
  lib.rs
  error.rs
  domain/
    mod.rs
    spec.rs
    zeros.rs
    frequency.rs
    topology.rs
  approx/
    mod.rs
    chebyshev.rs
    polynomial.rs
    normalization.rs
  matrix/
    mod.rs
    coupling_matrix.rs
    builder.rs
    topology.rs
    synthesis.rs
  response/
    mod.rs
    solver.rs
    sample.rs
  workflow/
    mod.rs
    chebyshev.rs
  adapters/
    mod.rs
    python.rs
    export.rs
```

This is not required immediately, but it is the shape to grow toward once
`src/*.rs` files become crowded.

## What Should Stay Stable

The following ideas are worth keeping even if implementation details change:

- typed `FilterSpec` rather than a mutable mega-object
- a separate frequency-planning concept
- polynomial artifacts as explicit outputs
- coupling matrix as a first-class domain artifact
- orchestration as a thin layer over lower-level stages
- a crate-wide error type

These fit Rust well and will scale better than a stateful object hierarchy.

## What Should Change As The Crate Grows

Several current concepts are useful placeholders but should evolve.

### 1. `spec` should become more expressive

Today `FilterSpec` mixes "family", "shape", and "band intent" in a compact
form. That is fine for now, but a better long-term model is:

- `FilterClass`: low-pass, band-pass, band-stop, high-pass
- `ApproximationFamily`: Chebyshev, Butterworth, elliptic, ...
- `PerformanceSpec`: return loss, ripple, attenuation targets
- `TransmissionZeros`: explicit typed set with domain-aware constructors

This keeps the semantic axes independent.

### 2. `FrequencyPlan` should move closer to domain types

Frequency mapping is not just a utility trait. It is part of the domain model.
In the long run, frequency plans should probably sit under `domain::frequency`
and expose:

- physical-domain validation
- normalized-domain mapping
- reverse mapping when mathematically valid
- domain-specific helpers for zero normalization

The important part is to avoid ambiguous "is this value normalized or in Hz?"
APIs.

### 3. polynomial work deserves its own home

`PolynomialSet` is central enough that it should not remain just a bag of
vectors forever. It should eventually gain:

- representation invariants
- degree checks
- root and symmetry helpers
- normalization conventions
- real vs. complex coefficient policy

That suggests a dedicated `approx::polynomial` module.

Current status:

- `PolynomialSet` now validates its base coefficient layout
- `generalized_chebyshev` now holds richer complex-polynomial helper artifacts
- the next step is deciding whether the generalized artifact should become the
  primary prototype representation rather than an optional attachment

### 4. matrix storage and matrix synthesis should be separated

Earlier bootstrap code mixed artifact storage with synthesis logic via
`CouplingMatrix::from_polynomials(...)`. The crate now uses a separate
matrix-synthesis boundary, and future work should continue in that direction.

Recommended split:

- `CouplingMatrix`: validated storage and accessors only
- `CouplingMatrixBuilder`: construction helper
- `MatrixSynthesizer` or `CouplingSynthesizer`: turns polynomial artifacts into
  a matrix
- `TopologyTransform`: converts one equivalent form into another

This keeps the matrix type lightweight and makes synthesis strategies swappable.

### 5. response solving should have a real numeric boundary

Response code will likely be the first place where backend choice matters.
When that happens, keep the dependency behind an internal boundary:

- public API returns plain Rust structs
- internal solver may use `nalgebra` or another backend
- backend-specific types should not leak into public signatures

That avoids locking the crate into a math library at the API level.

Current status:

- the solver now computes lossless response by explicitly inverting the
  coupling-matrix response operator
- `group_delay` is available on response samples
- `ResponseSettings` carries source and load resistance without leaking solver
  internals into the rest of the crate

## Recommended Public API Shape

The crate should expose two API levels.

### Low-level API

This is for advanced users, testing, and future bindings.

Examples:

- construct `FilterSpec`
- construct frequency plans
- run an approximation engine directly
- synthesize a coupling matrix explicitly
- evaluate a response solver explicitly
- inspect optional generalized Chebyshev artifacts when available

### High-level API

This is for the common "give me a filter flow" use case.

Examples:

- `synthesize_chebyshev(...)`
- `synthesize_and_evaluate_chebyshev(...)`
- orchestration structs such as `ChebyshevSynthesis`

Guideline:

The high-level API should be a thin convenience layer over the low-level API,
not a separate implementation path.

## Data Ownership Rules

For this project, the safest default is:

- inputs are immutable value types
- algorithms borrow inputs and return new outputs
- mutation is local to builders or solver internals
- public APIs do not require shared mutable state

This gives us reproducibility and easier testing.

Avoid introducing:

- hidden caches in public structs
- stageful objects that become invalid after partial use
- interior mutability unless profiling proves it is necessary

## Error Model

The current unified error type is the right idea. The next refinement should be
more structured categories, not more string payloads.

Recommended categories:

- input validation errors
- unsupported configuration errors
- synthesis failure errors
- numerical stability or convergence errors
- adapter or export errors

A good rule is:

if the caller can reasonably react differently, it deserves its own variant.

## Numeric Strategy

Do not choose a heavy math backend too early.

Recommended policy:

1. keep public artifacts backend-neutral
2. keep dense manual storage where logic is still moving
3. introduce a backend only when a real solver or transform needs it
4. keep backend use internal unless there is a compelling API reason

Most likely path:

- stay with plain `Vec<f64>` for storage while synthesis is immature
- adopt a mature dense linear algebra backend when response solving becomes
  real
- consider complex arithmetic entering only at the solver boundary

## Testing Strategy

The project should treat testing as part of architecture, not cleanup work.

Recommended test layers:

### 1. Invariant tests

Test constructors and basic object validity.

Examples:

- invalid order
- invalid return loss
- invalid frequency ranges
- invalid transmission-zero values
- matrix dimension mismatches

### 2. Formula tests

Test small deterministic mathematical pieces.

Examples:

- frequency mapping formulas
- ripple-factor computation
- polynomial generation helpers
- Cameron recursion
- generalized Chebyshev helper polynomials
- transform identities

### 3. Artifact tests

Test properties of outputs, not only exact numbers.

Examples:

- polynomial degree
- matrix symmetry
- source/load placement
- normalized center frequency mapping to zero

### 4. Regression tests

Freeze benchmark fixtures and compare with tolerances.

Examples:

- known low-order synthesis cases
- selected transmission-zero layouts
- representative response curves
- selected generalized Chebyshev helper outputs against Python reference cases

### 5. Cross-validation tests

Compare against trusted external references only after the internal model is
clear. The old Python code is useful here as a reference implementation, not as
an architectural template.

## Recommended Development Order

If we want the cleanest path with the least rework, the implementation order
should be:

1. finalize domain vocabulary and constructors
2. finish frequency mapping and zero-normalization semantics
3. build robust polynomial infrastructure
4. integrate generalized Chebyshev helper outputs into the main approximation
   path
5. split matrix storage from matrix synthesis
6. add topology transforms
7. implement a real response solver
8. add adapters such as Python bindings

This order keeps dependencies flowing one way and prevents backend details from
shaping the top-level API too early.

## Concrete Near-Term Refactors

These are the most worthwhile refactors before the crate becomes larger:

1. Move `TransmissionZero` normalization policy into a more explicit domain API.
2. Keep matrix synthesis outside `CouplingMatrix` itself and evolve the
   synthesizer boundary as algorithms become more realistic.
3. Let `response` expose only evaluated response types and keep backend details
   private.
4. Keep expanding the generalized Chebyshev helper layer until it can replace
   placeholder approximation coefficients in the main path.
5. Continue using directory modules such as `src/approx/` and `src/response/`
   once a subsystem has multiple responsibilities, and apply the same pattern
   to `matrix` or `spec` when it meaningfully reduces coupling.

## Coding Guidelines

Recommended style for this crate:

- prefer total, validated constructors
- keep domain types small and explicit
- favor pure functions for formulas
- keep orchestration methods short
- avoid trait abstraction unless there is a real second implementation coming
- add comments only where mathematical intent is not obvious

For naming:

- reserve `Spec` for user intent
- reserve `Plan` for domain mapping rules
- reserve `Set` or `Artifact` for algorithm outputs
- reserve `Solver` and `Synthesizer` for active computation types

## Definition Of Done For New Algorithm Work

A new synthesis stage is not complete when the formula merely compiles.
It is complete when:

- the input and output types are clear
- invariants are enforced
- there are direct unit tests for the math
- there is at least one benchmark or regression case
- the public API does not expose accidental implementation details

## Final Recommendation

The best architecture for this Rust library is not an object-oriented
"filter instance that owns everything". It is a typed pipeline of explicit
artifacts:

`Spec -> Normalization -> Approximation -> Matrix -> Response`

That model matches the problem domain, matches Rust's strengths, and leaves
plenty of room for future bindings and more serious numerical kernels.
