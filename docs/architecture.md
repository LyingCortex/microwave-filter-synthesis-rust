# MFS Rust Library Architecture

## Overview

`mfs` is the Rust rewrite target for a microwave filter synthesis library.
The goal is to move from a Python prototype into a strongly typed, testable,
and performance-oriented core library while preserving the domain workflow:

`Filter specification -> approximation -> coupling matrix -> response`

At the current stage, the crate provides a stable skeleton and a minimal
end-to-end pipeline. The numerical algorithms are intentionally incomplete;
the main purpose of this document is to lock down the architectural shape
before deeper synthesis logic is ported.

## Design Goals

- Build a reusable synthesis core instead of a script-oriented toolkit.
- Keep domain concepts explicit in types rather than implicit in object state.
- Prefer immutable value objects for specifications, polynomial sets, and
  coupling matrices.
- Separate mathematical synthesis, physical frequency mapping, and response
  evaluation into independent modules.
- Make the core library suitable for future Python bindings and CLI tooling.
- Keep errors structured and recoverable.

## Non-Goals

- Reproducing the full Python feature set in the first pass.
- Adding plotting, GUI, or notebook-specific helpers in the core crate.
- Binding to a specific matrix backend before the algorithmic interfaces are
  stable.

## Architectural Style

The crate follows a layered domain-core structure:

1. `spec`
   Encodes engineering intent and filter requirements.
2. `freq`
   Maps physical frequency values into normalized prototype space.
3. `approx`
   Produces low-pass prototype polynomial data from a specification.
4. `matrix`
   Holds coupling matrix data and transformation primitives.
5. `response`
   Evaluates electrical response from a coupling matrix over a frequency grid.
6. `error`
   Defines shared error types used across all layers.
7. `synthesis`
   Provides orchestration helpers that compose the lower-level stages.

The top-level crate root exposes a thin facade over these modules, including
`synthesize_chebyshev(...)` as a simple high-level entry point.

## Module Responsibilities

### `spec`

Purpose:
Capture the user-facing design intent with validation close to data creation.

Current types:

- `FilterSpec`
- `FilterType`
- `TransmissionZero`

Responsibilities:

- Represent order, return loss, filter family, and transmission zeros.
- Enforce basic construction invariants such as positive order and return loss.
- Remain independent from numerical backends and response solvers.

Future expansion:

- Ripple and passband specification models
- Cross-coupling constraints
- Topology requirements
- Source/load termination models

### `freq`

Purpose:
Separate physical frequency plans from normalized low-pass prototype math.

Current types:

- `FrequencyPlan` trait
- `LowPassPlan`
- `BandPassPlan`
- `FrequencyGrid`
- `NormalizedSample`

Responsibilities:

- Convert `Hz` inputs into normalized frequency coordinates.
- Hold domain-specific mapping formulas.
- Offer reusable grid generation for response evaluation.

Design note:
This layer is intentionally independent from approximation and matrix code so
multiple synthesis and evaluation engines can share the same mapping logic.

Future expansion:

- `HighPassPlan`
- `BandStopPlan`
- Multi-band mapping
- Reverse mapping from normalized domain back to physical domain
- Unit-safe wrappers if frequency unit handling becomes richer

### `approx`

Purpose:
Convert filter specifications into prototype polynomial objects.

Current types:

- `ApproximationEngine` trait
- `ChebyshevApproximation`
- `PolynomialSet`
- `PrototypePoint`

Responsibilities:

- Encapsulate approximation-family-specific synthesis logic.
- Return typed polynomial artifacts instead of mutating shared state.
- Keep approximation methods swappable through trait-based APIs.

Design note:
The initial implementation uses placeholder coefficients. That is deliberate:
the crate shape should stabilize before Cameron/Amari equations are ported in.

Future expansion:

- Generalized Chebyshev polynomial synthesis
- Transmission-zero normalization
- Polynomial utilities for root extraction and stability checks
- Other approximations such as Butterworth or elliptic variants

### `matrix`

Purpose:
Represent the coupling matrix as a stable domain value object.

Current types:

- `CouplingMatrix`
- `CouplingMatrixBuilder`
- `MatrixShape`

Responsibilities:

- Store coupling matrix data with dimension validation.
- Provide a foundation for future topology and rotation transforms.
- Decouple matrix ownership from solver implementation details.

Design note:
The current representation uses a flat `Vec<f64>`. This keeps the core crate
lightweight while we decide whether to remain with a manual dense layout or
adopt a backend such as `nalgebra` later.

Future expansion:

- Topology transforms: folded, arrow, wheel
- Sparse or structured matrix representations
- Source/load node tagging
- Matrix rotation and extraction operations

### `response`

Purpose:
Evaluate S-parameter responses from a coupling matrix over a frequency grid.

Current types:

- `ResponseSolver`
- `SParameterResponse`
- `ResponseSample`

Responsibilities:

- Accept coupling matrix plus frequency grid input
- Return typed response samples
- Isolate numerical solver logic from synthesis logic

Design note:
The current solver is a placeholder. Later it should evolve into a proper
complex linear algebra solve, likely with a dedicated numeric backend.

Future expansion:

- Complex-valued solver backend
- `S11`, `S21`, group delay, and phase
- Normalized-domain and physical-domain evaluation
- Export adapters for Touchstone and Python ecosystems

### `error`

Purpose:
Provide a consistent error model across the crate.

Current types:

- `MfsError`
- `Result<T>`

Responsibilities:

- Keep failures explicit and typed
- Avoid ad hoc `String` errors in public APIs
- Support ergonomic propagation with `Result`

Future expansion:

- Numerical convergence failures
- Unsupported topology errors
- Invalid synthesis-state transitions

## Data Flow

The intended happy-path flow is:

1. Construct a validated `FilterSpec`
2. Choose a `FrequencyPlan`
3. Run an approximation engine to obtain a `PolynomialSet`
4. Build or synthesize a `CouplingMatrix`
5. Evaluate response across a `FrequencyGrid`

In code, the facade currently looks like this:

```rust
let spec = FilterSpec::chebyshev(4, 20.0)?
    .with_transmission_zeros(vec![TransmissionZero::finite(-1.25)]);
let plan = BandPassPlan::new(6.75e9, 300.0e6)?;

let (polynomials, matrix) = synthesize_chebyshev(&spec, &plan)?;
let grid = FrequencyGrid::linspace(6.0e9, 7.0e9, 201)?;
let response = ResponseSolver::default().evaluate(&matrix, &grid)?;
```

This is intentionally functional in style: each stage returns data for the next
stage rather than mutating a long-lived filter object with hidden intermediate
state.

## Why This Shape Fits Rust Well

Compared with the Python prototype, Rust is a strong fit for this problem
because:

- domain invariants can be enforced at construction time
- value ownership is clear
- numerical kernels can be made efficient without changing the public API
- testing can focus on pure functions and deterministic value objects
- future Python bindings can wrap a stable native core rather than reimplement
  logic

The architecture therefore avoids a large mutable `BaseFilter` object and
instead favors small validated types plus explicit stage outputs.

## Public API Strategy

The crate should expose two API levels:

### Low-level API

For advanced users and future internal composition:

- construct `FilterSpec`
- choose `FrequencyPlan`
- call approximation engine directly
- build or transform a `CouplingMatrix`
- run `ResponseSolver`

### High-level API

For common workflows:

- `synthesize_chebyshev(...)`
- `synthesize_and_evaluate_chebyshev(...)`
- `ChebyshevSynthesis`
- future convenience functions for common prototype and topology patterns

This split keeps the crate ergonomic without hiding important engineering
artifacts.

## Current High-Level Entry Points

The crate currently supports two ergonomic orchestration styles:

```rust
let (polynomials, matrix) = mfs::synthesize_chebyshev(&spec, &plan)?;
```

and

```rust
let outcome = mfs::ChebyshevSynthesis::default()
    .synthesize_and_evaluate(&spec, &plan, &grid)?;
```

The second form is the better long-term home for orchestration because it keeps
the top-level crate root lighter and leaves room for additional synthesis
families later.

## Dependency Strategy

Current dependency policy:

- keep the crate dependency-free while interfaces are still settling

Recommended next step:

- evaluate `nalgebra` once the response solver and matrix transforms become
  real implementations

Reasoning:

- early over-commitment to a numeric backend can leak backend choices into the
  public API
- the current phase benefits more from architectural clarity than from backend
  optimization

## Mapping From the Python Prototype

The original Python prototype suggests the following conceptual mapping:

- `base.py` -> split between `spec`, facade helpers, and orchestration
- `physics/frequency_plan.py` -> `freq`
- `filters/chebyshev.py` -> `approx` plus high-level synthesis helpers
- `response/sparameters.py` -> `response`

This split is intentional. In the Python version, workflow and state were
starting to converge into a single object model. In Rust, it is cleaner to
turn those stages into explicit modules with typed boundaries.

## Planned Evolution

### Phase 1

- finalize type layout and module ownership
- keep placeholder algorithm implementations
- add more unit tests for invariants and API behavior

### Phase 2

- port frequency mapping formulas completely
- implement generalized Chebyshev approximation
- add validated polynomial and transmission-zero helpers

### Phase 3

- implement real coupling matrix synthesis
- add topology transforms
- choose and integrate a complex matrix backend

### Phase 4

- add Python bindings through `pyo3`
- add import/export compatibility helpers
- add reference-validation tests against the Python prototype

## Testing Strategy

The crate should eventually have three test layers:

1. Construction tests
   Validate specs, plans, and matrix dimensions
2. Numerical unit tests
   Validate mapping formulas, polynomial coefficients, and matrix transforms
3. Reference regression tests
   Compare Rust outputs against trusted benchmark data or the Python prototype

At the current stage, the crate includes only minimal smoke tests to ensure the
pipeline compiles and basic API contracts hold.

## Open Decisions

- Whether `CouplingMatrix` should remain backend-neutral or become a thin wrapper
  around `nalgebra`
- Whether complex arithmetic enters at the matrix layer or only inside response
  solving
- Whether topology transforms belong in `matrix` or deserve a dedicated
  `topology` module
- Whether unit-safe newtypes for Hz, GHz, and normalized frequency are worth the
  extra API surface

## Recommended Near-Term Work

- complete the frequency mapping layer first
- then port generalized Chebyshev synthesis
- then implement real response solving
- only after that decide the matrix backend

This order keeps the design grounded in the actual synthesis flow and avoids
locking in low-level implementation choices too early.
