# MFS Refactor Roadmap

## Purpose

This document turns the architectural guidance into an executable refactor
plan for the current Rust crate.

It is based on the code as it exists today, not on the old Python structure.
The focus is:

- reduce architectural coupling
- make future algorithm work easier
- avoid public API churn later

## Current Snapshot

The crate already has a healthy top-level split:

- `spec`
- `freq`
- `approx`
- `matrix`
- `response`
- `synthesis`
- `error`

That is a good start. The main issue is not missing modules, but mixed
responsibilities inside some of them.

## Progress Status

Completed in the current codebase:

- transmission-zero normalization was moved out of `spec` and into `freq`
- matrix synthesis was separated from `CouplingMatrix` into
  `CouplingMatrixSynthesizer`
- `approx` was reshaped into directory modules with a dedicated polynomial
  boundary
- `response` was reshaped into directory modules with a private backend helper
- the response backend now performs real lossless S-parameter solves and group
  delay extraction
- generalized Chebyshev helper routines from the Python core have been ported
  into `approx::generalized_chebyshev`

Still pending:

- full integration of generalized Chebyshev helper outputs into the main
  approximation coefficients
- stronger polynomial validation and normalization rules
- topology-specific matrix transforms
- coupling-matrix recovery and optimization logic from the Python core

## Main Architectural Tensions

### 1. `spec` knew too much about normalization

This was the first issue addressed. `TransmissionZero::to_normalized(...)` was
removed, and normalization now lives in the frequency-mapping layer.

Why this matters:

- domain objects should stay close to intent
- normalization policy may vary by filter class
- later zero semantics will likely become richer than a single `f64`

Current code points:

- [src/spec.rs](/c:/Users/eynulai/Downloads/mfs/src/spec.rs)
- [src/freq.rs](/c:/Users/eynulai/Downloads/mfs/src/freq.rs)

### 2. `matrix` used to mix storage with synthesis logic

This was also addressed. Matrix synthesis now goes through
`CouplingMatrixSynthesizer` instead of a constructor on `CouplingMatrix`.

Why this matters:

- `CouplingMatrix` should be a validated artifact, not an algorithm host
- later there may be multiple synthesis strategies
- topology transforms do not belong in the same abstraction as raw storage

Current code point:

- [src/matrix/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/matrix/mod.rs)
- [src/matrix/coupling_matrix.rs](/c:/Users/eynulai/Downloads/mfs/src/matrix/coupling_matrix.rs)
- [src/matrix/synthesis.rs](/c:/Users/eynulai/Downloads/mfs/src/matrix/synthesis.rs)

### 3. `approx` now has richer rules, but the main approximation path is still partial

`PolynomialSet` is now validated, and generalized Chebyshev helpers exist for
complex-polynomial stages, but the main approximation output still uses
placeholder `e/f` coefficient generation.

Why this matters:

- polynomial degree and coefficient conventions need explicit validation
- real and complex coefficient handling may diverge later
- the generalized Chebyshev implementation will need helper utilities

Current code point:

- [src/approx/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/mod.rs)
- [src/approx/polynomial.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/polynomial.rs)

### 4. `response` now has a real backend, but it is still an internal solver

The response API is now backed by a real complex matrix inversion path, which
is a major improvement. The remaining question is long-term numeric strategy,
not whether there is a real solver at all.

Why this matters:

- backend choice should stay internal
- response samples should remain stable even if the solver backend changes
- testing the solver gets easier if internal math is isolated

Current code point:

- [src/response/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/response/mod.rs)
- [src/response/backend.rs](/c:/Users/eynulai/Downloads/mfs/src/response/backend.rs)

### 5. `synthesis` is thin in a good way, but it depends on unstable seams

`ChebyshevSynthesis` is already acting like the right orchestration layer.
The problem is that it currently wires together pieces whose boundaries are not
settled yet.

Current code point:

- [src/synthesis.rs](/c:/Users/eynulai/Downloads/mfs/src/synthesis.rs)

## Recommended Phases

## Phase 1: Stabilize Domain Semantics

Goal:
clarify the meaning of current input types before deeper algorithms land

Tasks:

1. Separate filter shape from approximation family in the spec model.
2. Make transmission-zero semantics more explicit than `(value, domain)`.
3. Decide whether normalized zeros are always low-pass prototype values.
4. Move zero normalization policy out of `TransmissionZero` methods and into a
   dedicated helper or domain service.

Suggested target shape:

```rust
pub enum FilterClass {
    LowPass,
    HighPass,
    BandPass,
    BandStop,
}

pub enum ApproximationFamily {
    Chebyshev,
}

pub enum TransmissionZero {
    Normalized(f64),
    PhysicalHz(f64),
}
```

Exit criteria:

- input-domain types are semantically unambiguous
- no normalization algorithm lives inside plain spec value objects

## Phase 2: Strengthen Frequency Mapping Boundaries

Goal:
make `freq` the authoritative place for domain mapping rules

Tasks:

1. Keep `FrequencyMapping` as the only place that knows how physical frequencies
   map to normalized coordinates.
2. Add explicit helpers for zero normalization rather than calling back from
   `spec`.
3. Decide whether reverse mapping should always return a single value or expose
   branch-aware behavior for band-pass and band-stop transforms.
4. Add tests for invalid physical-domain edge cases.

Important design note:

For band-pass mappings, inverse transforms can be branch-sensitive. That should
be reflected in API naming before real solver work depends on it.

Exit criteria:

- all physical-to-normalized mapping lives in `freq`
- transmission-zero normalization depends on `freq`, not on `spec`

## Phase 3: Promote Polynomial Infrastructure

Goal:
prepare `approx` for real generalized Chebyshev logic

Tasks:

1. Extract polynomial helpers into a dedicated submodule.
2. Introduce validation helpers for degree and coefficient layout.
3. Clarify coefficient ordering conventions in docs and tests.
4. Keep `PolynomialSet` as an artifact, but stop treating it as an unvalidated
   bag of vectors.

Suggested split:

- `approx::chebyshev`
- `approx::polynomial`
- `approx::normalization`

Exit criteria:

- polynomial conventions are explicit
- Chebyshev implementation can grow without bloating one file

## Phase 4: Split Matrix Artifact From Matrix Synthesis

Goal:
make coupling-matrix generation pluggable and easier to test

Tasks:

1. Remove synthesis logic from `CouplingMatrix`.
2. Introduce `CouplingMatrixSynthesizer` or equivalent.
3. Keep `CouplingMatrixBuilder` focused on validated assembly only.
4. Reserve `matrix` for artifact storage, builders, and transforms.

Suggested direction:

```rust
pub struct CouplingMatrixSynthesizer;

impl CouplingMatrixSynthesizer {
    pub fn synthesize(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        // ...
    }
}
```

Exit criteria:

- `CouplingMatrix` no longer owns approximation-to-matrix logic
- synthesis strategies can vary without changing the matrix artifact type

## Phase 5: Prepare For Topology Work

Goal:
make room for equivalent-network transforms without overloading core storage

Tasks:

1. Introduce a topology-focused module once a second transform exists.
2. Represent transform operations separately from the base matrix type.
3. Keep topology equivalence testing at the response or invariant level.

Exit criteria:

- topology logic has a dedicated extension point
- matrix storage stays small and understandable

## Phase 6: Create A Real Response-Solver Boundary

Goal:
allow numerical backend upgrades without API churn

Tasks:

1. Keep `ResponseSample` and `SParameterResponse` as public output artifacts.
2. Move backend-specific numeric implementation behind private helpers or
   internal modules.
3. Decide when complex arithmetic enters the system.
4. Add regression tests once response values become physically meaningful.

Exit criteria:

- public response API is backend-neutral
- solver implementation can change without touching callers

## Phase 7: Reshape The Source Tree

Goal:
convert the current flat file layout into ownership-based directories only when
the module contents justify it

Recommended order:

1. split `approx.rs`
2. split `matrix`
3. split `response.rs`
4. split `spec.rs` and `freq.rs` into a future `domain/` area if the semantic
   model becomes richer

Do not restructure too early just for aesthetics. Split when it removes real
cognitive load.

## Concrete Backlog For The Current Codebase

These are the highest-value next tasks, in order.

### Backlog A

Move transmission-zero normalization out of
[src/spec.rs](/c:/Users/eynulai/Downloads/mfs/src/spec.rs) and into
[src/freq.rs](/c:/Users/eynulai/Downloads/mfs/src/freq.rs) or a dedicated
normalization helper.

Status:
completed

Why first:

- it improves domain boundaries immediately
- it is low risk
- it will help every later approximation algorithm

### Backlog B

Extract matrix synthesis from
[src/matrix/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/matrix/mod.rs) and update
[src/synthesis.rs](/c:/Users/eynulai/Downloads/mfs/src/synthesis.rs) to depend
on a synthesizer instead of `CouplingMatrix::from_polynomials(...)`.

Status:
completed

Why second:

- it breaks the biggest current responsibility leak
- it prepares the crate for multiple synthesis strategies

### Backlog C

Split polynomial helpers out of
[src/approx/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/mod.rs) before adding
real generalized Chebyshev math.

Status:
completed

Why third:

- current placeholder logic is still small
- moving now avoids a harder refactor later

### Backlog D

Refactor
[src/response/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/response/mod.rs) so the
public sample types remain stable while solver internals become private helper
functions or a backend module.

Status:
completed

### Backlog E

Strengthen the domain model so filter class, approximation family, and
performance targets are represented independently.

Suggested code points:

- [src/spec.rs](/c:/Users/eynulai/Downloads/mfs/src/spec.rs)
- [src/lib.rs](/c:/Users/eynulai/Downloads/mfs/src/lib.rs)

Why next:

- this is the biggest remaining semantic ambiguity
- future approximation families will otherwise force avoidable API churn

Status:
completed

### Backlog F

Add stronger polynomial validation and normalization rules around
[src/approx/polynomial.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/polynomial.rs).

Status:
started

Why next:

- it prepares the crate for real generalized Chebyshev equations
- it keeps future numerical bugs from hiding behind raw `Vec<f64>` fields

### Backlog G

Integrate
[src/approx/generalized_chebyshev.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/generalized_chebyshev.rs)
more deeply into
[src/approx/chebyshev.rs](/c:/Users/eynulai/Downloads/mfs/src/approx/chebyshev.rs)
so placeholder `e/f` coefficient generation can be replaced by real synthesis
artifacts.

Why next:

- the generalized helper chain now exists
- the remaining gap is in connecting helper outputs to the main approximation
  artifact model

### Backlog H

Port coupling-matrix recovery and optimization logic from the Python core,
including the stages around circuit-method matrix recovery and Amari-style
gradient work when the Rust artifact model is ready.

Why later:

- these depend on a more complete prototype artifact representation
- forcing them in too early would create churn in the matrix and approximation
  layers

## Suggested Deliverables Per Phase

For each phase, aim to produce:

1. one focused refactor PR or commit
2. one corresponding doc update
3. direct unit tests for the new boundary
4. no behavior regression in the high-level flow

## Recommended Immediate Next Refactor

If starting today, the single best next step is:

strengthen the domain model in `spec` so filter class, approximation family,
and performance intent are separate concepts has already been completed.

The highest-leverage remaining refactor is now integrating generalized
Chebyshev helper outputs into the primary approximation pipeline.
