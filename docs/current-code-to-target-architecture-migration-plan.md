# Current Code to Target Architecture Migration Plan

## Purpose

This document records how the codebase now maps onto the target
filter-synthesis library architecture.

Earlier versions of this note were written as a migration proposal. The crate
has since moved far enough that this document is now better used as:

- a snapshot of what has already been reshaped
- a guide to which boundaries are now considered real
- a short list of the next refactors worth doing

## Current Top-Level Structure

The crate now has the following top-level layout:

- `src/spec/`
- `src/approx/`
- `src/matrix/`
- `src/response/`
- `src/synthesis/`
- `src/transform/`
- `src/verify/`
- `src/freq.rs`
- `src/error.rs`
- `src/prelude.rs`

At a high level, the crate now expresses this workflow:

```text
FilterSpec -> approximation -> canonical matrix -> topology transform
           -> response / verification
```

That is much closer to the intended library architecture than the original
flat-file layout.

## What Has Already Landed

The following architecture changes are already implemented:

- `FilterParameter` has been replaced by `FilterSpec`
- `src/spec.rs` has been split into `src/spec/`
- `FilterSpec` now supports explicit `Chebyshev` and
  `GeneralizedChebyshev` family selection
- `src/synthesis.rs` has been split into `src/synthesis/`
- canonical synthesis now belongs to `CanonicalMatrixSynthesis`
- section-oriented synthesis now belongs to `SectionSynthesis`
- the low-level synthesis engine now lives in `MatrixSynthesisEngine`
- topology conversion is exposed through `transform`
- section extraction is now exposed through `transform` facades and reused by `synthesis`
- response comparison and topology-shape checks now live in `verify`
- section-verification helpers now live in `verify`
- `CouplingMatrix` now carries explicit topology metadata
- advanced trisection extraction now enforces Arrow preconditions
- transform outcomes now include a minimal `TransformReport`
- transform and section-transform flows can attach optional response comparisons
- orchestration outcomes can now report both approximation-stage kind and
  matrix-synthesis method
- `prelude` now exposes a coherent high-level workflow surface
- synthesis-specific tests have been moved out of `matrix`
- topology-shape tests have been moved out of `matrix`

In short, the project is no longer only planning the target architecture. It is
already operating inside most of it.

## Current Module Assessment

### `src/spec/`

Current role:

- defines `FilterSpec`
- defines `FilterClass`
- defines `ApproximationFamily`
- defines return-loss and transmission-zero input types
- provides `FilterSpecBuilder`

Assessment:

- this is now a real first-class `spec` subsystem
- the public language is already aligned with the design docs
- the remaining work is refinement, not relocation

Recommended next direction:

- keep `FilterSpec` as the dominant user-facing term
- continue tightening validation rules as more synthesis paths land
- add richer spec fixtures as approximation support grows

### `src/approx/`

Current role:

- approximation traits and engines
- Chebyshev implementation
- reusable complex-polynomial primitives and root-solving strategy objects
- generalized Chebyshev helpers
- shared generalized-path `w <-> s` transforms and recurrence helpers
- explicit strict generalized-Chebyshev orchestration path
- polynomial artifacts and utilities

Assessment:

- still one of the strongest matches to the target architecture
- already serves as a stable approximation layer
- now distinguishes between default classical flow and explicit
  generalized-Chebyshev flow
- future work is mostly about deeper approximation coverage and richer
  generalized fidelity

Recommended next direction:

- keep the current directory structure
- continue separating reusable math from orchestration where helpful
- add more literature-backed approximation fixtures
- keep making the generalized path more mathematically complete rather than
  adding new facade layers first

### `src/matrix/`

Current role:

- defines `CouplingMatrix` and related matrix artifacts
- provides builders and low-level matrix operations
- stores topology metadata used by higher-level precondition checks
- retains internal extraction primitives used by transform-layer facades

Assessment:

- this module has been substantially reduced in scope
- it no longer owns the synthesis subsystem
- advanced section extraction is no longer part of the intended public matrix surface
- it is now much closer to a true matrix-domain layer

Recommended next direction:

- keep `matrix` focused on data structures and low-level operations
- avoid moving new workflow orchestration back into `matrix`
- continue trimming any remaining topology-specific convenience logic that
  belongs more naturally in `transform`

### `src/response/`

Current role:

- response solver
- response samples and response settings
- solver backend

Assessment:

- still a strong analysis-oriented layer
- the current implementation is already reusable across synthesis and
  verification tests

Recommended next direction:

- keep functionality here for now
- consider a later public rename to `analysis` only if it improves API clarity

### `src/synthesis/`

Current role:

- high-level orchestration through `ChebyshevSynthesis`
- canonical synthesis through `CanonicalMatrixSynthesis`
- low-level synthesis execution through `MatrixSynthesisEngine`
- residue/admittance helpers
- placeholder construction
- section-oriented synthesis through `SectionSynthesis`
- stage reporting through `ApproximationStageKind` and
  `MatrixSynthesisMethod`

Assessment:

- this is now a real top-level synthesis subsystem
- the earlier problem of matrix-owned synthesis logic has been largely fixed
- the code structure now matches the intended architecture well

Recommended next direction:

- keep new synthesis work inside `synthesis`
- continue splitting internals only when it improves maintainability
- grow fixtures and tests alongside new algorithms

### `src/transform/`

Current role:

- topology conversion facade
- folded and arrow conversion entry points
- wheel conversion entry point
- section-extraction facades
- transform engine and transform outcome/report types

Assessment:

- this subsystem now exists as a first-class layer
- topology conversion has already been pulled out of `matrix` at the API level
- section-oriented transform execution is now starting to move here as well

Recommended next direction:

- keep new topology workflows here
- expand transform reporting as more advanced topologies are added
- move additional topology semantics here rather than reintroducing them into
  matrix-oriented types

### `src/verify/`

Current role:

- response comparison helpers
- response tolerance types
- topology-shape checks for folded and arrow forms
- section-verification reports and helpers

Assessment:

- this subsystem is now real rather than aspirational
- it already supports the most important invariance and structure checks

Recommended next direction:

- add richer verification helpers for advanced section extraction
- keep verification logic centralized instead of scattering assertions through
  tests

### `src/freq.rs`

Current role:

- frequency grids
- physical-to-normalized mapping
- low-pass and band-pass mapping helpers

Assessment:

- stable support module
- already fits the target architecture well enough

Recommended next direction:

- keep as-is unless a future workspace split makes a rename worthwhile

### `src/error.rs`

Current role:

- library error types

Assessment:

- stable shared foundation

Recommended next direction:

- keep centralized
- add more specific errors only where they improve user-facing diagnostics

### `src/prelude.rs`

Current role:

- ergonomic re-exports for the high-level library workflow

Assessment:

- important usability layer
- makes the new architecture visible to users without exposing every module
  detail

Recommended next direction:

- keep the prelude focused on the intended public workflow
- avoid re-exporting low-level legacy concepts just for convenience

## Mapping to the Target Architecture

The target architecture from the design documents includes these major areas:

- `spec`
- `approx`
- `matrix`
- `synthesis`
- `transform`
- `analysis`
- `verify`
- `export`

The current crate maps to that target as follows:

### Target `spec`

Current source:

- `src/spec/mod.rs`
- `src/spec/types.rs`
- `src/spec/builder.rs`

Status:

- landed

### Target `approx`

Current source:

- `src/approx/mod.rs`
- `src/approx/chebyshev.rs`
- `src/approx/complex_poly.rs`
- `src/approx/generalized_chebyshev.rs`
- `src/approx/generalized_ops.rs`
- `src/approx/polynomial.rs`

Status:

- landed

### Target `matrix`

Current source:

- `src/matrix/mod.rs`

Status:

- landed in a narrower and healthier form

### Target `synthesis`

Current source:

- `src/synthesis/mod.rs`
- `src/synthesis/engine.rs`
- `src/synthesis/canonical.rs`
- `src/synthesis/orchestration.rs`
- `src/synthesis/residues.rs`
- `src/synthesis/placeholder.rs`
- `src/synthesis/sections.rs`

Status:

- landed

### Target `transform`

Current source:

- `src/transform/mod.rs`
- `src/transform/folded.rs`
- `src/transform/arrow.rs`

Status:

- landed

### Target `analysis`

Current source:

- `src/response/mod.rs`
- `src/response/backend.rs`

Status:

- partially landed under the current `response` name

### Target `verify`

Current source:

- `src/verify/mod.rs`

Status:

- landed

### Target `export`

Current source:

- no dedicated export subsystem yet

Status:

- not started as a first-class module

## What Still Looks Incomplete

The highest-value unfinished areas are now narrower than before:

- richer generalized Chebyshev approximation fidelity inside the main
  user-facing workflow
- broader topology coverage beyond folded and arrow
- stronger verification helpers for advanced extracted sections
- a dedicated export/reporting layer
- more literature-backed regression fixtures and benchmarks

These are now feature gaps, not foundational architecture gaps.

## Recommended Next Steps

If implementation work continues from here, the most useful next moves are:

1. deepen approximation coverage rather than reshuffling modules again
2. add advanced topology workflows to `transform`
3. add advanced verification helpers alongside each new topology
4. introduce export/reporting only when the current artifacts stabilize
5. keep tests aligned with module ownership as new code lands

## Summary

The original migration problem was that `matrix` held too many unrelated
responsibilities and the public API language was still shaped by older
experimental code.

That is no longer the main story.

The current story is:

- `FilterSpec` is the real spec model
- `synthesis` is a real subsystem
- `transform` is a real subsystem
- `verify` is a real subsystem
- `matrix` has been reduced toward its proper scope

The codebase has not finished the library, but it has already crossed the line
from "planned architecture" into "implemented architecture."
