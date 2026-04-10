# Public API Sketch for a Filter Synthesis Library

## Purpose

This document sketches a public Rust API for an open-source microwave filter
synthesis library.

It complements the architecture and crate-layout documents by answering a
practical design question:

- what should users of the library actually write in code

The goal is not to freeze every function name up front. The goal is to define
the shape of an API that is:

- coherent
- discoverable
- hard to misuse
- flexible enough for research and engineering workflows

## API Design Goals

The public API should:

- model the filter-design workflow directly
- expose domain objects instead of raw matrices whenever possible
- make advanced transforms explicit about their preconditions
- preserve access to low-level operations for expert users
- return structured reports rather than opaque success flags
- support both script-like and library-style usage

The core design rule is:

- common workflows should read like a design pipeline, not like manual matrix
  bookkeeping

## Primary User Personas

The API should serve at least three kinds of users:

### 1. Application User

Wants to:

- define a filter specification
- synthesize a matrix
- convert it to a realizable topology
- analyze or export the result

Needs:

- simple builders
- stable defaults
- minimal exposure to internal algorithm details

### 2. Research User

Wants to:

- inspect intermediate results
- force a specific canonical form
- perform custom transform sequences
- compare alternate synthesis procedures

Needs:

- lower-level access to matrices and rotations
- deterministic reporting of intermediate operations

### 3. Contributor

Wants to:

- add a new approximation method
- add a new topology transform
- improve verification or export tooling

Needs:

- stable extension points
- clear type boundaries
- explicit preconditions and invariants

## API Surface Overview

At the highest level, the public API should revolve around a few central types:

- `FilterSpec`
- `Approximation`
- `CouplingMatrix`
- `TopologyKind`
- `NetworkResponse`
- `TransformReport`
- `VerificationReport`

These types should be easy to discover from the crate root and from a small
prelude module.

## Recommended Top-Level Exports

The crate root should expose the main user-facing types:

```rust
pub use crate::approx::Approximation;
pub use crate::analysis::NetworkResponse;
pub use crate::matrix::CouplingMatrix;
pub use crate::spec::FilterSpec;
pub use crate::topology::TopologyKind;
pub use crate::transform::TransformReport;
pub use crate::verify::VerificationReport;
```

For most users, these should be enough to get started.

## Builder-Oriented Specification API

The specification API should be builder-based so users can express intent
clearly and validation can happen once at `build()`.

Example:

```rust
use filter_synthesis::prelude::*;

let spec = FilterSpec::builder()
    .order(6)
    .center_frequency_ghz(2.45)
    .fractional_bandwidth(0.08)
    .return_loss_db(22.0)
    .transmission_zeros(vec![
        1.8,
        3.1,
    ])
    .symmetric(false)
    .build()?;
```

Suggested builder methods:

- `order(...)`
- `passband(...)`
- `center_frequency_ghz(...)`
- `fractional_bandwidth(...)`
- `return_loss_db(...)`
- `ripple_db(...)`
- `transmission_zeros(...)`
- `symmetric(...)`
- `prototype_kind(...)`
- `build()`

Design rule:

- mutually exclusive options should be validated in the builder rather than
  deferred to synthesis-time failure

## Approximation API

Approximation should feel like a derived object, not a manually assembled bag
of polynomials.

Example:

```rust
let approx = Approximation::from_spec(&spec)?;
```

If multiple methods are supported, use an explicit selection mechanism:

```rust
let approx = Approximation::builder()
    .method(ApproximationMethod::GeneralizedChebyshev)
    .from_spec(&spec)?;
```

Suggested capabilities:

- inspect characteristic polynomials
- inspect zeros and poles
- export a normalized summary

Example:

```rust
println!("E(s) = {}", approx.e_polynomial());
println!("Transmission zeros = {:?}", approx.transmission_zeros());
```

## Canonical Matrix Synthesis API

Canonical synthesis should be explicit about the requested form.

Example:

```rust
let matrix = CouplingMatrix::synthesize(&approx, CanonicalForm::Transversal)?;
```

Alternative style:

```rust
let matrix = SynthesisEngine::new()
    .canonical_form(CanonicalForm::Transversal)
    .build_matrix(&approx)?;
```

For first release ergonomics, the associated-function style is probably better.

Suggested capabilities on `CouplingMatrix`:

- query matrix order
- inspect entries
- inspect topology label
- inspect source/load indices
- export as dense table
- compute a graph view

Example:

```rust
assert_eq!(matrix.order(), 6);
println!("Topology: {:?}", matrix.topology());
println!("M[0, 1] = {}", matrix.entry(0, 1)?);
```

## High-Level Transform API

High-level topology transforms should be discoverable and safe.

Example:

```rust
let folded = matrix.transform_to(TopologyKind::Folded)?;
let arrow = matrix.transform_to(TopologyKind::Arrow)?;
```

This API should:

- enforce obvious preconditions
- attach transform reports
- preserve access to the original matrix unless the user chooses in-place
  mutation

Suggested result style:

```rust
let outcome = matrix.transform(TopologyKind::Arrow)?;

let arrow = outcome.matrix();
let report = outcome.report();
```

This is often better than returning only the transformed matrix because the
report is operationally important.

The current crate is converging on helper-oriented variants of this shape:

```rust
let outcome = transform_matrix_with_response_check(
    &matrix,
    TopologyKind::Arrow,
    &grid,
    ResponseTolerance::default(),
)?;

assert!(outcome.report.passes());
if let Some(comparison) = &outcome.report.response.comparison {
    println!("max |S21| deviation = {}", comparison.max_s21_magnitude_deviation);
}
```

The current crate also exposes detail-preserving orchestration helpers so
callers can see whether a run stayed on the classical path, attached
generalized helper data, or fell back to a placeholder matrix-construction
path:

```rust
let outcome = synthesize_and_evaluate_chebyshev_with_details(
    &spec,
    &mapping,
    &grid,
)?;

println!("{:?}", outcome.approximation_kind);
println!("{:?}", outcome.matrix_method);
```

Internally, the current crate has also started to separate the approximation
layer into smaller responsibilities:

- `approx::complex_poly` for reusable complex-polynomial storage and root solving
- `approx::generalized_ops` for generalized-path `w <-> s` transforms and
  recurrence-side helpers
- `approx::generalized_chebyshev` for Cameron/generalized pipeline stages

That split is still mostly an implementation boundary rather than a polished
top-level public surface, but it is already useful guidance for contributors.

## Advanced Section Extraction API

Advanced transforms should read as specialized operations, not generic topology
 requests.

Example:

```rust
let arrow = matrix.transform_to(TopologyKind::Arrow)?;

let trisections = arrow.extract_trisection_cascade()?;
let quartets = trisections.merge_adjacent_trisections()?;
let boxed = trisections.to_box_sections()?;
```

Important rule:

- APIs for advanced sections should encode the expected precursor form

That means `extract_trisection_cascade()` should fail clearly if called on a
non-arrow matrix.

Possible advanced APIs:

- `extract_trisection_cascade()`
- `shift_trisection_left(index)`
- `merge_adjacent_trisections()`
- `to_box_sections()`
- `to_extended_box_sections()`

These should return structured outcomes and not silently mutate the matrix
unless the user explicitly requests it.

The current crate is converging on dedicated section-transform helpers:

```rust
let triplet = extract_triplet_section_with_response_check(
    &matrix,
    -1.3,
    2,
    &grid,
    ResponseTolerance::default(),
)?;

assert!(triplet.passes());
```

For callers that want to stay at the synthesis layer rather than dropping to
transform-only helpers, the current crate is also converging on:

```rust
let sections = SectionSynthesis::default();
let verified = sections.synthesize_trisection_with_response_check(
    &polynomials,
    -1.25,
    (2, 4),
    &grid,
    ResponseTolerance::default(),
)?;

assert!(verified.passes());
```

## Analysis API

Users should be able to analyze any valid matrix with minimal ceremony.

Example:

```rust
let response = matrix.analyze(
    FrequencySweep::new()
        .start_ghz(1.0)
        .stop_ghz(4.0)
        .points(2001)
)?;

println!("Peak insertion loss: {}", response.max_insertion_loss_db());
println!("Minimum return loss: {}", response.min_return_loss_db());
```

Suggested analysis entry points:

- `analyze(sweep)`
- `s11(sweep)`
- `s21(sweep)`
- `group_delay(sweep)`
- `transmission_zero_locations()`

Suggested response methods:

- `samples()`
- `s11_db()`
- `s21_db()`
- `group_delay_ns()`
- `min_return_loss_db()`
- `max_insertion_loss_db()`

## Verification API

Verification should be explicit and composable.

Example:

```rust
let report = folded.verify_against(&matrix)
    .check_response_invariance()
    .check_topology_constraints()
    .run()?;

assert!(report.passed());
```

Alternative ergonomic style:

```rust
let report = Verifier::new()
    .reference(&matrix)
    .candidate(&folded)
    .check_response_invariance()
    .check_minimum_path_rule()
    .run()?;
```

Suggested checks:

- response invariance
- topology shape
- transform preconditions
- transmission-zero capacity
- symmetry and reciprocity

The output should not just be `bool`.

Suggested `VerificationReport` contents:

- overall pass/fail
- check-by-check results
- maximum deviation metrics
- tolerance used
- warnings

## Export API

Export should be convenient but not dominate the core interface.

Example:

```rust
let json = matrix.to_json()?;
let markdown = matrix.to_markdown_report()?;
let graph = matrix.to_graph()?;
```

If export variants grow, a formatter-based API may scale better:

```rust
let report = MarkdownReport::new()
    .include_spec(&spec)
    .include_approximation(&approx)
    .include_matrix(&matrix)
    .include_response(&response)
    .render()?;
```

This is especially useful for CLI tooling and batch workflows.

## Low-Level Expert API

The library should also support expert workflows that operate directly on the
matrix model.

Example:

```rust
let rotated = matrix.rotate((2, 4), theta)?;
let annihilated = matrix.annihilate((0, 4), using_pivots(2, 3))?;
let shortest_path = matrix.shortest_source_load_path_len()?;
```

Suggested low-level methods:

- `entry(i, j)`
- `set_entry(i, j, value)` if mutable forms are allowed
- `rotate((i, j), theta)`
- `annihilate(target, pivots)`
- `is_symmetric(tol)`
- `pattern(tol)`
- `shortest_source_load_path_len()`

These operations are valuable, but they should be clearly documented as expert
tools rather than the recommended starting point.

## In-Place vs Immutable API

The public API should prefer immutable transforms first.

Good default:

```rust
let arrow = matrix.transform_to(TopologyKind::Arrow)?;
```

Optional advanced variant:

```rust
matrix.transform_to_in_place(TopologyKind::Arrow)?;
```

Why:

- immutable APIs are easier to reason about
- transform comparisons become easier
- contributors are less likely to introduce accidental state corruption

In-place methods can still exist for performance-sensitive workflows.

## Error Design

Errors should be precise and domain-specific.

Examples:

- `InvalidSpecification`
- `UnsupportedApproximation`
- `SynthesisFailed`
- `InvalidTopologyPrecondition`
- `SingularTransform`
- `VerificationFailed`

Error messages should help the user recover.

Example:

```text
cannot extract trisection cascade from topology Folded; expected Arrow
```

This is much better than a generic "invalid state" error.

## Report-Oriented Operations

Many high-value operations should return reports.

Examples:

- synthesis can return a `SynthesisReport`
- transform can return a `TransformReport`
- verification can return a `VerificationReport`

Why this matters:

- users often need to inspect what happened, not just the final object
- reports improve debuggability
- reports make CLI and export features easier to build later

Suggested `TransformReport` fields:

- input topology
- output topology
- rotation sequence
- annihilated couplings
- created couplings
- maximum response deviation
- warnings about numerical conditioning

The current implementation already has a smaller but real version of this:

- input and output topology metadata
- pattern verification result
- a shared `ResponseCheckReport`
- implementation notes

The same reporting style is now starting to show up at the section-synthesis
layer too, via `VerifiedSectionSynthesis`.

## Prelude Sketch

A small prelude can make simple scripts pleasant.

Example:

```rust
use filter_synthesis::prelude::*;
```

Suggested prelude exports:

- `FilterSpec`
- `Approximation`
- `ApproximationMethod`
- `CouplingMatrix`
- `CanonicalForm`
- `TopologyKind`
- `FrequencySweep`
- crate `Result`

Keep it intentionally small.

## Example Workflow APIs

### Example 1: Basic End-to-End Use

```rust
use filter_synthesis::prelude::*;

let spec = FilterSpec::builder()
    .order(4)
    .return_loss_db(20.0)
    .build()?;

let approx = Approximation::from_spec(&spec)?;
let matrix = CouplingMatrix::synthesize(&approx, CanonicalForm::Transversal)?;
let folded = matrix.transform_to(TopologyKind::Folded)?;
let response = folded.analyze(FrequencySweep::normalized(-3.0, 3.0, 1601))?;
```

### Example 1b: Explicit Generalized Chebyshev Request

```rust
use filter_synthesis::prelude::*;

let spec = FilterSpec::generalized_chebyshev(4, 20.0)?
    .with_transmission_zeros(vec![
        TransmissionZero::normalized(-1.5),
        TransmissionZero::normalized(1.8),
    ]);

let mapping = LowPassMapping::new(1.0)?;
let outcome = synthesize_chebyshev_with_details(&spec, &mapping)?;

assert_eq!(
    outcome.approximation_kind,
    ApproximationStageKind::GeneralizedChebyshev,
);
println!("{:?}", outcome.matrix_method);
```

### Example 2: Transform with Audit Trail

```rust
let outcome = matrix.transform(TopologyKind::Arrow)?;

println!("Created couplings: {:?}", outcome.report().created_entries());
println!("Removed couplings: {:?}", outcome.report().annihilated_entries());
```

### Example 3: Advanced Expert Workflow

```rust
let theta = matrix.solve_annihilation_angle(
    target_entry(0, 4),
    pivot_pair(2, 3),
)?;

let candidate = matrix.rotate((2, 3), theta)?;

let verification = Verifier::new()
    .reference(&matrix)
    .candidate(&candidate)
    .check_response_invariance()
    .run()?;
```

## Trait and Extension Strategy

Traits should be used sparingly and only for stable extension points.

Good candidates:

- `ApproximationMethod`
- `TopologyTransform`
- `Analyzer`
- `Exporter`
- `Verifier`

Each trait should reflect a concept users may reasonably want to extend.

Avoid introducing traits that exist only to generalize one implementation.

## API Stability Guidelines

For open-source adoption, stability matters as much as algorithm quality.

Recommended rules:

- keep the number of top-level public types small
- prefer adding methods over reshaping core objects
- keep low-level modules public only when they are intended for external use
- document preconditions as part of the public contract
- treat report object shapes as public API once external tools depend on them

## Suggested MVP Public API

For an initial release, the public API can stay compact.

Recommended MVP surface:

- `FilterSpec::builder()`
- `Approximation::from_spec(...)`
- `CouplingMatrix::synthesize(...)`
- `CouplingMatrix::transform_to(...)`
- `CouplingMatrix::analyze(...)`
- `Verifier`
- `TopologyKind`
- `CanonicalForm`
- `FrequencySweep`

This is enough to support a complete workflow without overcommitting too early.

## Summary

The best public API for this library is one that makes the design process
legible in code:

```text
specification -> approximation -> canonical matrix -> transform -> analyze
```

The most important API design decisions are:

- center the API around domain objects
- prefer builder-based specification input
- keep canonical synthesis explicit
- make advanced transforms precondition-aware
- return structured reports for synthesis, transform, and verification
- support low-level matrix control without forcing all users into it

If the API follows these principles, the library can serve both practical
engineering workflows and deeper research-oriented experimentation.
