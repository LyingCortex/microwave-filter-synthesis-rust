# Testing Strategy for a Filter Synthesis Library

## Purpose

This document proposes a testing strategy for an open-source microwave filter
synthesis library.

The central challenge in this kind of software is that correctness has several
layers:

- mathematical correctness
- numerical robustness
- topology-structure correctness
- response invariance under transformation
- API stability for downstream users

Because of that, the test strategy must be broader than ordinary unit testing.

## Testing Goals

The test suite should ensure that the library:

- computes valid canonical matrices from supported specifications
- preserves electrical response during topology reconfiguration
- rejects invalid inputs and invalid algorithm preconditions clearly
- remains numerically stable within defined tolerances
- keeps output formats and public API behavior consistent across releases

The most important rule is:

- no topology algorithm should be trusted unless both its structure and its
  electrical behavior are tested

## Core Testing Principles

The library should follow a few guiding rules:

- test domain invariants close to the code that enforces them
- test workflow guarantees at integration level
- make tolerances explicit and shared
- prefer literature-backed fixtures where available
- use regression tests for every bug that reaches a user-visible surface

The test suite should be treated as part of the mathematical specification of
the project, not just as a safety net for refactoring.

## Test Categories

The full strategy should include these categories:

- unit tests
- integration tests
- fixture-based reference tests
- invariance tests
- property-oriented tests
- regression tests
- export and CLI smoke tests
- performance benchmarks

Each category serves a different purpose and should not be collapsed into a
single style.

## 1. Unit Tests

### Purpose

Validate small mathematical and structural behaviors locally.

### Typical Targets

- polynomial evaluation
- tolerance-aware comparisons
- matrix indexing rules
- symmetry checks
- rotation matrix construction
- annihilation-angle solver helpers
- error construction and formatting

### Placement

Place these near the implementation module using Rust module-local tests.

### Example Questions

- does source index always map to `0`
- does load index always map to `N + 1`
- does a two-pivot rotation remain orthogonal within tolerance
- does threshold logic classify near-zero entries consistently

These tests should be fast and numerous.

## 2. Integration Tests

### Purpose

Validate user-visible workflows across module boundaries.

### Typical Targets

- specification to approximation to canonical synthesis
- canonical synthesis to response analysis
- canonical synthesis to folded conversion
- canonical synthesis to arrow conversion
- arrow conversion to trisection extraction

### Placement

Use the repository-level `tests/` directory.

### Example Questions

- can a valid specification produce a canonical transversal matrix
- does folded conversion preserve the response of the precursor matrix
- does trisection extraction reject non-arrow input

These tests should model how the library is actually used.

## 3. Fixture-Based Reference Tests

### Purpose

Compare outputs against known or literature-derived reference cases.

### Typical Targets

- low-order canonical matrices
- folded and arrow sparsity patterns
- known transmission-zero assignments
- expected report structures for standard transforms

### Fixture Sources

Prefer:

- analytically simple hand-checked cases
- published literature examples
- previously verified internal golden cases

### Storage Strategy

Keep fixtures centralized and documented.

Suggested locations:

- `src/fixtures/` for small reusable code-based fixtures
- `tests/fixtures/` for JSON, CSV, or golden-output assets

Fixture names should explain what they represent, not just their order number.

## 4. Response Invariance Tests

### Purpose

Verify the core promise of similarity transforms:

- topology changes must not change the ideal response

### Typical Targets

- transversal to folded
- transversal to arrow
- arrow to trisection conditioning
- trisection shifting
- trisection merging into quartets

### Test Pattern

For each transform:

1. create or load a valid precursor matrix
2. analyze its response over a fixed sweep
3. perform the transform
4. analyze the transformed response over the same sweep
5. compare deviation against explicit tolerances

### Metrics

Recommended comparison metrics:

- max absolute deviation of `|S11|`
- max absolute deviation of `|S21|`
- deviation in transmission-zero locations if extracted
- optional group-delay deviation in passband

### Important Rule

- structural success is not enough if response invariance fails

This category is one of the most important in the whole project.

## 5. Topology Shape Tests

### Purpose

Verify that a matrix claiming a certain topology actually has the intended
non-zero pattern.

### Typical Targets

- canonical transversal sparsity
- folded sparsity
- arrow sparsity
- section-local patterns for trisections or box sections

### Test Pattern

Use tolerance-aware pattern extraction rather than exact zero checks.

Example questions:

- are all non-nearest-neighbor source couplings removed in the folded form
- are cross-couplings concentrated in the expected row and column in the arrow
  form
- does the extracted trisection section have the expected local pattern

These tests are necessary, but never sufficient on their own.

## 6. Precondition and Validation Tests

### Purpose

Ensure the library fails safely when users call algorithms in the wrong state.

### Typical Targets

- invalid filter specifications
- unsupported approximation settings
- trisection extraction on non-arrow input
- quartet merging on a matrix without adjacent trisections
- invalid pivot selections for transforms

### Expected Behavior

Errors should be:

- deterministic
- domain-specific
- easy to interpret

These tests are critical for API trustworthiness.

## 7. Property-Oriented Tests

### Purpose

Check general invariants across many generated cases, not just a few fixtures.

### Suitable Properties

- matrix symmetry is preserved under supported transforms
- similarity transforms preserve eigenvalues within tolerance
- rotating by zero angle leaves the matrix unchanged
- applying a rotation and its inverse returns the original matrix
- shortest-path calculations are invariant under value-preserving graph relabels
  where applicable

### Tooling

If property-based testing is added, use it carefully and keep the generated
inputs meaningful.

Bad generated inputs can waste time on invalid states that teach little.

### Recommendation

Use property-oriented tests for foundational modules first:

- `math`
- `matrix`
- `transform`

This gives high value with manageable complexity.

## 8. Regression Tests

### Purpose

Prevent known failures from returning.

### When to Add One

Add a regression test whenever:

- a numerical bug is fixed
- a topology transform silently changed response
- a fixture case exposed incorrect indexing
- an export format changed unexpectedly
- a user-reported edge case is resolved

### Naming

Regression test names should describe the bug, not the implementation detail.

Example:

- `arrow_transform_preserves_response_for_asymmetric_case`
- `trisection_extraction_rejects_folded_input`

Regression tests should live close to the workflow they protect.

## 9. Export and CLI Smoke Tests

### Purpose

Ensure the user-facing tooling remains usable.

### Typical Targets

- JSON export validity
- Markdown report generation
- CSV formatting
- CLI command success on standard examples

### Test Style

These do not need deep mathematical coverage in every case. They should
confirm:

- outputs are produced
- formats are parseable
- core fields exist
- command wiring is correct

Mathematical correctness should still be covered primarily in core tests.

## 10. Benchmarking and Performance Regression

### Purpose

Detect major performance regressions once correctness is stable.

### Typical Targets

- polynomial evaluation
- canonical synthesis
- repeated matrix rotations
- response analysis at large sweep counts
- batch verification over many cases

### Guidance

Benchmarking should not replace correctness testing.

Add it after the library has a stable enough API and representative workloads.

## Tolerance Strategy

Tolerance handling should be explicit and centralized.

### Recommended Policy

Define shared tolerance profiles for:

- unit-level matrix comparisons
- response invariance checks
- fixture comparisons
- export formatting where numeric rendering matters

### Why This Matters

Without a central policy:

- every test invents its own thresholds
- failures become hard to interpret
- contributors cannot tell whether a deviation is acceptable

### Suggested Implementation

Use named profiles such as:

- `Tolerance::strict()`
- `Tolerance::default()`
- `Tolerance::analysis()`

Tests should reference these profiles rather than embed unexplained constants.

## Test Data Strategy

The project should treat test data as a maintained asset.

### Good Test Data Includes

- low-order symmetric cases
- asymmetric cases
- cases with and without direct source-load coupling
- cases with multiple finite transmission zeros
- cases that stress transform edge conditions

### Coverage Goal

Do not rely only on easy symmetric examples. A library can look correct on
simple cases while failing badly on asymmetric or tightly constrained designs.

## Recommended Test Matrix by Phase

### Early Phase Coverage

Focus on:

- specification validation
- matrix invariants
- polynomial utilities
- canonical synthesis for low-order examples

### Mid-Phase Coverage

Add:

- response analysis regression
- folded and arrow invariance tests
- topology pattern checks
- transform report checks

### Advanced Phase Coverage

Add:

- trisection workflow tests
- section-merging tests
- minimum-path capacity checks
- advanced precondition enforcement tests

The suite should grow with features rather than waiting for a final testing
push.

## Suggested Test Directory Layout

```text
tests/
├─ fixtures/
│  ├─ canonical/
│  ├─ transforms/
│  ├─ responses/
│  └─ exports/
├─ spec_validation.rs
├─ approximation.rs
├─ canonical_transversal.rs
├─ response_analysis.rs
├─ transform_folded.rs
├─ transform_arrow.rs
├─ transform_trisection.rs
├─ topology_capacity.rs
├─ regression_invariance.rs
├─ export_json.rs
├─ export_markdown.rs
└─ cli_smoke.rs
```

This layout keeps tests aligned with user-visible capabilities.

## Example Test Scenarios

### Canonical Synthesis Scenario

Verify that:

- a valid order-4 specification synthesizes successfully
- the resulting matrix has transversal sparsity
- the response contains the expected number of poles and zeros

### Folded Conversion Scenario

Verify that:

- a canonical transversal matrix converts to folded form
- source and load side couplings are reduced to the expected pattern
- `S11` and `S21` match the precursor response within tolerance

### Arrow to Trisection Scenario

Verify that:

- a valid arrow matrix accepts trisection extraction
- the first extracted trisection satisfies the expected local pattern
- the transformed matrix preserves the original response

### Error Scenario

Verify that:

- calling trisection extraction on a folded matrix returns
  `InvalidTopologyPrecondition`

These scenarios should exist early and evolve as the implementation grows.

## Review Checklist for New Features

Every new algorithm should arrive with tests that answer at least these
questions:

- does it enforce its documented preconditions
- does it preserve core matrix invariants
- does it produce the intended topology pattern
- does it preserve or intentionally modify the electrical response
- does it work on at least one asymmetric or nontrivial case
- does it add a regression test if it fixes a previously known bug

This checklist should become part of code review culture.

## Anti-Patterns to Avoid

The test strategy should explicitly reject these weak practices:

- checking only that a transform ran without panicking
- checking only exact zeros in floating-point matrices
- using ad hoc tolerances in every file
- relying only on one golden example
- merging advanced topology features without invariance tests
- testing CLI output while leaving core math under-tested

Avoiding these patterns is essential for credibility.

## Continuous Integration Guidance

CI should run tests in layers.

### Fast Default CI

Run:

- unit tests
- core integration tests
- regression tests

### Extended CI

Run:

- slower fixture suites
- export snapshots
- CLI smoke tests
- benchmarks or performance comparisons where practical

This keeps iteration speed reasonable while still protecting the project.

## Summary

The recommended testing strategy is built around one central idea:

- a filter synthesis library must be tested as mathematics, software, and
  engineering workflow at the same time

The most important elements are:

- strong unit tests for math and matrix primitives
- integration tests for end-to-end workflows
- explicit response invariance tests for all topology transforms
- centralized tolerance policy
- fixture and regression discipline as the project grows

If the library follows this strategy, it will be much more likely to earn user
trust and remain maintainable as the feature set expands.
