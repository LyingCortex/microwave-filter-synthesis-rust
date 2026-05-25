# Literature Fixture Catalog

## Purpose

This document records the small literature-backed and literature-shaped
fixtures currently used by the crate.

The goal is not to claim that every numeric value here is copied from a
published table. The goal is to keep a stable set of low-order reference cases
whose workflow intent is anchored in the Cameron/ZTE literature and whose
provenance is explicit.

The generalized helper layer that feeds these fixtures now uses
`num-complex` for baseline complex arithmetic. The fixtures therefore anchor
domain outputs such as recurrence coefficients, `F(s)`, `P(s)`, `A(s)`,
`E(s)`, and ripple parameters rather than anchoring a project-local complex
number implementation.

The current generalized helper pipeline also exposes structured `A`-stage and
`E`-stage intermediates. That means fixtures can now stabilize not only final
`A(s)` / `E(s)` results, but also whether the helper preserved the expected
internal stage shape.

## Current Fixtures

### `cameron_generalized_order4_spec`

Location:

- [src/fixtures/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/fixtures/mod.rs)

Purpose:

- stable generalized-Chebyshev low-order synthesis case
- drives orchestration through the generalized approximation path
- expected to use residue-based matrix synthesis in the current implementation

Reference intent:

- Cameron 2003 generalized coupling-matrix workflow

Expected behavior:

- selects the generalized approximation path
- exposes helper-backed generalized polynomial metadata
- currently resolves through residue-based matrix synthesis

### `cameron_style_section_polynomials`

Location:

- [src/fixtures/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/fixtures/mod.rs)

Purpose:

- stable low-order section-synthesis seed
- reused for triplet and trisection extraction
- paired with normalized response-invariance checks

Reference intent:

- Cameron/ZTE-style progression from canonical matrix to advanced extracted
  sections

Expected behavior:

- supports reported triplet synthesis with response checks
- supports reported trisection synthesis with response checks
- acts as a stable regression anchor for section workflows

### `literature_reference_grid`

Location:

- [src/fixtures/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/fixtures/mod.rs)

Purpose:

- shared normalized sweep used by transform and section-transform reference
  checks

Expected behavior:

- keeps response-invariance checks on a stable normalized sweep
- provides repeatable regression coverage for literature-style workflows

### `cameron_single_zero_exact_case`

Location:

- [src/fixtures/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/fixtures/mod.rs)

Purpose:

- exact low-order numeric anchor for the generalized helper layer
- fixes one concrete Cameron recurrence case in code
- gives the crate one true exact-value literature-shaped reference today

Reference intent:

- Cameron 2003 recurrence and generalized-polynomial helper flow

Expected behavior:

- reproduces the first-order Cameron recurrence exactly
- reproduces the low-order `P(s)` coefficients exactly
- anchors the helper epsilon computation to a fixed numeric case

### `cameron_order3_generalized_pipeline_exact_case`

Location:

- [src/fixtures/mod.rs](/c:/Users/eynulai/Downloads/mfs/src/fixtures/mod.rs)

Purpose:

- exact order-3 generalized helper pipeline anchor
- fixes one complete low-order helper output in code
- protects `F(s)`, `P(s)`, `A(s)`, `E(s)`, and ripple parameters together

Reference intent:

- Cameron 2003 generalized helper pipeline for low-order single-zero cases

Expected behavior:

- reproduces the full low-order generalized helper pipeline
- stabilizes `F(s)`, `P(s)`, `A(s)`, and `E(s)` against accidental drift
- preserves observable `A`-stage and `E`-stage helper artifacts
- anchors generalized helper outputs beyond workflow-only checks

## Current Usage

The current fixture-backed regression coverage lives in:

- [tests/literature_fixtures.rs](/c:/Users/eynulai/Downloads/mfs/tests/literature_fixtures.rs)

These tests currently protect:

- generalized main-flow orchestration selection
- triplet section synthesis with response checking
- trisection section synthesis with response checking
- one exact low-order Cameron recurrence and epsilon case
- one exact order-3 generalized helper pipeline case
- presence of structured `A`-stage and `E`-stage generalized artifacts

## Scope Boundary

At the current stage, these fixtures should be interpreted as:

- literature-backed workflow fixtures
- stable regression anchors
- one exact low-order numeric anchor for the generalized helper layer
- one exact order-3 helper-pipeline anchor for generalized polynomial data

They should not yet be interpreted as:

- exact reproductions of a published coupling table
- complete validation against every polynomial or physical realization formula

That stronger claim should wait until the crate carries exact published
reference cases with source-to-output traceability.
