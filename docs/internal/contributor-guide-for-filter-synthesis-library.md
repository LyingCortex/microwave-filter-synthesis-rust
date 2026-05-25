# Contributor Guide for a Filter Synthesis Library

## Purpose

This document is a contributor guide for an open-source microwave filter
synthesis library.

Its purpose is to help contributors understand:

- how the project is organized
- what engineering standards new code should meet
- how to add or change algorithms without weakening the mathematical model

This guide is written for contributors working on:

- approximation methods
- canonical matrix synthesis
- topology transforms
- response analysis
- verification infrastructure
- export and tooling

## Contribution Philosophy

This project should be maintained as a mathematical software library, not as a
collection of loosely connected research scripts.

That means contributions should aim for:

- correctness
- explicit assumptions
- reproducibility
- inspectability
- stable interfaces

The most important project rule is:

- every algorithmic claim should be backed by tests that check both structure
  and electrical behavior where applicable

## Read Before Contributing

New contributors should read these documents first:

1. [Index for Filter Synthesis Library Design Documents](/c:/Users/eynulai/Downloads/mfs/docs/index-for-filter-synthesis-library-design.md)
2. [Filter Synthesis Open-Source Library Architecture](/c:/Users/eynulai/Downloads/mfs/docs/filter-synthesis-open-source-library-architecture.md)
3. [Suggested Rust Crate and Module Layout for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/rust-crate-and-module-layout-for-filter-synthesis.md)
4. [Public API Sketch for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/public-api-sketch-for-filter-synthesis-library.md)
5. [Testing Strategy for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/testing-strategy-for-filter-synthesis-library.md)

If you are changing transform logic, also read:

6. [Cameron 2003 Coupling-Matrix Algorithm Notes](/c:/Users/eynulai/Downloads/mfs/docs/cameron-coupling-matrix-algorithms.md)

If you are comparing proposed changes with the current implementation, also
read:

7. [Cameron 2003 vs Current Rust Matrix Implementation](/c:/Users/eynulai/Downloads/mfs/docs/cameron-vs-rust-matrix-implementation.md)

## Types of Contributions

Good contributions include:

- new approximation methods behind a clear interface
- improvements to canonical matrix synthesis
- new topology transforms with documented preconditions
- stronger verification and invariance checks
- better fixtures from literature-backed examples
- cleaner exports and reporting
- examples and contributor-facing documentation

Less helpful contributions include:

- adding a large feature without tests
- bypassing canonical forms to reach a desired topology quickly
- introducing placeholder heuristics without clear labeling
- exposing unstable internal helpers as public API without a strong reason

## Architectural Boundaries

Contributors should preserve the project's core module boundaries.

### `spec`

Owns:

- user input model
- validation of external specifications

Should not own:

- matrix rotation logic
- topology-specific assumptions

### `approx`

Owns:

- approximation generation
- filtering polynomials and approximation outputs

Should not own:

- topology reconfiguration logic

### `matrix`

Owns:

- `N+2` coupling-matrix representation
- indexing rules
- symmetry and invariant enforcement
- low-level rotation primitives

Should not own:

- high-level folded, arrow, or trisection workflows

### `synthesis`

Owns:

- construction of canonical coupling matrices

Should not own:

- advanced topology extraction

### `transform`

Owns:

- topology-preserving matrix reconfiguration
- transform reports

Should not own:

- approximation generation

### `analysis`

Owns:

- ideal electrical response derived from matrices

### `verify`

Owns:

- invariance checks
- topology constraint checks
- precondition checks
- transmission-zero capacity checks

### `export`

Owns:

- serialization and reporting outputs

Should not own:

- business logic that belongs in synthesis or transform code

Contributions that blur these boundaries should be treated cautiously even if
the code appears to work.

## How to Add a New Algorithm

When adding a new algorithm, contributors should follow a predictable process.

### Step 1: Define the Role of the Algorithm

Clarify:

- what problem it solves
- what its inputs are
- what canonical form it requires
- what output it guarantees
- whether it preserves response or intentionally changes it

If this cannot be stated clearly, the algorithm is not ready to merge.

### Step 2: Choose the Correct Module

Examples:

- approximation generation belongs in `approx`
- canonical transversal synthesis belongs in `synthesis`
- arrow conversion belongs in `transform`
- response comparison belongs in `verify`

Avoid placing code where it is convenient rather than where it belongs.

### Step 3: Encode Preconditions

If the algorithm expects:

- arrow form
- adjacent trisections
- source-load coupling
- asymmetric polynomial support

then the code should check that explicitly.

Do not rely on comments alone to express required state.

### Step 4: Add Tests

At minimum, add:

- one success-path test
- one invalid-precondition test
- one regression or invariance test if the algorithm transforms a matrix

### Step 5: Document the Algorithm

Document:

- its precursor assumptions
- its output form
- its tolerance expectations
- the literature or rationale behind it

This can be a module-level note, doc comments, or a companion design document.

## Mathematical and Numerical Standards

Contributors should treat numerical discipline as part of correctness.

### Required Practices

- use shared tolerance helpers
- avoid exact zero checks on floating-point matrix entries
- preserve symmetry intentionally
- compare responses with explicit tolerances
- make thresholding rules visible in code

### Avoid

- ad hoc epsilon values copied into individual functions
- silent normalization or sign changes without documentation
- introducing alternate indexing conventions in local code

Small numerical shortcuts can create large downstream confusion.

## Transform-Specific Contribution Rules

Any contribution that changes topology transform behavior must satisfy extra
requirements.

### Required

- preserve matrix invariants
- preserve the ideal response when the transform is meant to be similarity-
  based
- produce a topology pattern consistent with the documented target form
- emit or update transform reports where relevant

### Tests Required

- topology shape test
- response invariance test
- invalid-precondition test where applicable

### Documentation Required

- which precursor form is required
- which couplings are expected to be removed or created
- whether any sign normalization occurs

This is one of the most sensitive parts of the project.

## Public API Contribution Rules

Changes to the public API should be conservative.

### Good Reasons to Expand the Public API

- a workflow is impossible without awkward internal access
- a type is clearly becoming a stable extension point
- repeated contributor use shows that an internal helper should become public

### Bad Reasons

- a test needs easier access
- one implementation detail is temporarily inconvenient
- the internal design has not stabilized yet

When in doubt:

- prefer `pub(crate)` first

Public APIs are expensive to maintain once users adopt them.

## Documentation Expectations

Contributor-facing code should be explainable without reading the entire code
base.

Contributors should add or update documentation when they:

- introduce a new algorithm
- change preconditions
- add a new topology kind
- alter report behavior
- introduce a new fixture family

Good documentation includes:

- what the code does
- why it belongs in this module
- what assumptions it depends on

## Test Expectations

Before proposing a change, contributors should ensure that:

- existing relevant tests pass
- new logic has focused tests
- invariance is tested for topology-changing algorithms
- error behavior is tested for invalid states

If a bug was fixed, a regression test should be added whenever practical.

The goal is not just to prove the new code works once, but to keep the project
from drifting later.

## Fixture Contribution Guidelines

Fixtures are part of the project's technical knowledge base.

When adding a fixture:

- name it clearly
- describe where it came from
- state whether it is literature-derived, hand-constructed, or generated
- document what behavior it is meant to exercise

Avoid adding opaque fixtures that nobody can interpret later.

## Suggested Pull Request Checklist

Contributors should be able to answer yes to these questions before opening or
merging a pull request:

- does the change respect module boundaries
- are preconditions explicit in code
- are tolerances handled through shared utilities
- does the change include the right level of tests
- if a topology changes, is electrical invariance checked
- if a public API changed, is that change intentional and documented
- if a bug was fixed, was a regression test added
- is the algorithm source or rationale documented

This checklist is short on purpose. It should be usable in practice.

## Code Review Expectations

Reviewers should focus on:

- mathematical correctness
- precondition clarity
- test sufficiency
- API stability
- module placement

Reviewers should be cautious about approving code that:

- "works" only on one fixture
- silently changes numerical semantics
- adds a transform without invariance coverage
- solves a topology problem by bypassing canonical precursor forms

The standard should be supportive but technically firm.

## Good First Contribution Areas

For new contributors, safer starting points include:

- adding small fixtures
- improving error messages
- strengthening validation tests
- adding export formats
- improving examples and docs
- writing response-comparison helpers

These contributions help the project without requiring immediate changes to the
most delicate mathematical code.

## Advanced Contribution Areas

More advanced contributions include:

- new approximation families
- canonical synthesis improvements
- folded and arrow transform refinements
- trisection or quartet workflows
- numerical robustness improvements in analysis

These usually require deeper familiarity with both the literature and the
library architecture.

## Anti-Patterns to Avoid

Contributors should avoid:

- hard-coding special cases into general algorithms without documentation
- introducing new matrix conventions in one module only
- adding advanced topology features before their precursor workflows are stable
- modifying public APIs casually
- merging transform logic that is only structurally tested
- relying on comments instead of enforceable code checks

Avoiding these is often more important than adding new features quickly.

## Communication Guidance

When discussing a contribution, it helps to state changes in this format:

- problem being solved
- precursor assumptions
- module being changed
- tests added
- expected effect on public API

This makes review faster and reduces ambiguity.

## Summary

The project should welcome contributions, but with strong guardrails.

The key expectations are:

- respect module boundaries
- make assumptions explicit
- test both structure and electrical behavior
- treat numerical discipline as part of correctness
- expand the public API carefully

If contributors follow these rules, the library can grow without losing the
trustworthiness that mathematical engineering software depends on.
