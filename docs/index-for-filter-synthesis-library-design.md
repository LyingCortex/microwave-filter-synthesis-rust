# Index for Filter Synthesis Library Design Documents

## Purpose

This document is a navigation page for the filter synthesis design notes in
this repository.

The documents currently fall into two groups:

- paper-oriented Cameron notes
- software-architecture notes for a reusable open-source library

The first group is useful for understanding the algorithmic literature. The
second group is useful for turning those ideas into a maintainable Rust code
base.

## Current Implementation Status

The codebase has now adopted most of the library-facing structure described in
these design notes.

What is already reflected in `src/`:

- `FilterSpec` and `FilterSpecBuilder` are the active spec entry points
- `ApproximationFamily` now distinguishes explicit Chebyshev and generalized
  Chebyshev requests
- `approx` is now internally split into approximation engines, reusable
  complex-polynomial primitives, and generalized-domain helper modules
- `synthesis` is a directory-based subsystem with canonical, section, residue,
  engine, and orchestration modules
- `transform` is a dedicated subsystem for topology conversion and section extraction
- `verify` provides reusable response-invariance, section, and topology-shape checks
- `matrix` has topology metadata and has been reduced toward a lower-level matrix-domain layer
- transform workflows now emit minimal reports, with optional response-comparison summaries
- orchestration outcomes now surface `ApproximationStageKind` and
  `MatrixSynthesisMethod`
- `prelude` exposes a coherent high-level user workflow

The most up-to-date implementation snapshot is:

- [Current Code to Target Architecture Migration Plan](/c:/Users/eynulai/Downloads/mfs/docs/current-code-to-target-architecture-migration-plan.md)
  Summarizes which target boundaries have already landed in code and which
  feature gaps remain.

## Reading Paths

### If you want the algorithm background first

Start with:

1. [Cameron 2003 Coupling-Matrix Algorithm Notes](/c:/Users/eynulai/Downloads/mfs/docs/cameron-coupling-matrix-algorithms.md)
2. [Cameron 2003 vs Current Rust Matrix Implementation](/c:/Users/eynulai/Downloads/mfs/docs/cameron-vs-rust-matrix-implementation.md)

This path is best if your first question is:

- what does the literature say the synthesis and reconfiguration flow should
  look like

### If you want the software design view first

Start with:

1. [Filter Synthesis Open-Source Library Architecture](/c:/Users/eynulai/Downloads/mfs/docs/filter-synthesis-open-source-library-architecture.md)
2. [Suggested Rust Crate and Module Layout for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/rust-crate-and-module-layout-for-filter-synthesis.md)
3. [Public API Sketch for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/public-api-sketch-for-filter-synthesis-library.md)
4. [Implementation Roadmap for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/implementation-roadmap-for-filter-synthesis-library.md)
5. [Testing Strategy for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/testing-strategy-for-filter-synthesis-library.md)
6. [Contributor Guide for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/contributor-guide-for-filter-synthesis-library.md)

This path is best if your first question is:

- how should we build this into a real library from scratch

## Document Map

### Literature and Current-State Notes

- [Cameron 2003 Coupling-Matrix Algorithm Notes](/c:/Users/eynulai/Downloads/mfs/docs/cameron-coupling-matrix-algorithms.md)
  Summarizes the canonical `N+2` model, transversal synthesis, similarity
  transforms, folded and arrow forms, trisections, quartets, and box sections.

- [Cameron 2003 vs Current Rust Matrix Implementation](/c:/Users/eynulai/Downloads/mfs/docs/cameron-vs-rust-matrix-implementation.md)
  Compares the intended Cameron workflow with the current code structure and
  highlights alignment gaps.

### Open-Source Library Design Notes

- [Filter Synthesis Open-Source Library Architecture](/c:/Users/eynulai/Downloads/mfs/docs/filter-synthesis-open-source-library-architecture.md)
  Defines the overall domain model, pipeline, module boundaries, and design
  goals for a reusable library.

- [Suggested Rust Crate and Module Layout for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/rust-crate-and-module-layout-for-filter-synthesis.md)
  Proposes workspace structure, crate roles, module layout, testing layout,
  and visibility policy for a Rust implementation.

- [Public API Sketch for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/public-api-sketch-for-filter-synthesis-library.md)
  Sketches the user-facing Rust API, including builders, synthesis flow,
  topology transforms, analysis, verification, and export patterns.

- [Implementation Roadmap for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/implementation-roadmap-for-filter-synthesis-library.md)
  Breaks the project into phases from core domain types through advanced
  topology support, with suggested version milestones and exit criteria.

- [Testing Strategy for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/testing-strategy-for-filter-synthesis-library.md)
  Defines how to test mathematical correctness, topology structure, response
  invariance, regression behavior, exports, and CI layers.

- [Contributor Guide for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/contributor-guide-for-filter-synthesis-library.md)
  Explains contribution expectations, module boundaries, algorithm-addition
  rules, test requirements, and code review standards.

- [Current Code to Target Architecture Migration Plan](/c:/Users/eynulai/Downloads/mfs/docs/current-code-to-target-architecture-migration-plan.md)
  Maps the current implementation to the target architecture and records what
  has already been refactored into place.

- [Literature Fixture Catalog](/c:/Users/eynulai/Downloads/mfs/docs/literature-fixture-catalog.md)
  Records the current Cameron/ZTE-oriented reference fixtures used in tests.

## Suggested Order for Contributors

If you are about to implement or refactor code, the most practical reading
order is:

1. [Filter Synthesis Open-Source Library Architecture](/c:/Users/eynulai/Downloads/mfs/docs/filter-synthesis-open-source-library-architecture.md)
2. [Suggested Rust Crate and Module Layout for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/rust-crate-and-module-layout-for-filter-synthesis.md)
3. [Public API Sketch for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/public-api-sketch-for-filter-synthesis-library.md)
4. [Implementation Roadmap for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/implementation-roadmap-for-filter-synthesis-library.md)
5. [Testing Strategy for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/testing-strategy-for-filter-synthesis-library.md)
6. [Contributor Guide for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/contributor-guide-for-filter-synthesis-library.md)
7. [Current Code to Target Architecture Migration Plan](/c:/Users/eynulai/Downloads/mfs/docs/current-code-to-target-architecture-migration-plan.md)
8. [Cameron 2003 Coupling-Matrix Algorithm Notes](/c:/Users/eynulai/Downloads/mfs/docs/cameron-coupling-matrix-algorithms.md)

This order starts from software boundaries and then returns to the literature
details when implementing specific transforms.

## Suggested Order for Research-Oriented Readers

If you are validating the mathematics first, a better order is:

1. [Cameron 2003 Coupling-Matrix Algorithm Notes](/c:/Users/eynulai/Downloads/mfs/docs/cameron-coupling-matrix-algorithms.md)
2. [Cameron 2003 vs Current Rust Matrix Implementation](/c:/Users/eynulai/Downloads/mfs/docs/cameron-vs-rust-matrix-implementation.md)
3. [Filter Synthesis Open-Source Library Architecture](/c:/Users/eynulai/Downloads/mfs/docs/filter-synthesis-open-source-library-architecture.md)
4. [Public API Sketch for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/public-api-sketch-for-filter-synthesis-library.md)
5. [Current Code to Target Architecture Migration Plan](/c:/Users/eynulai/Downloads/mfs/docs/current-code-to-target-architecture-migration-plan.md)
6. [Testing Strategy for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/testing-strategy-for-filter-synthesis-library.md)
7. [Contributor Guide for a Filter Synthesis Library](/c:/Users/eynulai/Downloads/mfs/docs/contributor-guide-for-filter-synthesis-library.md)

This order starts with the algorithm model and then moves toward
implementation.

## Current Gaps

These documents now cover:

- algorithm notes
- architecture
- Rust layout
- public API
- implementation roadmap
- testing strategy
- contributor guidance

Natural next documents, if needed later, would be:

- coding conventions for numerical and matrix code
- glossary of coupling-matrix terminology
- fixture catalog for literature-derived reference cases
- physical realization adapter design

## Summary

This index is meant to keep the design documents readable as a set rather than
as isolated notes.

In short:

- read the Cameron notes for algorithm context
- read the architecture and API notes for software structure
- read the roadmap and testing notes before starting implementation work
