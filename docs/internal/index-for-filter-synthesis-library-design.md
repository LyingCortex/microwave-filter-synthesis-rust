# Index for Filter Synthesis Library Design Documents

## Purpose

This document is a navigation page for the filter synthesis design notes in
this repository.

The documents currently fall into two groups:

- paper-oriented Cameron notes
- software-architecture notes for a reusable open-source library

The first group is useful for understanding the algorithmic literature. The
second group is useful for turning those ideas into a reusable Rust library.

## Document status

Most of the "open-source library design" documents were written **before** the
current implementation (0.2–0.4 era) and describe *planned* boundaries, API
sketches, and phased roadmaps. They remain valuable as the rationale for the
module layout, but where they disagree with the code, the code wins.

Use this status column when reading:

| Status | Meaning |
|---|---|
| **Current** | Describes the shipped 0.6.0 code; verify against `src/` when in doubt |
| **Algorithm reference** | Literature notes (Cameron / Macchiarella); valid independent of code state |
| **Historical** | Design-time plan/sketch; superseded by the implementation, kept for context |

## Current Implementation Status

The codebase has adopted most of the library-facing structure described in
these design notes (as of 0.6.0):

- `FilterSpec` and `FilterSpecBuilder` are the active spec entry points
- the approximation layer exposes a generalized-only Chebyshev path
- `approx` is internally split into approximation engines, reusable
  complex-polynomial primitives, and generalized-domain helper modules
- `synthesis` is a directory-based subsystem with engine, residue, section,
  and orchestration modules (`build_filter_spec` is the single spec-builder
  shared by all entry points)
- `transform` is a dedicated subsystem for topology conversion and section
  extraction
- `verify` provides reusable response-invariance, section, and topology-shape
  checks
- `matrix` is a low-level matrix-domain module (core/builder/rotations/
  scaling/sections/cascade)
- `duplexer` / `multiplexer` / `spectral` / `optimize` / `touchstone` /
  `output` / `pipeline` / `fixtures` round out the crate
- the authoritative module map and data-flow diagrams live in
  [architecture.md](../architecture.md) (current)
- open work items are tracked in [todo.md](./todo.md) (current)
- the continuation plan for the active development areas (duplexer,
  multiplexer, multi-band, optimization) is
  [continuation-plan.md](../continuation-plan.md) (current)

The most recent assessment of the architecture is
[architecture-review-2026-09-25.md](./architecture-review-2026-09-25.md)
(current; its P1–P4 findings were remediated in 0.5.0).

The 0.6.0 whole-crate code review (2026-10-05) is
[review-findings-0.6.0.md](./review-findings-0.6.0.md) — 17 valid findings
(6 P1 / 5 P2 / 6 P3). The "silent errors" wave (P1-1/2/3/4/6, P2-1/2/3/5,
all six P3) was fixed and gate-anchored 2026-10-06; P1-5 + P2-4 (the CT/CQ
cascade re-implementation) remain open.

## Reading Paths

### If you want the algorithm background first

Start with:

1. [Cameron 2003 Coupling-Matrix Algorithm Notes](../cameron-coupling-matrix-algorithms.md) *(algorithm reference)*
2. [Cameron 2003 vs Current Rust Matrix Implementation](./cameron-vs-rust-matrix-implementation.md) *(historical gap analysis; re-check against code)*

This path is best if your first question is:

- what does the literature say the synthesis and reconfiguration flow should
  look like

### If you want the software design view first

Start with:

1. [architecture.md](../architecture.md) *(current — read this first)*
2. [Filter Synthesis Open-Source Library Architecture](./filter-synthesis-open-source-library-architecture.md) *(historical)*
3. [Suggested Rust Crate and Module Layout](./rust-crate-and-module-layout-for-filter-synthesis.md) *(historical)*
4. [Public API Sketch](./public-api-sketch-for-filter-synthesis-library.md) *(historical — the actual API is in `docs/api-reference.md` and rustdoc)*
5. [Contributor Guide](./contributor-guide-for-filter-synthesis-library.md) *(mostly current conventions)*

This path is best if your first question is:

- how should we build this into a real library from scratch

## Document Map

### Literature and Current-State Notes

- [Cameron 2003 Coupling-Matrix Algorithm Notes](../cameron-coupling-matrix-algorithms.md) — *(algorithm reference)*
  Summarizes the canonical `N+2` model, transversal synthesis, similarity
  transforms, folded and arrow forms, trisections, quartets, and box sections.
- [CT/CQ Topology Design](../ct-cq-topology-design.md) — *(algorithm reference, zh)*
  Chinese design notes for the cascaded-trisection / cascaded-quadruplet
  extraction implemented in `src/matrix/cascade.rs`.
- [Cameron 2003 vs Current Rust Matrix Implementation](./cameron-vs-rust-matrix-implementation.md) — *(historical)*
  Compares the intended Cameron workflow with the 0.3-era code structure.
  Several listed gaps have since been closed (e.g. `extract_trisection` now
  enforces the Arrow precondition).

### Open-Source Library Design Notes

- [Filter Synthesis Open-Source Library Architecture](./filter-synthesis-open-source-library-architecture.md) — *(historical)*
  Defines the overall domain model, pipeline, module boundaries, and design
  goals for a reusable library.
- [Suggested Rust Crate and Module Layout](./rust-crate-and-module-layout-for-filter-synthesis.md) — *(historical)*
  Proposes crate roles, module layout, testing layout, and visibility policy.
  The single-crate layout was kept; the per-module split was adopted.
- [Public API Sketch](./public-api-sketch-for-filter-synthesis-library.md) — *(historical)*
  API sketches from the design phase. Several proposed names
  (`extract_trisection_cascade`, `to_box_sections`, …) were not adopted;
  see `docs/api-reference.md` for the shipped API.
- [Implementation Roadmap](./implementation-roadmap-for-filter-synthesis-library.md) — *(historical)*
  Phases 1–7 from spec types through export/CLI; all phases landed.
- [Testing Strategy](./testing-strategy-for-filter-synthesis-library.md) — *(mostly current)*
  Test placement and invariance-check conventions still used in `tests/`.
- [Contributor Guide](./contributor-guide-for-filter-synthesis-library.md) — *(mostly current)*
  Contribution expectations, module boundaries, and review standards.
- [Migration Plan (Current Code → Target Architecture)](./current-code-to-target-architecture-migration-plan.md) — *(historical, 0.4 era)*
  Maps which target boundaries had landed; superseded by
  [architecture.md](../architecture.md) + the 2026-09-25 review.
- [Refactor Roadmap](./refactor-roadmap.md) — *(historical, input to the 0.5.0 remediation)*
- [Current Code Architecture Review 2026-09-25](./architecture-review-2026-09-25.md) — *(current, pre-0.6.0)*
- [Review Findings 0.6.0 (2026-10-05)](./review-findings-0.6.0.md) — *(current; 17 findings; "silent errors" wave remediated 2026-10-06, CT/CQ cascade P1-5/P2-4 open)*
- [Continuation Plan](../continuation-plan.md) — *(current, 0.6.0, 2026-10-05)*
  All Wave 1–5 items DONE: duplexer full Macchiarella path (D1–D4),
  multiplexer N-channel full path (B1–B4), true multi-passband synthesis
  (C0/C2/C3, dual-band route-α; K ≥ 3 refuted), and the LM optimizer
  (D1–D4-opt).
- [Development Guide](./development-guide.md) — *(mostly current)*
  Build/test/doc workflow for contributors.
- [Open-Source Release Checklist](./open-source-release-checklist.md) — *(current)*
- [Algorithm Roadmap](./algorithm-roadmap.md) — *(mostly current)*
  Status of approximation families (Chebyshev done; Butterworth/Elliptic open).
- [Literature Fixture Catalog](../literature-fixture-catalog.md) — *(current)*
  Cameron/Macchiarella reference fixtures used by the test suite.

### One-Off Records (keep, do not treat as living docs)

- [First GitHub Release Notes](./first-github-release-notes.md) — *(historical, 0.1 era)*
- [GitHub Profile Template](./github-profile-template.md) — *(template)*
- [api-sketch.md](./api-sketch.md) — *(historical superseded sketch)*

## Suggested Order for Contributors

If you are about to implement or refactor code:

1. [architecture.md](../architecture.md)
2. [architecture-review-2026-09-25.md](./architecture-review-2026-09-25.md)
3. [todo.md](./todo.md)
4. [continuation-plan.md](../continuation-plan.md) *(current — next steps for duplexer / multiplexer / multi-band / optimization)*
5. [Contributor Guide](./contributor-guide-for-filter-synthesis-library.md)
6. [Development Guide](./development-guide.md)
7. [Testing Strategy](./testing-strategy-for-filter-synthesis-library.md)
8. [Filter Synthesis Open-Source Library Architecture](./filter-synthesis-open-source-library-architecture.md)
9. [Cameron 2003 Coupling-Matrix Algorithm Notes](../cameron-coupling-matrix-algorithms.md)

## Suggested Order for Research-Oriented Readers

If you are validating the mathematics first:

1. [Cameron 2003 Coupling-Matrix Algorithm Notes](../cameron-coupling-matrix-algorithms.md)
2. [CT/CQ Topology Design](../ct-cq-topology-design.md)
3. [Cameron 2003 vs Current Rust Matrix Implementation](./cameron-vs-rust-matrix-implementation.md)
4. [Filter Synthesis Open-Source Library Architecture](./filter-synthesis-open-source-library-architecture.md)
5. [Public API Sketch](./public-api-sketch-for-filter-synthesis-library.md)
6. [Testing Strategy](./testing-strategy-for-filter-synthesis-library.md)
7. [Contributor Guide](./contributor-guide-for-filter-synthesis-library.md)

## Summary

This index is meant to keep the design documents readable as a set rather than
as isolated notes.

In short:

- read the Cameron notes for algorithm context
- read `docs/architecture.md` (current) and the 2026-09-25 review for the
  current software structure
- read the historical design notes for *why* the boundaries look the way they
  do, not for what the code should look like
