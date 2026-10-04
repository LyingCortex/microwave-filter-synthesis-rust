# MFS TODO List

**Status as of 0.5.0 (2026-09-27).** This file previously tracked work that has
since been completed or superseded; it is maintained as the single open-item
list. Closed items are listed with the version that closed them so history is
traceable. Authoritative architecture/limitation references:
`docs/architecture.md` (Known Limitations / Future Roadmap) and
`docs/internal/architecture-review-2026-09-25.md`.

## Current state (2026-09-28)

Working tree clean, `master` == `gitlab/master` (head `2d735c0`), all
CI-parity checks green on that head:

- `cargo test --all-targets` — green (lib 202 passed + 15 integration
  binaries, 0 failed)
- `cargo clippy --all-targets -- -D warnings` — zero warnings
- `cargo doc --no-deps` — zero warnings (a `rustdoc::private_intra_doc_links`
  warning in `FaerRootSolver`'s doc, introduced by the faer residual gate,
  was fixed by the 2026-09-28 docs-sync commit)

The uncommitted work found after the 2026-09-28 session interruption was
verified complete (not lost): functional code (faer backward-residual gate,
DK/Aberth NaN guards, `recover_channel_design`, `context_from_request`,
benchmark standard) was already committed in `3dab2ff`..`24e5880`; the
remaining 62 files were one whole-crate `cargo fmt` pass plus the doc/
changelog syncs, committed as `df8ee83` (style), `a59e055`
([Unreleased] changelog), `75c6dba` (docs sync) and pushed to
`gitlab/master`.

Next steps: the continuation plan for the four active areas (duplexer,
multiplexer, multi-band, optimization) is
[`docs/continuation-plan.md`](../continuation-plan.md) (2026-09-28). The Open
sections below remain the flat open-item list. The `[Unreleased]` changelog
section should be released as the next version tag when the maintainer cuts it.

## Open

### High priority (core functionality)

- **Passband transmission zeros** (`|zero| < 1`) — `generalized_ops.rs`
  rejects `|zero| < 1` with `NotImplemented` (`safe_sqrt_term` restriction).
  Needs the complex-zero handling for the stopband-only assumption.
- **Repeated poles in residue expansion** — `residues/expansion.rs` returns
  `NotImplemented` when a pole slope degenerates (no Laurent expansion).
- **Generalize the full Macchiarella spectral path to N channels** for
  mutual loading — the **duplexer** portion is done (2026-10-03):
  `full_iteration(true)` converges and is gate-anchored at Type-I orders 2–8
  and Type-II (incl. the §VI degree-40 GSM design). What remains is the
  N > 2 multiplexer (continuation plan B2). See `docs/multiplexer-code-review.md`.

### Medium priority (engineering quality)

- **Multi-band synthesis (scope-1, decided 2026-09-29)** - true
  multi-passband synthesis, not per-band composition. Mapping layer done
  (`MultiBandMapping` in `freq`, 2026-09-29); **C0 design doc done
  2026-10-05** (`docs/internal/multi-band-synthesis-design.md`, route-α
  committed); **E1 investigated 2026-10-05 and is NOT a blocker** (gap TZs
  are ordinary stopband zeros — the C0 doc's "passband" claim was wrong,
  corrected); next is **C2a/C2b/C2c** then C3. Sequencing and acceptance in
  `docs/continuation-plan.md` section C.
- **Fold the duplexer synthesis into the multiplexer (N=2)** — the
  duplicated `recover_filter_design` path is now a single
  `synthesis::recover_channel_design` (closed), but `DuplexerSynthesis`
  still carries its own 2-channel iteration loop; it could become an N=2
  `MultiplexerSynthesis` call so the Macchiarella-Tamiazzo iteration is
  maintained once.
- **Builder symmetry ergonomics** — `CouplingMatrixBuilder::set` writes a
  single entry (symmetry is the caller's job; `set_symmetric` exists).
  Consider a symmetric-by-default builder for the public API.

### Low priority (ecosystem)

- **Butterworth / Elliptic approximation families** in `approx`.
- **Real Wheel topology transform** — `transform/wheel.rs` is still an alias
  of the Arrow reduction path (deprecated enum variant).
- **Docs.rs polish** — more doc examples on the module APIs (crate docs and
  `docs/api-reference.md` already cover the main workflows).

## Closed

| Item | Closed in | Notes |
|---|---|---|
| Touchstone (.s2p) export + import | 0.5.0 | `src/touchstone.rs`, `build_network` / `cascade` on top of touchstone 0.14.5 |
| Python bindings (pyo3) | 0.5.0 | `src/python.rs` behind `python` feature, `pyproject.toml` version-synced |
| CLI tool | 0.5.0 | `src/bin/mfs_cli.rs` (design / multiplexer subcommands, stage/resume) |
| Performance benchmarks | 0.5.0 | `benches/synthesis.rs`, `benches/root_solvers.rs` (criterion) |
| `coupling_matrix.rs` split | 0.5.0 | now the `src/matrix/` module (core/builder/rotations/scaling/sections/cascade) |
| Unified error semantics | 0.5.0 | `MfsError` has 9 purpose-specific variants; no `Unsupported` catch-all |
| Placeholder fallback | 0.5.0 | removed; failures are typed `MfsError`s |
| Structured convergence errors | 0.5.0 | `NumericalFailure` carries solver name + residual; all three solvers residual/finite-gated |
| `unwrap_or_default` in matrix access | 0.5.0 | `CouplingMatrix::get` has `debug_assert` bounds checks; remaining usages are in tests / report formatting |
| Root-residual validation in `FaerRootSolver` | 0.5.0+ | backward-residual gate (`≤ 1e-3`) + non-finite check; fall-through to iterative methods; DK/Aberth in-loop NaN guards |
| Duplicate `recover_filter_design` path | 0.5.0+ | single `mfs::synthesis::recover_channel_design` used by both duplexer and multiplexer |
| `HighPassMapping` / `BandStopMapping` | 0.5.0+ | `freq` now covers all four prototype families; singular points (`f = 0` for HP, `f = f0` for BS) clamp to `-SINGULAR_CLAMP`, edges map to ±1 exactly under the geometric-mean convention (E3, continuation plan) |
| High-order integration tests (N=8) | 0.5.0+ | `tests/benchmark_standard.rs::mfs_synthesized_folded_order8_reproduces_the_printed_book_matrix` pins order 8 against the printed Cameron 2018 `EIGHT_THREE_FOLDED` anchor (E4, continuation plan) |
