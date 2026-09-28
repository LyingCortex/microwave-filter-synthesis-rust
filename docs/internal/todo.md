# MFS TODO List

**Status as of 0.5.0 (2026-09-27).** This file previously tracked work that has
since been completed or superseded; it is maintained as the single open-item
list. Closed items are listed with the version that closed them so history is
traceable. Authoritative architecture/limitation references:
`docs/architecture.md` (Known Limitations / Future Roadmap) and
`docs/internal/architecture-review-2026-09-25.md`.

## Open

### High priority (core functionality)

- **Passband transmission zeros** (`|zero| < 1`) — `generalized_ops.rs`
  rejects `|zero| < 1` with `NotImplemented` (`safe_sqrt_term` restriction).
  Needs the complex-zero handling for the stopband-only assumption.
- **Repeated poles in residue expansion** — `residues/expansion.rs` returns
  `NotImplemented` when a pole slope degenerates (no Laurent expansion).
- **Converge the full Macchiarella spectral path (order ≥ 4)** for
  N-channel mutual loading — `full_iteration(true)` exists on
  `DuplexerSynthesis` but does not converge at order ≥ 4 (typed error; order 2
  exhausts iterations). See `docs/multiplexer-code-review.md`.

### Medium priority (engineering quality)

- **HighPass / BandStop frequency mappings** — `freq` currently supports
  LowPass and BandPass only.
- **Fold the duplexer synthesis into the multiplexer (N=2)** — the
  duplicated `recover_filter_design` path is now a single
  `synthesis::recover_channel_design` (closed), but `DuplexerSynthesis`
  still carries its own 2-channel iteration loop; it could become an N=2
  `MultiplexerSynthesis` call so the Macchiarella-Tamiazzo iteration is
  maintained once.
- **High-order integration tests (N=8)** — root-solver unit tests cover
  degrees 10/20/30 and `examples/high_order_stability.rs` stresses up to 30,
  but the documented "order 2–8" fixture band has no N=8 regression entry in
  `tests/`.
- **Builder symmetry ergonomics** — `CouplingMatrixBuilder::set` writes a
  single entry (symmetry is the caller's job; `set_symmetric` exists).
  Consider a symmetric-by-default builder for the public API.

### Low priority (ecosystem)

- **Butterworth / Elliptic approximation families** in `approx`.
- **Multi-band frequency mapping**.
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
| CLI context builder duplication | 0.5.0+ | `mfs_cli.rs::build_context_from_request` replaced by `mfs::pipeline::context_from_request`, shared with `run_full_pipeline` |
