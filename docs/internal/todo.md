# MFS TODO List

**Status as of 0.6.1 (2026-10-07).** This file is maintained as the single
open-item list. Closed items are listed with the version that closed them so
history is traceable. Authoritative architecture/limitation references:
`docs/architecture.md` (Known Limitations / Future Roadmap) and
`docs/internal/architecture-review-2026-09-25.md`.

## Current state (0.6.1, 2026-10-07)

Working tree clean; `master` head `ec4482c` (Wave 7: cross-backend gate).
All CI-parity checks green on that head:

- `cargo test --all-targets` — green (lib 256 passed + 15 integration
  test binaries, 0 failed)
- `cargo clippy --all-targets -- -D warnings` — zero warnings
- `cargo doc --no-deps` — zero warnings

0.6.1 closes the 0.6.0 whole-crate review end-to-end: the "silent errors"
wave (P1-1/2/3/4/6 + P2-1/2/3/5 + all six P3, the "gate-green but wrong
result" class) and the CT/CQ cascade wave (P1-5 + P2-4), plus the follow-up
clean CQ quadruplet tiling (source-side non-overlapping, `order >= 4·pairs`,
order-6 2-pair rejected as topologically impossible) and a fresh-eyes review
that corrected the CQ load-port gate index and two doc/code contradictions.
All 17 valid review findings are now fixed and gate-anchored. Passband
transmission zeros remain refuted as unrealizable in the single-function
prototype (E1, 0.6.0).

The continuation plan for the four areas is
[`docs/continuation-plan.md`](../continuation-plan.md); all of its Wave 1–7
items are done (2026-09-29 → 2026-10-07). Wave 7 (2026-10-07, `ec4482c`)
carries the standard's `cross_backend_equivalence` gate in the rust-native
`synthesis` bench: a pinned python-backend eigenvalue reference
(`bench-data/golden/cross_backend_reference.json`, vendored from
`filter-bench/data/golden/`) is compared against mfs's folded
resonator-subblock eigenvalues at the contract tolerance (1e-6); measured
≤ 2.7e-12 across all 7 twoport cases, with a conformance test in
`tests/benchmark_standard.rs`. The Open sections below hold the
remaining low-priority ecosystem items and the refuted passband-TZ work.

## Open

### High priority (core functionality)

- **Passband transmission zeros** (`|zero| < 1`) — still rejected at the
  `approx` entry point (typed `MfsError::InvalidTransmissionZero`). The 0.6.0
  E1 investigation refuted the in-band-null prototype in the single-function
  form (Feldtkeller violation); realizing them needs the multifunction
  characteristic, which is genuine new approximation research, not a helper
  relaxation.


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
| Full Macchiarella spectral path (duplexer, Type-I 2–8 + Type-II) | 0.6.0 | `full_iteration(true)` converges and is gate-anchored (orders 2–8 + the §V order-7 golden + the §VI degree-40 GSM anchor); continuation plan D1–D4 |
| Repeated-pole (Laurent) residue expansion | 0.6.0 | E2 — `residues/expansion.rs` Laurent principal part + per-Y fast-path summation; transversal rejects a repeated pole typed; `tests/synthesis_properties.rs` round-trip |
| N > 2 multiplexer full path (B1–B4) | 0.6.0 | true (N+1)-port response (B1), shared Macchiarella full path generalized to N channels (B2), convergence + fidelity gates at N = 3/4 (B3), per-gate cm/star unitarity gates (B4); `tests/multiplexer.rs` |
| Multi-band synthesis (scope-1, dual-band route-α) | 0.6.0 | C0 design doc (route committed) + C2a/C2b/C2c (`multi_band_chebyshev_polynomials` → `MatrixSynthesisEngine`) + C3 (`tests/multi_band.rs` Fig. 25 anchor); K ≥ 3 refuted as a degrees-of-freedom wall (E1-style) |
| Passband transmission zeros (`|zero| < 1`) | 0.6.0 (refuted) | E1: no valid in-band-null prototype in the single-function real-rooted form (Feldtkeller violation); `approx` rejects typed. Realizing them needs the multifunction characteristic (new research) |
| 0.6.0 whole-crate review (17 findings) — "silent errors" wave + CT/CQ cascade + clean CQ tiling + fresh-eyes gate/doc fix | 0.6.1 | all 17 valid findings (6 P1 / 5 P2 / 6 P3) fixed and gate-anchored: `recovery.rs` rescale_d, `scaling.rs` denorm sign, `response/backend.rs` group-delay port nodes, `machiarella.rs` p0² rejection, `cascade.rs` CT (Fig. 21) + CQ (source-side tiling, `order ≥ 4·pairs`), `python.rs` length checks, CLI/schema/validator P2s, P3s; the CQ load-port gate index and two doc/code contradictions corrected by a fresh-eyes review (2026-10-07) |
| Fold duplexer into multiplexer (N=2) | 0.6.1 | already satisfied: the full Macchiarella–Tamiazzo iteration (the item's stated goal) is maintained once in `src/synthesis/machiarella` and shared by both `DuplexerSynthesis` (full path) and `MultiplexerSynthesis` (N≥2, B2). The only non-shared loop is the *simplified* perturbation heuristic (duplexer bespoke 2-channel vs multiplexer N-channel), which the parity test `multiplexer_diplexer_parity_with_duplexer` documents as legitimately divergent (asserts band-placement + validity only, not value-identity). Measured 2026-10-07 (throwaway probe): forcing the fold would change the duplexer default path's recovered matrices by O(1–5) (ideal Type-I max ΔM 2.44, Type-II 5.54); the multiplexer's simplified N=2 path exhausts its 40-iter budget (conv_err 2.7e-3, unconverged) and leaves its recovered centers at the input, versus the duplexer's perturbed 15-iter convergence — a behavior change, not a dedup. Not folded. |
| Builder symmetry ergonomics | 0.6.2 (unreleased) | `CouplingMatrixBuilder::build()` now rejects a one-sided off-diagonal `set` by default (typed `PreconditionViolation`), since every downstream consumer assumes a symmetric matrix; `allow_asymmetric(true)` opts out (used by the deliberately asymmetric folded-mask gate fixtures). 3 new lib tests pin the invariant. |
