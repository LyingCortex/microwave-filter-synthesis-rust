# MFS TODO List

**Status as of 0.6.3-unreleased (2026-10-13).** This file is maintained as the single
open-item list. Closed items are listed with the version that closed them so
history is traceable. Authoritative architecture/limitation references:
`docs/architecture.md` (Known Limitations / Future Roadmap),
`docs/algorithms.md` (implementation-grounded algorithm reference), and
`docs/internal/architecture-review-2026-09-25.md`.

## Current state (0.6.3-unreleased, 2026-10-13)

`master` head `7dd617a`; the working tree carries the round-2 and round-3
review remediation (uncommitted). Round-2 report:
`docs/reviews/mfs-review-round2.md`; round-3 report:
`docs/reviews/mfs-review-round3.md`. All CI-parity checks green on this tree
(local, WSL2 / Rust 1.95.0):

- `cargo test` — **457 passed, 0 failed** (default features; = lib 308 +
  integration 132 + doctests 17; does not build examples/benches as test
  targets)
- `cargo test --features cli,python` — **473 passed, 0 failed** (adds the
  `mfs_cli` binary + `cli_integration` / `cli_save_resume_e2e` suites, +16)
- `cargo test --all-targets` — **441 passed, 0 failed** (= lib 308 +
  integration 132 + 1 example `#[test]` in
  `examples/high_order_stability.rs`; runs benches/example unit tests but
  **not** doctests). The 457 vs 441 gap: default runs 17 doctests and no
  example test; all-targets runs the 1 example test and no doctests.
- `cargo clippy --all-targets -- -D warnings` — zero warnings
- `cargo clippy --features cli,python --all-targets -- -D warnings` — zero
  warnings

Since the 0.6.1 snapshot: the 0.6.1 whole-crate review wave (3 P1 + 14 P2 +
11 P3) is **fully closed and gate-anchored** (released as 0.6.2, 2026-10-11:
fast-path group delay + three-algorithm verification, §D.7 dB single source,
CI wheel/FFI gate). On top of 0.6.2, in 0.6.3-unreleased: `flip_sign` /
`flip_sign_inplace` deduplicated into one shared `pub(crate)` helper
(`matrix/mod.rs`), bit-exact A/B verified (9 prototype cases, orders 4–10,
byte-identical folded/CT/CQ matrices + 401-point response); the P2-1
`PoleExpansionData::from_matrix` transversality guard (structural
precondition + LU fallback, 5 regression tests); and the P1-1 all-finite
feasibility domain guard in `find_eps` (typed `PreconditionViolation` for
`eps <= 1`, 7 tests), with the `eps > 1` domain documented in
`docs/algorithms.md` §2.4 and `docs/architecture.md`. Round-3 (P2-2
port-diagonal guard + 3 tests, complex S21+S11 fast-path gate + 5 tests,
measured high-order boundary 28/24 + example regression test, `synthesize`
doc + `verify_against_spec` doctest, with-zeros group-delay FD anchor test)
is fully closed and test-anchored — see the 0.6.3-unreleased CHANGELOG entry.
The Open sections below hold the remaining low-priority ecosystem items and
the refuted passband-TZ work.

## Open

### High priority (core functionality)

- **Passband transmission zeros** (`|zero| < 1`) — still rejected at the
  `approx` entry point (typed `MfsError::InvalidTransmissionZero`). The 0.6.0
  E1 investigation refuted the in-band-null prototype in the single-function
  form (Feldtkeller violation); realizing them needs the multifunction
  characteristic, which is genuine new approximation research, not a helper
  relaxation.

- **Lossy (Q > 0) pole-expansion response path** — `design.response()` /
  `response_normalized()` / `response_bandpass()` have the lossless
  `PoleExpansionData` fast path (transversal partial fraction, O(N) per
  point) with a 5-point LU verification fallback
  (`src/design.rs:419-448`). `response_lossy()` / `response_lossy_normalized()`
  go straight to `ResponseSolver` per-point LU (`src/design.rs:636/58`) —
  no fast path. The lossy math supports the same expansion with a one-line
  change: finite-Q loss is a *uniform* resonator-diagonal shift — the
  prototype-domain response matrix is `A(ω) = M + ω·I_r − jR − (j/Qu)·I_r`
  (`src/response/backend.rs:129-138`, `dissipation = 1/Qu` added identically
  to every resonator diagonal, real-ω convention, `ω` carried as
  `Complex64::new(omega, 0.0)` at `pole_expansion.rs:177`) — so every pole
  shifts from `−B_k` to `−B_k + j/Qu` while the residues are unchanged
  (they are pure coupling products `M_Sk² / M_Sk·M_kL / M_kL²`,
  `src/response/pole_expansion.rs:96-107`). Concretely: add a
  `q: Option<f64>` (or `shift: f64`) parameter to
  `PoleExpansionData::from_matrix` and `evaluate_at` so the pole offset
  becomes `−B_k + j/Qu`; `evaluate_at` sums `r_k / (ω − (−B_k + j/Qu))`.
  The existing 5-point LU verification gate is reusable unchanged
  (evaluate both sides with the lossy `ResponseSettings`). This matters
  because the GUI's Python binding (`fuqtgui:
  src/domain/synthesis_engine.py:246-249`) takes `response_lossy` whenever
  `q > 0` (GUI default `Q = 1000`), so the lossy path is the one the
  default GUI actually exercises. Scale: ~30–60 lines in Rust + tests.
  Measured on the same machine: lossless 201 points 0.043–0.085 ms (fast
  path) vs lossless-LU 0.23–0.55 ms; the lossy 1801-point sweep was ~13 ms
  in the fuqtgui serve path (`fuqtgui: docs/roadmap/maintenance-backlog.md`
  §9.0 table), which is where this optimization pays off.
  (Logged 2026-10-10. Math **verified** same day, throwaway script — no
  production code touched: numpy re-derivation of the 2×2 Schur
  partial-fraction `evaluate_at` formula, lossless and
  `−B_k + j/Qu`-shifted, against the real `mfs` crate's own
  `response` / `response_lossy` outputs (the exact slow path this fast
  path would replace): prototype-domain golden cases 4-4 / 6-3 / 8-3,
  Q=1000, 201 points → max |ΔS21| ≤ 6.5e-13 dB, |ΔS11| ≤ 1.3e-11 dB
  (machine precision, incl. the all-pole 4-4 case with complex M_SL).
  Physical-domain BPF (order 8, 3.5 GHz / 200 MHz, 401 Hz-grid points,
  zeros mapped via the standard L→BP transform, physical Q = 1000) →
  1.2e-8 / 6.1e-9 dB (LU-vs-partial-fraction floating-point rounding only).
  **Pitfall found during verification:** the bandpass physical domain
  normalizes the physical Q at response time,
  `Qu_norm = Q·bw/center` (`src/design.rs:629`) — the lossy fast path's
  pole offset must use `Qu_norm`, not the raw physical Q (using Q=1000
  directly where Qu≈57.1 is correct produces 24–40 dB errors at the
  transmission zeros). Same-basis numpy magnitude: per-point LU (n=10)
  3–5 ms vs partial-fraction 0.17–0.22 ms over 401 points (~20×).)

- **Pole-expansion fast path vs per-point LU — scope, not a second solver.**
  Two distinct things; don't conflate them:
  - **Fast path is the sole source of the Rust speedup.** 1801-pt sweep:
    ~0.05 ms (fast) vs ~13 ms (per-point LU) vs ~100 ms (Python skrf). The
    bench `response` stage (`benches/synthesis.rs:187`, the 12×–39× headline
    over Python) times per-point LU, *not* the fast path — the Rust response
    advantage only materializes when the fast path is taken; forcing LU
    shrinks it to a language-only gap.
  - **LU is the correctness fallback + the fast path's verification anchor.**
    Only path that covers non-transversal (folded/arrow, complex coupling)
    and lossy matrices, and the ground truth the fast path verifies against
    on a 5-point grid before being trusted (`src/design.rs:419-448`).
  - **Fast-path group delay — implemented (2026-10-11).**
    `PoleExpansionData::evaluate_at` now returns `(s11, s21, group_delay)` and
    emits a finite group delay at O(N)/point: group delay = `−Im(d/dω ln S21)`,
    and S21 is a known rational function of ω in the Schur 2×2 form
    (`src/response/pole_expansion.rs:150-254`), so its log-derivative is a
    per-pole partial-fraction sum (one extra O(N) sum alongside the existing
    S-parameter sum). It is `None` exactly at a transmission zero
    (`|S21| → 0`, `inv_10` degenerate), mirroring the LU backend's
    `denominator.norm_sqr() <= 1e-18` guard — explicit undefined, not a fake 0.
    Both the normalized and bandpass fast-path callers
    (`evaluate_normalized` / `evaluate_bandpass`) forward it, so
    `design.response()` / `response_bandpass()` now deliver engine-native
    group delay. **Verification (three independent algorithms agree):**
    (1) fast vs LU on a 3-zero filter, full band, to 1e-7 (the cross-algorithm
    guarantee — `group_delay_matches_across_fast_and_lu_backends`, Check 1);
    (2) fast vs an all-pole-ladder central difference of the S21 phase, the
    value anchor (Check 2); (3) cross-language, the new mfs fast-path value
    converted to seconds via the exact BP chain rule
    `dω/df = (1/FB)(1/fc + fc/f²)` matches **both** skrf `group_delay` (the
    same third-party engine the GUI uses for its group-delay plots) and an
    independent central difference of the S21 phase to **16.6 ps** over 268
    points of a 2 GHz / 200 MHz band. The 16.6 ps is the 2nd-order
    central-difference's own O(h²) truncation (Δf = 0.9 MHz) — skrf uses the
    identical `−dφ/dω` central difference on the same S21, so it agrees with
    the FD reference to machine precision (0.000 ps) and sits 16.6 ps from the
    analytic fast path on both sides. Measured O(N): 1801-pt sweep
    **2.0 µs/pt** (prototype) / **2.2 µs/pt** (bandpass) with the derivative
    active, release build.
    **Check 3 is reproducible** from the consumer repo at
    `fuqtgui/src/scripts/verify_fastpath_gd.py` (needs mfs ≥ 0.6.2 + skrf in
    the venv; run `uv run python src/scripts/verify_fastpath_gd.py`, exit 0 on
    pass). Checks 1–2 are the permanent CI legs (Rust test, no third-party dep).
  - **Mistake recorded (2026-10-11) — 2π units slip in the *verification* FD
    reference, misattributed to "FD truncation near the transmission zeros".**
    Root cause: the cross-check built a finite-difference reference as
    `−d(phase)/df` (per **hertz**), but group delay is `−dφ/dω` with
    `ω = 2πf` — the reference was `2π`× too large, so `fast−FD` and
    `skrf−FD` both read `(2π−1)·gd ≈ 85 ns`. I first (wrongly) wrote that up as
    a genuine "85 ns FD truncation near the zeros, not a fast-path error" and
    shipped it in two docs. Why it looked plausible: the offset was a smooth
    function of |S21| and a *coarser* central difference, so the reflex
    diagnosis was "truncation + steep phase." **How to catch it next time** —
    the two tell-tale symptoms that rule out truncation and point at a
    units/convention slip instead: (1) the error *does not shrink* when the
    grid is refined (truncation would collapse ~16× at 16× finer Δf); (2) it
    *grows with |S21|* toward the passband, opposite of where a zero-driven
    steep-phase/ill-conditioning error would peak. A `2π`-family offset that
    scales as `(2π−1)·max(gd)` is the signature of a per-hertz-vs-per-radian
    (or missing `1/(2π)`) slip, not an algorithm bug. **Fix:** differentiate
    w.r.t. `2πf` (i.e. `central_diff(2πf, S21)`), which makes the FD a genuine
    ~20 ps independent anchor and skrf-vs-FD drop to machine precision —
    confirming the fast path was correct all along. The units + conversion +
    the two symptoms are now pinned on `ResponseSample::group_delay`
    (`src/response/mod.rs`) and `response_full` (`src/python.rs`) so the next
    reader hits the convention before reaching for a per-hertz derivative.


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
| Builder symmetry ergonomics | 0.6.2 (unreleased) | The symmetry invariant now lives in the primitive `CouplingMatrix::new_with_topology` (tolerance-based, scale-aware: `1e-12*max(1,|a|,|b|)` per pair) rather than `CouplingMatrixBuilder::build()`, so every construction path — including the CLI's user-supplied JSON matrix (`parse_matrix_json`) — is guarded, not just the builder. The builder routes to the strict constructor by default and to the explicit `new_with_topology_asymmetric` bypass via `allow_asymmetric(true)` (the folded-mask fixtures). 9 tests pin it (6 constructor — incl. two boundary tests that pin the tolerance from both sides, relative ~1e-13 accepted / ~1e-6 rejected — plus the 3 builder tests from the prior commit); a CLI unit test + e2e pin the user-input rejection. |
