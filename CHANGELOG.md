# Changelog

## [Unreleased]

## [0.6.0] — 2026-10-05

### Added

- **Multi-band synthesis design doc (C0) — route committed.**
  `docs/internal/multi-band-synthesis-design.md` fixes the multi-passband
  synthesis plan (`docs/continuation-plan.md` §C, scope-1): **route-α**
  (multi-band Chebyshev via frequency transformation), anchored on
  Macchiarella–Tamiazzo, *Design techniques for dual-passband filters*,
  IEEE TMTT 53(11) 2005 (DOI `10.1109/TMTT.2005.855749`); route-β deferred
  with rationale. The doc specifies the characteristic function, the
  `MultiBandSpec`/`BandSpec`/`RejectionSpec` surface, the transmission-zero
  placement law (outer zeros stopband; **inter-band zeros also stopband** —
  corrected 2026-10-05, see the Fixed entry below), the target $N+2$
  transversal topology (inter-band TZs become trisection/extracted-zero
  sections), and the first acceptance fixture
  `Cameron2011_10_5_dual_band_symmetric` (Cameron 2011 Fig. 25), already
  in-repo. (C2/C3 subsequently shipped 2026-10-05 — see below.)

- **True multi-passband (multi-band) synthesis — C2, route-α (2026-10-05).**
  A dual-band filter whose response is equiripple in **two** passbands with
  inter-band and outer stopband rejection — one structure, not two
  single-band designs wired together. Per the C0 design doc, route-α
  (multi-band Chebyshev via frequency transformation, Macchiarella–Tamiazzo
  2005, DOI `10.1109/TMTT.2005.855749`):
  - **C2b — approximation** (`src/approx/multi_band.rs`):
    `multi_band_chebyshev_polynomials(bands, order, rl, tz)` → `PolynomialSet`.
    The positive-axis reflection zeros are solved by a damped
    Levenberg–Marquardt **equiripple** step (a square system: the `n_rz + 1`
    local maxima of the characteristic function `|F/P|` in one positive
    passband interval all share the requested ripple level), then the
    Feldtkeller `E(s)` spectral factor and `PolynomialSet` assembly **reuse the
    single-band helpers verbatim** — no new machinery.
  - **C2a — spec surface** (`src/spec/multi_band.rs`): `MultiBandSpec` /
    `BandSpec` / `RejectionSpec`, a separate type from `FilterSpec` (a K-band
    filter is not a K-channel multiplexer), with `MultiBandMapping`'s
    invariants (≥ 2 disjoint center-sorted bands).
  - **C2c — synthesis** reuses `MatrixSynthesisEngine::synthesize` with
    **no changes** (the output is an ordinary `PolynomialSet`).
  - **Dual-band only (K = 2):** the equiripple solve is *square* only for a
    single positive-axis passband interval (true for a symmetric design like
    `[-1,-0.35]` ∪ `[0.35,1]`). K ≥ 3 — a center band spanning the origin, or
    two or more positive-axis passbands — **overdetermines the real-rooted
    characteristic by M − 1** (M = number of positive-axis intervals) and has
    no equiripple solution; this is a fundamental degrees-of-freedom limit of
    the single-function characteristic, **measured** 2026-10-05 (M=1 control
    residual 1.7e-9 vs. 5.7e+2 / 4.2e-2 for the K=3 cases, zeros pinned to the
    band edges). The solve rejects it with a typed `MfsError::NotImplemented`
    stating the reason; a genuine K ≥ 3 equiripple needs a complex-rooted
    (Zolotarev) or multifunction characteristic (out of scope). See
    `docs/internal/multi-band-synthesis-design.md` §"K ≥ 3: the DOF Wall".

- **Multi-band verification surface — C3 (2026-10-05).** `tests/multi_band.rs`
  (`dual_band_response_reproduces_cameron_fig25_fixture`) anchors the design to
  the in-repo literature fixture `Cameron2011_10_5_dual_band_symmetric` (Cameron
  2011 Fig. 25: order 10, 22 dB RL, dual-band `±[0.35,1]`, 5 TZs) on a
  4001-point grid over `-2..2`. The synthesized **transversal** matrix and the
  fixture's printed `trisection_cascade_with_extracted_zero` matrix are
  *different topologies of the same network*, so "reproduces the fixture" is
  asserted at the **response** level (topologically invariant), plus absolute
  spec-level checks: (1) response invariance — max|ΔS11| 4.8e-4, max|ΔS21|
  1.7e-4 (both < 1e-3); (2) per-band equiripple — worst in-band |S11| = −22.00
  dB (the RL spec, ±1); (3) unitarity — max‖|S11|²+|S21|²−1‖ = 2.2e-15 (< 1e-10);
  (4) transmission zeros — |S21| = −261.8 / −60.4 / −70.5 dB at w = 0, ±0.2542,
  ±1.1573 (all < −40); (5) outer stopband — worst |S21| = −16.6 dB over
  |w|∈[1.10, 2.0] (< −15). `examples/multi_band_response.rs` renders the
  response (per-band / gap / outer table + a Touchstone .s2p).

- **`benches/optimizer.rs` — a Criterion bench for the LM optimizer (D4-opt).**
  Times `optimize_matrix` (`tune` scenario: order-4/6/8 folded, 4 inter-resonator
  couplings perturbed 7 %, re-optimized from the perturbed good guess) and
  `extract_matrix` (`extract` scenario: order-4/6, cold inline-chain guess),
  printing the iteration count, final cost, and physical quality next to each
  timing (each run 2.8–17 ms warm). **Deliberately outside the cross-language
  benchmark standard** — `gates.json` defines no `optimizer` family, so case ids
  and tolerances would have no provenance (same stance as `synthesis.rs`'s
  program-generated sweep). Order-8 `extract` is not timed: from the cold chain
  guess it converges to a local minimum far from the target (measured: cost 39.2,
  max|ΔS| ≈ 1.0, 435 iterations / ~85 s), so a timing would measure noise.

- **`tests/optimizer_regression.rs` — extraction from a noisy response
  (D4-opt).** Pins the optimizer's *extraction* contract beyond the well-posed
  cases: order-6, deterministic 0.1 dB magnitude noise (LCG). Asserts
  convergence, fidelity to the underlying filter (max|ΔS| < 5e-3; measured
  1.22e-3 against a 1.0e-2 noise floor), and that the fit does **not** chase the
  noise (residual against the noisy data < 4× the noise floor).

- **Optimizer is now reachable from the CLI and Python (D3-opt).**
  `mfs_cli` gains two subcommands, both thin wrappers over the shared
  `optimize_matrix` core:
  - `optimize` — input JSON `{order, topology, data, samples}` (the matrix as a
    flat row-major `data`, the target response as `samples` with
    `frequency_hz`/`s11_re`/`s11_im`/`s21_re`/`s21_im`; a Touchstone `.s2p`
    target is also accepted via `--target`). Knobs: `--order`,
    `--max-iterations`, `--tolerance`, `--lambda`, `--q`, `--output`.
  - `extract` — input JSON `{samples:[...]}` (or `--target` `.s2p`) + `--order`,
    calls `extract_matrix` directly; `--output` optional.
  Both print `{order, topology, matrix, data, cost, iterations, converged,
  recovered:{samples}}` (or a `--format table` view).
  `src/python.rs` gains a module-level `extract(freqs, s21_re, s21_im,
  s11_re, s11_im, order) -> {matrix, cost, iterations, converged}` (a
  `PyFilterDesign.tune` already wrapped `optimize_matrix` since D1-opt).
  Two end-to-end black-box regression tests in `tests/cli_integration.rs`
  (gated `#[cfg(feature = "cli")]`) drive the real binary and assert physical
  quantities: `test_optimize_subcommand_recovers_perturbed_matrix` (order-4
  folded, 3 couplings ×1.05 → recovered to < 1e-4 relative) and
  `test_extract_subcommand_recovers_response` (order-4 transversal → max|ΔS|
  < 1e-3). The bin also gains a `#[cfg(test)] optimizer_tests` module covering
  the flat/nested matrix-JSON and response-JSON parsers.

- **Benchmarks are now driven by the cross-language benchmark standard.**
  `benches/synthesis.rs` no longer hardcodes its case matrix: cases come from
  `filter-bench/data/cases/twoport.json`, timing policy from `harness.json`, and
  stage names from `gates.json`, all through the vendored copy in `bench-data/`
  (`tests/bench_support/`). A program-generated order sweep is kept, and
  documented as deliberately outside the standard (no provenance, no invariants).
- **`benches/multiplexer.rs`** benchmarks the standard's `multiplexer` cases that
  mfs can actually run (tag `type1`/`synthetic`, non-resonant junction), reporting
  the `poly` (per-channel iterative synthesis) and `response` stages. The
  contract's `cm`/`star` stages have no mfs counterpart and are omitted rather
  than fabricated; every gate requiring them is `SKIP`.
- **`tests/benchmark_standard.rs`** enforces the standard from mfs's side:
  contract self-check, vendored-data drift guard (byte-for-byte against the
  sibling repo when present, `MANIFEST.sha256` otherwise), the `folded_mask` gate
  on every `twoport` case, and publication-anchored response checks that replay
  the book's printed Cameron matrices through mfs's response engine.
- **`tests/folded_response_invariance.rs` consumes the standard** instead of
  hand-copied constants: the six Cameron cases are now looked up by global id
  (`twoport_*`), so a change to the standard cannot silently desynchronize them.
- `MultiplexerSynthesisResult::max_iterations` — the iteration budget in force,
  so a caller holding only the result can distinguish convergence from an
  exhausted budget.
- **`MultiBandMapping`** — a piecewise multi-passband frequency mapping in
  `freq` (one band-pass prototype per passband, geometric-mean center
  convention, bands validated disjoint and center-sorted). The forward map
  assigns each physical frequency to the nearest band by prototype-transform
  distance; the inverse is band-local. It does not implement
  `FrequencyMapping` because the trait's inverse has no band-index argument.
- **`HighPassMapping` and `BandStopMapping`** — `freq` now covers all four
  prototype families. Both follow the cross-language benchmark standard's
  sign conventions (`ω = -fc/f`, `ω = -1/Ω_BPF(f)`), clamp to
  `±SINGULAR_CLAMP` at their singular points (`f = 0` for high-pass,
  `f = f0` for band-stop), and implement exact forward/inverse round trips —
  band-stop edges hit ±1 exactly under the geometric-mean convention and the
  inverse is defined for both signs of `ω`.
- **N=8 regression entry in the documented order 2–8 fixture band**
  (continuation plan E4). `tests/benchmark_standard.rs` gains
  `mfs_synthesized_folded_order8_reproduces_the_printed_book_matrix`:
  the standard's only order-8 case (`twoport_8_3`, Cameron 2018 Sec 9.6 /
  Fig 9.22a) is synthesized and folded, and the folded matrix + response are
  pinned against the printed `EIGHT_THREE_FOLDED` anchor at the 4-decimal
  printing quantum.
- **`tests/optimizer_integration.rs` — the optimizer's RED→GREEN gate**
  (continuation plan D1-opt). Two order-6 integration tests, both
  `#[ignore]`d with their current failure documented: (a) synthesize →
  perturb four inter-resonator couplings by 7 % → re-optimize from the
  perturbed good-guess matrix and assert 1e-4 relative coupling recovery
  (currently stops at ~8.3e-3); (b) `extract_matrix` from its crude uniform
  guess (currently wanders into a singular response solve). They pin the
  fixed-δ finite-difference Jacobian's scale sensitivity — the defect
  D2-opt (scale-aware step) is meant to fix — and will be un-`#[ignore]`d
  when that lands.
- **`DuplexerSynthesisResult::diplexer_response`** — junction-level 3-port
  response of the full Macchiarella–Tamiazzo path (continuation plan D1).
  When `full_iteration(true)` converges, the result captures the common-domain
  diplexer polynomials (N, D, Pt, Pr) and the p₀ scalings; this method
  evaluates the paper's physical 3-port model `S11 = N/D`,
  `S21 = √(|p0t|²)·Pt/D`, `S31 = √(|p0r|²)·Pr/D` at any frequency sweep.
  It is `None` for the simplified path or a non-converged run. This is the
  observable that carries the junction's TX/RX power split; a recovered
  channel taken alone (the standalone `tx`/`rx` matrices) carries only its
  share (`p0t²` / `p0r²`) and cannot reach 0 dB.
- **Duplexer full-path D1 gate family** (continuation plan D1/D3), now
  **GREEN** and no longer `#[ignore]`d:
  `full_iteration_order{2,3,4,5,6,8}_quasi_equiripple` plus the
  Section-V printed-design golden benchmark
  `full_iteration_order7_section5_golden`. Each asserts the junction-level
  diplexer response: a captured fixed point, both channel-port passbands
  within 0.5 dB of unity, and quasi-equiripple antenna reflection at the four
  band edges (the 22/20 dB RL spec is pinned there by eq. 13, not at the band
  centers, where odd orders place a Chebyshev reflection zero). The golden
  test additionally matches the paper's printed order-7 N/D/Pt/Pr/p0
  coefficients to print precision.
- **Type-II (resonating-node) junction through the full spectral path**
  (continuation plan D4), with two gates: `full_iteration_type2_order4_quasi_equiripple`
  (symmetric order 4, same quasi-equiripple assertions as the Type-I family)
  and `full_iteration_type2_section6_anchor` (the paper's §VI GSM design: RX
  10 / TX 9 poles, 22 dB, transmission zeros, real reflection zero 1.5 ⟹
  `c0 = 2/3`, degree-40 R(s)). The §VI anchor is *response-based*, not
  coefficient-based, because §VI prints coupling matrices (Tables I/II)
  rather than the N/D/Pt/Pr polynomials: it asserts 0 dB channel passbands,
  quasi-equiripple antenna |S11| (local peaks within ~2 dB of the 22 dB spec,
  peak-to-peak ripple < 2.5 dB — the paper reports < 1.5 dB, Fig. 9), and the
  lossless para-identity `D·D̃ = R` on the passband axis. The paper's eq. (11)
  is implemented exactly: N/D degree `n_tx+n_rx+1` with the extra pure-real
  reflection zero at `s = 1/c0`, and `S_TX·S_RX = (c0/2)(D − N)` with D
  normalized exactly monic so the leading terms cancel to machine precision.

- **True (N+1)-port multiplexer response (B1).** `Multiplexer::response`
  now composes the full complex (N+1)×(N+1) S-matrix at the common node:
  each channel filter is a lossless reciprocal 2-port whose input port sits
  at the node and whose output port is the external channel port. For
  lossless channels the result is unitary and reciprocal, and N=1 reduces
  exactly to the channel's own 2-port S-matrix (asserted by gate tests).
  New API: `Multiplexer::response_matrix` (per-frequency complex S-matrix)
  and `Multiplexer::response_with(.., model)`; the new
  `ResponseModel::{Admittance, GeometricMean}` names the two models. The
  legacy geometric-mean magnitude model is retained as
  `ResponseModel::GeometricMean` for one release cycle — it approximates
  the true model in the deep stopband but underestimates the common-port
  reflection by up to ~150 dB inside a channel's passband.
  `MultiplexerResponse` field semantics are unchanged in shape; the dB
  values are now the true (N+1)-port S-magnitudes
  (`isolation[i][j]` = |Sᵢ₊₁,ⱼ₊₁|, channel j driven).
- **Shared Macchiarella–Tamiazzo full spectral path for N channels (B2).**
  The 2-channel full-spectral iteration is now in
  `src/synthesis/machiarella.rs` and is called by both
  `DuplexerSynthesis::full_iteration(true)` and the new
  `MultiplexerSynthesis::full_iteration(true)`, so the N = 2 and N ≥ 3
  multiplexer/diplexer full paths run the same loop. The N ≥ 3
  generalization implements the N×N p₀² system (`nalgebra` LU with partial
  pivoting), per-channel outer-band return-loss edge selection, and
  contiguous N-block root assignment; the N = 2 branch keeps the
  historical Cramer's-rule 2×2 and ±j edge selection with
  bit-identical operation order (the 10 gate-anchored duplexer tests,
  including the paper §V/§VI printed-design goldens, pass unchanged).
  `MultiplexerSynthesisResult` gains `junction_common` (common-domain N, D,
  Ptᵢ, p0ᵢ², band constants) when the full path converges;
  `mfs::synthesis::junction_response` evaluates the (N+1)-port junction
  response `S_ant = N/D`, `Sᵢ = √p0ᵢ·Ptᵢ/D`. For N > 2 only the ideal
  Type-I junction (n = 1, b0 = 0) and Type-II are supported — a
  non-ideal Type-I is rejected with `PreconditionViolation` (no paper
  basis for that generalization). New gates:
  `multiplexer_full_iteration_triplexer_converges` /
  `_quadruplexer_converges` (convergence, 0 dB per-channel passbands,
  antenna RL at the 22 dB spec at each outer band edge, para-identity
  D·D̃ = R) and `_rejects_nonideal_type_i_for_n_gt_2`.
- **Convergence + fidelity gates for the full path; vendored convergence gate
  flips FAIL→PASS (B3).** `MultiplexerSynthesisResult` gains a `converged`
  flag — `true` when the iteration stops on its own convergence criterion
  (strict root-displacement tolerance, or the high-degree stagnation accept
  for the full spectral path) rather than exhausting the budget. New
  `tests/multiplexer.rs` gates assert, on the full path:
  `multiplexer_full_iteration_well_separated_matches_machine_tol` (far-apart
  channels reach `convergence_error < 1e-6`, in fact ~0 — the contrast against
  the simplified path's documented bounded residual),
  `_adjacent_bounded_loading` (overlapping channels show bounded loading with
  a quasi-equiripple common port), and `_vendored_zhao_n4_converges` (the
  vendored cross-language standard's N=4 case, with transmission zeros). Each
  asserts convergence, 0 dB per-channel passbands, the antenna RL pinned
  within 2 dB of the per-channel spec at each outer band edge, and the
  losslessness para-identity D·D̃ = R (regime-aware bound: strict → 1e-4,
  high-degree stagnation → 5e-2, the duplexer order-7/8 pattern). The
  `benches/multiplexer.rs` `poly` stage now runs the full spectral path and
  `convergence_gate` reads the `converged` flag (an exhausted budget is a
  non-convergence whatever the residual), so the vendored-standard
  `convergence` gate is **PASS** on every runnable case: `multiplexer_a`
  (Zhao N=4, 9 iters, conv_err 1.7e-4 via stagnation) and
  `multiplexer_diplexer` (12 iters, 2.6e-9 strict); `type2`/`contiguous`
  cases remain capability-SKIP.
- **Bench-standard `cm`/`star` gates reported per gate (B4).** The blanket
  `SKIP` for the multiplexer's `cm`/`star`-requiring gates is replaced with
  per-gate, reason-specific reporting. `unitarity` is now **PASS** (measured,
  not skipped): mfs computes the full (N+1)×(N+1) common-junction S-matrix
  directly at the node, so the contract's star-junction power-conservation
  invariant `max|1 − |S11|² − Σᵢ|S1i|²|` is evaluated on that response
  (`gates.rs::multiplexer_unitarity_residual` + `unitarity_tol`); measured at
  machine precision (3.8e-15 / 1.3e-15 on the two runnable standard cases,
  tol 1e-8). The runnable cases' aggregate status moves from SKIP to **PASS**.
  `matrix_size`, `equiripple`, `s11_match`, and `rz_distribution` remain
  `SKIP`, each with the specific quantity mfs's algorithm does not produce
  (the single multiport star-junction CM, the per-channel determinants |Δᵢ|,
  the multiport rational-function S11, and the rational-function pole/zero
  set, respectively) — the quasi-equiripple antenna RL that `equiripple`
  would assert is instead covered by the B2/B3 full-path fidelity gates.


- **Repeated-pole (Laurent) residue expansion (E2, 2026-10-05).** The
  residue expansion now supports repeated poles, closing the plan's last
  open Wave-5 numeric enabler (needed by B3's fidelity work for degenerate
  root sets):
  - `ResiduePole` (`src/synthesis/residues/expansion.rs`) gains
    `laurent: Option<Vec<ComplexCoefficient>>` — the principal part
    `b_0..b_{m-1}` of the repeated pole, with `residue` = `b_0` (the
    leading, most-singular coefficient). Simple poles keep `laurent: None`.
  - `residue_expansion` clusters coalesced roots (single-linkage at 1e-5;
    the solver splits a multiplicity-m root into m roots within ~1e-7–1e-5,
    while distinct filter poles are ≥ ~1.6e-3 apart), and takes the
    **multiplicity `m` as the cluster size** (not re-derived from the
    pole-shifted coefficient vector, which is ill-conditioned for high-degree
    denominators — a scale-relative "vanishing" threshold miscounts small
    genuine `q_1` as zero and inflated `m` on order-16 simple poles). The
    principal part is computed by power-series division of the shifted
    numerator by `x^m·R(x)`, then leading near-zero coefficients are trimmed
    (a numerator–denominator cancellation lowers the effective pole order).
  - `PoleExpansionData::evaluate_at` (`src/response/pole_expansion.rs`) sums
    the full principal part per repeated pole in the fast path, and each Y
    parameter is summed with **its own** principal part — a shared pole's
    effective order can differ between the three after cancellation.
  - `build_transversal_from_residues` rejects a repeated pole with a typed
    `MfsError::NotImplemented`: a transversal (inline-chain) matrix has one
    simple pole per resonator and cannot represent `1/(s-p)^m`. (This was
    already in place from B3; the Laurent data now flows to it.)
  - `tests/synthesis_properties.rs` `property_2` round-trip now reconstructs
    the Laurent form, so the data-preserving gate covers repeated poles.
  Valid filter specs have distinct poles, so the Laurent path is exercised
  only by degenerate/synthetic inputs; the new unit tests
  (`double_pole_laurent_principal_part` with hand-computed coefficients,
  `repeated_pole_rejected_with_not_implemented`,
  `laurent_principal_part_contributes_to_sum`) pin the exact coefficients,
  the typed rejection, and the fast-path summation.


### Fixed

- **Multi-band `K ≥ 3` equiripple — refuted as infeasible for the
  single-function characteristic (measured, 2026-10-05).** The C2
  multi-band approximation (`src/approx/multi_band.rs`) was documented as
  "dual-band only, K ≥ 3 a Low roadmap item". Measuring that extension
  (a throwaway probe generalizing the equiripple solve to `M` disjoint
  positive-axis intervals) showed it is **not** a bounded feature: the
  characteristic `K(w) = F(w)/P(w)` is real-rooted, so equiripple on `M`
  positive-axis passband intervals imposes `N/2 + M − 1` independent
  constraints on the `N/2` reflection zeros — square only for `M = 1`
  (the symmetric dual-band). `M = 2` / `M = 3` (a center band spanning the
  origin, or 2+ positive-axis bands) are overdetermined by `M − 1`; the
  probe measured residuals `5.7e+2` / `4.2e-2` (zeros pinned to band edges)
  versus the `M = 1` control's `1.7e-9`. The rejection is now a typed
  `MfsError::NotImplemented` whose message states the degrees-of-freedom
  reason instead of "currently requires a single contiguous interval"; the
  design doc gains a "K ≥ 3: the DOF Wall" section with the proof and the
  measurement. A genuine `K ≥ 3` equiripple design needs a complex-rooted
  (Zolotarev) or multifunction characteristic — out of scope for this route.

- **`PoleExpansionData::evaluate_at` — per-Y principal-part length (review fix, 2026-10-05).**
  The fast path summed all three Y parameters with a single shared length
  (`m = l11.len()`, then indexing `l12[k]`/`l22[k]` over it), assuming the
  three principal parts are always equal-length. They are not: each Y
  parameter is a distinct rational function, and `trim_leading_zeros`
  (numerator–denominator cancellation) can shorten one of them
  independently. With y11 shortest, extra y12/y22 terms were silently
  dropped; with y11 longest, the indexing panicked. Unreachable through
  `from_matrix` (all-`None`) and by the hand-built unit tests (equal
  lengths), but reachable through `from_polynomials` on a degenerate root
  set with per-parameter cancellation — exactly the E2 input class. Fixed
  by summing each Y parameter with its own principal part via a
  `laurent_sum` helper; new regression test
  `repeated_pole_with_cancellation_uses_each_y_own_order` pins a
  hand-computed case where y11/y12 keep the double-pole part and y22
  collapses to a simple residue (the old code silently dropped the y11/y12
  principal parts there).

- **`PoleExpansionData::evaluate_at` — Laurent principal-part power indexing.**
  The fast path summed a repeated pole's principal part with the wrong
  powers: `b_k` was multiplied by `1/(s-p)^(k+1)` instead of
  `1/(s-p)^(m-k)`, so the most-singular coefficient `b_0` landed on a
  constant (order-0) term and the result was wrong for every multiplicity ≥ 2
  pole. This bug was latent (repeated poles never occur for valid filter
  specs) and was caught by the new
  `laurent_principal_part_contributes_to_sum` test, which hand-computed
  `s11 = 7/25`, `s21 = -24j/25` and observed the buggy `s11 = 475/493`
  before the fix. Fixed by accumulating the inverse
  powers `1/denom, 1/denom², …` and reversing so index `k` holds
  `1/denom^(m-k)`.

- **E1 (passband transmission zeros) — investigated, refuted, and the
  design-doc premise corrected.** The plan (`docs/continuation-plan.md` §E
  item 1) and the C0 design doc claimed (a) that relaxing `safe_sqrt_term`
  (`src/approx/generalized_ops.rs`) would enable in-band transmission zeros
  and (b) that multi-band route-α's inter-band zeros are passband zeros that
  make E1 a hard prerequisite for C2. Both are wrong:
  - A literal relaxation (imaginary `√(1 − 1/z²)`) produces a complex
    **non-Hermitian** `F(s)` (`F(jω) ∈ ℂ`), so the Feldtkeller identity
    `|E|² = |F|² + |P|²/ε²` fails at O(1) for `|z| < 1` (measured: residual
    ~1e-1…1e1), versus ≤ 1e-13 for `|z| ≥ 1`. The relaxed recurrence is not a
    valid prototype; a physical in-band null needs the *multifunction*
    characteristic, which is new approximation research, not a helper change.
  - A `BandPassMapping`'s `|Ω| = 1` boundary *is* the band's own edge, so
    inter-band gap frequencies map to `|Ω| ≥ 1` on **both** adjacent bands
    (worked example, bands 10–11/12–13 GHz: the 11–12 GHz gap is
    `Ω ∈ [+1.00, +2.83]` on the lower band and `[−3.18, −1.00]` on the upper).
    Route-α's gap zeros are ordinary stopband zeros → **E1 is not a blocker
    for C2**.
  - `generalized_chebyshev_polynomials` now rejects `|zero| ≤ 1` at the
    **entry point** with a typed `MfsError::InvalidTransmissionZero` (was an
    untyped `NotImplemented` leaking from the helper); `safe_sqrt_term` keeps
    its guard, reworded. The spec layer stays permissive.
  - New integration test `tests/passband_tz.rs` pins the typed rejection, the
    unchanged stopband path, and the corrected gap→stopband TZ law.
  - Docs corrected: plan §E item 1, §C status/route-α rationale/table/Wave 5
    line, and the C0 design doc (TZ law, acceptance criteria, risks).

- **The LM optimizer (`src/optimize.rs`) — four correctness fixes; the two
  order-6 integration gates go RED→GREEN (D2-opt).**
  (a) `optimize_matrix` declared convergence by comparing the cost against the
  running *best*, which equals the current cost on the iteration after any
  accepted step — a false "no change" that stopped the loop after 2 iterations
  (~8e-3 relative, far from the 1e-4 gate). Convergence is now judged on the
  *successive-iteration* cost change, and only after an accepted step (a
  rejected step leaves the cost trivially unchanged — that is stagnation, not
  convergence).
  (b) `extract_matrix`'s crude initial guess was star-like — every resonator
  coupled to both ports with the inter-resonator chain couplings frozen at
  zero — a template that cannot represent a transversal (inline-chain) target;
  the LM steps then wandered into singular response solves. The guess is now a
  crude inline chain (source `M[0,1]`, load `M[N,N+1]`, chain `M[k,k+1]`).
  (c) A trial step whose response evaluation failed (a singular matrix)
  propagated the error out of `optimize_matrix` instead of being a
  rejection; it is now λ×2.
  (d) The finite-difference Jacobian step is scale-aware:
  `δ = 1e-7·max(1,|x_i|)` with up to 4 halvings when a perturbed evaluation
  returns `Err` (the plan's D2-opt spec — implemented, but not the primary
  fix; (a)–(c) were).
  `tests/optimizer_integration.rs`: both `#[ignore]`d RED gates are now run and
  green — order-6 folded perturb→recover converges in 6 iterations to
  ~3e-10 relative (gate: 1e-4); order-6 transversal extract-from-crude-guess
  converges in 9 iterations to max|ΔS| ~5e-10 (gate: 1e-3). Known limitation
  (documented): extraction recovers within the guess's topology — a folded or
  arrow target is not representable from the inline-chain start.

- **The D1 full-iteration gate is green at orders 2, 3, 4, 5, 6, and 8.**
  The Tx/Rx diplexers (`full_iteration(true)`, ideal Type-I junction, 22 dB)
  all converge to a physically valid quasi-equiripple diplexer at the
  junction level: 0 dB passbands at both ports, and band-edge antenna return
  loss ≈ −22 dB (within ~2 dB of the spec) at the four passband limits where
  eq. 13 pins the RL. Orders 2/3/4/5/6 converge to the strict 1e-8 tolerance
  (conv 6e-9 … 2.5e-9); orders 7/8 capture their fixed point via the
  stagnation criterion (degree-14/16 root-extraction floor exceeds 1e-8). The
  original failure was the gate asserting on the *standalone* recovered
  channels (each carries only its share of the junction power split, capping
  at −0.6 dB) rather than the diplexer's port response; the gate now asserts
  the junction-level `diplexer_response` (Added). The deep odd-order antenna
  null at a band *center* is the channel's Chebyshev reflection zero — a
  physical feature, not a defect (an earlier note mislabeled it "degenerate").
- **Order-7 full-iteration convergence via a stagnation criterion, plus the
  Section-V printed-design anchor.** At degree 14 the per-iteration
  root-extraction noise floor (~5e-6) exceeds the strict `convergence_tol`
  (1e-8), so the order-7 loop previously exhausted its budget without
  capturing the (physically converged) fixed point. The iteration now accepts
  the solution via a stagnation criterion — accept once the root displacement
  stops improving and is already small (`best_delta < 1e-3` after 3 flat
  iterations, min 5 iters) — the standard remedy for root-iteration, with the
  strict fast path unchanged for orders 2/4/6. New gate
  `full_iteration_order7_section5_golden` is a **golden benchmark against the
  paper's printed Section-V design** (7/7, 20 dB, Type-I `n = 1.47`,
  `b0 = −0.171`): the converged common-domain polynomials are compared
  coefficient-by-coefficient against the paper's printed N, D, Pt, Pr and the
  p0 scalings. Measured agreement (2026-10-03): |p0t|² 1.2e-4, |p0r|² 1.6e-4,
  D 1.2e-2, Pt 9.4e-3, Pr 5.9e-3, N 1.5e-2 — all at the paper's print-rounding
  level (p0 to 4 sig figs, polynomials to 2–3 sig figs), i.e. the code
  reproduces the printed design, not a near-miss. The gate also asserts the
  para-identity `D(s)·D̃(−s) = R(s)` on-axis (the exact losslessness invariant
  of the spectral factorization) and the physical §V response (0 dB
  passbands, band-edge antenna RL within ~2 dB of the 20 dB spec).
- **The vendored standard exposed a folded-convention mismatch in the contract,
  not in mfs.** mfs's even-order `to_folded` leaves the reduction pivot at index
  sum `N+2`, like the book's own printed folded matrices, and outside the
  contract's original canonical-only `folded_mask`. The standard now accepts both
  legal folded conventions (canonical ∪ mirrored), which is what the anchor data
  demanded; mfs's folding is unchanged.

- **Faer root solver now validates its eigenvalues instead of silently
  shipping garbage.** `FaerRootSolver` (the default top-of-chain solver in
  `AdaptiveRootSolver`) can return companion-matrix eigenvalues that are not
  roots of the polynomial (nor of any nearby one) when the eigenvalue
  problem is ill-conditioned — e.g. roots spanning many decades of
  magnitude — without the faer call reporting an error. It now gates its
  output on a backward-residual check
  (`max_backward_residual ≤ 1e-3`, a threshold calibrated against production
  F_s/E_s/P_s polynomials up to order 30, worst case ~6e-6, i.e. ~170×
  margin) and a non-finiteness check, and returns
  `MfsError::NumericalFailure` when it fires, so the adaptive chain falls
  through to Durand-Kerner / Aberth-Ehrlich.
- **`FaerRootSolver` balances its companion matrix before the eigenvalue
  solve** (Parlett–Reinsch diagonal-similarity, LAPACK `xGEbal` style), so the
  backward-residual gate above no longer fires on well-behaved-but-ill-scaled
  inputs. A companion matrix whose coefficients span many decades (the §VI
  type-II diplexer R(s) spans ~19, from 1.6e-18 to 5.6e1) is severely
  unbalanced — one row/column has a 1-norm ~1e-18 while another ~1e2 — which
  makes the raw eigenvalue problem ill-conditioned; the gate measured a
  backward residual of ~0.75 (garbage) on it and forced a fall-through to the
  iterative solvers, which also failed at degree 40. Balancing applies a
  sequence of eigenvalue-preserving similarity transforms (row i scaled by
  `d_i`, column i by `d_i⁻¹`, `d_i = sqrt(col_norm/row_norm)`), bringing the
  §VI R(s) backward residual to ~1e-9 and recovering the degree-40 factor the
  unbalanced path could not. The gate remains as a safety net; two existing
  regression tests that pinned the *pre-balancing* premise (faer must be
  gated on a 21-decade geometric polynomial) now assert the recovered root
  quality instead, and `faer_solver_recovers_ill_conditioned_geometric_polynomial`
  locks the improvement.
- **Durand-Kerner and Aberth-Ehrlich no longer return NaN root sets as
  "converged".** Their in-loop convergence check compared `max_delta`
  (accumulated with `f64::max`, which swallows NaN) against the tolerance,
  so an iteration that overflowed to NaN read as converged with
  `max_delta = 0.0`. Both solvers now also require all roots to be finite
  before the early return.
- **Type-II junction rejects non-positive `c0`** with a typed
  `PreconditionViolation`. A zero/negative resonating-node capacitance is
  physically meaningless and would otherwise silently corrupt the Type-II
  reconstruction (`S_TX·S_RX = (c0/2)(D − N)` → zero or negative product).
  Regression: `type2_junction_rejects_nonpositive_capacitance`.
- **The faer backward-residual gate's accept predicate is unit-tested
  directly** (`faer_gate_accepts_predicate`). After the balancing fix, no
  production input reaches the gate's rejection path, so the wiring that
  turns a non-finite / over-threshold root set into a `NumericalFailure`
  (and the adaptive chain's fall-through on it) needed its own lock.
- **Stale "full path not converged at order ≥ 4" documentation reconciled**
  across `docs/architecture.md` (Known Limitations #4 + roadmap),
  `docs/internal/todo.md`, and `docs/continuation-plan.md` (baseline, known
  limitations, B2).

### Changed (refactor)

- **Duplicated channel-recovery path consolidated.** The duplexer and
  multiplexer each carried a byte-identical copy of the post-iteration
  E/F reconstruction (`recover_filter_design`); both now call the single
  `mfs::synthesis::recover_channel_design` (`src/synthesis/recovery.rs`).
  Behavior unchanged.
- **CLI context construction moved into `pipeline`.** The CLI's local
  `build_context_from_request` (a duplicate of the `run_full_pipeline`
  prologue) is replaced by the exported
  `mfs::pipeline::context_from_request`, now the single canonical
  request→context prologue for both the full pipeline and the CLI's
  incremental `--stage` execution.

## [0.5.0] — 2026-09-25

### Changed (architecture remediation)

- **`PolynomialSet::new` takes a `PolynomialSetParams` struct** (8 positional
  args → 1). New `pub struct PolynomialSetParams { order, ripple_factor, eps,
  eps_r, transmission_zeros_normalized, e, f, p }` (`mfs::approx`). **Breaking**
  (tests-only in-tree callers updated). Same `QuadrupletExtraction`
  parameter struct for the 8/9-arg quadruplet section-extraction APIs.

- **N-channel mutual loading now propagates into output coupling matrices (P1).**
  The duplexer/multiplexer `recover_filter_design` rebuilt E/F from the
  perturbed S-roots but wrote only `PolynomialSet.e/f`; the synthesis engine
  reads `generalized.e_s/f_s`, so the loading correction was silently dropped
  and output matrices were bit-identical to independent per-channel synthesis.
  Both recovery paths now write the perturbed E/F back to `generalized.e_s/f_s`.
  (`src/multiplexer.rs`, `src/duplexer.rs`)
- **`ResponseSample.group_delay` is now `Option<f64>` (P9).** `None` when the
  solver cannot compute it (pole-expansion fast path, degenerate LU case)
  instead of a silent `0.0`. Terminal output renders "n/a"; response
  verification compares group delay only where both samples have it.
  **Breaking.**
- **Pipeline `mapping` now drives physical-domain evaluation (P3).** A
  `bandpass`/`lowpass` mapping maps a physical-Hz grid through the frequency
  mapping before response evaluation; previously it was silently ignored (a
  physical grid was misread as normalized ω, e.g. 6.75 GHz → |S21| ≈ −490 dB).
  Unknown or incomplete mappings are rejected with `MfsError::InvalidFrequency`.
  **Breaking** for requests that relied on the old (ignored) behavior.
- **Legacy crate-root free functions are deprecated (P2).** `filter_spec`,
  `lowpass`, `bandpass`, `normalize_transmission_zeros_hz`,
  `generalized_chebyshev`, `generalized_chebyshev_polynomials`,
  `generalized_chebyshev_with_response`, `synthesize_canonical_matrix`,
  `synthesize_canonical_matrix_with_details`, and `synthesize_matrix_with_topology`
  are `#[deprecated]` (pointing at the module APIs / `FilterDesign`). Still
  functional; removal deferred to a future breaking release.
- **Band-center convention unified to the geometric mean (P2/P9).** The CLI
  `design` subcommand and duplexer `tx()`/`rx()` used the arithmetic mean while
  `FilterDesign::bandpass_edges` used the geometric mean. All edge-based entry
  points now share `freq::edges_to_center_bandwidth`, so `Ω(f1)=-1`/`Ω(f2)=+1`
  holds consistently.

### Added

- **Single-sourced spec construction (P2).** `synthesis::build_filter_spec` is
  the one FilterSpec-assembly helper used by `FilterDesign`, the JSON pipeline,
  the CLI, and the duplexer/multiplexer — replacing the copy-pasted
  order/RL/zeros/unloaded-Q blocks.
- **Duplexer full-Macchiarella iteration is opt-in.** `DuplexerSynthesis::
  full_iteration(bool)` (default `false`) replaces a hardcoded flag, making the
  full spectral path reachable and testable without changing default behavior.
- **Root-solver accuracy tests + criterion bench (P4).** `benches/root_solvers.rs`
  times the four solvers at degrees {10, 20, 30}; new lib tests pin accuracy and
  document a genuine defect: the nalgebra `CompanionMatrixRootSolver`
  complex-coefficient branch (2n×2n real embedding) can return conjugate-flipped
  non-roots, so it stays out of the production `AdaptiveRootSolver` chain
  (faer companion → Durand-Kerner → Aberth) and is used only as a real-
  coefficient test oracle.
- **Engine seam test (P1).** `matrix_synthesizer_is_sensitive_to_generalized_e_s`
  proves the engine responds to `generalized.e_s` — the field it actually
  consumes — guarding the P1 no-op class at the seam.
- **Propagation regression tests (P1).** `adjacent_channels_show_loading`
  (loading must alter the recovered matrix) and `recovered_channels_remain_
  passband` (propagation must not degenerate the channel) in
  `tests/multiplexer.rs`.

### Removed

- **Dead code (P7).** `response::pole_expansion::evaluate_with_pole_expansion_
  verified` (an unused dB-tolerance dispatcher), `spectral`'s
  `balanced_real_companion_roots`/`balanced_complex_companion_roots` + no-op
  `balance_matrix`, and `duplexer`'s `wilson_spectral_factor` + its exclusive
  `fit_polynomial_ls`/`polynomial_max_diff` helpers. ~340 lines.

### Dependencies

- **MSRV raised to 1.86** (`rust-version` + CI test matrix). The locked
  dependency floor (faer 0.20 → 1.81, clap 4.6 → 1.85) made the old 1.75 CI
  leg unbuildable; 1.86 is a safe bound above it. Downstream consumers that
  build this crate's Python bindings (e.g. `maturin develop --features
  python`) must use a Rust toolchain ≥ 1.86.
- **touchstone 0.12.1 → 0.14.5.** Only breaking change touched:
  `Network::new` now returns `Result` (parse errors are mapped to
  `MfsError::PreconditionViolation` instead of the old panic). The upgrade
  exposes `touchstone::NetworkBuilder` (crate root, not feature-gated) for
  programmatically generating S-parameter networks, plus in-memory
  `Network::from_str`/`from_bytes`.

### Fixed

- **`touchstone::cascade` rejected per-port reference impedance as a panic.**
  Networks parsed from a 2-port v2 file with a per-port `[Reference] 50 75`
  keyword line carry `ReferenceImpedance::PerPort` while their scalar `z0`
  field holds the first port's value, so they passed the `z0` pre-validation
  and hit an upstream panic. They are now rejected with
  `MfsError::PreconditionViolation` (regression test added).

### Documentation

- `docs/architecture.md` module list and dependency table now match the code;
  the root-finder fallback order is stated once, canonically (it was previously
  documented four contradictory ways). See
  `docs/internal/architecture-review-2026-09-25.md` for the full assessment.
- Version aligned across `Cargo.toml`, `pyproject.toml`, and `README.md`.
- Post-review doc hygiene: the 0.4.1 release notes (below) were separated out
  from the 0.5.0 heading; the multiplexer "Known limitation" and the two
  multiplexer design docs no longer assert the pre-P1 bit-identical no-op;
  the numerical-optimization report's fallback chain now reads faer →
  Durand-Kerner → Aberth (matching the code and README).

## [0.4.1] — 2026-09-25

### Fixed

- **Odd-order generalized-Chebyshev synthesis produced spurious passband
  transmission zeros and destroyed return loss.** Two compounding defects:
  (a) `build_e_polynomial_stage`
  (`src/approx/generalized_chebyshev_helpers.rs`) did not normalize F(w)/P(w)
  to monic, omitted the `j` factor on the P(w)/eps term of Cameron's
  `E(w) = F(w)/eps_r + j·P(w)/eps`, and did not re-normalize E(s) to monic;
  (b) `synthesize_admittance_polynomials`
  (`src/synthesis/residues.rs`) built the `E(-s)*`/`F(-s)*` terms without the
  global `(-1)^N` factor, zeroing the admittance denominator's leading
  coefficient for odd N and pushing the poles off the imaginary axis. Even
  orders were only unaffected because `j^N` is real. After the fix, orders 3–8
  (all-pole and finite-zero, including the 7p2z B42 case) match the fuqtgui
  Python engine to 0.000 dB. Analysis:
  `docs/odd-order-e-polynomial-fix.md`.
- **Folded reduction corrupted self-couplings for even-order asymmetric filters**
  (`CouplingMatrix::flip_sign` in `src/matrix/core.rs`). The diagonal element
  was negated a third time after the row/column negation, flipping the
  self-coupling sign of any flipped node and silently dropping the transmission
  zeros from the folded response. Removed the redundant negation so `flip_sign`
  is a proper similarity transform (diagonal unchanged), matching the Cameron
  reference. After the fix, folding preserves the response to `≤ 1e-12` for all
  tested orders/parities and reproduces the published book matrices
  entry-for-entry.
- **Fully canonical synthesis (`nfz == order`) was rejected** in
  `PolynomialSet::validate` (`src/approx/polynomial.rs`). The guard required
  `P(s)` to have at most `order` coefficients (degree `≤ order-1`), which
  excludes the fully canonical case where `P(s)` has degree `order`. Corrected
  the bound to `≤ order + 1` (reject only `nfz > order`). The 4-4 fully
  canonical case now synthesizes and folds correctly, reproducing the book's
  Fig 8.21 matrix including the source-load coupling `MSL = 0.0151`.

### Added

- `benches/synthesis.rs` — aligned native Criterion benchmark with
  `fuqtgui/src/scripts/bench_synthesis.py`. Uses the same case matrix (Cameron
  golden + order/TZ sweep) and the same three-stage split
  (`synthesis` / `folding` / `response`), so `cargo bench` native numbers are
  directly comparable to the Python (numpy) and Python-via-FFI (`mfs` module)
  numbers. Methodology in
  `fuqtgui/docs/analysis/benchmark-python-vs-rust-2026-09.md`.
- `tests/folded_response_invariance.rs` — real synthesized-filter regression
  coverage for folded-response fidelity, using Cameron book cases (4-2, 6-3,
  8-3, 6-2, 7-1-2, 4-4 & 6-6 fully canonical, 10-4, all-pole). Fails on the
  even-order asymmetric cases with the old `flip_sign`; passes with the fix.
- `docs/folded-fidelity-flip-sign-fix.md` — analysis of the folding bug and the
  fully canonical synthesis rejection: root causes, fixes, impact, regression
  coverage, and the reconfigured-layout equivalence note.
- `docs/odd-order-e-polynomial-fix.md` — analysis of the odd-order
  generalized-Chebyshev synthesis breakage: root causes in `E(s)` construction
  and admittance parity, the two fixes, and the 11-config cross-validation
  against the fuqtgui Python engine.
- **N-channel multiplexer synthesis and analysis** (`src/multiplexer.rs`):
  - `MultiplexerSynthesis` — builder for iterative N-channel synthesis over a
    common (antenna) junction, generalizing the 2-channel diplexer's simplified
    perturbation path (`channel`/`channels`/`junction`/`max_iterations`/
    `convergence_tol`/`synthesize`).
  - `MultiplexerSynthesisResult` — per-channel `FilterDesign`s, iteration count,
    convergence error, junction, and per-channel S-polynomial roots.
  - `Multiplexer` — N+1-port analysis (`response`, `channel`, `channel_folded`,
    `channel_matrices`) with a geometric-mean response model (max-S11 common
    reflection, product-of-S21 isolation).
  - `MultiplexerResponse` — frequency grid, antenna S11 (dB), per-channel
    S21/S11 (dB), and pair-wise isolation (dB).
  - Re-exported from `crate::prelude`; wired into the CLI (`multiplexer`
    subcommand, `--features cli`) and Python bindings (`PyMultiplexer` /
    `multiplexer()`, `--features python`).
- `docs/multiplexer-synthesis.md` — design, math generalization, API surface,
  phase breakdown, and as-built implementation deviations.
- `tests/multiplexer.rs` — 6 integration tests (N=2 diplexer parity, triplexer /
  quadruplexer synthesis & convergence, N-port response shape, common-port and
  isolation structure).
- **Known limitation** (shared with `DuplexerSynthesis`): mutual loading uses
  the simplified perturbation path, not the full Macchiarella N×N p₀ spectral
  solve. As of 0.5.0 (P1) the perturbed E/F propagates into the output coupling
  matrices, but the model is approximate: it does not converge to machine
  tolerance for well-separated channels and adjacent channels pay a passband-peak
  cost. (The earlier "fully-reflective recovered channel" symptom was the
  pre-0.4.1 odd-order residue-extraction bug, since fixed — see `### Fixed`
  above — and the diagnostic test `multiplexer_recovered_channel_is_passband`
  is un-`#[ignore]`d and passing.)

## [0.2.0] — 2025-05-25

### Added

- **Touchstone (.s2p) file export** — `design.save_touchstone("filter.s2p")`
  - Supports v1.0 and v2.0 formats
  - Supports RI, MA, DB data formats
  - Supports Hz, kHz, MHz, GHz frequency units
  - Configurable reference impedance
  - Auto-generated comment headers with design metadata
  - Round-trip verified (export → import → compare)
- Touchstone file import via `touchstone::read_touchstone(path)`
- `FilterDesign::to_touchstone()` builder for configured export
- Python `save_touchstone()` and `to_touchstone()` methods
- `TouchstoneVersion::V2` support with `[Network Data]`/`[End]` markers

### Dependencies

- Added `touchstone` v0.12.1 (for file parsing only)

## [0.1.0] — 2025-05-25

### Added

- `FilterDesign` high-level API with `bandpass()` and `prototype()` entry points
- Generalized Chebyshev (Cameron) polynomial synthesis for arbitrary order
- Coupling matrix synthesis via residue expansion
- Topology transforms: Transversal → Folded, Arrow
- S-parameter response evaluation with two backends:
  - LU-based solver (exact, O(N³) per frequency point)
  - Pole-expansion solver (fast, O(N) per frequency point, auto-verified)
- Adaptive root-finding with three-level fallback:
  - Durand-Kerner (fastest, degree ≤ ~28)
  - Aberth-Ehrlich (cubic convergence, degree 25-40)
  - Companion Matrix eigenvalues (most robust, any degree)
- In-place Givens similarity rotations (O(N) per rotation vs O(N³))
- Band-pass scaling and external Q computation
- Section extraction: triplet, quadruplet, trisection
- Python bindings via PyO3 (optional `python` feature)
- 164 unit and integration tests
- High-order stability verified up to order 30

### Performance

- Matrix topology transforms: ~500× faster than naive implementation
- Frequency response (pole expansion): ~100-500× faster than LU for large sweeps
- Numerical precision: machine epsilon (2.2e-16) symmetry error at order 30
