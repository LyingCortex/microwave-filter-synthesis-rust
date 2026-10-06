# MFS Continuation Plan

**Date:** 2026-09-28 · **Baseline head:** `2d735c0` (`master`, pushed to gitlab)
**Status:** all CI-parity checks green on baseline (lib 202 + 15 integration
binaries, clippy `-D warnings`, rustdoc — zero warnings).
**Current (0.6.0, 2026-10-05):** all Wave 1–5 items DONE; `master` head
`2106465`; lib 250 + 19 integration suites, 0 failed; clippy `-D warnings`;
doc — zero warnings. Tag `v0.6.0`.

This is the single authoritative continuation plan for the four active
development areas — **duplexer, multiplexer, multi-band, optimization** — plus
the cross-cutting numeric enablers they depend on. `docs/internal/todo.md`
remains the open-item list; this document gives each item its context,
sequence, and acceptance criteria. Per-area design docs
(`docs/duplexer-synthesis.md`, `docs/multiplexer-synthesis.md`,
`docs/numerical-optimization-report.md`) carry the math; this document
carries *what to do next and how to know it is done*.

---

## 0. Baseline: what is shipped and how good it is

| Area | Shipped state | Fidelity class |
|---|---|---|
| Two-port synthesis | Full generalized-Chebyshev pipeline, orders 2–8 fixture-anchored, 10–30 numerically stable | High (publication-anchored) |
| Duplexer | Simplified frequency-offset perturbation (default); full Macchiarella–Tamiazzo spectral path behind `DuplexerSynthesis::full_iteration(true)` **converges at Type-I orders 2–8 and Type-II (incl. the §VI degree-40 GSM design), all gate-anchored** (2026-10-03) | High (publication-anchored, full path) / Medium (simplified default) |
| Multiplexer (N-channel) | Independent per-channel synthesis + simplified perturbation (default); **full Macchiarella–Tamiazzo spectral path for N channels behind `MultiplexerSynthesis::full_iteration(true)`** (shared with the duplexer, `src/synthesis/machiarella.rs`; B2, 2026-10-04) converges at N = 3/4 order-4 (gate-anchored); response is the **true (N+1)-port admittance model** (B1) | High (publication-anchored, full path) / Medium (simplified default) |
| Multi-band | **True multi-passband synthesis shipped (C2/C3, 2026-10-05)** — route-α multi-band Chebyshev: `MultiBandMapping` (C1) + `MultiBandSpec` (C2a) + `multi_band_chebyshev_polynomials` (C2b) → `MatrixSynthesisEngine` (C2c, unchanged) → Fig. 25 dual-band anchor; response invariance, per-band RL, TZ nulls, unitarity, outer-stopband all gate-anchored (C3, `tests/multi_band.rs`); **K ≥ 3 infeasible** for this characteristic (DOF wall: equiripple is square only for one positive-axis interval; M ≥ 2 overdetermined — measured, design doc "K ≥ 3: the DOF Wall") | High (publication-anchored, dual-band) |
| Optimization | LM `tune_matrix` / `extract_matrix` (hand-rolled Gaussian elimination, **scale-aware FD Jacobian `δ = 1e-7·max(1,|x_i|)`** with halving retry, topology frozen to template; convergence on successive-iteration cost after accepted steps; singular trials rejected, not errors; **`extract_matrix` starts from a crude inline-chain guess**); 3 lib unit tests + **2 order-6 integration gates (synthesize→perturb→recover, extract-from-crude-guess) GREEN (D2-opt, 2026-10-04)** + **CLI `optimize`/`extract` subcommands and Python `extract`/`tune` (D3-opt, 2026-10-04), e2e-gated** + **`benches/optimizer.rs` & 0.1 dB-noise extraction regression (D4-opt, 2026-10-04)** | Medium (order ≤ 6 extraction/recovery gate-anchored; CLI/Python e2e-gated; benched; order-8 cold-start extraction and topology-bound recovery documented as limits, report §6) |

Known limitations that gate the areas below (full list in
`docs/architecture.md` "Known Limitations"):

1. **Passband transmission zeros** (`|zero| < 1`) — rejected with a typed
   `MfsError::InvalidTransmissionZero` at the `approx` entry point; E1's
   2026-10-05 investigation showed the single-function prototype cannot
   realize them (Feldtkeller violation) — see §E item 1.
2. ~~**Repeated poles**~~ — handled by the Laurent residue expansion (E2, 2026-10-05): `residues/expansion.rs` returns the principal part for a repeated pole and the pole-expansion fast path sums it; the transversal (inline-chain) builder rejects a repeated pole with a typed `NotImplemented` (a chain cannot represent `1/(s-p)^m`) — see §E item 2.
3. ~~**HighPass / BandStop mappings**~~ — implemented in `freq` (E3, 2026-09-29: `HighPassMapping`, `BandStopMapping` with full forward + inverse and singular-point clamps).
4. **N-channel mutual-loading fidelity — RESOLVED for the full path (B2, 2026-10-04):**
   the full spectral path is now generalized to N channels (shared
   `src/synthesis/machiarella.rs`; N×N p₀² LU solve, per-channel outer-edge
   return loss, N-block root assignment; N = 2 branch bit-identical to the
   historical loop). The multiplexer's *default* still uses the simplified
   under-relaxed perturbation; B3 adds the convergence/fidelity gates that
   flip the vendored-standard convergence gate from FAIL to PASS.

---

## A. Duplexer — converge the full Macchiarella–Tamiazzo path

**Goal:** make the full spectral iteration (`full_iteration(true)`) converge
for order ≥ 4, so the duplexer can deliver true quasi-equiripple antenna
return loss instead of the simplified perturbation's approximation.

### Status (evidence)

- Simplified path (default): `src/duplexer.rs` — frequency-offset perturbation
  loop; channels converge to valid passband filters; antenna S11 is the
  geometric-mean model `√(|S11_tx|·|S11_rx|)`.
- Full path: fixed reflection-zero `N` (channel F zeros + the Type-II real
  zero) → `solve_scaling_constants` → `spectral_factor` (cepstrum +
  Parlett–Reinsch-balanced faer companion fallback in `src/spectral.rs`) →
  `compute_s_product` → `assign_roots`. Wired behind
  `DuplexerSynthesis::full_iteration(true)` (default `false`).
  **Converged and gate-anchored at Type-I orders 2–8 and Type-II
  (including the §VI degree-40 GSM design)** — see the D1–D4 items below and
  `duplexer::tests::full_iteration_*`.

### Work items (ordered)

1. **D1 — Make the full path testable and failing-visible.** Add lib tests
   that call `full_iteration(true)` at order 2 (TX/RX of order 2 each):
   assert convergence, quasi-equiripple antenna RL within ~2 dB of the
   independent-filter bound, and passband peaks within 0.5 dB of the
   per-channel spec. **Orders 2–8 are all GREEN** (2026-10-03): the gate
   family `full_iteration_order{2,3,4,5,6,8}_quasi_equiripple` (no longer
   `#[ignore]`d) plus the Section-V printed-design golden benchmark
   `full_iteration_order7_section5_golden` all pass. Each asserts the
   junction-level diplexer response (`result.diplexer_response`, added
   2026-10-03) rather than the standalone recovered channels, which each
   carry only their share of the junction power split and cannot reach 0 dB.
   See the 2026-10-03 note in D2 for the full rationale (including the
   retraction of the "degenerate odd orders" claim).
2. **D2 — Diagnose the order ≥ 4 non-convergence.** Instrument the
   iteration (residual per iteration, root-displacement via
   `compute_root_displacement`, spectral-factor relative error from
   `spectral.rs`) in a diagnostic example. Root-cause candidates to check,
   in order: (a) cepstrum spectral factorization mis-assigns LHP roots at
   higher degree (verify `select_lhp_roots` against a known factorization),
   (b) `solve_scaling_constants` ill-conditioning as the p₀ system grows,
   (c) under-relaxation `alpha` schedule too aggressive for the N(s) root
   movement. Fix the root cause, not the tolerance.
   **D2 findings (2026-09-29, verified via instrumented trajectories):**
   - **(a) spectral factorization ruled out.** With p₀ held fixed, the
     iteration is a clean contraction (root displacement 0.55 → 0.02 over
     37 iterations, monotone); `select_lhp_roots` / cepstrum is not the
     breaker.
   - **(b) The 2×2 edge-RL p₀ solve, verified against the code.**
     `solve_scaling_constants` (`duplexer.rs:720`) implements the
     unitary-condition form with the antenna RL spec substituted in
     (`|D|² = |N|²/ε` at each edge) and `n₀` fixed at 1 — coherent with
     the unitary condition, and `n₀ = 1` is **confirmed by the paper**
     (eq. 12, verified by OCR: “|n₀| = 1 for both junction types”). The
     doc's Step-3 line (`|N|² + |p₀|²|P|² = 10^(−RL/10)`) is the
     unphysical transcription (it predicts 0 dB antenna RL at the edge).
     Measured: the 2×2 is singular for the D1 gate's *symmetric* setup
     (equal order/RL/adjacent bands): `det → 0` at iter 0, degenerate
     fallback fires (p₀² ≈ 707), p₀² then runs to ~1e8 and the
     iteration plateaus at max_delta ≈ 0.07.
   - **Fixed p₀ = 1 is not the answer.** Unity passband gain removes the
     divergence but converges to a *physically wrong* diplexer: recovered
     channel passband peaks −4.9 / −5.4 dB (spec ≈ 0), antenna ≈ −1 dB
     (spec ≈ −22). p₀ = 1 is a false lead — it is a contraction to the
     wrong fixed point.
   - **Asymmetric probe resolved (2026-09-29).** An asymmetric 2/4
     diplexer (different order/RL, adjacent bands,
     `full_iteration(true)`, 40 iters) diverged:
     convergence_error → 2.89, recovered channels degenerate (Tx S21
     ≈ −46 dB, Rx ≈ −128 dB). `n₀ = 1` is confirmed correct by the paper
     (eq. 12), so the divergence is a construction mismatch in N/Pt/Pr/D,
     not the n₀ normalization — and not a symmetric-degeneracy artifact.
   - **D3 unblocked (2026-09-29, paper OCR).** The Macchiarella–Tamiazzo
     paper is vendored in the repo (`tests/Macchiarella and Tamiazzo -
     2006 - Novel Approach to the Synthesis of Microwave Diplexers.pdf`);
     its equations are image-rendered but were recovered by OCR. Verified
     from the paper: `|n₀| = 1` (eq. 12); the p₀ equation the code
     implements (eq. 13 reduced via the unitary condition + edge RL spec);
     loop steps 1–5. **Section V prints a complete worked Type-I design**
     (RX 14.9–15.1 / TX 15.15–15.35 GHz, 20 dB, 7 poles each,
     `n = 1.47`, `b₀ = −0.171`) with all coefficients: `n₀N`, `D`, `Pt`,
     `Pr`, `p₀t = 4.717e−4+1.743e−4i`, `p₀r = 4.336e−4+1.6e−4i` (complex
     p₀; the code currently keeps p₀ real). That worked design is the
     regression fixture: run the loop on it and match the printed
     polynomials.
   - **Corrected root cause (2026-10-01, decisive probe against the printed
     Section-V design).** The earlier "construction mismatch / 2×2
     near-singular" conclusions are **wrong**. The printed design is
     self-consistent to print precision (unitary identity
     `D·D* ≈ |n₀N|² + |p₀tPt|² + |p₀rPr|²` holds to ~1e-3 at high order,
     degrading only in the tiny low-order coefficients, consistent with the
     4–5-sig-fig print). Decisive tests:
     - The eq.-13 2×2 is **well-conditioned** (cond ≈ 1.09 at iter 0;
       `pt_tx≈pr_rx≈117.8`, `pr_tx≈pt_rx≈5.3`), **not** near-singular.
     - The eq.-13 2×2 **formula is correct**: at the printed design's two
       outer passband edges (`|S11|=0.1`, `s≈±j0.74`) it reproduces the
       printed `|p₀t|²=2.5e-7` / `|p₀r|²=2.1e-7` (got 4.0e-7 / 1.7e-7).
     - The breaker is a **domain + normalization mismatch**: the code builds
       a narrow-band-linearized "common domain"
       (`Δ = 2·Δcenter/avg_bw ≈ 2.5`) and leaves **N un-normalized**, while
       the paper works in the true bandpass w-domain (RX/TX centered
       ≈ ∓j0.55/+j0.56, `Δ≈1.11`) with **monic N**. In the stretched domain
       `|p₀|²` explodes (order 2: 1.6e3 → 1.4e8 over 5 iters) and the root
       set goes non-finite at **every** order 2/4/6/7.
     - **D3 is blocked on a domain redesign, not a p₀-formula fix**: the full
       path must be re-based on the paper's bandpass w-domain (or the
       narrow-band linearization replaced by the exact transform) and N
       normalized monic before the eq.-13 loop can converge. Evidence
       scripts: `examples/paper_p0_edge_test.py`,
       `examples/paper_domain_probe.py`, `examples/diag_paper_section5.rs`.
   - **(2026-10-02, order-2 re-probe on `wip/d3-duplexer-full-iteration`.)**
     With the loop left as-is (common domain, un-normalized N) the order-2
     D1 gate *does* converge (`err 4.5e-10`, 16 iters) and the recovered
     pair is **lossless in the common domain** (`E·E~ − F·F~ = p₀t²·P·P~`,
     round-off residues ~1e-8; the converged `p₀t² = 0.872`, `p₀r² = 0.813`).
     The common-domain recovered F has off-axis (non-minimum-phase) zeros,
     so the *standalone*-recovered channel's peak S21 = `√p₀t² = −0.564 dB`
     (not 0 dB). Recomputing eps/eps_r from the recovered E/F does **not**
     restore 0 dB (F's zeros are off-axis, so the S11 null is absent) and was
     reverted. **The gate conflict is not a p₀/implements bug**: the recovered
     channel's own S21 is `p₀TX·P_TXn/E_TX` by eq. model, and `p₀TX²` is the
     junction's TX/RX power split — *physically impossible to remove from the
     isolated 2-port*. The D1 gate measures `Duplexer::response.tx_s21_db`,
     which samples each recovered channel as an isolated 2-port (its own
     bandpass mapping, own matrix), **not** the junction-combined TX-port
     transmission. The diplexer-level response (computed from the recovered
     E/F/P and p₀ in the common domain) is the correct physical observable
     and **passes**: TX-band `|S21| = 0.962 (−0.34 dB)`, antenna `|S11| =
     0.068 (−23 dB)`, unitarity = 1.0000. **Decision required from the
     maintainer**: either (a) re-specify the D1 gate to assert the
     diplexer-level TX transmission (requires `Multiplexer::response` to
     combine the recovered matrices through the junction), or (b) document
     that the standalone-channel peak −0.6 dB is the correct full-path
     result and the gate was mis-specified. Do NOT weaken the gate to hide a
     defect — the diplexer-level evidence shows it is a test-target bug, not
     a synthesis defect.
   - **(2026-10-03, D1 gate re-specified — GREEN.)** Adopted option (a): the
     gate now asserts the *junction-level* diplexer response rather than the
     standalone recovered channels. `DuplexerSynthesisResult` gained a
     `diplexer_common` snapshot (converged N, D, Pt, Pr, p0t², p0r², f0, b)
     and a `pub fn diplexer_response(start, stop, points) -> Option<
     DuplexerResponse>` that evaluates `S11 = N/D`, `S21_tx = √p0t²·Pt/D`,
     `S23 = √p0r²·Pr/D` in the common bandpass variable `s = j·Ω(f)` (|n0| = 1
     for the ideal Type-I junction). The D1 test
     (`full_iteration_order2_quasi_equiripple`) now calls
     `result.diplexer_response(...)` and asserts: convergence < 1e-8 (14
     iters), Tx-port S21 > −0.5 (got −0.10 dB), Rx-port S21 > −0.5 (got
     −0.11 dB), antenna RL ≤ −20 dB at both band centers (got −23.4 dB), and
     no degenerate antenna port. The gate is **no longer `#[ignore]`d** and
     passes. The standalone-channel `Duplexer::response` (in `duplexer.rs`)
     is left as-is: it is still the correct observable for the *simplified*
     path (no full spectral iteration, no junction-level model); extending
     it to a true junction-level 3-port is out of scope here.
   - **Clippy parity restored on the WIP branch.** Removed the now-dead
     Steps 16-17 `rederive_channel_d` / `rederive_channel_d_at` path and the
     `d_tx_common`/`d_rx_common` loop state (recovery uses the
     normalization-independent `D = S − F⁰`), the `rescale_d` argument to
     `recover_channel_design` (restored master's always-rescale 7-arg
     signature — a no-op for the monic `S − F⁰` D), the `/tmp` debug writes,
     and the throwaway `examples/diag_*.rs`. `cargo clippy --all-targets --
     -D warnings` is clean.
   - **(2026-10-03, order 4/6 verified GREEN; order 7 resolved via stagnation.)** The
     earlier claim that "order ≥ 4 fails in `recover_channel_design` with a
     typed pole-placement error" was for the *old* narrow-band common domain.
     The exact bandpass common-domain transform (in the WIP rewrite) already
     converges order 4 (14 iters, conv 8.9e-9) and order 6 (9 iters, conv
     7.7e-9) to a physically valid quasi-equiripple diplexer: Tx/Rx-port S21
     ≈ −0.03 dB, antenna RL ≈ −22 dB at all four band edges and both band
     centers. The D1 gate was refactored into
     `assert_diplexer_quasi_equiripple(order)` and extended with
     `full_iteration_order4` and `full_iteration_order6` tests — all green.
     **Order 7 (and 5, 8) — resolved (2026-10-03).** The ~1e-6 floor is
     NOT a divergence: p0t²/p0r² and min|D| freeze to 5 sig figs by iter 3,
     and the product-roots freeze to 6 digits by iter 2 — the iteration
     converges to a stable fixed point, but the degree-14 root-extraction
     noise floor (~5e-6) sits above the strict 1e-8 tolerance, so
     `max_delta` never triggers the capture. Fix: a **stagnation criterion**
     (accept once `best_delta < 1e-3` after 3 flat iterations, min 5 iters)
     — the standard remedy for root-iteration, preserving the strict fast
     path for orders 2/4/6.
     - **Odd orders are valid, not degenerate (corrected 2026-10-03).** An
       earlier note called the symmetric ideal-junction order-7 "degenerate"
       (spurious center null, near-singular 2×2). Both parts were wrong:
       (i) the deep antenna null at an odd-order band *center* is the channel's
       Chebyshev reflection zero — a physical feature, and the §V golden
       (order 7, which matches the paper's printed coefficients) shows the
       same null; (ii) the 2×2 p0 system is well-conditioned in the exact
       bandpass domain (det/scale² ≈ 0.6–0.92 at both symmetric and §V
       junctions) — the near-singular 2×2 was an artifact of the *old*
       narrow-band linearization.
     - The **Section-V printed design** (`n = 1.47`, `b0 = −0.171`, 20 dB,
       RX 14.9–15.1 / TX 15.15–15.35 GHz, **all-pole inline Chebyshev**) is a
       *golden benchmark* and is now matched by
       `full_iteration_order7_section5_golden`. The §V spec is the same
       standard-Chebyshev geometry the code synthesizes, so the converged
       common-domain polynomials are compared coefficient-by-coefficient
       against the paper's printed N, D, Pt, Pr and p0. Measured (2026-10-03):
       |p0t|² 1.2e-4, |p0r|² 1.6e-4, D 1.2e-2, Pt 9.4e-3, Pr 5.9e-3, N 1.5e-2
       (relative, after leading normalization) — all at the paper's
       *print-rounding* level (p0 to 4 sig figs, polynomials to 2–3 sig figs),
       i.e. the code reproduces the printed design, not a near-miss. The gate
       also asserts the physical §V response (0 dB passbands, band-edge
       antenna RL within ~2 dB of the 20 dB spec).
     - **Note for future work:** an earlier note (2026-10-03) wrongly said the
       printed coefficients "cannot be matched" because the paper imposes
       specific transmission zeros (method [11]). That confusion came from §VI
       (the Type-II GSM design, which does use imposed zeros); §V is all-pole
       Chebyshev and matches to print precision.

3. **D3 — Rebase the full path onto the paper's bandpass w-domain, then
   converge order ≥ 4 and pin it. (Core deliverable DONE, 2026-10-03.)**
   The rebase onto the exact bandpass common-domain transform (Δ from the
   true `Ω(f)=(f/f₀−f₀/f)/B`, N from the fixed reflection zeros) is in the
   WIP rewrite; it converges order 2 (14 iters, 6.1e-9), 4 (14 iters, 8.9e-9),
   6 (9 iters, 7.7e-9) to the strict tolerance, and orders 3/5/7/8 (strict
   1e-8 at 3/5; the stagnation criterion at 7/8 where the degree-14/16
   root-extraction floor exceeds 1e-8). D1 gates
   `full_iteration_order{2,3,4,5,6,8}` and the **§V golden benchmark**
   `full_iteration_order7_section5_golden` are all green. The golden test
   matches the paper's printed order-7 coefficients (N, D, Pt, Pr, p0) to
   print precision — p0t² 1.2e-4, p0r² 1.6e-4, D 1.2e-2 (relative, after
   leading normalization) — i.e. the implementation reproduces the published
   design. The earlier claim that the symmetric ideal-junction odd orders land
   on a "degenerate" fixed point is retracted: a per-order scan (2026-10-03)
   shows all orders 2-8 give a valid quasi-equiripple diplexer (0 dB passbands,
   band-edge antenna RL ~ -22 dB; the deep odd-order *center* nulls are the
   Chebyshev reflection zeros, a physical feature the §V golden also shows),
   and the 2×2 p0 system is well-conditioned (det/scale² ≈ 0.6-0.92) at both
   symmetric and §V junctions. **Type-II (resonating-node) junctions are now
   also run through the full path (2026-10-03):** the paper's eq. (11) is
   implemented exactly — the extra pure-real reflection zero at `s = 1/c0`,
   D degree `n_tx+n_rx+1`, and `S_TX·S_RX = (c0/2)(D − N)` with D normalized
   exactly monic (the `D−N` leading-term cancellation). Gated by
   `full_iteration_type2_order4_quasi_equiripple` and the §VI GSM anchor
   `full_iteration_type2_section6_anchor` (RX 10 / TX 9 poles, degree-40
   R(s)); the anchor reproduces the paper's quasi-equiripple |S11| (peaks at
   the 22 dB spec, peak-to-peak < 1.5 dB per Fig. 9). Enabling this required
   **Parlett–Reinsch diagonal-similarity balancing of the faer companion
   matrix** (`FaerRootSolver`): an unbalanced solve of the ~19-decade-spread
   §VI R(s) gave backward residual ~0.75 (garbage), balancing brings it to
   ~1e-9. Changing the *default* to `full_iteration(true)` is still a
   separate, deliberate behavior change with its own changelog entry.
4. **D4 — Document the as-built fidelity. (DONE, 2026-10-03.)**
   `docs/duplexer-synthesis.md` updated: the Type-II polynomial forms now
   match the paper's eq. (11) (the N/D cross term, the extra real zero at
   `s=1/c0`, the degree `n_tx+n_rx+1`, and the `S_TX·S_RX = (c0/2)(D−N)`
   Step-4 reconstruction with D normalized exactly monic); the as-built
   status states the full path converges at both junction types (Type-I
   orders 2–8 + §V golden; Type-II §VI GSM anchor) and names the
   Parlett–Reinsch-balanced companion-matrix spectral factorization; the
   stale `wilson_spectral_factor` / `fit_polynomial_ls` references and the
   "not numerically converged" status are removed.

**Acceptance (section done):** `cargo test --lib -- duplexer` green with
full-path tests at Type-I orders 2/3/4/5/6/8 + §V golden and Type-II
order-4 + §VI anchor; `docs/duplexer-synthesis.md` accurate; no regression
on the simplified path (its tests unchanged).

---

## B. Multiplexer (N-channel) — fidelity, then unification

**Goal:** (1) a true (N+1)-port response model, (2) high-fidelity mutual
loading via the full spectral path generalized to N, (3) one
Macchiarella–Tamiazzo iteration loop shared by duplexer and multiplexer.

### Status (evidence)

- `src/multiplexer.rs` `synthesize`: independent per-channel synthesis,
  junction coupling (`TypeI`/`TypeII`), under-relaxed frequency-offset
  perturbation, recovery through the shared
  `mfs::synthesis::recover_channel_design` (single source since 0.5.0+;
  writes perturbed E/F back to `generalized.e_s/f_s`, so the correction
  reaches the output matrices — P1).
- **Model caveat (documented, `docs/multiplexer-synthesis.md` deviation 5):**
  the simplified perturbation does not converge to independent synthesis
  even for well-separated channels; the cross-channel term keeps a non-zero
  residual, so far-apart channels pay a small passband-peak cost.
- **Response is now a true (N+1)-port model (B1, 2026-10-04):**
  `Multiplexer::response` defaults to the admittance node model — each channel
  is a lossless reciprocal 2-port (port A at the common node) and the full
  complex (N+1)×(N+1) S-matrix is composed at the node (unitary + reciprocal;
  N=1 reduces exactly to the channel's 2-port). `response_matrix` exposes the
  complex matrix; `response_with(.., ResponseModel::GeometricMean)` keeps the
  legacy magnitude model for one release cycle.
- **Full-spectral path now generalizes to N channels (B2, 2026-10-04):** the
  shared `src/synthesis/machiarella.rs` loop powers both
  `DuplexerSynthesis::full_iteration` and `MultiplexerSynthesis::
  full_iteration` (opt-in; the default stays the simplified perturbation).
- **Benchmark-standard gates:** the vendored standard's `cm`, `star`,
  `s11_match`, `rz_distribution` gates have no mfs counterpart and are
  `SKIP`. The convergence gate is now **PASS** on every runnable case
  (B3, 2026-10-04): `benches/multiplexer.rs` runs the `poly` stage through
  the full spectral path, and `convergence_gate` reads the synthesis's own
  `converged` flag — case `a` (Zhao N=4) converges in 9 iters (conv_err
  1.7e-4, stagnation), `diplexer` in 12 (2.6e-9, strict); `type2`/
  `contiguous` cases remain capability-SKIP.

### Work items (ordered; depends on A)

1. **B1 — True (N+1)-port response (magnitude model out). [DONE 2026-10-04]**
   Common-node model: each channel is a lossless reciprocal 2-port (port A at
   the common node, port B external); the full complex (N+1)×(N+1) S-matrix is
   composed at the node via the wave-variable node equation.
   `Multiplexer::response_matrix` exposes the complex matrix;
   `Multiplexer::response` defaults to this model (dB output); the legacy
   geometric-mean model is kept as `ResponseModel::GeometricMean` via
   `response_with` for one release cycle. Gate tests: unitarity + reciprocity
   (1e-9, lossless inputs), N=1 exact reduction to the channel 2-port, stopband
   agreement between the two models. Note: the models do *not* agree in the
   passband (geomean underestimates the common-port reflection by up to ~150
   dB) — the original "agree to ~0.1 dB where accurate" criterion only holds in
   the all-reflective stopband, which is what the gate test asserts.
2. **B2 — Generalize the full spectral path to N channels. [DONE 2026-10-04]**
   Ported `duplexer.rs`'s full-path loop into `src/synthesis/machiarella.rs`
   (shared by `DuplexerSynthesis::full_iteration` and
   `MultiplexerSynthesis::full_iteration`), generalized to N channels: the
   N×N p₀² system (Cramer 2×2 verbatim for N = 2, LU for N ≥ 3), per-channel
   outer-edge return-loss selection (±j for N = 2 exactly as before; each
   channel's own outer band edge mapped to the common axis for N ≥ 3), and
   contiguous N-block root assignment in ascending channel-center order. The
   N = 2 branch is bit-identical: all 10 gate-anchored duplexer tests
   (orders 2–8 quasi-equiripple, §V/§VI printed-design goldens) pass
   unchanged. Non-ideal Type-I at N > 2 is rejected with an explicit
   `PreconditionViolation` (no paper basis for that generalization).
3. **B3 — Convergence + fidelity gates for N = 3, 4. [DONE 2026-10-04]**
   Added full-path gates in `tests/multiplexer.rs`:
   `multiplexer_full_iteration_well_separated_matches_machine_tol` (two
   far-apart order-5 channels — the iteration reaches `convergence_error <
   1e-6`, in fact ~0, the contrast against the simplified path's documented
   bounded residual), `_adjacent_bounded_loading` (overlapping channels —
   bounded loading, quasi-equiripple common port), and
   `_vendored_zhao_n4_converges` (the vendored standard's N=4 case with
   transmission zeros). Each asserts convergence, 0 dB per-channel
   passbands, the antenna RL pinned within 2 dB of the per-channel spec at
   each outer band edge, and the losslessness para-identity `D·D̃ = R`
   (regime-aware bound: strict-converged → 1e-4, high-degree
   stagnation-accepted → 5e-2, the duplexer order-7/8 pattern). **Flips the
   vendored-standard convergence gate FAIL→PASS:** `benches/multiplexer.rs`
   now runs the `poly` stage through the full spectral path and
   `convergence_gate` reads the synthesis's own `converged` flag (an
   exhausted budget is a non-convergence whatever the residual) — case `a`
   (Zhao N=4) converges in 9 iters (conv_err 1.7e-4, stagnation), `diplexer`
   in 12 (2.6e-9, strict); `type2`/`contiguous` remain capability-SKIP.
   Note: the original "well-separated channels converge to *independent*
   synthesis to 1e-6" wording is re-scoped — the full path converges to the
   *loaded* fixed point (displacement from independent ~1e-2, bounded), not
   to independent synthesis; the machine-tolerance claim holds for the
   *iteration's* convergence error, which is what the gate asserts.
4. **B4 — Close the bench-standard gaps. [DONE 2026-10-04]** The
   `cm`/`star` gates are now reported per gate with specific reasons instead
   of a blanket `SKIP`:
   - `unitarity` → **PASS** (measured, not skipped). mfs computes the full
     complex (N+1)×(N+1) common-junction S-matrix directly at the node
     (`Multiplexer::response_matrix`, the B1 wave-variable node equation on
     lossless per-channel 2-ports), so the contract's star-junction unitarity
     invariant — antenna-row power conservation
     `max|1 − |S11|² − Σᵢ|S1i|²| ≤ 1e-8` — is measured on that response.
     `gates.rs` gains `multiplexer_unitarity_residual` + `unitarity_tol`;
     `benches/multiplexer.rs` sets the gate in the `response` stage.
     Measured on both runnable standard cases: `multiplexer_a` (Zhao N=4)
     residual 3.8e-15, `multiplexer_diplexer` 1.3e-15 (machine precision —
     the KCL junction of lossless 2-ports conserves power by construction).
     The case aggregate status moves from SKIP to **PASS**.
   - `matrix_size` → SKIP: mfs does not assemble the contract's single
     multiport star-junction CM (n_total+n_ports[+1]); it synthesizes each
     channel to its own (nᵢ+2) CM and composes them at the node, so
     `cm.shape == invariants.cm_size` has no mfs quantity to measure.
   - `equiripple` → SKIP: the contract metric compares |S11| against the
     per-channel determinants |Δᵢ| at band extrema; mfs does not expose |Δᵢ|
     (a property of the multiport CM star transform). The quasi-equiripple
     antenna RL is instead asserted by the B2/B3 full-path fidelity gates.
   - `s11_match` → SKIP: cross-checks the multiport CM response against the
     `poly`-stage rational-function S11; mfs produces neither a multiport CM
     nor a rational-function S11, so there is no second S11 to match.
   - `rz_distribution` → SKIP (unchanged): mfs exposes coupling matrices, not
     the contract's rational-function pole/zero set.
   Capability-blocked cases (type2/contiguous) still SKIP all gates with the
   deciding tag in `detail`.

**Acceptance (section done):** B1–B4 green; `docs/multiplexer-synthesis.md`
deviation list 2 and 3 resolved or explicitly re-scoped; multiplexer N=2
output matches duplexer output (single loop proves parity).

---

## C. Multi-band — true multi-passband synthesis (scope-1)

**Goal:** filters with more than one passband in a single structure —
equiripple response in ≥ 2 passbands, rejection in between and outside.
**Scope decision (2026-09-29, maintainer):** scope-1 confirmed — a real
multi-passband synthesis, not per-band design composition. This is the
largest planning gap of the four areas: no code, no design doc, one
roadmap line — and the one that most needs new mathematics.

### Status

- `src/freq.rs` implements all four prototype families
  (`LowPassMapping`, `HighPassMapping`, `BandPassMapping`,
  `BandStopMapping`) plus `MultiBandMapping` (E3/C1, 2026-09-29).
- **C0 design doc done (2026-10-05):** route committed (**C-route-α**),
  math/spec/TZ-law/topology/fixture fixed in
  `docs/internal/multi-band-synthesis-design.md`.
- **C2 done (2026-10-05):** `src/approx/multi_band.rs` (route-α multi-band
  Chebyshev, `PolynomialSet` output, LM equiripple solve) + `MultiBandSpec`
  (C2a, `src/spec/multi_band.rs`) + `MatrixSynthesisEngine` reuse (C2c, no
  new engine). **K ≥ 3 infeasible** (measured 2026-10-05) — the equiripple
  solve is square only for a single positive-axis passband interval; a center
  band spanning the origin (or 2+ positive-axis bands) overdetermines the
  real-rooted characteristic by M − 1 (DOF wall; see the design doc's
  "K ≥ 3: the DOF Wall" for the proof and the measurement).
- **C3 done (2026-10-05):** `tests/multi_band.rs` (Fig. 25 anchor — response
  invariance, per-band RL, TZ nulls, unitarity, outer stopband) +
  `examples/multi_band_response.rs`.
- The vendored benchmark standard has **no** multi-band case family —
  there is no external gate; the acceptance bar is the in-repo literature
  fixture `Cameron2011_10_5_dual_band_symmetric` (Cameron 2011 Fig. 25).

### Why generalized Chebyshev cannot do this

The shipped approximation family represents one lowpass prototype (one
passband, stopband zeros) mapped to a band. A K-passband spec needs an
approximation with **K equiripple passbands and 2K−1+ stopbands**, i.e. a
different E(s)/F(s) structure. Two candidate routes (the design doc
C0 must decide between them before any code):

- **C-route-α — multi-band Chebyshev approximation.** Equiripple in each
  band via a product of per-band band-Chebyshev factors
  (∏ over bands of the band mapping's (1+s²)-type term), with the
  characteristic polynomial's zeros placed in every stopband. Yields a
  Cameron-style transversal ladder: each stopband gap contributes
  transmission zeros, so a 2-band filter already carries inter-band TZs.
  Pro: reuses `MatrixSynthesisEngine` if the output is expressed as a
  `PolynomialSet`; consistent with the crate's model. The inter-band TZs
  are physical gap stopband zeros; relative to each band's prototype they
  map to `|Ω| ≥ 1` (E1's 2026-10-05 investigation corrected the earlier
  "passband zeros" claim), so route-α does **not** need passband-TZ
  support. (The C0 design doc's "E1 must be designed together with C"
  premise is withdrawn.)
- **C-route-β — step-impedance-resonator / comb-line ladder.** Classical
  multi-band topology; per-band impedance levels set by the band
  Chebyshev ladder; response is a rational function with known zeros.
  Pro: literature is deep and each band's behavior is independently
  verifiable. Con: the synthesis route is not the residue-expansion
  engine (it is impedance-level extraction), so it is a genuinely new
  `synthesis` path with its own verification surface.

### Work items (ordered)

1. **C0 — Design doc — DONE (2026-10-05).**
   `docs/internal/multi-band-synthesis-design.md` (written the way
   `docs/multiplexer-synthesis.md` was). **Route committed: C-route-α**
   (multi-band Chebyshev via frequency transformation), anchored on
   Macchiarella–Tamiazzo, *Design techniques for dual-passband filters*,
   IEEE TMTT 53(11), 2005, pp. 3265–3271 (DOI `10.1109/TMTT.2005.855749`);
   route-β deferred with rationale (it needs a new synthesis engine).
   The doc fixes: the characteristic function (product of per-band
   band-Chebyshev factors on the `MultiBandMapping` transform), the
   `MultiBandSpec`/`BandSpec`/`RejectionSpec` surface (invariants inherited
   `MultiBandMapping`), the **TZ placement law** (outer zeros are
   stopband — and, as E1's investigation established, so are inter-band
   zeros: every gap frequency maps to `|Ω| ≥ 1` on both adjacent bands,
   so route-α does **not** require passband TZs; the C0 doc's original
   "inter-band TZs are passband zeros" claim was wrong), the target
   topology ($N+2$ transversal; inter-band TZs become
   trisection/extracted-zero sections — no new `MatrixTopology`), and the
   **first acceptance fixture: `Cameron2011_10_5_dual_band_symmetric`
   (Cameron 2011 Fig. 25, order 10, 22 dB RL, 2 bands ±[0.35,1.0], 5 TZs)**,
   already in-repo at `tests/cameron2011_fixtures.json`.
2. **C1 — `MultiBandMapping` in `freq` — closed (2026-09-29).** Disjoint,
   center-sorted passbands; piecewise forward map (nearest band by
   prototype-transform distance `|Ω_b(f)|`), band-local inverse (the type
   does not implement `FrequencyMapping` because the trait's inverse has
   no band-index argument). Per-band round-trip and S4 edge-consistency
   tests green.
3. **C2 — Approximation + synthesis — DONE (2026-10-05).** route-α:
   `src/approx/multi_band.rs` produces a `PolynomialSet` from
   `multi_band_chebyshev_polynomials(bands, order, rl, tz)`; the
   positive-axis reflection zeros are solved by a damped
   Levenberg–Marquardt equiripple step (square system: `n_rz + 1`
   local maxima of `|K|` share one ripple level per positive interval),
   then the Feldtkeller `E(s)` / `PolynomialSet` assembly reuses the
   single-band helpers verbatim. `MatrixSynthesisEngine::synthesize`
   (C2c) consumes the output with **no changes** — no new engine. The
   C2a spec surface `MultiBandSpec`/`BandSpec`/`RejectionSpec` is in
   `src/spec/multi_band.rs`. 2-band first (the Fig. 25 fixture); **K ≥ 3 is
   infeasible for this characteristic (measured 2026-10-05)** — the
   equiripple solve is square only for a single positive-axis passband
   interval (`n_rz + 1` maxima sharing one level); a center band spanning
   the origin, or 2+ positive-axis bands, overdetermines the real-rooted
   characteristic by M − 1 (M = number of positive-axis intervals). The
   probe measured residuals 5.7e+2 / 4.2e-2 (M = 2/3) vs. the M = 1 control's
   1.7e-9, with the zeros pinned to the band edges. The solve rejects K ≥ 3
   with a typed `MfsError::NotImplemented` stating the DOF reason; a genuine
   K ≥ 3 equiripple needs a complex-rooted (Zolotarev) or multifunction
   characteristic (see the design doc's "K ≥ 3: the DOF Wall").
4. **C3 — Verification surface — DONE (2026-10-05).** `tests/multi_band.rs`
   (`dual_band_response_reproduces_cameron_fig25_fixture`) asserts, on a
   4001-point grid over `-2..2`: (1) **response invariance** — max|ΔS11|,
   max|ΔS21| < 1e-3 vs the fixture's printed matrix (the synthesized
   transversal and the fixture's `trisection_cascade_with_extracted_zero`
   are different topologies of the same network, so the gate is at the
   response level; measured 4.8e-4 / 1.7e-4); (2) **per-band equiripple** —
   worst in-band |S11| = −22 dB (the RL spec) ±1; (3) **unitarity** —
   max‖|S11|²+|S21|²−1‖ < 1e-10; (4) **transmission zeros** — |S21| < −40
   dB in a ±0.005 window at each spec TZ (measured −261.8, −60.4, −70.5
   dB). `examples/multi_band_response.rs` renders the response (per-band /
   gap / outer table + Touchstone .s2p).

**Acceptance (section done):** C0 doc reviewed and route committed;
C1 round-trip tests green; a 2-band design reproduces the C0 literature
fixture at the **response** level (max|ΔS11|, max|ΔS21| < 1e-3; per-band
ripple at the spec RL; TZ nulls < −40 dB; unitarity < 1e-10) — the
synthesized transversal and the fixture's printed trisection matrix are
different topologies of the same network, so coupling equality is **not**
the gate; `docs/architecture.md` module list and Known Limitations
updated (item "multi-band" removed from the roadmap table).

---

## D. Optimization — from "works at order 3" to a supported feature

**Goal:** make `tune_matrix` / `extract_matrix` trustworthy and reachable
from the CLI/Python surfaces, with integration-level evidence.

### Status (evidence)

- `src/optimize.rs`: `OptimizeConfig` (200 iters, tol 1e-8, LM λ=1e-3),
  `extract_parameters`/`rebuild_matrix` freeze the topology to the
  template's non-zero pattern (>1e-15), **scale-aware finite-difference
  Jacobian (`delta = 1e-7 * max(1, |x_i|)`, halving retry on singular
  perturbed evaluations)**, hand-rolled Gaussian elimination in `lm_step`,
  convergence on the *successive-iteration* cost (judged only after an
  accepted step — the old running-`best` comparison falsely terminated on
  iteration 2), a singular trial step is a **rejection** (λ×2), and
  `extract_matrix` starts from a **crude inline-chain** guess (the old
  star-like guess was unrepresentable for a transversal target).
- 3 lib unit tests, all order 3–4, well-posed. **2 `tests/optimizer_integration.rs`
  gates now GREEN (D2-opt, 2026-10-04): order-6 folded perturb→recover (6 iters
  → ~3e-10 relative) and order-6 transversal extract-from-crude-guess (9 iters
  → max|ΔS| ~5e-10).** **Exposed in the CLI (`optimize`/`extract` subcommands, D3-opt, 2026-10-04) and Python (module-level `extract`, `tune` since D1-opt), e2e-gated; benched (`benches/optimizer.rs`, tune o4/6/8 + extract o4/6) and noise-regression-gated (`tests/optimizer_regression.rs`, D4-opt, 2026-10-04).** Known limits (order-8 cold-start extraction, topology-bound recovery) documented in §6 of the report.
- Design/rationale: `docs/numerical-optimization-report.md` (the three
  numeric optimizations: O(n) Givens, selective linear solve, adaptive
  root finding — all shipped; **§6 now covers the LM optimizer's status and
  limits**).

### Work items (ordered)

1. **D1 (opt) — Integration test: synthesize → perturb → extract — RED gate committed (2026-09-29), GREEN (2026-10-04).**
   `tests/optimizer_integration.rs` (two order-6 gates, both now run and pass):
   (a) perturb four inter-resonator couplings by 7 %, re-optimize from the
   perturbed good-guess matrix, assert 1e-4 relative recovery — converges in
   6 iterations to ~3e-10; (b) `extract_matrix` from its crude guess —
   converges in 9 iterations to max|ΔS| ~5e-10.
2. **D2 (opt) — Scale-aware Jacobian step — DONE (2026-10-04).** The fixed
   1e-7 became `delta = 1e-7 * max(1, |x_i|)` with up to 4 halving retries when
   a perturbed evaluation returns `Err` (singular). **Honest attribution:
   the scale-aware δ was *not* the primary fix.** Scenario (a) was broken by
   the convergence check comparing cost against the running *best* (a false
   "no change" on iteration 2 — it stopped at ~8.3e-3 regardless of δ, and
   fixed/scale/central differences all failed identically); it is now judged
   on the successive-iteration cost after an accepted step. Scenario (b) was
   broken by `extract_matrix`'s star-like initial guess (every resonator
   coupled to both ports, chain couplings frozen at zero — a template that
   cannot represent a transversal inline-chain target), now a crude inline
   chain; and a trial step whose response evaluation failed (singular)
   propagated out as an error instead of being a λ×2 rejection. Both gates
   are un-`#[ignore]`d and green; the hand-rolled elimination is kept (not a
   bottleneck — 2 gates in ~1 s of test time).
3. **D3 (opt) — CLI + Python surface — DONE (2026-10-04).** `mfs_cli.rs` gains
   two subcommands, both delegating to the shared optimizer core (`one
   canonical entry, thin wrappers`):
   - **`optimize`** — input JSON `{order, topology, data, samples}` (the
     matrix as a flat row-major `data`, the target as `samples` with
     `frequency_hz`/`s11_re/im`/`s21_re/im`) or a Touchstone `.s2p` via
     `--target` (the matrix must then still come from the JSON); knobs
     `--order`, `--max-iterations`, `--tolerance`, `--lambda`, `--q`,
     `--output`. Prints `{order, topology, matrix, data, cost, iterations,
     converged, recovered:{samples}}` (or a table).
   - **`extract`** — input JSON `{samples:[...]}` (or a `.s2p` via `--target`),
     `--order`, `--output`; calls `extract_matrix` directly.
   - `src/python.rs` gains a module-level **`extract(freqs, s21_re, s21_im,
     s11_re, s11_im, order)`** returning `{matrix, cost, iterations,
     converged}`; it delegates to `extract_matrix`. (`PyFilterDesign.tune`
     already wrapped `optimize_matrix` — D1-opt — so no new tune surface was
     needed.)
   - End-to-end black-box regression in `tests/cli_integration.rs`
     (`#[cfg(feature = "cli")]`): `test_optimize_subcommand_recovers_perturbed_matrix`
     (order-4 folded, 3 couplings ×1.05 → recovered to < 1e-4 relative) and
     `test_extract_subcommand_recovers_response` (order-4 transversal →
     max|ΔS| < 1e-3), both asserting physical quantities, plus a `#[cfg(test)]`
     `optimizer_tests` module in the bin covering the JSON parsers (flat +
     nested) and the two optimizer cores.
4. **D4 (opt) — Bench + regression — DONE (2026-10-04).**
   - **Bench:** `benches/optimizer.rs` (Criterion, registered in `Cargo.toml`,
     `harness = false`). Times `tune` at orders 4/6/8 and `extract` at orders
     4/6 (each run 2.8–17 ms warm; measured cost reduction and iteration count
     printed alongside: tune 6/6/6 iters → worst_rel ≈ 3e-10/1e-11/7e-12;
     extract 12/9 iters → max|ΔS| ≈ 2e-12/5e-10). **Deliberately outside the
     cross-language standard** (no `optimizer` family in `gates.json` —
     inventing case ids/tolerances without provenance is forbidden, same
     stance as synthesis.rs's `SWEEP_CASES`). **Order-8 `extract` is not
     timed**: measured from the cold chain guess it converges to a local
     minimum far from the target (cost 39.2, max|ΔS| ≈ 1.0 after 435 iters /
     ~85 s) — timing it would measure noise. See §6.3 of the report.
   - **Regression:** `tests/optimizer_regression.rs` — order-6 extraction from
     a 0.1 dB-noise-floor response (deterministic LCG). Asserts convergence,
     fidelity to the **truth** (max|ΔS| < 5e-3; measured 1.22e-3 against a
     1.0e-2 noise floor) and that it does **not** chase the noise (residual vs.
     noisy data < 4× the noise floor). Couplings are not asserted: the
     transversal parameterisation is redundant, so the physical contract is
     the response, not the entries.
   - Report: `docs/numerical-optimization-report.md` gains **第六部分** —
     "耦合矩阵 LM 优化器 —— 状态与边界" (algorithm, gate evidence, and the four
     documented limits: topology-bound recovery, order-8 cold-start, response
     vs. matrix uniqueness, FD-Jacobian cost).

**Acceptance (section done):** D1-opt green; D2-opt shows no worse (and
typically better) convergence than fixed-δ on the D1-opt suite; CLI and
Python entry points covered by `tests/cli_integration.rs`
(`test_optimize_subcommand_recovers_perturbed_matrix`,
`test_extract_subcommand_recovers_response`) and `tests/optimizer_regression.rs`;
`docs/numerical-optimization-report.md` §6 "耦合矩阵 LM 优化器 —— 状态与边界".

---

## E. Cross-cutting numeric enablers (blockers for A/B, and C via route-α)

Ordered by how much they block the areas above:

1. **Passband transmission zeros (`|zero| < 1`) — investigated 2026-10-05;
   the planned work item is mathematically impossible as specified.**
   The plan's route ("relax `safe_sqrt_term`" +
   `generalized_ops.rs:96-105`) was checked end-to-end and fails:
   - `P(w) = ∏(w − z_k)` always has **real** roots, so a zero at `|z| < 1`
     *does* impose a finite real-frequency null (`|S21| = 0` at `w = z`).
     That part of the premise is fine.
   - But the literal relaxation — making `√(1 − 1/z²)` imaginary — makes
     `F(s)` complex **non-Hermitian** (`F(jω) ∈ ℂ`), whereas the physical
     characteristic function `K(jω) = F/P` must be **real** for the
     Feldtkeller identity `|E|² = |F|² + |P|²/ε²` to hold. Measured
     (order-4, RL 20 dB, `z = 0.5` versus `z = 1.5`): the identity
     residual is O(1) at every ω for `|z| < 1`, versus ≤ 1e-13 for
     `|z| ≥ 1`. The relaxed recurrence is not a valid prototype.
   - The design doc's claim that inter-band gap frequencies are `|Ω| < 1`
     per band is **false**. For adjacent passbands `[f1,f2]`, `[f3,f4]`, a
     band's `|Ω| = 1` boundary *is* its own edge, so every gap frequency
     maps to `|Ω| ≥ 1` on both bands: for 10–11 GHz / 12–13 GHz, the
     11–12 GHz gap maps to `Ω ∈ [+1.00, +2.83]` on the lower band and
     `[−3.18, −1.00]` on the upper — ordinary **stopband** zeros.
   **Consequence:** an equiripple prototype with an in-band null requires
   the *multifunction* characteristic (a ∏-form `K(s)` whose in-band
   extrema are not the roots of `F`), not the single-function generalized
   Chebyshev `F/P`. That is genuine new approximation research, not a
   helper relaxation, and it is **not** the mechanism multi-band route-α
   needs (gap zeros are ordinary stopband zeros). Spec-layer decision:
   `FilterSpec` stays permissive, and `approx` now rejects `|zero| < 1`
   with a typed error (`MfsError::InvalidTransmissionZero`) at the
   entry point instead of leaking `NotImplemented` from a helper.
   `|zero| = 1` (imaginary-axis ripple zero) is likewise rejected.
2. **Repeated-pole residue expansion (Laurent) — closed (2026-10-05).**
   `ResiduePole` gained `laurent: Option<Vec<ComplexCoefficient>>` (the
   principal part `b_0..b_{m-1}`, `residue` = `b_0`); `residue_expansion`
   now clusters coalesced roots (single-linkage, 1e-5 — the solver splits a
   multiplicity-m root into m roots within ~1e-7–1e-5, while distinct filter
   poles are ≥ ~1.6e-3 apart) and takes the **multiplicity `m` as the cluster
   size** (not re-derived from the pole-shifted coefficient vector, which is
   ill-conditioned for high-degree denominators and miscounted order-16
   simple poles as repeated); the principal part is computed by power-series
   division of the shifted numerator, then leading near-zero coefficients are
   trimmed (numerator–denominator cancellation lowers the effective pole
   order). The pole-expansion fast path sums the full principal part
   (`PoleExpansionData.evaluate_at`), summing each Y parameter with its own
   principal part — a shared pole's effective order can differ between
   y11/y12/y22 after cancellation, so the three need not be equal-length.
   `build_transversal_from_residues` rejects a repeated pole with a typed
   `MfsError::NotImplemented` — a transversal chain has one simple pole per
   resonator and cannot represent `1/(s-p)^m`. Gates:
   `expansion::tests::double_pole_laurent_principal_part` (hand-computed
   coefficients `-2/3, 7/9, 2/9`),
   `transversal::tests::repeated_pole_rejected_with_not_implemented`,
   `pole_expansion::tests::laurent_principal_part_contributes_to_sum`,
   `pole_expansion::tests::repeated_pole_with_cancellation_uses_each_y_own_order`
   (per-Y principal-part length), plus the `property_2` residue round-trip
   (now covers the Laurent form).
   Valid filter specs have distinct poles, so this path is only exercised by
   degenerate/synthetic inputs.
3. **HighPass / BandStop mappings — closed (2026-09-29).** `HighPassMapping`
   and `BandStopMapping` are implemented in `freq` with the shared
   `SINGULAR_CLAMP` convention at their singular points. Unblocks
   stop-band-anchored multi-band specs (C1).
4. **High-order integration tests (N=8 fixture band) — closed (2026-09-29).**
   `tests/benchmark_standard.rs` gained
   `mfs_synthesized_folded_order8_reproduces_the_printed_book_matrix`:
   the standard's only order-8 case (`twoport_8_3`, Cameron 2018 Sec 9.6 /
   Fig 9.22a) is synthesized and folded, and its matrix + response are
   pinned against the printed `EIGHT_THREE_FOLDED` anchor at the 4-decimal
   printing quantum.

---

## Suggested execution order

```
Wave 1 (parallel):  E4 (N=8 regression) · D1-opt (RED)   [E3 + C1 + E4 closed; D1-opt RED gate committed 2026-09-29, GREEN 2026-10-04]
Wave 2 (serial):    D1–D4 (duplexer full path)   [DONE 2026-10-03 — Type-I orders 2–8 + §V golden + Type-II §VI anchor, all gate-anchored]
Wave 3 (parallel):  D2-opt [DONE 2026-10-04 — scale-δ + early-stop + chain-guess, both order-6 gates green] · D3-opt [DONE 2026-10-04 — CLI optimize/extract + Python extract, e2e-gated] · D4-opt [DONE 2026-10-04 — `benches/optimizer.rs` + noise regression + report §6] · C0 [DONE 2026-10-05 — design doc, route-α committed]
Wave 4 (serial):    B1 [DONE 2026-10-04] → B2 [DONE 2026-10-04] → B3 [DONE 2026-10-04] → B4 [DONE 2026-10-04] (multiplexer; Wave 2 done — full wave complete)
Wave 5:             E1 [INVESTIGATED 2026-10-05 — refuted: no valid passband-TZ prototype in the single-function form; approx now rejects |z|<1 typed; route-α does NOT need it] · C2 [DONE 2026-10-05 — route-α approx + spec + engine reuse] → C3 [DONE 2026-10-05 — Fig. 25 anchor GREEN + example] · E2 (Laurent) [DONE 2026-10-05 — Laurent principal part + fast-path summation; transversal rejects repeated pole typed; gates green]
Wave 6 (review):    0.6.0 whole-crate review "silent errors" wave (review-findings-0.6.0.md §F.1) [DONE 2026-10-06 — P1-1/2/3/4/6 + P2-1/2/3/5 + all six P3 fixed and gate-anchored] · CT/CQ cascade (§F.2) [DONE 2026-10-06 — P1-5 + P2-4 fixed: `to_cascaded_trisections` via the edge-capable `extract_triplet` (reproduces Cameron 2011 Fig. 21), `to_cascaded_quadruplets` via the verified `extract_cascade_quadruplets` with a **source-side non-overlapping tiling** (pair *i* spans resonators `4i-3..=4i`, `order >= 4·num_pairs` enforced) — clean order-6 single quadruplet, reproduces the verified `CascadeQuartets` `[[1,3],[1,4],[5,7],[5,8]]` at order 8, over-constrained order-6 2-pair rejected; gated by Fig-21 trident + value anchor, order-6 mainline, order-4/order-6/order-8 quadruplet patterns, CQ≠CT, order-rejection]. Note: an earlier same-day RE-SCOPED note claimed the fix was REFUTED because `extract_trisection` (internal-center-only) cannot produce edge trisections — that tested the *wrong* primitive; the fix uses the distinct `extract_triplet` (center=1 allowed), which does reproduce Fig. 21 (see review-findings-0.6.0.md P2-4 Status, both bullets). The earlier "order-6 two-pair ~1e-2 residual" follow-up is resolved as **topologically impossible** (two 4-resonator quadruplets cannot tile 6 resonators non-overlappingly), not a fixable bug — now rejected with a typed error.
```

Rationale: the duplexer full path (Wave 2) was the only item on the critical
path and is now complete (2026-10-03); the multiplexer's high-fidelity work
(B2/B3) reuses it, and porting a *converged* loop keeps the debug surface
tractable. Everything in Waves 1, 3, and 4 (beyond B2's porting) is
independent and safe to parallelize.

---

## Doc hygiene — fixed 2026-09-28 (with this plan)

- `docs/multiplexer-synthesis.md` acceptance item 6 claimed "3 pre-existing,
  unrelated lib failures remain" — stale; lib suite is fully green
  (202 passed, 0 failed).
- `docs/duplexer-synthesis.md` "Full Algorithm Infrastructure" named
  `wilson_spectral_factor` and `fit_polynomial_ls` (deleted in 0.5.0 dead-code
  removal) and "Future Work" listed Type-II junction (implemented,
  `JunctionType::TypeII`).
- `docs/architecture.md` claimed `spec` enforces "Transmission zero
  magnitudes ≥ 1.0" — it does not; the rejection happens at the
  approximation layer.
- `docs/internal/todo.md` baseline head updated to `2d735c0`.
- `docs/internal/index-for-filter-synthesis-library-design.md` now points
  at this document as the continuation plan.
