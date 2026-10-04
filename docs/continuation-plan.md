# MFS Continuation Plan

**Date:** 2026-09-28 · **Baseline head:** `2d735c0` (`master`, pushed to gitlab)
**Status:** all CI-parity checks green on baseline (lib 202 + 15 integration
binaries, clippy `-D warnings`, rustdoc — zero warnings).

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
| Multiplexer (N-channel) | Independent per-channel synthesis + under-relaxed perturbation; shared `synthesis::recover_channel_design`; response is a **geometric-mean magnitude model**, not an (N+1)-port solve | Medium (channels are valid passband filters; isolation/antenna values are model approximations) |
| Multi-band | **Mapping layer in place** — `freq` supports all four prototype families plus `MultiBandMapping` (E3/C1, 2026-09-29); no design doc, no multi-band approximation, no synthesis path | — |
| Optimization | LM `tune_matrix` / `extract_matrix` (hand-rolled Gaussian elimination, fixed 1e-7 FD Jacobian, topology frozen to template); 3 lib unit tests at order 3–4 | Low (works on small well-posed cases; no integration coverage, no CLI/Python surface) |

Known limitations that gate the areas below (full list in
`docs/architecture.md` "Known Limitations"):

1. **Passband transmission zeros** (`|zero| < 1`) — `generalized_ops.rs`
   rejects with `NotImplemented`.
2. **Repeated poles** — `residues/expansion.rs` returns `NotImplemented`.
3. ~~**HighPass / BandStop mappings**~~ — implemented in `freq` (E3, 2026-09-29: `HighPassMapping`, `BandStopMapping` with full forward + inverse and singular-point clamps).
4. **N-channel mutual-loading fidelity** — the duplexer's full spectral path
   converges (Type-I orders 2–8, Type-II, both gate-anchored; 2026-10-03), but
   it is not yet generalized to N channels (multiplexer B2), and the
   multiplexer's default still uses the simplified under-relaxed perturbation.

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
- **Response is a magnitude model:** antenna S11 = geometric mean
  `(∏ᵢ|S11ᵢ|)^(1/N)`, isolation[i][j] = product of S21 magnitudes
  (`multiplexer.rs` N+1 analysis). No complex (N+1)-port solve, no
  common-port admittance sum `Y_in = Σ 1/Z_in,i`.
- **No full-spectral path** for N-channel (the full path exists only in
  `duplexer.rs`; it is converged and gate-anchored there — section A, done
  2026-10-03 — but not yet generalized to N > 2).
- **Benchmark-standard gates:** the vendored standard's `cm`, `star`,
  `s11_match`, `rz_distribution` gates have no mfs counterpart and are
  `SKIP`; the convergence gate is a documented FAIL
  (`benches/multiplexer.rs`, `bench-data/gates.json`). These are the
  quantitative definition of "not high-fidelity yet".

### Work items (ordered; depends on A)

1. **B1 — True (N+1)-port response (magnitude model out).** Implement the
   common-port admittance model from `docs/multiplexer-synthesis.md`
   §"New API Surface": `Y_in = Σᵢ 1/Z_in,i`, complex-valued antenna S11 and
   channel-to-channel S-parameters. Keep the geometric-mean model behind a
   named flag for one release cycle; tests assert the two agree at N=2 to
   within ~0.1 dB where the model is known to be accurate. Acceptance:
   unitary-check test — for lossless inputs the (N+1)×(N+1) S-matrix is
   power-conserving to 1e-9 at sampled frequencies.
2. **B2 — Generalize the full spectral path to N channels.** Port
   `duplexer.rs`'s full-path helpers (`solve_scaling_constants`, `compute_s_product`,
   etc.) into a shared module (candidate: `src/synthesis/macchiarella.rs`) so the
   N=2 duplexer and the N-channel multiplexer share one iteration loop.
   This subsumes the todo.md item "Fold the duplexer synthesis into the
   multiplexer (N=2)". Must happen *after* D3/D4 (done) — porting the now
   converged path keeps the N > 2 debug surface tractable.
3. **B3 — Convergence + fidelity gates for N = 3, 4.** Extend the
   `separation_limit` / `adjacent_channels_show_loading` test pair
   (`tests/multiplexer.rs`) to the full path: well-separated channels
   converge to independent synthesis to 1e-6 (machine tolerance, unlike the
   simplified model's documented residual), adjacent channels show bounded
   loading with antenna RL quasi-equiripple within ~2 dB. This flips the
   vendored-standard convergence gate from FAIL to PASS.
4. **B4 — Close the bench-standard gaps.** Implement what mfs can honestly
   implement of the `cm`/`star` gates (per-channel coupling-matrix checks,
   star-domain response) and mark remaining gaps with reasons in
   `benches/multiplexer.rs` instead of blanket `SKIP`.

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
- No multi-band approximation, no multi-band synthesis path.
- The vendored benchmark standard has **no** multi-band case family —
  there is no external gate; literature fixtures will have to define the
  acceptance bar (the fixture tool is `tests/extract_from_pdf.py`).

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
  `PolynomialSet`; consistent with the crate's model. Con: the
  inter-band TZs are *passband* zeros relative to the per-band
  prototypes — this is the direct, large-scale application of E1
  (passband TZ handling), which must therefore be designed together
  with C, not after it.
- **C-route-β — step-impedance-resonator / comb-line ladder.** Classical
  multi-band topology; per-band impedance levels set by the band
  Chebyshev ladder; response is a rational function with known zeros.
  Pro: literature is deep and each band's behavior is independently
  verifiable. Con: the synthesis route is not the residue-expansion
  engine (it is impedance-level extraction), so it is a genuinely new
  `synthesis` path with its own verification surface.

### Work items (ordered)

1. **C0 — Design doc (gate for everything else).** Write
   `docs/internal/multi-band-synthesis-design.md` *before any code*, the
   way `docs/multiplexer-synthesis.md` was: the chosen route (α or β) with
   the math, the spec structure (`MultiBandSpec`: per-band
   `{center_hz, bandwidth_hz, return_loss_db}` + out-of-band rejection
   requirements), the TZ placement law (which zeros in which gap, and the
   per-band passband-zero analysis that E1 must cover), the target matrix
   topology, and a 2-band reference design from the literature as the
   first acceptance fixture.
2. **C1 — `MultiBandMapping` in `freq` — closed (2026-09-29).** Disjoint,
   center-sorted passbands; piecewise forward map (nearest band by
   prototype-transform distance `|Ω_b(f)|`), band-local inverse (the type
   does not implement `FrequencyMapping` because the trait's inverse has
   no band-index argument). Per-band round-trip and S4 edge-consistency
   tests green.
3. **C2 — Approximation + synthesis (route-dependent).** Per C0:
   route-α → new family in `approx` (candidate `approx/multi_band.rs`)
   producing a `PolynomialSet`-shaped output so `MatrixSynthesisEngine`
   is reused; route-β → new impedance-level extraction path in
   `synthesis`. Either way: 2-band first (3-band is the extension
   check), with the literature fixture from C0 as the RED→GREEN anchor.
4. **C3 — Verification surface.** Response invariance checks per band
   (each band's response must match its standalone equiripple spec to
   within the ripple tolerance), stopband rejection assertions in every
   gap, and an example in `examples/` rendering the multi-band response.

**Acceptance (section done):** C0 doc reviewed and route committed;
C1 round-trip tests green; a 2-band design reproduces the C0 literature
fixture (all couplings within 1e-4 relative, per-band ripple within the
spec); `docs/architecture.md` module list and Known Limitations updated
(item "multi-band" removed from the roadmap table).

---

## D. Optimization — from "works at order 3" to a supported feature

**Goal:** make `tune_matrix` / `extract_matrix` trustworthy and reachable
from the CLI/Python surfaces, with integration-level evidence.

### Status (evidence)

- `src/optimize.rs`: `OptimizeConfig` (200 iters, tol 1e-8, LM λ=1e-3),
  `extract_parameters`/`rebuild_matrix` freeze the topology to the
  template's non-zero pattern (>1e-15), finite-difference Jacobian with a
  **fixed `delta = 1e-7`**, hand-rolled Gaussian elimination in `lm_step`,
  `extract_matrix` starts from a crude uniform-transversal guess.
- 3 lib unit tests, all order 3–4, well-posed. **No `tests/` integration
  coverage, no bench, not exposed in the CLI or Python bindings.**
- Design/rationale: `docs/numerical-optimization-report.md` (the three
  numeric optimizations: O(n) Givens, selective linear solve, adaptive
  root finding — all shipped; that doc covers performance, not the LM
  optimizer's roadmap).

### Work items (ordered)

1. **D1 (opt) — Integration test: synthesize → perturb → extract — RED gate committed (2026-09-29).**
   `tests/optimizer_integration.rs` (two `#[ignore]`d RED tests, order 6):
   (a) perturb four inter-resonator couplings by 7 %, re-optimize from the
   perturbed good-guess matrix, assert 1e-4 relative recovery — currently
   stops at ~8.3e-3; (b) `extract_matrix` from its crude uniform guess —
   currently wanders into a singular response solve. Both `#[ignore]`d with
   the failure documented; they are the RED→GREEN gate for the rest of this
   section and the scale-sensitivity exposure of the fixed-δ Jacobian that
   D2-opt is meant to fix.
2. **D2 (opt) — Scale-aware Jacobian step.** Replace the fixed 1e-7 with
   `delta = 1e-7 * max(1, |x_i|)` (or central differences if the FD cost
   stays acceptable at our problem sizes — measure first, don't assume).
   Keep the hand-rolled elimination only if the D1-opt integration tests
   show it is not the bottleneck; otherwise move the normal equations to
   `nalgebra` (already a dependency) — one line, removes a bug surface.
3. **D3 (opt) — CLI + Python surface.** `mfs_cli.rs` `optimize`
   subcommand (input matrix JSON/Touchstone + target response JSON or
   Touchstone → optimized matrix, same `--format` conventions as the
   existing subcommands) and a thin `PyFilterDesign.optimize(...)` /
   module-level `extract` in `src/python.rs`. Both delegate to the same
   `optimize_matrix` core (the `pipeline::context_from_request` pattern —
   one canonical entry, thin wrappers).
4. **D4 (opt) — Bench + regression.** Criterion bench (tune and extract at
   orders 4/6/8, measuring iterations and cost reduction); a regression
   entry in `tests/synthesis_regression.rs`-style file for an
   extract-from-noise case (measured-style S-data with 0.1 dB noise floor)
   so the "extraction" contract is pinned beyond the well-posed cases.

**Acceptance (section done):** D1-opt green; D2-opt shows no worse (and
typically better) convergence than fixed-δ on the D1-opt suite; CLI and
Python entry points covered by `tests/cli_integration.rs` and the python
doctests; `docs/numerical-optimization-report.md` gains a short
"Optimizer status & limits" section.

---

## E. Cross-cutting numeric enablers (blockers for A/B, and C via route-α)

Ordered by how much they block the areas above:

1. **Passband transmission zeros (`|zero| < 1`).** Blocks high-fidelity
   duplexer/multiplexer designs where TZs land in-band after the N-channel
   perturbation shifts the roots — and is a hard dependency of multi-band
   route-α (section C), whose inter-band TZs are passband zeros relative to
   the per-band prototypes. Work: complex-zero handling in
   `generalized_ops.rs` (relax `safe_sqrt_term`), plus the spec-layer
   decision — today `FilterSpec` *accepts* `|z| < 1` (see
   `spec/types.rs` test storing 0.75) and the rejection is deferred to the
   approximation layer; the plan is to keep spec permissive and fail at
   `approx` with a typed error, then extend.
2. **Repeated-pole residue expansion (Laurent).** Blocks any design whose
   perturbed root set degenerates — a realistic risk once the full
   spectral path actually moves roots (A/B). Work in
   `residues/expansion.rs`.
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
Wave 1 (parallel):  E4 (N=8 regression) · D1-opt (RED)   [E3 + C1 + E4 closed; D1-opt RED gate committed 2026-09-29]
Wave 2 (serial):    D1–D4 (duplexer full path)   [DONE 2026-10-03 — Type-I orders 2–8 + §V golden + Type-II §VI anchor, all gate-anchored]
Wave 3 (parallel):  D2–D4-opt (optimizer) · C0 (multi-band design doc + route choice)
Wave 4 (serial):    B1 → B2 → B3 → B4 (multiplexer; Wave 2 done — unblocked)
Wave 5:             E1 (passband TZ) · E2 (Laurent) — needed by B3's fidelity work and
                    by C2 when route-α is chosen (inter-band TZs); then C2 → C3
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
