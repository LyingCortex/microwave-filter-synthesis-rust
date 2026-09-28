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
| Duplexer | Macchiarella–Tamiazzo **simplified** frequency-offset perturbation; full spectral path exists behind `DuplexerSynthesis::full_iteration(true)` but is **untested and non-convergent at order ≥ 4** | Medium (model is approximate; antenna-RL quasi-equiripple not asserted) |
| Multiplexer (N-channel) | Independent per-channel synthesis + under-relaxed perturbation; shared `synthesis::recover_channel_design`; response is a **geometric-mean magnitude model**, not an (N+1)-port solve | Medium (channels are valid passband filters; isolation/antenna values are model approximations) |
| Multi-band | **Nothing** — `freq` supports LowPass + BandPass only; no design doc, no mapping, no synthesis path | — |
| Optimization | LM `tune_matrix` / `extract_matrix` (hand-rolled Gaussian elimination, fixed 1e-7 FD Jacobian, topology frozen to template); 3 lib unit tests at order 3–4 | Low (works on small well-posed cases; no integration coverage, no CLI/Python surface) |

Known limitations that gate the areas below (full list in
`docs/architecture.md` "Known Limitations"):

1. **Passband transmission zeros** (`|zero| < 1`) — `generalized_ops.rs`
   rejects with `NotImplemented`.
2. **Repeated poles** — `residues/expansion.rs` returns `NotImplemented`.
3. **HighPass / BandStop mappings** — not implemented in `freq`.
4. **N-channel mutual-loading fidelity** — simplified perturbation does not
   converge to independent synthesis for well-separated channels; the full
   spectral path is not converged at order ≥ 4.

---

## A. Duplexer — converge the full Macchiarella–Tamiazzo path

**Goal:** make the full spectral iteration (`full_iteration(true)`) converge
for order ≥ 4, so the duplexer can deliver true quasi-equiripple antenna
return loss instead of the simplified perturbation's approximation.

### Status (evidence)

- Simplified path (default): `src/duplexer.rs` — frequency-offset perturbation
  loop; channels converge to valid passband filters; antenna S11 is the
  geometric-mean model `√(|S11_tx|·|S11_rx|)`.
- Full path: `compute_n_polynomial` → `solve_scaling_constants` →
  `spectral_factor` (cepstrum + balanced-companion fallback in
  `src/spectral.rs`) → `compute_s_product` → `assign_roots`. Wired behind
  `DuplexerSynthesis::full_iteration(true)` (default `false`).
  **No test, example, or bench anywhere invokes it.** It returns a typed
  error at order ≥ 4 and exhausts iterations at order 2
  (`docs/architecture.md` Known Limitations #5).

### Work items (ordered)

1. **D1 — Make the full path testable and failing-visible.** Add lib tests
   that call `full_iteration(true)` at order 2 (TX/RX of order 2 each):
   assert convergence (`converged == true`), quasi-equiripple antenna RL
   within ~2 dB of the independent-filter bound, and passband peaks within
   0.5 dB of the per-channel spec. Expected RED at the current state — this
   is the RED gate for everything else in this section.
2. **D2 — Diagnose the order ≥ 4 non-convergence.** Instrument the
   iteration (residual per iteration, root-displacement via
   `compute_root_displacement`, spectral-factor relative error from
   `spectral.rs`) in a diagnostic example. Root-cause candidates to check,
   in order: (a) cepstrum spectral factorization mis-assigns LHP roots at
   higher degree (verify `select_lhp_roots` against a known factorization),
   (b) `solve_scaling_constants` ill-conditioning as the p₀ system grows,
   (c) under-relaxation `alpha` schedule too aggressive for the N(s) root
   movement. Fix the root cause, not the tolerance.
3. **D3 — Converge order ≥ 4 and pin it.** Extend the D1 tests to order 4
   and order 6 (per channel); assert the same fidelity band. Only after D3,
   consider changing the default to `full_iteration(true)` — that is a
   separate, deliberate behavior change with its own changelog entry.
4. **D4 — Document the as-built fidelity.** Update
   `docs/duplexer-synthesis.md` Future Work (currently lists Type-II
   junction as future — implemented at `src/duplexer.rs` JunctionType;
   and names `wilson_spectral_factor` / `fit_polynomial_ls`, deleted in
   0.5.0) to state what the full path achieves and at what cost.

**Acceptance (section done):** `cargo test --lib -- duplexer` green with
full-path tests at orders 2/4/6; `docs/duplexer-synthesis.md` accurate;
no regression on the simplified path (its tests unchanged).

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
  `duplexer.rs`, and is itself not yet converged — see section A).
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
   `duplexer.rs`'s full-path helpers (`compute_n_polynomial` et al.) into a
   shared module (candidate: `src/synthesis/macchiarella.rs`) so the
   N=2 duplexer and the N-channel multiplexer share one iteration loop.
   This subsumes the todo.md item "Fold the duplexer synthesis into the
   multiplexer (N=2)". Must happen *after* D3 — porting a non-convergent
   path multiplies the debug surface.
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

## C. Multi-band — design first, then scope

**Goal:** filters with more than one passband. This is the largest planning
gap of the four areas: no code, no design doc, one roadmap line.

### Status

- `src/freq.rs` implements `LowPassMapping` and `BandPassMapping` only.
  `BandPassMapping`'s inverse takes the positive branch
  (`freq.rs` `map_normalized_to_hz`).
- No multi-band mapping, no multi-band approximation, no tests.
- The vendored benchmark standard has **no** multi-band case family, so
  there is no external gate to aim at yet.

### Scope decision (needs maintainer input before code)

Multi-band has two fundamentally different targets, and the plan branches
on which one is wanted:

- **C-scope-1 — "multi-band" = several passbands in one filter** (true
  multi-band, e.g. comb-line or step-impedance resonator response). This
  requires a new approximation family — generalized Chebyshev cannot
  represent multi-passband specs. This is a research-grade item and
  belongs in `approx` as a new family (the todo.md
  "Butterworth/Elliptic families" slot is the right home for a second
  family too).
- **C-scope-2 — "multi-band" = design convenience for several unrelated
  bands** (design N separate filters, cascade/combine their responses).
  This is pure composition of existing machinery: per-band
  `FilterDesign` + a `freq` multi-grid + response combination. Small.

**Recommended sequence (defaulting to C-scope-2, C-scope-1 deferred):**

1. **C1 — Multi-band frequency mapping.** `MultiBandMapping` in `freq`:
   a list of `(center_hz, bandwidth_hz)` bands; forward map is piecewise
   (each physical band edge maps to ±1 in its own normalized domain),
   inverse is band-local. Tests: round-trip per band, band-edge
   consistency with the existing `edges_to_center_bandwidth` geometric-mean
   convention (S4 invariant).
2. **C2 — Multi-band response composition.** A
   `MultiBandResponse` helper (in `response` or a new `src/multiband.rs`)
   that evaluates per-band `FilterDesign`s on a shared physical grid and
   composes them; document the composition law (independent filters →
   sum of through-couplings at the common port, mirroring B1's model).
3. **C3 — Design doc for C-scope-1** (only if the maintainer wants true
   multi-band synthesis): a `docs/internal/multi-band-synthesis-design.md`
   covering the resonator topology, the approximation family, and the
   synthesis route — written before any code, the way
   `docs/multiplexer-synthesis.md` was.

**Acceptance (C1+C2 done):** `MultiBandMapping` round-trip tests green;
a 2-band example in `examples/` produces per-band specs that match their
standalone `FilterDesign` counterparts; `docs/architecture.md` module list
updated.

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

1. **D1 (opt) — Integration test: synthesize → perturb → extract.**
   Generate an order-6 folded matrix, evaluate its response, perturb 3–5
   entries by 5–10 %, run `extract_matrix` with a *good* initial guess (the
   unperturbed matrix) and `tune_matrix` from the crude guess; assert
   recovered couplings within 1e-4 relative. This is the RED→GREEN gate for
   the rest of the section and exposes the fixed-δ Jacobian's scale
   sensitivity (1e-7 is only safe when couplings are O(1)).
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

## E. Cross-cutting numeric enablers (blockers for A/B)

Ordered by how much they block the areas above:

1. **Passband transmission zeros (`|zero| < 1`).** Blocks high-fidelity
   duplexer/multiplexer designs where TZs land in-band after the N-channel
   perturbation shifts the roots. Work: complex-zero handling in
   `generalized_ops.rs` (relax `safe_sqrt_term`), plus the spec-layer
   decision — today `FilterSpec` *accepts* `|z| < 1` (see
   `spec/types.rs` test storing 0.75) and the rejection is deferred to the
   approximation layer; the plan is to keep spec permissive and fail at
   `approx` with a typed error, then extend.
2. **Repeated-pole residue expansion (Laurent).** Blocks any design whose
   perturbed root set degenerates — a realistic risk once the full
   spectral path actually moves roots (A/B). Work in
   `residues/expansion.rs`.
3. **HighPass / BandStop mappings.** Independent of A/B but cheap;
   `freq` only. Needed before multi-band (C1) if stop-band-anchored specs
   are wanted.
4. **High-order integration tests (N=8 fixture band).** The documented
   fixture band is order 2–8 but `tests/` has no N=8 regression entry;
   add one (root-solver unit tests already cover 10/20/30).

---

## Suggested execution order

```
Wave 1 (parallel):  E3 (HP/BS maps) · E4 (N=8 regression) · D1-opt (RED)
Wave 2 (serial):    D1–D4 (duplexer full path)  ← critical path for B
Wave 3 (parallel):  D2–D4-opt (optimizer) · C1–C2 (multi-band scope-2)
Wave 4 (serial):    B1 → B2 → B3 → B4 (multiplexer; needs Wave 2)
Wave 5:             E1 (passband TZ) · E2 (Laurent) — as B3's fidelity work
                    starts to need them; C3 (scope-1 design doc) if chosen
```

Rationale: the duplexer full path (Wave 2) is the only item on the critical
path — the multiplexer's high-fidelity work (B2/B3) reuses it, and porting
a non-convergent loop (Wave 4 before Wave 2) would multiply debug surface.
Everything in Waves 1 and 3 is independent of it and safe to parallelize.

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
