# Changelog

## [Unreleased]

### Added

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

### Fixed

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
