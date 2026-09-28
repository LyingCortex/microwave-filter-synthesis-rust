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

### Fixed

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
- **Durand-Kerner and Aberth-Ehrlich no longer return NaN root sets as
  "converged".** Their in-loop convergence check compared `max_delta`
  (accumulated with `f64::max`, which swallows NaN) against the tolerance,
  so an iteration that overflowed to NaN read as converged with
  `max_delta = 0.0`. Both solvers now also require all roots to be finite
  before the early return.

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
