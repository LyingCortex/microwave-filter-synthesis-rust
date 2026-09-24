# Changelog

## [Unreleased]

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
- **Known limitation** (shared with `DuplexerSynthesis`): the perturbation
  path's `recover_filter_design` reuses the original `D` polynomial with the
  converged `S` roots, yielding fully-reflective channels (peak S21 = `-inf`).
  The design doc and an `#[ignore]`d diagnostic test document this; a real fix
  requires re-deriving `D` by interpolation at the converged roots via the full
  Macchiarella spectral path (deferred as polish).

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
