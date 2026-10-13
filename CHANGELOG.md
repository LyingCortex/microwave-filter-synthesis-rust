# Changelog

## [Unreleased]

### Changed

- **Round-3 review remediation: port-diagonal guard, complex fast-path
  gate, measured high-order boundary, and with-zeros group-delay anchor.**
  Closes the three round-3 findings (`docs/reviews/mfs-review-round3.md`).
  (1) **`from_matrix` port-diagonal guard (P2-2):** the Schur-2×2 fast path
  assumes pure series terminations (`M[0,0] = M[N+1,N+1] = 0`), but
  `assert_transversal_structure` only scanned the resonator block, so a
  legitimate `CouplingMatrix` carrying a port shunt was accepted and its
  shunt **silently dropped** (measured ~1e-1 in |S21|). The structure check
  now rejects a nonzero port diagonal above the same scale-aware tolerance
  as the resonator block; the `from_matrix` contract doc now states the port
  requirement. Three new `pole_expansion` tests pin it (source reject, load
  reject, tolerance boundary: 1e-12 accepted / 1e-8 rejected).
  (2) **Fast-path verification gate now compares complex S21 + S11 (P2-2 /
  §六):** the 5-point cross-check behind `response_normalized` /
  `response_bandpass` compared only `|S21|`, so a phase-rotated or S11-only
  distortion would slip. It is now a shared `fast_path_gate` helper comparing
  the **complex** S21 (magnitude and phase) and S11, treating any non-finite
  compared value as an automatic LU fallback. Five `design` tests pin it,
  including the regression case (equal |S21|, rotated phase) the old
  magnitude-only gate would have passed.
  (3) **High-order boundary made honest and test-anchored (P2-1):** README
  and `architecture.md` no longer claim "stable to 30+"; the measured
  boundary (all-pole OK through order 28, finite-TZ through order 24 at
  RL 20–22 dB; order 30+ typed-rejected by the adaptive pole-axis tolerance)
  is now the documented claim. `examples/high_order_stability.rs` splits into
  a **verified set that must succeed** (exit non-zero on any regression) and
  an **explicitly-labeled exploratory scan** (order 30+, rejection expected),
  and a `#[test]` in the example locks the boundary (order 28 succeeds,
  order 30 rejected).
  (4) **`synthesize` doc + with-zeros group-delay anchor (P3-1 / P3-2):** the
  `synthesize` doc now states a successful result does not certify the design
  meets its RL spec (use `verify_against_spec`), with a failure-handling
  doctest on `verify_against_spec`. A new `response` test adds the
  previously-missing *independent* group-delay value anchor on a finite-TZ
  design (central difference of S21 phase, where the phase is smooth), plus
  a step-convergence guard (n=401 vs n=801) and a zero-region check that
  `group_delay` is `None` (not a fake 0.0) at a transmission zero on both
  backends.

- **O(N) pole-expansion fast path now enforces its transversal precondition
  instead of silently dropping resonator-resonator couplings.**
  `PoleExpansionData::from_matrix` (the O(N) partial-fraction fast path behind
  `FilterDesign::response_normalized` / `response_bandpass`) only models a
  *transversal* coupling matrix: the resonator block `M[1..=N][1..=N]` is
  diagonal and each resonator couples only to the source (`M[0,k]`) and the
  load (`M[k,N+1]`). Folded / arrow / wheel matrices carry inter-resonator
  couplings `M[k,k']` that this path *cannot* represent — but it read only
  `M[0,k]`, `M[k,k]`, `M[k,N+1]` and `M[0,N+1]`, so any real `M[k,k']` was
  **silently discarded** and a physically wrong (but silently computed)
  response was returned. The only existing safeguard was the 5-point LU
  cross-check in `design.rs`, which compares |S21| and therefore can *miss*
  the error when the dropped coupling does not move |S21| (e.g. a complex
  phase-only or source/load-only perturbation), and it only fires when the
  `from_matrix` fast path is even attempted.
  `from_matrix` now runs a **structural** precondition
  (`assert_transversal_structure`): it scans every off-diagonal resonator
  entry and rejects any whose magnitude exceeds a scale-aware tolerance
  `1e-9 · max(1, max|M|)` — far above f64 round-off (~1e-16 on a normalized
  matrix) and far below any meaningful O(1) coupling — returning
  `MfsError::PreconditionViolation` (the same structural variant
  `assert_symmetric` uses; no new public variant, so no serde /
  exhaustive-match break). The check inspects the *entries*, not the
  topology label, so a matrix labelled `Transversal` that still holds an
  inter-resonator coupling is caught too. On rejection the main-flow
  `response_normalized` / `response_bandpass` fall through to the LU backend
  (`ResponseSolver`), which handles arbitrary symmetric matrices, so a
  folded/arrow design still returns the *correct* response — just via LU.
  Five regression tests in `pole_expansion::tests` pin the contract: a
  genuinely transversal matrix is accepted; a ~1e-12 round-off residual is
  accepted while an O(1) coupling in the same matrix is rejected (proving
  the check is not a no-op); folded and arrow matrices from
  `transform_matrix` are rejected with a typed `PreconditionViolation`; and
  `response_normalized` on a folded-matrix design matches the LU reference to
  1e-6 (the fallback-protection behaviour, end to end).

- **Infeasible all-finite transmission-zero specs now return a
  domain-specific `PreconditionViolation` instead of a misleading generic
  "coefficients must be finite" error.**
  The all-finite formulation uses `eps_r = eps / sqrt(eps^2 - 1)`, which is
  real-valued only for `eps > 1`. When every transmission zero is finite and
  lies near the band edge (`|z| ~ 1`) or the return loss is high, `eps <= 1`
  and the square root is `0` (at `eps = 1`) or `NaN` (below), so `eps_r` is
  singular/complex. Previously this propagated a `NaN 1/eps_r` into the
  E-polynomial stage and surfaced later, in
  `complex_poly::ComplexPolynomial::new`, as the unrelated
  `"coefficients must be finite"` finiteness error with no hint of the actual
  cause. `find_eps` (shared by the single-band
  `synthesize_generalized_chebyshev_data` and the multi-band
  `multi_band_chebyshev_polynomials` paths) now rejects `eps <= 1` up front
  with a `PreconditionViolation` that names the value, the domain, and the
  feasible alternative (a mixed zero set — keep at least one zero at
  infinity, which forces `eps_r = 1` and is always synthesizable — or move
  the finite zeros off the band edge). No new public error variant, so no
  serde / exhaustive-match break. The feasible all-finite domain (`eps > 1`)
  is unchanged: well-separated zeros still synthesize and hit the RL spec.
  Seven new tests pin the contract (5 unit tests in
  `generalized_chebyshev_helpers` — mixed-zero control, `eps < 1` reject, the
  exact `eps = 1` singular reject, and two pipeline-level infeasible-spec
  rejects — plus 2 integration tests in `synthesis_regression` through the
  public `generalized_chebyshev_polynomials` API). The all-finite feasibility
  domain (`eps > 1`) is now documented in `docs/algorithms.md` §2.4 and
  `docs/architecture.md` Known Limitations.

- **`flip_sign` / `flip_sign_inplace` deduplicated into one shared helper.**
  Two copies of the similarity sign flip (`R·M·Rᵀ`, `R[i][i] = -1`: negate
  resonator `i`'s row and column, leaving the diagonal unchanged) lived in
  `matrix::core` (`flip_sign`, the immutable `&self -> Self` form used by
  `to_folded`) and `matrix::sections` (`flip_sign_inplace`, the `&mut self`
  form used by the python-feature CT/CQ `extract_*` topology methods),
  cross-referenced by two "consider deduplicating" TODOs. Both now delegate
  to a single `pub(crate) fn flip_sign_inplace(&mut self, index)` in
  `src/matrix/mod.rs`; `core::flip_sign` is a thin `clone + call` wrapper.
  The row/column-negation logic is byte-for-byte the old `inplace` body, so
  the behaviour is unchanged. **Bit-exact A/B verified:** the same 9
  prototype cases (orders 4–10, including even-order asymmetric filters
  that exercise the sign flip) were synthesized and run on the pre-refactor
  commit (`353f506~1`) and this commit; the `to_folded` matrices, the `to_ct` /
  `to_cq` matrices (the `flip_sign_inplace` paths), and the 401-point folded
  S-parameter responses came out **byte-identical** (5364-line dump, `diff`
  empty). Regression test
  `matrix::core::flip_sign_negates_row_and_column_only_and_preserves_diagonal`
  pins the three semantics a naive edit breaks — the flipped node's diagonal
  is negated twice and comes out **unchanged** (the historical third-negation
  bug dropped transmission zeros for even-order asymmetric filters), only
  the flipped row/column change sign, and the immutable/in-place forms agree
  entry-for-entry. All 289 `cargo test --lib` tests green; clippy clean
  (default + `cli,python` features).

- **Docs synced to current state + new implementation-grounded algorithm
  reference.** Added `docs/algorithms.md`: a per-stage writeup of what the
  code *actually* computes — the Cameron w-domain recurrence (with the
  corrected V-update, no spurious `w` factor) and the `E(w) = F/eps_r +
  j·P/eps` upper-half-plane-monic stage; admittance Y-polynomial assembly
  (N-mod-4 leading-coefficient guard + odd-order parity sign); residue
  expansion (1e-5 single-linkage clustering, Laurent principal parts, per-Y
  independent summation); transversal assembly (adaptive pole-axis tolerance,
  M_SL sources, repeated-pole `NotImplemented`); the `to_folded` / `to_arrow`
  Givens sequences + the shared `flip_sign_inplace`; CT/CQ cascade tiling;
  LU + pole-expansion response (5-point verification gate, group-delay units
  + 2π trap); the Macchiarella–Tamiazzo N-channel iteration; dual-band
  route-α + the K≥3 DOF wall; cepstrum + balanced-companion spectral
  factorization; LM optimization/extraction — each with `file:line` anchors
  and a consolidated numerical-threshold table. Also corrected the stale
  "group delay is `None` on the fast path" claims in `README.md` and
  `docs/api-reference.md` (both solver paths now compute it; `None` only at a
  transmission zero), updated the `matrix/`/`synthesis/`/`bin/` module tree
  and test counts in `docs/architecture.md`, and refreshed the
  `docs/internal/todo.md` status header to the 0.6.2-unreleased head
  (`cargo test` 435/0 across 22 binaries).

- **§D.7 dB single source — library side closed.** A crate-wide shared
  helper `response::magnitude_to_db(mag) -> f64` (−300 dB floor at
  `mag <= 1e-15`, else `20·log10(mag)`) is now the single source of truth for
  every magnitude→dB conversion. The seven library hand-rolled lambdas
  (`duplexer.rs` antenna-S11 + isolation + the common-domain `to_db`;
  `multiplexer.rs` geometric-mean S11 + pairwise isolation + the
  `response_admittance`/`response_matrix` dB; `synthesis/machiarella.rs`
  `junction_response` `to_db`) and the public `ResponseSample::s11_db()`/
  `s22_db()` all route through it. The `s*_db` threshold tightens from
  `mag <= 0.0` to `mag <= 1e-15`, which clamps the machine-epsilon null
  region (a true transmission zero lands at `mag ~ 1e-17`) to a *stable*
  −300.0 — previously `20·log10(1e-17)` amplified run-to-run residual noise
  into a large, unstable dB swing. This is what makes the `cli_save_resume_
  e2e` deep-null gate and the CLI −300-floor tests reproducible. No test pins
  the old `(0, 1e-15]` window; the full suite is green. The `multiplexer.rs`
  `peak_db` `NEG_INFINITY` site is a *degeneracy sentinel* (a fully-reflective
  channel has no passband), not a serialization floor, and is intentionally
  left distinct.

### Added

- **CI wheel/pyo3 build + FFI gate** (`.github/workflows/ci.yml` `wheel`
  job, Linux). The Python wheel was previously built only by the manual
  `scripts/build_wheel.sh` (needs uv + docker), so the cdylib/FFI path had no
  CI gate: a breaking change in the `python` feature or the pyo3 surface
  would surface only at wheel-build time, downstream. The job builds the
  wheel with the same `uv build --wheel` command as `build_wheel.sh`,
  installs it into a clean venv, and asserts a pure-Rust-computation value —
  the k12 coupling of a 4th-order RL=20 dB filter, stable across glibc
  2.35→2.39 (libm has no algorithm change there) — against the Rust-side
  reference to 1e-6 (the same cross-glibc tolerance `build_wheel.sh` uses).
  The glibc-2.36 *runtime* load check still needs docker + credentials (out
  of scope for CI; see review-findings §F).


## [0.6.2] — 2026-10-11

### Added

- **Fast-path group delay (pole-expansion).** `PoleExpansionData::evaluate_at`
  now returns `(s11, s21, group_delay)` instead of `(s11, s21)` and emits a
  finite group delay at O(N)/point. The fast path previously returned
  `group_delay: None` (explicit-undefined, T4.3) because the group-delay
  consumer didn't exist yet; it is now computed as `−Im(d/dω ln S21)` where
  S21 is a known rational function of ω in the Schur 2×2 form
  (`src/response/pole_expansion.rs`), so its log-derivative is a per-pole
  partial-fraction sum — one extra O(N) sum alongside the existing
  S-parameter sum. It is `None` exactly at a transmission zero
  (`|S21| → 0`), mirroring the LU backend's degenerate guard (explicit
  undefined, not a fake 0.0). Both the normalized and bandpass fast-path
  callers (`evaluate_normalized` / `evaluate_bandpass`) forward it, so
  `design.response()` / `response_bandpass()` now deliver engine-native group
  delay (the Python binding's `response_full` exposes it via `group_delay`).
  **Verified against three independent algorithms:** (1) fast vs LU on a
  3-zero filter, full band, to 1e-7 — `group_delay_matches_across_fast_and_lu_
  backends` (Check 1, the cross-algorithm guarantee); (2) fast vs an
  all-pole-ladder central difference of the S21 phase (Check 2, value anchor);
  (3) cross-language, the new value converted to seconds via the exact BP
  chain rule `dω/df = (1/FB)(1/fc + fc/f²)` matches **both** skrf `group_delay`
  — the same third-party engine the GUI uses for its group-delay plots — and an
  independent central difference of the S21 phase to **16.6 ps** over 268
  points of a 2 GHz / 200 MHz band. The 16.6 ps is the 2nd-order
  central-difference's own O(h²) truncation at Δf = 0.9 MHz: skrf computes
  `group_delay` with the identical `−dφ/dω` central difference on the same S21
  values, so it agrees with the FD reference to machine precision (0.000 ps)
  and sits 16.6 ps from the analytic fast path on both sides. Measured O(N):
  1801-pt sweep **2.0 µs/pt** (prototype) / **2.2 µs/pt** (bandpass) with the
  derivative active. **Check 3 is reproducible** from the consumer repo at
  `fuqtgui/src/scripts/verify_fastpath_gd.py` (needs this version + skrf in the
  venv; `uv run python src/scripts/verify_fastpath_gd.py`, exit 0 on pass).
  Checks 1–2 are the permanent CI legs (Rust test, no third-party dep).

- **Whole-crate code review of 0.6.1 (2026-10-09).**
  `docs/internal/review-findings-0.6.1.md`: six parallel read-only review
  slices over all of `src/` (design/spec/freq API, matrix/transform/output,
  approx/response/spectral, synthesis/verify, duplexer/multiplexer/optimize,
  pipeline/CLI/Python/touchstone) plus an independent adjudication pass that
  re-derived every confirmed finding against source and **empirically
  reproduced the two highest-stakes claims** with a temporary throwaway
  example (deleted after use). Tally:
  **3 P1 / 13 P2 / 10 P3**. The three P1s are the "gate-green but wrong
  result" class on the *default* paths: (1) the simplified duplexer/
  multiplexer recovery path (`recover_channel_design`, `rescale_d=true`)
  forces D to monic and silently destroys the passband return-loss spec —
  the "D's leading coefficient cancels to ~0" premise in the code comment is
  wrong, because `sub` trims the cancelled degree-n term and the stored D is
  genuinely degree n−1 with an O(1) leading coefficient (reproduced: worst
  case passband RL −20.52 dB → −2.53 dB, spec 22 dB); (2)
  `Multiplexer::response()` (default Admittance model) reads wrong
  column-major S-matrix indices for `channel_s21_db[i]` (i≥1) and
  `isolation_db[i][j]` (all pairs) (reproduced: true −72.31 dB
  antenna→channel-1 reported as −0.38 dB); (3) the response solver exploits
  symmetry without a consumer-side guard, so the public opt-out asymmetric
  constructors yield silently wrong group delay. Notable P2s:
  `to_cascaded_trisections` order=1 `usize` underflow, `f_start_hz`/
  `f_stop_hz` returning arithmetic-mean edges instead of the library's
  geometric-mean `band_edge_pair` (breaks `bandpass_edges` round-trip),
  `default_response` negative window for fbw>2/3, unvalidated `unloaded_q`
  in `response_lossy*`, `PoleExpansionData::from_polynomials`
  pole-domain inconsistency, pipeline `unloaded_q` dead parameter,
  `--resume` evaluating the untransformed matrix, and missing resource
  upper bounds on `order`/`grid.points`. Verification at assessment:
  `cargo test --all-targets` 395 passed / 0 failed, clippy `-D warnings`
  clean, docs clean. **The assessment's findings were closed out in the same
  session — see the "Fixed" section below.**

- **`cross_backend_equivalence` gate carried by the rust-native synthesis bench**
  (`benches/synthesis.rs`). The contract's fifth twoport gate — "the two
  backends build the same filter" — now runs in CI without the sibling python
  backend present. The reference is the **python backend's** folded
  resonator-subblock eigenvalues, pinned in
  `bench-data/golden/cross_backend_reference.json` (vendored from
  `filter-bench/data/golden/`, the authoritative home). The mfs side computes
  the same eigenvalues from its own folded matrix (real Schur, nalgebra/faer)
  and compares sorted-paired at the contract tolerance (1e-6). Measured
  deviation across all 7 standard twoport cases: **≤ 2.7e-12** (machine
  precision), so the two independently-implemented backends produce the same
  filter. The eigenvalue multiset is a similarity invariant — it absorbs the
  folding convention, plane rotation, and global sign that make element-wise
  comparison impossible. New typed accessor
  `bench_support::data::CrossBackendReference` and metric
  `bench_support::gates::cross_backend_equivalence_metric`; new conformance
  test `cross_backend_equivalence_matches_python_reference` (pins the 1e-6
  bound and the reference's real-valuedness) in `tests/benchmark_standard.rs`.
  No `contract_version` bump: the gate, its metric, tolerance, and report
  schema are unchanged — this adds authoritative data and a rust-native
  carrier for an existing gate.

### Changed
- **The coupling-matrix symmetry invariant now lives in the primitive
  constructor, not the builder.** A one-sided off-diagonal entry
  (`M[i][j]` set without `M[j][i]`) silently produced an asymmetric matrix
  that every consumer assumes is symmetric — and `response/backend.rs`
  *exploits* that symmetry to cut a solve (it reads the last row of A⁻¹ as the
  transpose of the last column), so an asymmetric matrix yields a silently
  wrong S21/S12, not a loud failure (the same silent-corruption class the
  `flip_sign` diagonal bug was). `CouplingMatrix::new_with_topology` (and
  `new`) now rejects such a matrix with a typed `PreconditionViolation` naming
  the offending pair and the fix; `allow_asymmetric(true)` on the builder, or
  the explicit `CouplingMatrix::new_with_topology_asymmetric` constructor,
  routes around the check for the deliberately asymmetric folded-mask
  fixtures in `tests/benchmark_standard.rs`.

  Sinking it to the constructor — rather than the `build()`-only check this
  work first shipped as — also closes the **CLI user-input bypass**:
  `mfs_cli parse_matrix_json` builds from `new_with_topology` directly, so
  a user-supplied asymmetric JSON matrix is now rejected before it reaches
  the solver (pinned by the new e2e `test_optimize_rejects_asymmetric_user_matrix`).

  The check is **tolerance-based and scale-aware** (`1e-12 * max(1, |a|, |b|)`
  per pair): a folded (Givens-rotated) matrix is symmetric only to within
  rounding (~1 uLP), so exact `a != b` would false-reject legitimate rotated
  matrices (the default gate caught one in `optimizer_integration` at
  `M[4][5]` vs `M[5][4]`, 1 uLP apart). The threshold still rejects a gross
  one-sided entry while permitting rotation/rounding noise.

  Six new lib tests pin the constructor invariant (symmetric accept,
  one-sided rejection + error-naming, asymmetric bypass, diagonal-only
  accept, plus two boundary tests that pin the tolerance from both sides:
  a relative ~1e-13 pair is accepted, a relative ~1e-6 pair is rejected).
  The two `transpose()` tests that use deliberately asymmetric data now route
  through `new_with_topology_asymmetric`; the builder keeps its
  three tests (symmetric build, one-sided rejection, opt-out). Closes the open
  "builder symmetry ergonomics" item.


- **The D3-opt optimizer e2e gate is now live in CI** (it was dormant).
  `tests/cli_integration.rs` is `#[cfg(feature = "cli")]` and the `mfs_cli`
  bin has `required-features = ["cli"]`, so the default `cargo test
  --all-targets` / `clippy --all-targets` CI steps never built the binary —
  the e2e test compiled to zero tests and the `cli`/`python` feature code was
  never linted. `.github/workflows/ci.yml` gains a `cli` job (3-OS matrix:
  `cargo build --features cli --bin mfs_cli` then `cargo test --features cli
  --test cli_integration`, which builds the binary via `CARGO_BIN_EXE_mfs_cli`
  and runs the 2 optimizer e2e gates — `optimize` order-4 folded recovery and
  `extract` order-4 transversal — plus the 4 pipeline gates) and a
  `feature-lints` job (`cargo clippy --features cli,python --all-targets
  -- -D warnings`). Wiring the features surfaced lints the default gate had
  been hiding, now fixed: `run_design` 14 args → grouped `DesignArgs` struct
  and `run_remaining_stages` collapsible `match` arm (`mfs_cli.rs`); range
  index loops → `iter().enumerate()` (`cli_integration.rs`); two redundant
  `_` wildcard match arms and pyo3 `tune`/`duplexer` signatures with
  `#[allow(clippy::too_many_arguments)]` (`python.rs`). The
  `cli_integration.rs` binary-missing guard changed from a silent skip to a
  hard `assert!` — a gate that cannot fail is no gate.

### Fixed
- **0.6.1 review fix wave (2026-10-09/10)** — the three P1s and all thirteen
  P2s from the [0.6.1 review](docs/internal/review-findings-0.6.1.md)
  (documented in "Added" above as an assessment) are fixed in the working tree.
  The 72 dB→−0.38 dB class of "gate-green but wrong result" defects on the
  *default* paths is closed. 395 → **406 passed** (11 new regression tests),
  clippy/docs/bench/conformance all green (see the findings doc §F).
  - **P1-1 simplified-path recovery rescale** (`src/synthesis/recovery.rs`):
    the `rescale_d` flag and the whole rescale branch are **deleted**; every
    caller now supplies D already in S's normalization and it is used as-is.
    The "D's leading coefficient cancels to ~0" premise was wrong — `sub`
    trims the cancelled degree-n term so the stored D is genuinely degree n−1
    with an O(1) leading coefficient, and the monic rescale moved the
    reflection zeros (worst-case passband RL −20.52 → −2.53 dB, spec 22 dB).
    New regression `multiplexer_recovered_channel_preserves_return_loss`
    (zero-perturbation limit, E/F coefficient match within 1e-8).
  - **P1-2 `Multiplexer::response()` S-index** (`src/multiplexer.rs`): the
    column-major extraction read the wrong elements —
    `channel_s21_db[i]` (i≥1) and `isolation_db[i][j]` (all pairs). Now
    `channel_s21_db[i] = |S[i+1][0]|`, `isolation_db[i][j] = |S[i+1][j+1]|`;
    the struct doc rewritten to state the column-major storage. New regression
    `multiplexer_response_db_fields_match_s_matrix` (asserts each dB field
    against the raw `response_matrix()` element).
  - **P1-3 solver symmetry guard** (`src/response/backend.rs`):
    `evaluate_response` now calls `matrix.assert_symmetric()?` before the LU
    solve (the consumer-side half of the c0902ba invariant), so a
    deliberately asymmetric matrix returns a typed `PreconditionViolation`
    instead of a silently wrong group delay.
  - **P2-1** `to_cascaded_trisections` order=1 `usize` underflow →
    `order < zeros.len() + 2` + `saturating_sub(2)`; test
    `ct_zero_guard_order1_returns_error_not_panic`.
  - **P2-2** `f_start_hz`/`f_stop_hz` return the geometric-mean
    `band_edge_pair(c, bw)` edges (was arithmetic-mean `c ∓ bw/2`); the
    `bandpass_edges` round-trip test re-anchored to the correct convention.
  - **P2-3** `default_response` window derived from `band_edge_pair` edges and
    floored positive (was `center − 1.5·bw`, which went < 0 Hz for fbw > 2/3);
    test `response_default_wideband_positive_window`.
  - **P2-4** `response_lossy*` validate `unloaded_q` (non-finite/non-positive
    rejected; ∞ = lossless kept) — was silently dropped to a lossless
    response.
  - **P2-5** `compare_responses` `group_delay_compared` flag: set when
    *at least one* sample pair is measured (any-pair semantics — group delay
    is legitimately `None` at transmission zeros), so an entirely-unmeasured
    comparison (the pole-expansion fast path) fails a finite group-delay
    tolerance instead of passing on a silent 0.0; test
    `compare_responses_uncompared_group_delay_fails_finite_tolerance`.
  - **P2-6** `MultiplexerSynthesis::synthesize` validates the TypeII junction
    `c0` (finite, positive) — was NaN/Inf-poisoned; test
    `multiplexer_type_ii_rejects_nonpositive_capacitance`.
  - **P2-7** the duplexer simplified (default) synthesis path only implements
    the ideal Type-I junction, but silently discarded a requested non-ideal
    Type-I `b0` and remapped a Type-II junction to ideal. The simplified path
    now rejects a non-ideal Type-I (`n≠1` or `b0≠0`) or any Type-II junction
    with a typed `PreconditionViolation` pointing at `.full_iteration(true)`
    (which consumes the parameters exactly); the ideal Type-I default is
    unchanged; `junction()` doc states the path constraint; regression
    `simplified_path_rejects_non_ideal_junction` (4 parameter sets) and
    `macchiarella_type_ii_junction` re-anchored to the full path.
  - **P2-8** pipeline `unloaded_q` is no longer dead: the response stage
    consumes it via `lossy_settings` (normalized `Q·bw/center`) — was a
    lossless response with no warning; test
    `response_evaluation_stage_applies_lossy_unloaded_q`.
  - **P2-9** `--resume` no longer drops the transform/response artifacts:
    `ContextSnapshot::into_context` rejects snapshots whose metadata
    disagrees with the restored artifacts and clears the non-restorable
    stages (`topology_transform`, `response_evaluation`) from
    `stages_executed` (trimming the parallel `stage_timings_ms`), so a resume
    re-runs them instead of skipping by stale metadata and leaving a context
    that claims the stages ran but has no artifacts. `run_remaining_stages`
    moved from the CLI binary into `mfs::pipeline` so the real `--resume`
    driver is covered by the library suite; regression
    `resume_of_completed_context_reproduces_full_pipeline` (a fully-completed
    context, saved, loaded, and resumed, reproduces the full pipeline's
    transform topology and S-parameters).
  - **P2-10** context serialization: infinite transmission zeros → `null`
    (finite-only Q), response dB through the clamped `s*_db()` (−300 dB
    floor) — was a `serde_json` `NonFiniteValue` failure.
  - **P2-11** public request validator enforces `MAX_ORDER=256` /
    `MAX_GRID_POINTS=100_000` (was unbounded — a `{"order": 100000}` request
    OOMs). The bounds now live at the choke points every public entry funnels
    through — `FilterSpec::new` (`spec::MAX_ORDER`, order > 256 →
    `InvalidOrder`) and `FrequencyGrid::linspace` (`freq::MAX_GRID_POINTS`,
    points > 100 000 → `InvalidGridSize`) — so the CLI `design`/`multiplexer`
    subcommands and the Python `bandpass`/`prototype` bindings (which
    construct the spec/grid directly, bypassing `SynthesisRequest::validate`)
    are bounded as well; `pipeline::request` re-exports the constants.
    `InvalidOrder`/`InvalidGridSize` messages now state the supported range.
    Regressions `spec_rejects_order_above_max`, `linspace_rejects_points_above_max`.
  - **P2-12** `PoleExpansionData::from_polynomials` (s-plane poles, wrong
    S-parameters) deleted; only the real ω-domain `from_matrix` remains, with
    the convention documented.
  - **P2-13** `run_design` `--topology` has explicit `folded`/`transversal`/
    `arrow` arms; unknown values error instead of silently degrading to
    Folded.
  - **P2-14 (newly surfaced this pass) broken `mfs_cli` unit test + CI gate
    gap** (`src/bin/mfs_cli.rs`, `.github/workflows/ci.yml`): a pre-existing
    test broke at HEAD and was invisible because the `cli` job only *built*
    the binary and ran the black-box `cli_integration` test, never the
    binary's own unit tests. `parse_matrix_json_flat_and_nested` fed the
    asymmetric order-1 fixture `[0,0.7,0.7,0,0.7,0,0,0,0]` (`M[0][1]=0.7`
    vs `M[1][0]=0`), which fails now that the symmetry invariant lives in
    `new_with_topology` (commit `c0902ba`). Fixed by re-anchoring the fixture
    to a symmetric order-1 ladder and adding `cargo test --features cli
    --bin mfs_cli` to the CI `cli` job so the binary's unit tests are gated
    on every OS.
  - **P3 #11 response-stage matrix source** (`src/pipeline/execution.rs`):
    `run_response_evaluation_stage` (the `--resume` / `--stage` driver) always
    evaluated `context.matrix` (the untransformed transversal synthesis
    matrix), while `run_full_pipeline` evaluated the transformed (folded)
    `transform.matrix` when a transform exists. The two are similar
    representations of the same network, so S-parameters agreed to machine
    precision everywhere except at transmission-zero nulls, where the `20·log10`
    transform amplified a ~1e-16 residual into an ~11 dB spread in the
    serialized dB column. The stage now selects `transform.matrix` when a
    transform is present, exactly mirroring `run_full_pipeline`; a resumed
    context reproduces the full pipeline's S-parameters bit-for-bit (found
    while verifying the P2-9 fix).
  - **P3 #1 + #9 pole-expansion fast path** (`src/response/pole_expansion.rs`,
    `src/design.rs`): a near-singular Schur determinant (`|det|² ≤ 1e-30`, a
    network pole) previously returned a fake `(S11, S21) = (1, 0)` total-
    reflection sample, and the verification loop skipped every check point in
    the deep stopband (|S21| ≤ 1e-10) so the fast path could be accepted
    **unverified**. `evaluate_at` now returns `Result` and surfaces a typed
    `NumericalFailure` at the near-singular point; `evaluate_normalized` /
    `evaluate_bandpass` propagate it; and `response_normalized` /
    `response_bandpass` fall back to LU on that error and accept the fast path
    only when at least one check point was actually compared (an all-skipped
    check grid now falls back to LU). Regressions
    `near_singular_schur_determinant_is_an_error`, plus the fast-path LU-
    fallback behavior in the design tests.
  - **P3 #2 three-state check outcome** (`src/verify/mod.rs`):
    `ResponseCheckReport::passes()` returned `true` for skipped checks,
    conflating "not checked" with "passed". New `CheckOutcome` enum
    (Skipped/Pass/Fail) with `ResponseCheckReport::outcome()`; the terminal
    and markdown summaries now render the outcome (`skipped`/`pass`/`fail`).
    `CheckOutcome` is re-exported from the crate root and the prelude.
  - **P3 #3 arrow-pattern load column — investigated, REVERTED as a false
    positive** (`src/verify/mod.rs`): the review claimed `matches_arrow_pattern`
    should reject a nonzero load column (index `side−1`). Empirically that is
    wrong: the real `to_arrow` reduction **requires** nonzero load-column
    cross-couplings (they encode the transmission zeros) — e.g. order 4 has
    `M(2,5)`, `M(3,5)` nonzero alongside the load coupling `M(4,5)`. Enforcing
    a structural-zero load column rejected every legitimate arrow/wheel matrix
    (5 existing tests failed). The band loop `(row+2)..=order` is the correct,
    intended check; the only genuinely-unchecked residue is a hypothetical
    direct source→load path, which the reduction never creates. No code change.
  - **P3 #4 prototype `freq_range` Hz misinterpretation**
    (`src/design.rs`, `TouchstoneBuilder`): an explicit physical-Hz `freq_range`
    on a prototype design (no center/bw) is silently reinterpreted as
    **normalized ω** while the Touchstone frequency column still shows the Hz
    values (e.g. 6.5e9–7.0e9 → deep stopband data). `freq_range` on a
    prototype design now returns a typed `InvalidFrequency`; the normalized
    form is available via `freq_range_normalized`.
  - **P3 #5 multi-band order / repeated-zero validation**
    (`src/spec/multi_band.rs`, `src/approx/multi_band.rs`): the "order ≥ number
    of finite transmission zeros" doc invariant is now enforced early with a
    typed `DimensionMismatch` (was a misleading equiripple non-convergence
    deep in the approximation stage); repeated (multiplicity) zeros are no
    longer silently deduplicated to first order when building P.
  - **P3 #6 `RejectionSpec` index validation** (`src/spec/multi_band.rs`):
    `InterBand { index }` is now validated against `band_count − 1` at
    construction (was an unenforced documented invariant).
  - **P3 #7 error-message / condition mismatch** (`src/freq.rs`,
    `src/spec/multi_band.rs`): the "0 < f1 <= f2" message is aligned to the
    enforced condition `f1 < f2` (f1 == f2 is rejected), so the message no
    longer promises equality is allowed.
  - **P3 #8 legacy `filter_spec` `unloaded_q` default** (`src/lib.rs`): the
    deprecated `filter_spec` shim defaulted `unloaded_q` to `Some(2000.0)`
    while the new builder APIs default to `None`; it now delegates to
    `synthesis::build_filter_spec` and matches the `None` (lossless)
    convention, removing the dead `DEFAULT_UNLOADED_Q` const. Metadata-only
    field (not consumed by synthesis); legacy_api_still_works test still
    passes.
  - **P3 #10 pipeline/CLI/Python hygiene batch**:
    (a) `pipeline/schema.rs` output schema now matches the real
    `SynthesisContext` serialization — `transform.response_invariant` (a
    phantom field) replaced by `notes: Vec<String>`, `matrix` gains `side`,
    `spec` gains `unloaded_q`, and both topology enums gain `Wheel`;
    (b) `scripts/build_wheel.sh` `repair()` now retags **only** the newly
    built wheel (`NEW_WHEEL`) instead of every non-target wheel in `dist/`
    (which also rewrote prior build-original archives and made `NEW_REPAIRED`
    point at the last glob file);
    (c) `mfs_cli` `main()`'s JSON error fallback no longer
    `format!`-interpolates the message (which could emit invalid JSON) — it
    now builds the object via `serde_json::json!` so the message is escaped
    correctly;
    (d) `write_touchstone` is now atomic (write to a temp file in the same
    directory, then rename), so a mid-write failure cannot leave a truncated
    `.s2p`;
    (e) Python `tune()` now rejects empty `freqs` (matching `extract()`).
    The dead duplicate `pipeline::schema::validate_request` (a `pub`
    `&Value` validator with zero callers, shadowing the canonical typed
    `SynthesisRequest::validate`) and its 8 tests are removed.
  - **Open (not fixed):** none of the 13 P2s remain. All P3 items are now
    closed: #1, #2, #9, #11 fixed; #3 investigated and reverted (false
    positive); #4-#8 and #10 fixed this session (see items above).
    **Verification (this session):** `cargo test --features cli
    --all-targets` **414 passed / 0 failed**; `cargo clippy --features
    cli,python --all-targets -- -D warnings` clean. (The 8 now-deleted
    `pipeline::schema::validate_request` tests are the only test removal;
    the net change adds the P3 #5 regression and P3 #4/#6 regression tests.)

- **0.6.1 working-tree second-pass review fix wave (2026-10-10)** — a
  read-only re-review of the post-round-1 working tree found and fixed a set
  of "gate-green but wrong result" and silent-degradation defects that the
  round-1 wave had left open (the round-1 fixes touched the same files, so a
  second pass re-derived each known risk against the current code and hunted
  for new ones). All findings are fixed in the working tree.
  - **C2 spectral-factor companion fallback unverified + RHP fill**
    (`src/spectral.rs`): the cepstrum path verified the para-product identity
    `D(jω)·D̃(jω) = R(jω)` (= `|D(jω)|² = R(jω)`) at 8 axis points, but the
    companion-matrix fallback returned its factor **unverified** — and
    `select_lhp_roots`' last-resort fill has no `re < 0` check, so it can
    admit a right-half-plane root when the LHP count is short. A factor whose
    magnitude response is wrong would silently corrupt the downstream
    Machiarella coupling matrices. The 8-point check is now a shared
    `para_product_verified` helper applied to **both** paths (cepstrum and
    companion); a factor that fails it returns a typed `NumericalFailure`.
    This is a **magnitude** backstop, not a minimum-phase enforcer: for a real
    R the min- and max-phase factors share `|D(jω)|²` (both pass), so
    minimum-phase *selection* remains `select_lhp_roots`' LHP bias. New
    regression `para_product_verification_checks_magnitude_not_min_phase`
    (pins both the reject behavior and the real-R limitation).
  - **Type-II junction degenerate `d_lead`** (`src/synthesis/machiarella.rs`):
    `compute_s_product`'s Type-II branch monic-normalizes the spectral factor
    by dividing by `d_poly.leading_coefficient()`; `spectral.rs::monic` leaves
    a degenerate factor (|lead| < 1e-15) unnormalized, so the division would
    scale the polynomial across ~200 decades and produce catastrophic roots.
    Now a typed `NumericalFailure` instead. New regression
    `type_ii_s_product_rejects_degenerate_leading_coefficient`.
  - **F1/F2/F3 pipeline snapshot round-trip** (`src/pipeline/execution.rs`):
    infinite transmission zeros and the lossless-`∞` Q were stored as
    non-finite `f64` and broke the `serde_json` round trip; now stored as
    `Option` (finite) with `None` ↔ `∞` on restore. `into_context` now
    rejects a snapshot whose metadata claims a stage ran whose artifact is
    absent (was an unbounded `stage_timings_ms[i]` index — a hand-crafted /
    truncated file would panic), and the dead `has_transform`/`has_response`
    metadata fields are removed (written, never read). Regressions
    `save_and_load_context_round_trip_with_infinite_zero` and
    `loaded_context_rejects_crafted_stage_metadata`.
  - **F4/F5/F6 pipeline schema + lossy Q** (`src/pipeline/schema.rs`,
    `src/pipeline/stages.rs`): the input `topology` enum was missing `Wheel`;
    the output `response` schema was an object-with-`samples` but the real
    serializer emits a **direct array** (now matches). A "lossy" request with
    a non-positive / non-finite `unloaded_q` (other than the lossless `∞`
    sentinel) silently degraded to a lossless response (the backend only
    applies the term when `is_finite() && > 0.0`) — now rejected with a typed
    error; a finite Q with a low-pass mapping (incommensurate grid) is also
    rejected. `ResponseEvaluationStage`'s doc now states the lossless
    normalized path. Regressions `lossy_settings_rejects_invalid_q`,
    `lossy_settings_rejects_finite_q_with_lowpass`,
    `response_stage_rejects_invalid_q_on_normalized_path`.
  - **MultiBandSpec order bound** (`src/spec/multi_band.rs`): `new` enforced
    `order > 0` but not the `MAX_ORDER = 256` cap the single-band
    `FilterSpec::new` applies (an OOM vector: an unbounded order sizes an
    `(order+2)²` matrix). Now bounded; new regression
    `multi_band_spec_order_is_bounded_by_max_order`.
  - **Multiplexer default response window** (`src/bin/mfs_cli.rs`,
    `src/python.rs`, `src/design.rs`): `resolve_multiplexer_range` and the
    CLI/Python/`cascade` docs still used the arithmetic `center ∓ 1.5·bw`
    window, which goes **< 0 Hz** for wideband channels (fbw > 2/3) and fails
    in `map_hz_to_normalized`. Now derived from the geometric-mean
    `band_edge_pair` edges padded by 0.5·bw and floored positive — the same
    formula as `FilterDesign::default_response` (the P2-3 fix). New regression
    `multiplexer_default_range_is_edge_based_and_positive`.
  - **Polynomial-summary zero nullability** (`src/pipeline/context.rs`,
    `src/pipeline/schema.rs`): the public `SynthesisContext` `Serialize` impl
    serialized `PolynomialsSummary.transmission_zeros_normalized` as a bare
    `Vec<f64>`, but a zero at infinity (legal — see
    `approximation_accepts_explicit_infinite_transmission_zeros`) lands in
    that list. The P2-10 fix had applied the finite→`null` mapping to the
    *spec* summary but missed the *polynomial* summary, leaving the output
    (which serde_json emits as `null` for `inf`) inconsistent with the schema
    (which declared `items: number`). Now mapped to `Vec<Option<f64>>` like
    the spec summary; schema items retyped to `["number", "null"]`. New
    regression `context_serialize_emits_null_for_infinite_normalized_zero`.
  **Verification (this wave):** `cargo test --all-targets` **411 passed / 0
  failed** (the round-1 406 + 5 new second-pass regressions);
  `cargo clippy --features cli,python --all-targets -- -D warnings` clean;
  `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` clean;
  `cargo bench --no-run` compiles.
- **0.6.1 review third-pass closure (2026-10-11)** — the last open follow-ups
  from the 0.6.1 review (the §C "risks needing verification" and the §D
  architecture recommendations) are now closed, plus the four §F "not
  executed" verification legs. All changes build on the round 1+2 wave
  (`8f3191a`); details in `review-findings-0.6.1.md` §I.
  - **§C.2 `select_lhp_roots` RHP fill** (`src/spectral.rs`): the last-resort
    fill (reaching `target_count` when the LHP pass is short) had no `re < 0`
    check and could admit a right-half-plane or pure-imaginary root — a silent
    non-minimum-phase spectral factor feeding the full Machiarella path. The
    fill now admits only strictly-LHP roots (`re < -1e-9`). New regression
    `select_lhp_roots_fill_rejects_rhp_when_lhp_short`. (Closes the §C.2 risk
    that the round-2 C2 `para_product_verified` backstop did **not** cover —
    the backstop is a magnitude check, the fill was a min-phase *selection*
    gap.)
  - **§C.1 faer backward-residual gate** (`src/approx/complex_poly.rs`):
    **verified, no issue.** The gate accepts `max_backward ≤ 1e-3` and
    **rejects** residuals in (1e-3, 1e-2), falling through to Durand-Kerner /
    Aberth (the safe direction). Probed on production E_s/F_s/P_s polynomials
    at orders 4..40: backward residual ≤ ~1e-5, none in the band. New
    regression `faer_gate_does_not_fire_on_production_polynomials_up_to_order_40`;
    the gate comment now states the accept/reject semantics explicitly.
  - **§C.3 high-order normal-equation conditioning** (`src/spectral.rs`):
    **verified, perf-only.** At order 40 the cepstrum path falls through to
    the (correct) companion fallback — the O(N log N) intent is lost but there
    is no correctness or timeout problem. Measurement only, no code change.
  - **§C.4 external-Q convention** (`src/matrix/scaling.rs`): **verified,
    anchor added.** `Qe = 1/(M01²·Δ)` (vs the textbook `1/(M01·Δ)`) is
    self-consistent under this code's coupling convention where `M01 = √rate`;
    a physical order-1 RL=20 bandpass (`Qse = 2·f0/BW3dB`) matches `1/(M01²·Δ)`
    to <1% and the `1/(M01·Δ)` form is off by a factor `M01`. New anchor
    `external_q_convention_matches_physical_order1_bandpass`.
  - **§C.5 `FullIterationState` one-iteration lag** (`src/synthesis/
    machiarella.rs`): **verified, documented.** The lag is real and bounded by
    the convergence error (< `STAGNATION_ACCEPT` = 1e-3); it now carries a
    comment at the capture site. Fixing it would change convergence behavior;
    not worth the risk.
  - **§C.6 `to_proto` misleading error** (`src/synthesis/machiarella.rs`):
    the "non-finite root set" error now names the actual cause — an out-of-range
    *finite* root (|r| ≥ 50, the historical stability filter), which is distinct
    from a genuinely non-finite root. The pre-existing message-assertion test
    was aligned to the new message.
  - **§C.7 CLI `optimize`/`extract` grid contract** (`src/bin/mfs_cli/
    optimize.rs`, `main.rs`): the normalized contract is documented, and a
    magnitude guard rejects a physical-Hz grid (span > 1e6 → typed
    `PreconditionViolation`), so feeding a real `.s2p` no longer "converges"
    the optimizer at ±6e9 normalized frequency. New regression
    `response_grid_rejects_physical_hz_and_accepts_normalized`.
  - **§D.3 absolute spec checks** (`src/verify/mod.rs`): new
    `verify_against_spec(spec, response)` checks a response **against the
    design spec** (which `compare_responses` — two responses agreeing — never
    did): passband return loss (worst-case min RL ≥ spec − 0.1 dB), passband
    insertion loss (|S21| in [−1, +0.5] dB), transmission-zero depth (≤ −40 dB
    within ±2 grid steps; zeros at infinity → Skipped), and stopband rejection
    (always Skipped — `FilterSpec` carries no normalized stopband requirement,
    so there is nothing to gate on without inventing a requirement). Missing
    data is `Skipped`, never a silent `Pass`. Exposed in `lib.rs`/`prelude.rs`.
  - **§D.5 split `mfs_cli.rs`** (`src/bin/mfs_cli/`): the 1911-line
    `src/bin/mfs_cli.rs` is now `main.rs` (entry/dispatch) + `design.rs` +
    `multiplexer.rs` + `optimize.rs` + `target.rs` (shared
    `load_target`/`parse_matrix_json`) + `io.rs`. Line-by-line move, no CLI
    flag/help/output change.
  - **§D.7 dB single source** (`src/bin/mfs_cli/io.rs`): the five hand-rolled
    `10·log10(re²+im²)` sites in the CLI now route through
    `ResponseSample::s11_db()`/`s21_db()` (−300 dB floor at zero magnitude — a
    true transmission zero no longer serializes as `null`); the new
    `verify_against_spec` reuses the same methods (a private `db_from_mag` it
    first shipped was deduplicated onto them in the review follow-up). The
    `multiplexer.rs` `min_db`/`max_db` operate on library-precomputed `Vec<f64>`
    dB (no `ResponseSample`) and are left as-is, cited for §D.7. **Scope note:**
    pre-existing library-side hand-rolled `mag → dB` lambdas
    (`duplexer.rs:101/113/276`, `multiplexer.rs:462/958/979/1028`,
    `machiarella.rs:754`) predate this commit and were not unified — a full
    §D.7 closure is open and non-blocking.
  - **§D.6 save/resume closed-loop e2e in CI** (`tests/cli_save_resume_e2e.rs`,
    `.github/workflows/ci.yml`): a permanent, CLI-side successor to the deleted
    P3 #11 harness. `mfs_cli --input` (full pipeline) vs `mfs::pipeline::
    save_context` + `mfs_cli --resume` must reproduce the same response —
    non-zero samples within 1e-9 dB, transmission-zero samples as deep S21
    nulls in both runs (deep-null assertion, not bit-equal dB: 20·log10
    amplifies the ~1e-17 residual into a large dB swing). Wired into the CI
    `cli` job's exact test command.
  - **§F not-executed legs:** bench-data conformance tiers **executed green**
    (`benchmark_standard` 20/0, both the sibling byte-for-byte and the manifest
    sha256 tiers ran — the sibling is present here, correcting the "only tier
    available" framing; `folded_response_invariance` 10/0). The Python **wheel
    build** now produces an installable artifact (k12 = 0.414467071,
    bit-identical to the Rust probe) and `scripts/build_wheel.sh` is patched to
    recognize and repair maturin's untagged `linux_x86_64` output (no
    `patchelf` on the build host) rather than aborting. The **python runtime**
    was exercised through pyo3 (`import mfs` + `bandpass`/`duplexer`/
    `prototype`; `extract`/`bandpass_edges` present) — no longer "statically
    reviewed only". Non-Linux platforms remain **not executable here** (CI's
    3-OS × 2-rust matrix is the coverage; a CI job for the wheel/pyo3 runtime
    is a noted future gap).
  **Verification (this wave):** `cargo test --all-targets` all green;
  `cargo clippy --all-targets -- -D warnings` clean;
  `cargo clippy --features cli,python --all-targets -- -D warnings` clean;
  `cargo doc --no-deps` clean; the exact CI command
  `cargo test --features cli --test cli_integration --test cli_save_resume_e2e`
  green.

## [0.6.1] — 2026-10-07

### Added

- **Whole-crate code review of 0.6.0 (2026-10-05).**
  `docs/internal/review-findings-0.6.0.md`: five parallel review slices over
  all of `src/` plus an independent adjudication pass. 18 reported findings →
  **17 valid (6 P1 / 5 P2 / 6 P3) + 1 false positive** (the `response_lossy`
  Q formula, re-derived and confirmed correct). The six P1s are the
  "gate-green but wrong result" class: full-path channel recovery rescaling
  (`recovery.rs`), sign-mirrored band-pass denormalization (`scaling.rs`),
  LU-backend group-delay port nodes (`response/backend.rs`), silently clamped
  overlapping-channel p0² (`machiarella.rs`), no-op CQ merge (`cascade.rs`),
  and a Python `tune()` length-check gap. Assessment only at review time;
  the "silent errors" wave (P1-1/2/3/4/6, P2-1/2/3/5, all six P3) is fixed
  2026-10-06 (see Fixed below); the CT/CQ cascade wave (P1-5 + P2-4, doc §F.2)
  is fixed the same day (see Fixed below).

### Fixed

- **0.6.0 "silent errors" review wave (2026-10-06)** — the five
  "gate-green but wrong result" P1s plus P2-1/2/3/5 and all six P3
  (review-findings-0.6.0.md §F.1):
  - **P1-1 full-path channel recovery** (`src/synthesis/recovery.rs`):
    `recover_channel_design` gained a `rescale_d` flag; the full path passes
    `false` (D = S − F⁰ is already in S's monic normalization), so the a
    priori reflection zeros are no longer destroyed (was ~1e30 root
    corruption on the §V golden spec). Simplified path keeps the rescale.
    New value-level gate `multiplexer_full_iteration_preserves_a_priori_reflection_zeros`.
  - **P1-2 band-pass denormalization sign** (`src/matrix/scaling.rs`):
    `denormalize_resonator_frequency` now uses the positive-frequency branch
    `(t + √(t²+4c²))/2` (was the sign-mirrored `c(√(…) − t)`), and
    `normalize_resonator_frequency` routes through
    `BandPassMapping::map_hz_to_normalized` so the two cannot drift. New
    cross-check test vs `freq.rs`.
  - **P1-3 LU-backend group delay** (`src/response/backend.rs`): the numerator
    fold is `(1..side−1)` (resonator block only), not `(0..side)` — the spurious
    port-node terms are gone (was 2.679 vs 0.905 at ω=0). New value-level test
    vs the numerical S21 phase derivative.
  - **P1-4 overlapping-channel p0²** (`src/synthesis/machiarella.rs`): a
    genuinely negative p0² (`< −1e-12`, calibrated from the measured
    legitimate floor ~4e-12 and the non-physical negatives O(1e-2)) is now
    rejected with a typed `PreconditionViolation` instead of silently clamped;
    the `1e-15` floor is kept only for ~0 round-off.
  - **P1-6 Python `tune()`** (`src/python.rs`): all five companion arrays are
    length-checked (was only `target_s21_re`), raising `PyValueError` instead
    of crashing the interpreter.
  - **P2-1 CLI `optimize`** reads the matrix from `--input`/stdin even when the
    target is a `.s2p` (was "invalid matrix JSON").
  - **P2-2 CLI `design`** reports an `InputParse` error for an incomplete
    band-pass parameter set (was a silently-ignored center, normalizing to a
    prototype).
  - **P2-3 pipeline schema/validator** accept the canonical PascalCase
    `MatrixTopology` (`Transversal`/`Folded`/`Arrow`), matching the serde
    request path (was lowercase-only).
  - **P2-5 `DuplexerSynthesisResult`** now carries `pub converged: bool`
    (was dropped; multiplexer had it).
  - **P3** `FilterDesignBuilder` band-pass validation at the choke point (P3-1),
    `tune()` uniform-grid contract documented (P3-2), `bandpass_edges` validates
    `0 < f_start < f_stop` (P3-3), `from_polynomials_and_matrix` propagates the
    spec error instead of a silent RL=20 fallback (P3-4),
    `PoleExpansionData::from_polynomials` y22 length check (P3-5), Givens step-2
    comment corrected (P3-6).
- **0.6.0 "CT/CQ cascade" review wave (2026-10-06)** — P1-5 + P2-4
  (review-findings-0.6.0.md §F.2), the deferred second wave:
  - **P2-4 `to_cascaded_trisections` severs the mainline** (`src/matrix/cascade.rs`):
    the old per-zero loop re-conditioned the tail on each pass and rotated on the
    fixed load column, producing a dense folded/wheel matrix with a broken mainline.
    Now routes through the edge-capable `extract_triplet` primitive
    (`Arrow → extract_triplet(zeros[0], 1)`: condition-at-load-tail + peel-to-source),
    which produces the source/load-connected trisections. Reproduces Cameron 2011
    Fig. 21 (order 4, 22 dB, TZ 1.8082/1.3217): trident sparsity, max|ΔM| ~ 1e-4 vs
    the 4-decimal printed fixture, max|ΔS21| ~ 2e-15; order-6 mainline connected.
    Note: the earlier same-day note that "no path reproduces Fig. 21" had tested the
    *distinct* internal-center-only `extract_trisection` primitive — `extract_triplet`
    (center=1 allowed) is the correct edge-capable one. New gates:
    `ct_synthesis_reproduces_fig21_cameron2011` (trident + value anchor) and the
    order-6 mainline-connectivity assert in `ct_synthesis_produces_valid_matrix`.
  - **P1-5 `to_cascaded_quadruplets` merge was a no-op** (`src/matrix/cascade.rs`):
    the old `diagonal_rotation_angle` merge vanished for symmetric ±ω₀ pairs, leaving
    the CQ byte-identical to the CT. Now routes through the **verified
    `extract_cascade_quadruplets`** primitive (the same path as the Python `cq`
    binding) with a **source-side non-overlapping tiling**: pair *i* spans resonators
    `4i-3..=4i` (cross coupling `M[4i-3,4i]`), so `order >= 4·num_pairs` is required —
    otherwise a typed `PreconditionViolation` (use `to_cascaded_trisections` for denser
    configs). Measured: order-4 single pair → clean quadruplet (source→res 1 only,
    load→res 4 only, cross M[1,4] present, diagonals zeroed, max|ΔS21| 3e-15);
    **order-6 single pair (2 zeros) → clean single quadruplet** (cross M[1,4] present,
    mainline 1..=6 connected, source→res 1 / load→res 6, no other 2-skip,
    max|ΔS21| 2e-14); **order-8 two pairs (4 zeros) → reproduces the verified
    `CascadeQuartets` layout** `[[1,3],[1,4],[5,7],[5,8]]` (max|ΔS21| 1.4e-13,
    structurally distinct from the CT, max|ΔM| > 0.05). The over-constrained
    order-6 2-pair case (the earlier "~1e-2 residual" follow-up) is **topologically
    impossible** to tile cleanly and is now rejected. Gates:
    `cq_synthesis_single_quadruplet_structure`, `cq_synthesis_order6_single_quadruplet`,
    `cq_synthesis_order8_two_quadruplets` (verified layout + CQ≠CT + mainline),
    `cq_synthesis_rejects_too_small_order`.
- **Fresh-eyes review of the CT/CQ work (2026-10-07)** — no wrong-result
  (library) bug, but three gate/doc defects, all measured and fixed
  (`src/matrix/cascade.rs`, `src/design.rs`):
  - **CQ load-port gate never asserted the load column** (off-by-one). The load
    port is row `side-1` (`side = order+2`), so order-6 load = row 7, order-8
    load = row 9. The order-6 gate asserted `get(5,6)`/`get(4,6)` (the res5-res6
    mainline and res4-res6 2-skip, both already gated below) mislabelled as
    "load-res6/5"; the order-8 gate asserted `get(7,8)`. Since CQ is a
    similarity transform, S21-invariance cannot catch a load-port regression, so
    a wrong load coupling would have passed CI green. Corrected to the actual
    load column: order-6 `get(7,6)` > 0.05 (res6-load, measured 1.06) +
    `get(7,5)` < 1e-3 (res5-load, ~0); order-8 `get(9,8)` (res8-load, measured
    1.0443).
  - **CQ doc `common_resonator` corrected**: the source-side tiling pairs each
    quadruplet with a trisection starting at the same resonator, so
    `extract_cascade_quadruplets` selects `common_resonator = 1` on every
    quadruplet (the `= 4` load-end branch is never reached) — the doc had said 4.
  - **`to_ct` doc reconciled** with the single-`zeros[0]` implementation (the
    remaining zeros are already encoded in the arrow form; only the count is
    validated). Also fixed a `\\varepsilon` LaTeX typo in the multi-band DOF-wall
    doc and a stale "load-res7 ~0.023" comment in the order-8 gate (the measured
    inherent load couplings are load→res 2/4/6 = 0.010/0.035/0.020).

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
