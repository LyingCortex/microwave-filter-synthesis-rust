# MFS v0.6.1 Whole-Crate Code Review — Findings

**Scope:** Whole-crate defect hunt at the 0.6.1 boundary (`master`), all of `src/` plus
integration tests, CLI, Python bindings, and build scripts.
**Method:** Six parallel read-only review slices (S1 design/spec/freq API, S2
matrix/transform/output, S3 approx/response/spectral, S4 synthesis/verify, S5
duplexer/multiplexer/optimize, S6 pipeline/CLI/Python/touchstone) with
evidence-backed findings, then an independent adjudication pass: every confirmed
finding re-derived against source at `file:line`, and the two highest-stakes
numeric claims (simplified-path `rescale_d` recovery corruption, multiplexer
`response()` S-index) **empirically reproduced** with a temporary throwaway
example (deleted after use). This is an **assessment**, not a code change.
**Companion docs:** `review-findings-0.6.0.md` (prior wave, fixed & gated),
`architecture-review-2026-09-25.md`, `CHANGELOG.md` [Unreleased].

**Checks executed (2026-10-09, Linux x86_64, Rust 1.95.0):**

| Command | Result (assessment) | Result (after fix wave) |
|---|---|---|
| `cargo test --all-targets` | **395 passed, 0 failed, 0 ignored** | **406 passed, 0 failed, 0 ignored** (+11 regression tests) |
| `cargo clippy --all-targets -- -D warnings` | **0 warnings** | **0 warnings** |
| `cargo doc --no-deps` (CI gate `-D warnings`) | clean, no warnings | clean, no warnings |

**Fix status (2026-10-09/10, working tree, uncommitted).** The assessment's
findings were closed out in the same session. **All 3 P1s are fixed**
(P1-1 rescale branch deleted + simplified-path regression test; P1-2 S-index
corrected + N=2 S-matrix-match test; P1-3 consumer-side `assert_symmetric`
guard) and **all 13 P2s are fixed** (P2-1/2/3/4/5/6/7/8/9/10/11/12/13; P2-7
resolved by rejecting the non-ideal junction on the simplified path with a
typed error instead of documenting the silent ignore). One P3 was also
closed: #11 (the two response-evaluation drivers now select the same matrix),
found while verifying the P2-9 fix. Each finding below carries a `**Status:**` line.
---

## A. Verdict

**The core numerical architecture is sound and the CI gate is fully green — but the
green gate does not anchor the physical values of the default synthesis/recovery and
default analysis paths. Three P1 defects are "gate-green but wrong result" on the
default (non-opt-in) code paths:

1. **P1-1** — the default (simplified) duplexer/multiplexer recovery path
   (`rescale_d=true`) silently destroys the passband return-loss spec by forcing D
   to monic; the code comment's premise ("D's leading coefficient cancels to ~0")
   is a misconception about what `sub`/trim stores. Empirically reproduced:
   worst-case passband RL degraded from −20.52 dB to −2.53 dB (spec 22 dB).
2. **P1-2** — `Multiplexer::response()` (the default Admittance-model analysis)
   reads the wrong S-matrix elements for `channel_s21_db[i]` (i≥1) and
   `isolation_db[i][j]` (all pairs) due to a column-major vs row-major index
   mismatch. Empirically reproduced: a true −72.31 dB antenna→channel-1
   transmission is reported as −0.38 dB.
3. **P1-3** — the response solver exploits symmetry (last row of A⁻¹ = last column)
   without validating it on the consumer side; the public opt-out asymmetric
   constructors (`new_with_topology_asymmetric`, `builder().allow_asymmetric(true)`)
   can feed the solver and yield **silently wrong group delay** while S11/S21
   magnitudes stay correct.

These are the same defect class as the 0.6.0 wave (§A there), now on the
duplexer/multiplexer and analysis-API surfaces that the 0.6.0 wave did not gate.


---

## B. Confirmed defects

### P1-1 Simplified-path `recover_channel_design` rescale corrupts passband RL (default path)

- **Location:** `src/synthesis/recovery.rs:55-66` (`rescale_d` branch).
  Callers: `src/multiplexer.rs:334-343` (simplified path, `rescale_d=true`);
  `src/duplexer.rs:731-754` (`rescale_d = !self.full_iteration`). Both full-path
  call sites pass `false` (unaffected).
- **Description:** Recovery rebuilds each channel from converged S-roots via
  `E = S + D`, `F = S − D`, where `S = (E_i + F_i)/2` is **monic** (E and F are
  both explicitly normalized monic, `generalized_chebyshev_helpers.rs:425,452`)
  and `D_i = (E_i − F_i)/2`. In the unperturbed limit the identity holds
  **exactly without any rescale**: `S + D_i = E_i`, `S − D_i = F_i`. The rescale
  exists under the premise (recovery.rs:46-54 comment) that `D`'s leading
  coefficient "cancels to ~0" and is only rescaled "when numerically
  non-negligible" (historical behavior). That premise misreads the storage:
  `ComplexPolynomial::sub` calls `trim_trailing_complex_zeros`
  (`complex_poly.rs:198,833-844`), so the exactly-cancelled degree-n term of
  E−F is **trimmed away** and the stored D is genuinely degree n−1 with an O(1)
  leading coefficient — `leading_coefficient()` returns it, not ~0. The guard
  `d_leading.norm() > 1e-15` therefore **always fires**, rescales D to monic
  (factor `1/d_leading ≠ 1`), and `E_new = S + D′ ≠ E_i` / `F_new = S − D′ ≠ F_i`:
  the reflection zeros (roots of F) are changed and the passband-edge return loss
  is systematically destroyed. The full path's own comment (recovery.rs:50-54)
  states the correct convention — D is "degree n − 1, non-monic leading
  coefficient" and "rescaling it corrupts the exact reconstruction F = S − D" —
  yet the simplified path does exactly that.
- **Trigger:** any default (non-`full_iteration`) duplexer or multiplexer
  synthesis. Fires even for a zero-displacement iteration (recovery alone
  corrupts the channel).
- **Evidence (executed, §F):** order-4, RL 22 dB, 6.75 GHz/300 MHz:
  `E deg 4 lead 1.0000`, `F deg 4 lead 1.0000`, stored `(E−F) deg 3 lead
  2.3421` → D lead 1.1711 (O(1), not ~0); guard fires, factor 0.8539.
  Standalone worst-case passband RL **−20.52 dB** → recovered (simplified path)
  **−2.53 dB**.
- **Impact:** default-synthesized duplexer/multiplexer channels fail their return
  loss spec with no warning; existing tests assert only `peak S21 > −10 dB`
  (loose), so the gate stays green.
- **Minimal fix:** pass `rescale_d=false` at the simplified-path call sites
  (align with the full path) and delete/condition the rescale branch.
- **Regression:** assert the simplified-path recovered channel's passband-edge
  |S11| equals the standalone synthesis within 0.1 dB.
- **Confidence:** High (algebra + executed reproduction).
- **Status:** **FIXED.** The `rescale_d` flag and the entire rescale branch are
  deleted from `recover_channel_design` (recovery.rs); every caller now supplies
  D already in S's normalization and it is used as-is (full path: `D = S − F⁰`;
  simplified path: `D = (Eᵢ − Fᵢ)/2`, which `sub` trims to a genuinely degree n−1
  O(1)-leading polynomial — exactly the convention the full path's comment always
  stated). New regression `multiplexer_recovered_channel_preserves_return_loss`
  (multiplexer.rs) runs the zero-perturbation limit (`max_iterations(0)`) and
  asserts the recovered channel's monic-normalized E/F coefficients equal the
  standalone synthesis's within 1e-8 (the bug moved F's degree-(n−1) coefficient
  by ~0.24, five orders above the bound). The throwaway `examples/mux_diag.rs`
  instrument that confirmed the fix is deleted; the permanent test supersedes it.

### P1-2 `Multiplexer::response()` S-index off-by-one: wrong `channel_s21_db` / `isolation_db` (default analysis path)

- **Location:** `src/multiplexer.rs:700-720` (extraction in
  `response_admittance` from the flattened S-matrix of `compose_node_s`).
- **Description:** `compose_node_s` returns a nalgebra `DMatrix`; its `iter()`
  flattens **column-major** (`flat[k] = S[k%p][k/p]`) — verified from vendored
  nalgebra 0.33.3 source (`~/.cargo/registry/.../nalgebra-0.33.3/src/base/iter.rs`:
  inner loop walks the row dimension with `strides.0 = 1`, then jumps columns by
  `strides.1 = nrows`) and by the N=1 anchor test (`row[1]=S21`, `row[2]=S12`,
  multiplexer.rs:1072-1073). The extraction then indexes as if the flat order
  were the physical antenna→channel row:
  - `channel_s21_db[i] = row[i*p + 1]` — reads S[1][i] instead of the intended
    S[i+1][0] (antenna → channel i). Correct only for i=0 (coincidence).
  - `isolation_db[i][j] = row[i*p + (j+1)]` — reads S[j+1][i] instead of S[i+1][j+1].
    Wrong for **all** (i, j).
  - `channel_s11_db[i] = row[(i+1)*p + (i+1)]` — correct (diagonal).
  The struct docs (`multiplexer.rs:947` "row-major"; `:957-963` intended
  semantics) contradict both the actual storage and the code's indices.
- **Trigger:** any N≥2 multiplexer using the default `response()` (Admittance
  model) — the path the CLI `multiplexer` command and the Python bindings call.
- **Evidence (executed, §F):** N=2, p=3, at f≈ch0 center:
  `row[1]=S[1][0]=−10.82 dB` (correct ch_s21[0]), `row[2]=S[2][0]=−72.31 dB`
  (**correct** ch_s21[1]), `row[4]=S[1][1]=−0.38 dB` (**buggy** value); reported
  `channel_s21_db[1] = −0.38 dB`, equality with S[1][1] asserted true. A 72 dB
  isolation is reported as −0.38 dB.
- **Impact:** the default analysis API reports physically wrong quantities for
  every N≥2 design, in CLI tables/JSON and Python dicts; no test asserts these
  dB values (tests check unitarity/reciprocity/dimensions only).
- **Minimal fix:** `channel_s21_db[i] = row[i+1]` (S[i+1][0]);
  `isolation_db[i][j] = row[(j+1)*p + (i+1)]` (S[i+1][j+1]); fix the
  row-major/col-major doc contradiction.
- **Regression:** N=2 — assert `channel_s21_db[1]` ≈ 0 dB at channel 1's own
  center (passband peak) and deep-stopband at channel 0's center.
- **Confidence:** High (nalgebra source + index derivation + executed
  reproduction).
- **Status:** **FIXED.** The extraction now reads `channel_s21_db[i] = |S[i+1][0]|`
  (`row[i+1]`) and `isolation_db[i][j] = |S[i+1][j+1]|` (`row[(i+1)*p + (j+1)]`),
  matching the column-major flattening; `channel_s11_db[i]` (diagonal) was already
  correct and is unchanged. The struct doc (`isolation_db` field) was rewritten to
  state the column-major storage and the exact element mapping (S-index written in
  a code span so rustdoc does not parse it as an intra-doc link). New regression
  `multiplexer_response_db_fields_match_s_matrix` (multiplexer.rs) evaluates an
  N=2 mux and asserts every `channel_s21_db` / `channel_s11_db` / `isolation_db`
  entry equals the corresponding element of the raw `response_matrix()` S-matrix
  (the element the old index mis-selected is now asserted against directly).

### P1-3 Response solver does not validate symmetry — opt-out asymmetric matrices yield silent group-delay errors

- **Location:** `src/response/backend.rs:54-79` (group-delay numerator uses
  `row_last` from `A·y = e_{N−1}` — the last **column** of A⁻¹ — as the last
  **row**, valid only when `A = Aᵀ`, per the comment at :150-163);
  `validate_settings:97-110` checks resistances only, never matrix symmetry.
- **Description:** c0902ba sank the symmetry invariant into the **constructors**
  (`new_with_topology`), but the **consumer** (solver) has no corresponding
  guard. The public opt-outs `CouplingMatrix::new_with_topology_asymmetric`
  (core.rs:120) and `CouplingMatrixBuilder::allow_asymmetric(true).build()`
  (builder.rs:72) are documented "never feed the response solver" by text only —
  no mechanism enforces it. With an asymmetric M, `A ≠ Aᵀ`, `y` (last column) ≠
  last row, and **group delay is silently wrong** while S11/S21 magnitudes
  (which use `col_first` only) stay correct. The LU-fallback cross-check in
  `response_normalized` (design.rs:417-421) compares only |S21|, so it never
  fires on this.
- **Impact:** the silent-corruption class named in c0902ba's message remains
  open on the consumer side; wrong group delay with no error and no fallback.
- **Minimal fix:** in `evaluate_response`, validate symmetry with the
  `new_with_topology` tolerance before solving; return
  `PreconditionViolation` on violation.
- **Regression:** asymmetric matrix → `evaluate_normalized` returns `Err`;
  symmetric matrix → group delay matches the full-inverse result.
- **Confidence:** High.
- **Status:** **FIXED.** `evaluate_response` now calls `matrix.assert_symmetric()?`
  before the LU solve (backend.rs:31) using the same `1e-12·max(1,|a|,|b|)`
  tolerance the `new_with_topology` constructor uses, so the consumer enforces the
  symmetry invariant the constructors already establish. A deliberately asymmetric
  matrix (the public opt-out constructors) now returns a typed
  `PreconditionViolation` naming the offending pair instead of silently computing a
  wrong group delay; a symmetric matrix is unaffected. The opt-out's "never feed the
  solver" contract is now enforced by mechanism, not by comment. The constructor-side
  asymmetry rejection (core.rs) and its boundary tests already existed from the c0902ba
  wave; this closes the consumer-side gap.

### P2 confirmed (local defects / reliability)

| # | Location | Defect | Fix sketch | Status |
|---|---|---|---|---|
| P2-1 | `src/matrix/cascade.rs:36` (`to_cascaded_trisections`) | `zeros.len() > order - 2` is a `usize` subtraction: **underflows at order=1** (debug panic; release wraps to `usize::MAX`, guard becomes dead). Sibling `to_cascaded_quadruplets:77` uses the safe `order < 4*n` form. | `if order < zeros.len() + 2` with `order.saturating_sub(2)` in the message | **FIXED** — guard rewritten `order < zeros.len() + 2` with `saturating_sub(2)`; new test `ct_zero_guard_order1_returns_error_not_panic` |
| P2-2 | `src/design.rs:337-350` (`f_start_hz`/`f_stop_hz`); Python `f_start/f_stop` (`python.rs:76-84`) | Accessors return the **arithmetic-mean** edges `c ∓ bw/2`, but the library's band-pass convention is a geometric-mean center, so the Ω=±1 edges are `band_edge_pair` (`freq.rs:575-578`). `bandpass_edges` input edges cannot round-trip. | Return `band_edge_pair(c, bw)` (make it `pub(crate)`); fix the `bandpass_edges_equivalent` test | **FIXED** — accessors return `band_edge_pair(c, bw)`; test re-anchored to the geometric-mean round-trip (design.rs:947-950) |
| P2-3 | `src/design.rs:748-753` (`default_response`); `src/bin/mfs_cli.rs:530` (CLI default window) | Default window `start = center − 1.5·bw` goes **below 0 Hz when fbw > 2/3**; first sample fails in `map_hz_to_normalized`. | Clamp/derive the window from `band_edge_pair` with a pad, floor at a positive value | **FIXED** — window derived from `band_edge_pair` edges, floored at `1e-9·center`; new test `response_default_wideband_positive_window` |
| P2-4 | `src/design.rs:569-622` (`response_lossy*`); `src/bin/mfs_cli.rs:532-535` (CLI `--q`) | Negative or NaN `unloaded_q` never validated; backend silently drops the dissipation term → **lossless response with no warning**. | Validate at the `response_lossy*` entries (keep ∞ = lossless) | **FIXED** — `validate_unloaded_q` called at both `response_lossy`/`response_lossy_normalized` entries (∞ = lossless kept); finite-only Q documented |
| P2-5 | `src/verify/mod.rs:176-181` (`compare_responses`) | Group-delay comparison only runs when **both** samples have `group_delay = Some`; the pole-expansion fast path leaves it `None`, so the comparison is a **silent no-op**. | Add a `group_delay_compared` flag; treat "required but not compared" as failure | **FIXED** — `group_delay_compared` set when *at least one* sample pair is measured (any-pair semantics; `None` at transmission zeros is legitimate), so an entirely-unmeasured comparison fails a finite group-delay tolerance; new test `compare_responses_uncompared_group_delay_fails_finite_tolerance` |
| P2-6 | `src/multiplexer.rs:139-155` (entry validation) | Checks order/center/bw but **not** the junction: `TypeII { c0 }` with `c0 = 0/NaN/negative` → NaN/Inf junction coupling silently poisons all roots. | Validate the junction at the top of `synthesize` | **FIXED** — TypeII `c0` validated finite+positive at the top of `synthesize`; new test `multiplexer_type_ii_rejects_nonpositive_capacitance` |
| P2-7 | `src/duplexer.rs:625-629` (simplified loop) | `let (n_sq, _) = …` discards `b0`; `TypeII` remapped to `(1.0, 0.0)`. The junction parameters the user requested are silently ignored. | Either consume the parameters or document/reject them on the simplified path | **FIXED** — the simplified path now rejects a non-ideal Type-I (`n≠1` or `b0≠0`) or any Type-II junction with a typed `PreconditionViolation` telling the user to opt into `.full_iteration(true)` (which consumes the parameters exactly); the ideal Type-I remains the default; regression `simplified_path_rejects_non_ideal_junction` (4 parameter sets), `macchiarella_type_ii_junction` re-anchored to the full path |
| P2-8 | `src/pipeline/request.rs:84-85`; `src/pipeline/execution.rs:124-143` | `SynthesisRequest.unloaded_q` is a **dead parameter**: stored and echoed but the response stage uses `ResponseSettings::default()` (unloaded_q = INFINITY). | Pass the (normalized) Q into `ResponseSettings`, or remove the field | **FIXED** — `run_response_evaluation_stage` passes `context.spec.unloaded_q` into the response stage; `lossy_settings` normalizes `Q·bw/center`; new test `response_evaluation_stage_applies_lossy_unloaded_q` |
| P2-9 | `src/pipeline/execution.rs:253-275` (`into_context`); `src/bin/mfs_cli.rs:1216-1220` (`run_remaining_stages`) | `--resume` restores metadata (incl. `stages_executed` with `topology_transform`) but **not** the transform artifact; the remaining-stages loop skips the transform by metadata, then evaluates the **untransformed** `context.matrix`. | Clear non-restorable stages from `stages_executed` in `into_context`, or make the evaluator require the transform artifact | **FIXED** — `into_context` now rejects snapshots whose metadata disagrees with the restored artifacts (also closes risk R2) and clears `topology_transform`/`response_evaluation` from `stages_executed` (trimming the parallel `stage_timings_ms`), so `--resume` re-runs both instead of skipping by stale metadata; `run_remaining_stages` moved into `mfs::pipeline` and the regression `resume_of_completed_context_reproduces_full_pipeline` proves a saved-then-resumed completed context reproduces the full pipeline's transform topology and S-parameters |
| P2-10 | `src/pipeline/context.rs:219-224`; `:269-270` | `TransmissionZero.value` documents `f64::INFINITY` support, but the final `serde_json::to_string` fails with `NonFiniteValue`; same for `log10(0) = −inf` at a sample on a transmission zero. | Serialize infinite zeros as `null`; clamp dB via the `s*_db()` helpers | **FIXED** — infinite TZ serialized as `null` (`is_finite().then_some`) and infinite Q as absent; response dB routed through the clamped `s11_db()`/`s21_db()` (−300 dB floor) |
| P2-11 | `src/pipeline/request.rs:99-131`; `src/freq.rs:589-601` | No **upper** bounds on `order` or `grid.points`: a single `{"order": 100000}` allocates ~80 GB and dies in the OOM killer. | Add constant upper bounds (order ≤ 256, points ≤ 100k) with typed errors | **FIXED** — `MAX_ORDER=256` / `MAX_GRID_POINTS=100_000` enforced at the choke points every public entry funnels through: `FilterSpec::new` (`spec::MAX_ORDER`) and `FrequencyGrid::linspace` (`freq::MAX_GRID_POINTS`), so the CLI/Python direct-construction paths (which bypass `validate()`) are bounded as well; `request` re-exports the constants; `InvalidOrder`/`InvalidGridSize` messages state the range; regressions `spec_rejects_order_above_max`, `linspace_rejects_points_above_max` |
| P2-12 | `src/response/pole_expansion.rs:91` (`from_polynomials` stores s-plane poles) vs `:141` (`from_matrix` stores real ω-domain poles) | `evaluate_at` computes `denom = ω − pole`, only correct for the real ω-domain convention; `from_polynomials` produces **wrong S-parameters** and the two constructors are domain-inconsistent. | Delete `from_polynomials` or convert s-poles to ω-domain; document the convention | **FIXED** — `from_polynomials` deleted (only the real ω-domain `from_matrix` remains); the ω-domain convention documented on `from_matrix` |
| P2-13 | `src/bin/mfs_cli.rs:520-524` (`run_design`) | `--topology` match `_ => Folded`: any unrecognized/misspelled value **silently degrades to Folded** — contradicting the no-silent-degradation policy. | Explicit `"folded"` arm; error on unknown | **FIXED** — explicit `"folded"`/`"transversal"`/`"arrow"` arms; unknown values error with "unknown topology" |
| P2-14 | `src/bin/mfs_cli.rs` `optimizer_tests::parse_matrix_json_flat_and_nested`; `.github/workflows/ci.yml` `cli` job | **Broken test at HEAD** (pre-existing, surfaced this pass): the order-1 fixture `[0,0.7,0.7,0,0.7,0,0,0,0]` is asymmetric (`M[0][1]=0.7` vs `M[1][0]=0`), so it fails once `c0902ba` sunk the symmetry invariant into `new_with_topology`. It went uncaught because the `cli` job only *built* the binary and ran the black-box `cli_integration` test — never the binary's own unit tests. | Fix the stale fixture to a symmetric ladder; add the binary's unit tests to the CI `cli` gate | **FIXED** — fixture re-anchored to a symmetric order-1 ladder (`[0,0.7,0,0.7,0,0.7,0,0.7,0]`); CI `cli` job now runs `cargo test --features cli --bin mfs_cli` |

### P3 confirmed (non-blocking)

**Status:** #11 closed this wave (the two response-evaluation drivers now
select the same matrix, per the fix noted below). A hardening pass closed
**#1, #2, and #9**: #1 the pole-expansion near-singular Schur determinant now
surfaces a typed error that falls back to LU instead of a fake `(1, 0)`
sample; #2 `ResponseCheckReport::outcome()` introduces a three-state
`CheckOutcome` (Skipped/Pass/Fail) so a skipped check is no longer conflated
with a pass; #9 the pole-expansion fast path is only accepted when at least
one check point was actually compared. **#3 was investigated and REVERTED as a
false positive** — the real arrow form has required nonzero load-column
cross-couplings (see item 3), so the "fix" rejected every legitimate arrow
matrix. **All remaining items are now CLOSED** (#4-#8 and #10 fixed
2026-10-10, see the per-item `Status: FIXED` notes below); the dead duplicate
`pipeline::schema::validate_request` was also removed as part of the #10
hygiene batch (the canonical validator is `request.rs::SynthesisRequest::validate`).


1. `src/response/pole_expansion.rs:260-263` — near-singular Schur determinant
   (`|det| ≤ 1e-15`) silently returns (S11, S21) = (1, 0) "total reflection" at a
   resonance instead of the finite peak; explicit failure or regularization
   preferred (consistent with the `group_delay: None`-rather-than-fake-0.0
   philosophy, T4.3). **Status: FIXED** — `evaluate_at` now returns
   `Result` and a near-singular Schur determinant (`det.norm_sqr() ≤ 1e-30`)
   surfaces a typed `MfsError::NumericalFailure` (the network pole) instead of
   a fake `(1, 0)` sample; `evaluate_normalized`/`evaluate_bandpass` propagate
   it and the `design.rs` fast path falls back to LU on that error. Regression
   `near_singular_schur_determinant_is_an_error` pins the typed failure.
2. `src/verify/mod.rs:107` — `ResponseCheckReport::passes()` returns `true` for
   skipped checks ("not checked" conflated with "passed"); introduce a
   three-state outcome. **Status: FIXED** — new `CheckOutcome` enum
   (Skipped/Pass/Fail) with `ResponseCheckReport::outcome()`; the terminal and
   markdown summaries now render the outcome (`skipped`/`pass`/`fail`) so a
   skipped check is no longer conflated with a pass. `CheckOutcome` is
   re-exported from the crate root and prelude.
3. `src/verify/mod.rs:210-217` — `matches_arrow_pattern` row loop
   `0..order().saturating_sub(1)` + column loop `(row+2)..=order()` never checks
   the load column (index `side−1`): nonzero residues at (0, side−1) and
   (N−1, side−1) pass the structural check. **Status: FALSE POSITIVE / REVERTED.**
   The premise is wrong: the real `to_arrow` output **requires** nonzero
   load-column cross-couplings — they encode the transmission zeros. Measured
   (order 4, side 6): `M(2,5)=−0.348`, `M(3,5)=−0.401` nonzero alongside
   `M(4,5)=0.946` (the load coupling); order 6: `M(4,7)=0.263`, `M(5,7)=0.340`,
   `M(6,7)=−0.967` nonzero. The load column is therefore **not** a
   structural-zero column, and enforcing it rejected every legitimate arrow /
   wheel matrix (5 existing tests failed: `transform_matrix_wraps_arrow_
   conversion`, `wheel_transform_matches_current_wheel_shape_rules`,
   `arrow_transform_matches_arrow_shape_rules`, `topology_transform_stage_
   `converts_to_arrow`, and the wheel variant). The band loop's column range
   `(row+2)..=order` is the intended arrow structural check and correctly stops
   one column before the load (`load = order+1 = side−1`), so a hypothetical
   direct source→load coupling `(0, side−1)` is the only genuinely-unchecked
   residue — and the real reduction never creates one. No code change; the
   load-column check was reverted.
4. `src/design.rs:847-856` (`TouchstoneBuilder::get_response`) — for a
   **prototype** design (no center/bw), explicit `freq_range(Hz)` values are
   silently reinterpreted as **normalized ω** (e.g. 6.5e9–7.0e9 → deep stopband
   data) with the file's frequency column still showing the Hz values: silent
   wrong data. **Status: FIXED** — `freq_range` on a prototype design now
   returns a typed `MfsError::InvalidFrequency` naming the issue; the
   normalized form is available via `freq_range_normalized`.
5. `src/spec/multi_band.rs:130-144` — doc promises "order ≥ number of finite
   transmission zeros" but `new` doesn't check; failure surfaces late as a
   misleading `DimensionMismatch` in the approximation stage (single-band path
   reports a clear error in `pad_transmission_zeros`). Also: `s_domain_p` builds
   P from **deduplicated |z|** — repeated (multiplicity) zeros are silently
   reduced to first order. **Status: FIXED** — the invariant is enforced early
   with a typed `DimensionMismatch` in the approximation stage
   (`more_finite_zeros_than_order_is_a_dimension_error` regression); repeated
   zeros preserve their multiplicity when building P.
6. `src/spec/multi_band.rs:102-110` — `RejectionSpec::new` doesn't validate
   `InterBand { index } < band_count − 1`; no consumer today, but the
   documented invariant is unenforced at construction. **Status: FIXED** —
   `RejectionSpec` validates the inter-band index against the band count at
   construction.
7. `src/freq.rs:462-465` (and `spec/multi_band.rs:49-52`) — error message says
   "0 < f1 <= f2" but the enforced condition is `f1 < f2` (f1==f2 rejected).
   **Status: FIXED** — the message now reads "0 < f1 < f2", matching the
   enforced condition.
8. `src/lib.rs:114-116` — legacy `filter_spec` defaults `unloaded_q` to
   `Some(2000.0)` while the new builder APIs default to `None`; metadata-only
   inconsistency (field not consumed by synthesis), deprecated shim.
   **Status: FIXED** — the shim now delegates to
   `synthesis::build_filter_spec`, matching the `None` (lossless) convention;
   the dead `DEFAULT_UNLOADED_Q` const is removed.
9. `src/design.rs:505-526` (and `:412-429`) — pole-expansion verification skips
   all check points with |S21| ≤ 1e-10; if all 5 check points fall in the deep
   stopband, the fast path is returned **unverified** (verification spins idle
   exactly where cross-checking matters most; absolute data risk low).
   **Status: FIXED** — `response_normalized`/`response_bandpass` now count how
   many check points were actually compared (deep-stopband points, |S21| ≤
   1e-10, are skipped) and accept the fast path only when `compared > 0`; an
   all-skipped check grid falls back to LU instead of returning an unverified
   response. This is fixed together with #1 (same fast path).
10. Pipeline/CLI/Python hygiene batch: `schema.rs:196-201` describes a
    `transform.response_invariant` field that the real `TransformSummary`
    doesn't have (+ `spec` missing `unloaded_q`, topology enum missing `Wheel`);
    `scripts/build_wheel.sh` `repair()` retags **all** stale wheels in dist/ and
    `NEW_REPAIRED` ends up pointing at the last processed file; `main()`'s JSON
    error fallback `format!` interpolates unescaped text (can emit invalid
    JSON); `write_touchstone` is not atomic (truncated .s2p on mid-write
    failure); Python `tune()` doesn't check `freqs` non-empty (unlike
    `extract()`). **Status: FIXED** (2026-10-10) — (a) the output schema now
    matches the real `SynthesisContext` serialization (`response_invariant`
    → `notes: Vec<String>`, `matrix.side` added, `spec.unloaded_q` added,
    `Wheel` in both topology enums); (b) `repair()` retags only the new
    `NEW_WHEEL`; (c) the JSON error fallback is built via
    `serde_json::json!` (escaped); (d) `write_touchstone` is atomic
    (temp file + rename); (e) `tune()` rejects empty `freqs`. The dead
    duplicate `pipeline::schema::validate_request` and its 8 tests were also
    removed (canonical validator: `request.rs::SynthesisRequest::validate`).
11. **Two response-evaluation drivers selected different matrices** (found
    while verifying the P2-9 fix): `run_full_pipeline` (`src/pipeline/request.rs:255-266`)
    evaluates the **transformed (folded)** matrix when a transform exists, but
    `run_response_evaluation_stage` (`src/pipeline/execution.rs`, the
    `--resume` / `--stage` driver) always evaluated `context.matrix`
    (the untransformed **transversal** synthesis matrix). The two are
    similar representations of the same network, so the S-parameters agree to
    machine precision everywhere *except* at transmission-zero nulls
    (`|S| ~ 1e-15`), where the `20·log10` transform amplifies the ~1e-16
    residual into an ~11 dB spread in the serialized dB column. **FIXED** —
    `run_response_evaluation_stage` now selects `transform.matrix` when a
    transform is present (exactly mirroring `run_full_pipeline`'s
    `eval_matrix`); `into_context` never restores `transform` (so `load_context`
    + `run_remaining_stages` re-runs the transform before the response stage),
    making the two drivers bit-identical. Verified: full-pipeline vs
    `mfs_cli --resume` produce 8/8 identical S-samples including the
    transmission-zero points (S11 exactly 0.0, no 1e-16 residual).

---

## C. Risks needing verification (suspicious, not closed in-repo)

1. **faer backward-residual gate width** (`complex_poly.rs:860`, `FAER_BACKWARD
   RESIDUAL_TOL = 1e-3`): roots with residual in 1e-3..1e-2 are accepted
   (no fall-through to Durand-Kerner/Aberth which could be better). No
   production polynomial known to land in that band.
2. **`spectral_factor` companion fallback has no para-product validation**
   (`spectral.rs:354-384`; the cepstrum path validates at 8 points) and
   `select_lhp_roots` last-resort fill (:305-318) can admit RHP roots when the
   LHP count is short → a **silent non-minimum-phase spectral factor** feeding
   the full Machiarella path (machiarella.rs:580).
3. **`fit_polynomial_from_freq_samples` solves an ill-conditioned normal
   equation** (spectral.rs:126-184): at high order (N≈40) the cepstrum path
   most likely falls through to the companion fallback — correct via the
   fallback, but the O(N log N) intent is lost (performance, not correctness).
4. **External-Q convention** (`matrix/scaling.rs:144,154`): `Qe = 1/(g²·Δ)` vs
   the standard `Qe = 1/(g·Δ)` round-trips (tests green) but the correct answer
   depends on whether the matrix element is g or √g; no physical order-1 anchor
   test exists to pin it.
5. **Full-path `FullIterationState` one-iteration lag**
   (`machiarella.rs:596-622`): the captured `common` (from iteration k's R) and
   `s_roots` (iteration k+1) are off by one iteration; negligible at strict
   convergence, up to ~STAGNATION_ACCEPT (1e-3) on stagnation-accept.
6. **`to_proto` |r|≥50 root filter** (`machiarella.rs`): out-of-range *finite*
   roots produce the "non-finite root set" error — misleading message.
7. **CLI `optimize`/`extract` grid contract** (mfs_cli.rs:1128-1141, 1021-1029):
   the target response's physical Hz samples are used directly as the
   **normalized** evaluation grid (documented contract); feeding a real .s2p
   silently "converges" the optimizer at ±6e9 normalized frequency.
8. **`--resume` snapshot JSON has no internal consistency validation**
   (R2): metadata can claim `matrix_synthesis` done while `matrix: null`;
   stages then skip per metadata on untrusted input. — **Addressed by the
   P2-9 fix:** `into_context` now rejects any snapshot whose metadata claims a
   stage whose artifact is absent (and length-mismatched `stages_executed` /
   `stage_timings_ms`), instead of skipping on untrusted metadata.

## D. Architecture / maintenance recommendations

1. **Single source of truth for band edges.** `freq.rs::band_edge_pair` should
   be the only producer of "passband edges" (used by `f_start_hz/f_stop_hz`,
   `default_response`, CLI default window). The arithmetic-vs-geometric-mean
   drift behind P2-2/P2-3 will recur otherwise.
2. **Consumer-side symmetry preconditions.** Symmetry is now a library-level
   invariant (c0902ba); the solver, `multiply`, and `to_complex_dense` should
   enforce "symmetric input" once, instead of each consumer assuming (P1-3).
3. **verify needs absolute spec checks.** `compare_responses` is
   "two responses agree" only; there is no function that checks a response
   *against the spec* (RL gate, TZ depth gate, stopband depth gate). "Verify
   passes" currently never means "meets the design target".
4. **Unify `PoleExpansionData` pole-domain convention** (P2-12) — one
   constructor, or an explicit s→ω conversion + documented convention.
5. **`mfs_cli.rs` (1854 lines) split by subcommand** (design / optimize+extract
   sharing `load_target`+`parse_matrix_json` / multiplexer / io) — pure move.
6. **`save_context`/`load_context` ↔ CLI** either closed-loop (`--save` + e2e
   test) or documented, with a consistency check after `into_context` (P2-9).
7. **dB computation single source** — `context.rs`, the CLI table, and
   `ResponseSample::s*_db()` three hand-rolled variants with different
   zero-magnitude behavior; unify (P2-10).
8. **One authoritative request validator** — `schema.rs::validate_request`
   (test-only) vs `SynthesisRequest::validate` (production) have drifted; make
   `describe_schema`+validator part of the real entry or delete one.

## E. Checked, no issues found

- **Symmetry invariant boundary** (S2 audit): every constructor path that
  reaches the solver (new / new_with_topology / builder default / synthesis
  engine / CLI JSON parse / snapshot restore / scaling rebuild) validates
  symmetry or is trivially symmetric; opt-outs are public, documented, and
  deliberate. Tolerance `1e-12·max(|a|,|b|)` (core.rs:103) sits above
  rotation roundoff (~1e-15) and far below O(1) one-sided bugs, with boundary
  tests pinned at 1e-13 accept / 1e-6 reject. (Gap: consumer side, P1-3.)
- **Root-solver fallback chain** (S3): faer gate conditions (leading
  `norm_sqr ≤ 1e-24`, eigenvalue count, backward-residual + non-finiteness),
  Durand-Kerner and Aberth-Ehrlich convergence gates, and the adaptive
  fall-through are all precisely defined; NaN/Inf roots are caught at three
  checkpoints; no silent NaN propagation path found in production code.
- **Cameron recurrence boundedness** (S3): for |zero| ≥ 1 inputs the
  recurrence coefficients stay O(1) (reciprocals ≤ 1, safe-sqrt guards), no
  overflow at high order; eps/eps_r edge (RL → extreme) guarded.
- **Multi-band LM pipeline** (S3): K≥3 rejected with `NotImplemented`,
  infeasible RZ placement escalates λ until a typed error; band-edge handling
  correct.
- **Matrix transforms** (S2): Givens similarity, folded/arrow reduction,
  `flip_sign`, scaling round-trips, cascade tiling/extract primitives
  (except P2-1) verified numerically at orders 1–6.
- **CLI/Python input surface** (S6): 13-entry validation audit; `parse_matrix_json`
  (square/rect/side≥3/order match/symmetry/unknown-topology rejection,
  anchored by `test_optimize_rejects_asymmetric_user_matrix`),
  `optimize` q-filter, `extract` grid checks, `multiplexer` string spec
  parsing, Python FFI surface (no unsafe/raw pointers; `MfsError` →
  `PyValueError` complete; 0.6.0 P1-6 out-of-bounds fixed and length-checked).
  Residual gaps listed above (P2-13, P3 #10, risk #7/#8).
- **Version sync**: Cargo.toml 0.6.1 == pyproject.toml 0.6.1; CI
  `version-sync` covers both.
- **Concurrency**: entire crate is synchronous numeric code — no rayon,
  async, threads, or shared mutable state; **no data-race surface** (all six
  slices independently confirmed).

## F. Executed verifications (this review)

Environment: Linux x86_64 (WSL2), Rust 1.95.0 (CI matrix minimum 1.86),
`target/` warm (20 GB).

1. `cargo test --all-targets` → **395 passed, 0 failed** (exit 0).
2. `cargo clippy --all-targets -- -D warnings` → **0 warnings** (exit 0).
3. `cargo doc --no-deps` → clean (exit 0).
4. **P1-1 reproduction** (temporary `examples/` program, deleted):
   standalone order-4/RL-22 BP: `E deg 4 lead 1.0000`, `F deg 4 lead 1.0000`,
   stored `(E−F) deg 3 lead 2.3421` (D lead 1.1711, O(1) — refutes the
   "~0" premise), rescale guard fires with factor 0.8539; worst-case
   passband RL standalone **−20.52 dB** vs recovered-simplified **−2.53 dB**.
5. **P1-2 reproduction** (same temporary program): N=2 mux, raw
   `response_matrix` at f≈ch0 center: `row[2]=S[2][0]=−72.31 dB` (correct
   antenna→ch1) vs `row[4]=S[1][1]=−0.38 dB`; reported
   `channel_s21_db[1]=−0.38 dB` — equality with the wrong element asserted.

**Post-fix verification (2026-10-09/10, after applying the fix wave):**

1. `cargo test --all-targets` → **406 passed, 0 failed, 0 ignored** (exit 0);
   the +11 over the 395 assessment baseline are the new regression tests
   (P1-1/P1-2/P1-3 rescale+S-index+symmetry guards, P2-5 group-delay flag,
   P2-6 mux TypeII c0, P2-7 simplified-junction rejection, P2-8 lossy Q,
   P2-9 `resume_of_completed_context_reproduces_full_pipeline`, P2-11
   `spec_rejects_order_above_max`/`linspace_rejects_points_above_max`).
2. `cargo clippy --all-targets -- -D warnings` → **0 warnings** (exit 0).
3. `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` → **clean** (exit 0) — the CI
   docs gate; the `multiplexer.rs` S-index doc and the `pole_expansion.rs`
   `evaluate_at` reference were written to avoid rustdoc intra-doc-link errors.
4. `cargo bench --no-run` → compiles (exit 0); `cargo test --test
   benchmark_standard` → **20 passed, 0 failed** (conformance gates green).
5. **P3 #11 E2E** (temporary `examples/` harness, deleted): order-4/RL-20
   Folded design, `run_full_pipeline` + `save_context` vs `mfs_cli --resume` —
   8/8 S-samples bit-identical including the transmission-zero points
   (pre-fix: the two nulls differed ~11 dB in the serialized dB column).

**Not executed:** benchmark conformance tiers (bench-data), actual Python
wheel build, non-Linux platforms, `python`-feature runtime behavior (FFI
surface statically reviewed only).

## G. Top three priorities

**At review time** the three priorities were the three P1s (all now fixed,
§B): P1-1 (rescale, default synthesis path — broadest blast radius), P1-2
(response() S-index, default analysis path — wrong physics for every N≥2
design), P1-3 (solver symmetry guard — the c0902ba consumer-side gap).

**Remaining priorities (post-fix):** with all 3 P1s and all 13 P2s closed,
the open items are the non-blocking P3s and §D/§C follow-ups:

1. **P3 hardening pass** — the non-blocking items that are cheap and
   high-value: the Schur-determinant near-singular guard (#1), the
   `passes()` three-state outcome (#2, which subsumes the P2-5 flag semantics
   once fully realized), and the arrow-pattern load-column gap (#3).
2. **One authoritative request validator** (§D.8) — `schema.rs::validate_request`
   (test-only) vs `SynthesisRequest::validate` (production) have drifted; the
   P2-11 choke-point bounds (`FilterSpec::new` / `FrequencyGrid::linspace`)
   now backstop both, but the two validators should still be unified or one
   deleted.
3. **§C verification items** — the risks that need empirical confirmation
   (faer backward-residual gate width, CLI optimize/extract grid contract,
   etc.).

## H. Second-pass working-tree review (2026-10-10)

After the round-1 fix wave, a second read-only pass re-derived the §C known
risks against the **current** working tree (the round-1 fixes touched the same
files) and hunted for new defects the first wave had left open. Six findings,
all fixed in the working tree this session (the `C*`/`F*` labels below are the
second-pass slice findings, distinct from the round-1 P1/P2/P3):

| # | Severity | Location | Defect | Fix | Regression |
|---|---|---|---|---|---|
| C2 | P2 (highest) | `src/spectral.rs` (`companion_spectral_factor`, `select_lhp_roots`, `spectral_factor`) | The cepstrum path verified the para-product identity `|D(jω)|² = R(jω)` at 8 jω points, but the companion-matrix fallback returned its factor **unverified**; `select_lhp_roots`' last-resort fill has no `re < 0` check, so it can admit an RHP root when the LHP count is short. A wrong-magnitude factor then corrupts the downstream Machiarella coupling matrices with no warning. **Scope note:** the check is a *magnitude* backstop, not a min-phase enforcer — for a real R, the min- and max-phase factors share `|D(jω)|²` (both pass); min-phase *selection* is `select_lhp_roots`' LHP bias. | Extracted the 8-point check into a shared `para_product_verified` helper applied to **both** paths; a failing factor returns a typed `NumericalFailure` instead of a silent result. | `para_product_verification_checks_magnitude_not_min_phase` (pins the reject behavior **and** the real-R limitation) |
| C3 | P1 | `src/synthesis/machiarella.rs` (`compute_s_product`, Type-II arm) | The Type-II branch monic-normalizes D by dividing by `leading_coefficient()`; `spectral.rs::monic` leaves a degenerate factor (|lead| < 1e-15) unnormalized, so the division scales the polynomial across ~200 decades and produces catastrophic roots. | Guard: a near-zero `d_lead` returns a typed `NumericalFailure`. | `type_ii_s_product_rejects_degenerate_leading_coefficient` |
| F1 | P1 | `src/pipeline/execution.rs` (`SpecSnapshot`) | Infinite transmission zeros and the lossless-`∞` Q were stored as non-finite `f64`, breaking the `serde_json` round trip (`NonFiniteValue`). | Stored as `Option<f64>` (finite); `None` ↔ `∞` on restore. | `save_and_load_context_round_trip_with_infinite_zero` |
| F2 | P1 | `src/pipeline/execution.rs` (`into_context`) | A hand-crafted/truncated snapshot could claim more stages than it has per-stage timings, panicking on an out-of-bounds `stage_timings_ms[i]` index. | Reject up front when `stage_timings_ms.len() < stages_executed.len()`; the loop also rejects a stage whose artifact is missing. | `into_context_rejects_mismatched_stage_timings` |
| F3 | P2 | `src/pipeline/execution.rs` (`ContextSnapshot`) | Dead `has_transform`/`has_response` metadata fields: written but never read by `into_context`. | Removed. Old snapshots that still carry the fields still deserialize — `ContextSnapshot` uses plain `Deserialize` (no `deny_unknown_fields`), so unknown fields are ignored. | — |
| F4 | P2 | `src/pipeline/schema.rs` | Input `topology` enum was missing `Wheel`; output `response` schema was an object-with-`samples` but the real serializer emits a **direct array**. | Input enum gains `Wheel`; output `response` retyped to a direct array. | — |
| F4′ | P3 | `src/pipeline/context.rs` (`PolynomialsSummary`), `src/pipeline/schema.rs` | The public `SynthesisContext` `Serialize` impl serialized `PolynomialsSummary.transmission_zeros_normalized` as a bare `Vec<f64>`; a zero at infinity (legal) lands there. P2-10 fixed the *spec* summary's finite→`null` mapping but missed the *polynomial* summary, so the output (serde_json emits `null` for `inf`) was inconsistent with the schema (`items: number`). | Mapped to `Vec<Option<f64>>` like the spec summary; schema items retyped to `["number", "null"]`. (Caught in the round-2 deep review; the earlier numeric probe confirmed `null` is emitted and the field is reachable.) | `context_serialize_emits_null_for_infinite_normalized_zero` |
| F5 | P1 | `src/pipeline/stages.rs` (`lossy_settings`) | A "lossy" request with a non-positive / non-finite `unloaded_q` (other than lossless `∞`) silently degraded to a lossless response (the backend only applies the term when `is_finite() && > 0.0`). | `validate_response_q` rejects such Q; a finite Q with a low-pass mapping (incommensurate grid) is also rejected. | `lossy_settings_rejects_invalid_q`, `lossy_settings_rejects_finite_q_with_lowpass` |
| F6 | P2 | `src/pipeline/stages.rs` (`ResponseEvaluationStage`) | The normalized (no-mapping) path evaluated lossless but silently dropped a supplied Q; the doc claimed a lossy capability the code didn't have. | The normalized path validates the Q (rejects invalid) and the doc now states the lossless normalized path. | `response_stage_rejects_invalid_q_on_normalized_path` |
| — | P2 | `src/spec/multi_band.rs` (`MultiBandSpec::new`) | Enforced `order > 0` but not the `MAX_ORDER = 256` cap the single-band `FilterSpec::new` applies (an unbounded order sizes an `(order+2)²` matrix → OOM vector). | Bounded by `MAX_ORDER`. | `multi_band_spec_order_is_bounded_by_max_order` |
| — | P2 | `src/bin/mfs_cli.rs` (`resolve_multiplexer_range`), `src/python.rs`, `src/design.rs` | The multiplexer default response window and the CLI/Python/`cascade` docs still used the arithmetic `center ∓ 1.5·bw` window, which goes **< 0 Hz** for wideband channels (fbw > 2/3) and fails in `map_hz_to_normalized`. | Derived from the geometric-mean `band_edge_pair` edges padded by 0.5·bw and floored positive — the same formula as `FilterDesign::default_response` (the round-1 P2-3 fix); docs/help updated. | `multiplexer_default_range_is_edge_based_and_positive` |

**Second-pass verification (2026-10-10):**

1. `cargo test --all-targets` → **411 passed, 0 failed, 0 ignored** (exit 0)
   (410 lib+integration+conformance + the CLI/Python-gated suites; +5 new
   second-pass regressions over the 406 round-1 baseline, net of the 1
   `spectral` test that pins the shared verification).
2. `cargo clippy --features cli,python --all-targets -- -D warnings` →
   **0 warnings** (exit 0).
3. `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` → **clean** (exit 0).
4. `cargo bench --no-run` → compiles (exit 0).
5. Non-finding — `group_delay_matches_across_fast_and_lu_backends` appeared to
   fail flakily during the round-2 fix wave (one `cargo test --all-targets` run:
   `fast group delay at ω=-1.32: 3.2338 vs S21 phase derivative 3.2160`, tol
   0.0132). Ruled out as a numerical flake on 2026-10-10 by: (a) a margin scan
   over all 799 Check-2 grid points — worst gap is 49% of tolerance (ω=−1.105),
   Check 1's worst gap is 1.3e-12 against a 1e-7-class tolerance; (b) 20
   bit-identical re-runs of the full `synthesize → LU → pole-expansion` chain
   (no run-to-run drift); (c) 10/10 passes on the same binary under 12-process
   CPU load. Root cause: the `target/` debug test binary had become a **patch
   build mixing `.rcgu` objects from three source states** (cargo reuses one
   output filename per crate and re-links only changed objects), so the failing
   run executed stale `response`/`synthesis` code — the reported 3.23 values
   are not what the current source produces at that point (4.10). A forced
   clean rebuild (`rm` the binary + its objects, `cargo test --no-run`) restored
   a consistent binary; the test then passed 10/10 and under load. No code
   change; lesson: after reverting/re-applying edits, verify against a fresh
   binary before diagnosing a "flake".

**Third-pass closure (2026-10-11, §I):** §C risks 1–7 and §D items 3/5/6/7
were the last open follow-ups; all are now closed (§I). The remaining open
items after this session: §D.1/.2/.4/.8 (already satisfied at second pass)
are unchanged; no confirmed defect remains open.

**Status at 2026-10-10 (superseded by §I):** no confirmed defect remains open

- **§C risks 1–7: open** (not closed in-repo; no known production trigger):
  faer backward-residual gate width (1e-3), companion-fallback para-product
  validation + `select_lhp_roots` RHP fill (C2 closed the *fallback*
  para-product check, but the `select_lhp_roots` last-resort fill still has no
  `re < 0` guard — verify whether C2's shared `para_product_verified` helper
  should also gate it), high-order normal-equation conditioning (perf only),
  external-Q convention (g vs √g, no order-1 physical anchor), full-path
  `FullIterationState` one-iteration lag (≤ STAGNATION_ACCEPT), `to_proto`
  |r|≥50 misleading error message, CLI `optimize`/`extract` Hz-as-normalized
  grid contract (documented, silent by design). **§C risk 8: closed** by the
  P2-9 fix.
- **§D recommendations:** .1 band-edge single source — effectively satisfied
  (P2-2/P2-3 route through `band_edge_pair`); .2 consumer-side symmetry —
  satisfied (P1-3, backend `assert_symmetric`); .4 pole-domain convention —
  satisfied (P2-12, `from_polynomials` deleted); .8 request validator —
  satisfied (dead `schema::validate_request` deleted, P3 #10 batch). Still
  open: .3 (verify has no absolute spec check — "verify passes" ≠ "meets
  design target"), .5 (split 1854-line `mfs_cli.rs` by subcommand), .6
  (`save_context`/`load_context` closed-loop e2e in CI), .7 (three hand-rolled
  dB variants partially unified by P2-10/F4′ — context.rs and the CLI table
  still hand-roll vs `ResponseSample::s*_db()`; confirm no residual
  `log10(0)` path remains).
- **Delivery risk (largest open item):** the entire fix wave — all of the
  above plus F4′ — is **uncommitted** (34 modified files, +2414/−606 vs
  `5fa5c68`), and `review-findings-0.6.1.md` + `scripts/build_wheel.sh` are
  untracked. A working-tree loss (or a `git clean`/`checkout`) would drop the
  whole wave. Recommend committing in logical units (round-1 fixes, round-2
  fixes, F4′) before further work.
- **Not executed (carried from §F):** benchmark conformance tiers (bench-data),
  actual Python wheel build, non-Linux platforms, `python`-feature runtime
  (FFI statically reviewed only).


# I. Third-pass closure (2026-10-11)

Five parallel slices closed the last open items (§C.1–C.7, §D.3/D.5/D.6/D.7,
and the four §F "not executed" legs). All changes are on top of the committed
round 1+2 wave (`8f3191a`).

## I.1 §C risks 1–7 — verified closed

| # | Location | Outcome | Evidence |
|---|---|---|---|
| C.1 | `complex_poly.rs` faer backward-residual gate | **VERIFIED-NO-ISSUE** — the gate is `max_backward ≤ 1e-3` accept; residuals in (1e-3, 1e-2) are **rejected** → fall through to Durand-Kerner/Aberth (safe direction). Production polynomials probed orders 4..40: backward residual ≤ ~1e-5, none in the band. | `faer_gate_does_not_fire_on_production_polynomials_up_to_order_40` (PASS); comment at `complex_poly.rs:770` |
| C.2 | `spectral.rs` `select_lhp_roots` last-resort fill | **FIXED** — the fill now admits only strictly-LHP roots (`re < -1e-9`); a short LHP pass can no longer pull in an RHP or pure-imaginary root. The companion fallback's para-product backstop (C2, round 2) stays as the magnitude check. | `select_lhp_roots_fill_rejects_rhp_when_lhp_short` (PASS) |
| C.3 | `spectral.rs` high-order normal-equation conditioning | **VERIFIED-NO-ISSUE (perf only)** — at order 40 the cepstrum path falls through to the companion fallback; the fallback is correct (O(N log N) intent lost, no correctness/timeout issue). Measured 50.8 ms at order 40. | Measurement only; no code change |
| C.4 | `matrix/scaling.rs` external-Q convention | **VERIFIED-NO-ISSUE (anchor added)** — the `g²` convention is self-consistent: order-1 RL=20 BP `Qse = 2·f0/BW3dB` matches `1/(M01²·Δ)` to <1%; the "standard" `1/(M01·Δ)` is off by `M01` (rel diff 1.23). The g² convention is now pinned by a physical order-1 test. | `external_q_convention_matches_physical_order1_bandpass` (PASS) |
| C.5 | `machiarella.rs` `FullIterationState` one-iteration lag | **VERIFIED-NO-ISSUE (documented)** — the lag is real and bounded by the convergence error (< `STAGNATION_ACCEPT` 1e-3); it is now documented at the capture site. Fixing it would change convergence behavior and is not worth the risk. | Comment at `machiarella.rs:656` |
| C.6 | `machiarella.rs` `to_proto` misleading error | **FIXED** — the "non-finite root set" error now names the actual cause (out-of-range finite root |r| ≥ 50). The pre-existing message-assertion test was aligned to the new message. | `machiarella.rs:608`; existing test aligned |
| C.7 | `mfs_cli` optimize/extract grid contract | **FIXED (guard + docs)** — the normalized contract is documented, and a magnitude guard rejects physical-Hz grids (span > 1e6 → `PreconditionViolation("NORMALIZED…")`). Feeding a real .s2p no longer "converges" at ±6e9 normalized frequency. | `response_grid_rejects_physical_hz_and_accepts_normalized` (PASS); `optimize.rs:122-150`, `main.rs:176-119` |

## I.2 §D architecture items

- **§D.3 (verify absolute spec checks) — FIXED.** `verify::verify_against_spec`
  now checks a response **against a `FilterSpec`** (not just two responses
  against each other): passband return loss (min passband RL ≥ spec − 0.1 dB
  margin), passband insertion loss (|S21| in [−1, +0.5] dB), transmission-zero
  depth (≤ −40 dB within ±2 grid steps; infinite zeros → Skipped), and
  stopband rejection (always Skipped — `FilterSpec` carries no normalized
  stopband requirement; documented). Missing data → `Skipped`, never a silent
  `Pass` (P3 #2 semantics). Exposed in `lib.rs`/`prelude.rs` alongside
  `CheckOutcome`.
- **§D.5 (split `mfs_cli.rs`) — FIXED.** The 1911-line `src/bin/mfs_cli.rs`
  is split into `src/bin/mfs_cli/` — `main.rs` (entry/dispatch),
  `design.rs`, `multiplexer.rs`, `optimize.rs`, `target.rs` (shared
  `parse_response_json`/`parse_matrix_json`/`load_target`), `io.rs`
  (error/tables/formats). Line-by-line move, no CLI flag/help/output change.
- **§D.7 (dB single source) — FIXED (CLI + new verify code).** All five
  hand-rolled `10·log10(re²+im²)` sites in the CLI route through
  `ResponseSample::s11_db()`/`s21_db()` (−300 dB floor at zero magnitude —
  no more `null` at a true transmission zero). The new `verify_against_spec`
  code reuses the same `s*_db()` methods (a private `db_from_mag` it first
  shipped was deduplicated onto them in the review follow-up). **Scope
  note:** the raw-array `min_db`/`max_db` in `multiplexer.rs` operates on
  library-precomputed `Vec<f64>` dB (no `ResponseSample`), left as-is with a
  §D.7 citation. Pre-existing library-side hand-rolled `mag → dB` lambdas
  (`duplexer.rs:101/113/276`, `multiplexer.rs:462/958/979/1028`,
  `machiarella.rs:754`) are **outside this commit's scope** — they were there
  before and were not unified; a full §D.7 closure would route them through
  the same `s*_db()`/shared helper too (open, non-blocking).
- **§D.6 (save/resume closed-loop e2e in CI) — FIXED.** New
  `tests/cli_save_resume_e2e.rs` (gated `cli` feature) runs the permanent
  successor to the deleted P3 #11 harness: `mfs_cli --input` (full pipeline)
  vs `mfs::pipeline::save_context` + `mfs_cli --resume` must reproduce the
  same response — non-zero samples within 1e-9 dB, transmission-zero samples
  as deep S21 nulls in both runs. Zero-point is a deep-null assertion
  (≤ −60 dB S21 + near-total S11), not bit-equal dB (20·log10 amplifies the
  ~1e-17 residual into a large dB swing). The exact CI command
  (`cargo test --features cli --test cli_integration --test cli_save_resume_e2e`)
  is wired into `.github/workflows/ci.yml`.

## I.3 §F "not executed" legs

- **bench-data conformance tiers — EXECUTED (green).**
  `cargo test --test benchmark_standard` → 20/0 (both the sibling byte-for-byte
  tier and the manifest sha256 tier ran; the sibling is present in this
  environment, correcting the "only tier available" framing).
  `folded_response_invariance` → 10/0.
- **Python wheel build — PARTIAL → FIXED (script repaired).** The maturin
  wheel builds and installs (k12 = 0.414467071, bit-identical to the Rust
  probe). The script's wheel-ID glob aborted on maturin's `linux_x86_64`
  output (no `patchelf` on this WSL2 host); `scripts/build_wheel.sh` is now
  patched to accept and repair the untagged wheel. **Docker-based glibc-2.36
  verification was blocked by a dead containerd daemon in this environment** —
  the repaired wheel's max GLIBC requirement is statically verified as
  `GLIBC_2.35` (verneed nodes), so it is installable/importable on glibc ≥ 2.35.
- **`python`-feature runtime — EXECUTED (no longer "statically reviewed only").**
  `import mfs` + `bandpass`/`duplexer`/`prototype` all run through pyo3 at
  runtime; `extract`/`bandpass_edges` present and callable.
- **non-Linux platforms — NOT-EXECUTABLE-HERE (CI coverage).** No
  macOS/Windows toolchain in WSL2; `ci.yml` runs the 3-OS × 2-rust matrix for
  test and cli jobs. **Gap noted:** `ci.yml` has no job that builds the Python
  wheel or runs the pyo3 runtime (a delivery risk worth a future CI job).

## I.4 Third-pass verification

Recorded in the commit message and re-run at the end of this session:
`cargo test --all-targets` (all green),
`cargo clippy --all-targets -- -D warnings` (0),
`cargo clippy --features cli,python --all-targets -- -D warnings` (0),
`cargo doc --no-deps` (clean), and the exact CI command for the save/resume
gate.
