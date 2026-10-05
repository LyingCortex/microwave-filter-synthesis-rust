# MFS v0.6.0 Whole-Crate Code Review — Findings

**Scope:** Whole-crate defect hunt at the 0.6.0 release boundary (`master` @ `62e5dce`), all of `src/` plus integration tests.
**Method:** Five parallel read-only review slices (S1 approx/frequency, S2 matrix/transform, S3 synthesis/duplexer/multiplexer, S4 response/spectral, S5 API surfaces: CLI/Python/pipeline/output/spec/design) with evidence-backed findings, then an independent adjudication pass: every finding re-derived against source at `file:line`, and the single highest-stakes numeric claim (normalized lossy Q) independently re-derived analytically. S1's agent failed on a model length-stop; S1 was re-reviewed manually (see §E).
**Companion docs:** `architecture-review-2026-09-25.md` (prior wave), `continuation-plan.md` (roadmap), `CHANGELOG.md` [Unreleased].

This is an **assessment**, not a code change. **Status (2026-10-06):** the
"silent errors" wave (§F.1) — P1-1/2/3/4/6, P2-1/2/3/5, all six P3 — is **fixed
and gate-anchored** on `master`; each P1 fix ships with the value-level gate
that anchors the previously-unchecked observable (per-channel F zeros,
physical detuning vs `freq.rs`, group-delay value vs the numerical
derivative, p0² sign, array lengths). **P1-5 + P2-4 (the CT/CQ cascade
re-implementation, §F.2) remain open.**

---

## A. Verdict

**The core numerical architecture is sound; the release is fit for its gated scenarios. The residual risk is a consistent class of *silent* defects: paths where the shipped gates anchor on junction-level observables (N/D/Pt/Pr, p0², junction S21) or on response *magnitude*, so per-channel products and secondary observables can be wrong while every gate stays green.**

Findings tally after adjudication: **18 reported → 17 valid (6×P1, 5×P2, 6×P3) + 1 false positive.**

The six P1s are all "gate-green but wrong result":

1. Full-iteration per-channel recovery corrupts the a priori reflection zeros (`recovery.rs`).
2. Band-pass resonator-frequency denormalization is sign-mirrored vs the canonical mapping (`scaling.rs`).
3. LU-backend group delay includes port nodes (`response/backend.rs`).
4. Overlapping-channel p0² negatives are silently clamped instead of rejected (`machiarella.rs`).
5. `to_cascaded_quadruplets` does not merge anything (`cascade.rs`).
6. Python `tune()` can panic the interpreter on short companion arrays (`python.rs`).

Clean areas (explicitly verified, no findings): `src/approx/` (root-finding fallback chain, generalized Chebyshev helpers, multi-band LM solve), `src/freq.rs` (band-pass/high-pass/band-stop mappings are self-consistent round-trips with documented singular clamping), `src/error.rs` taxonomy, `src/touchstone.rs` (per-port reference impedance, verified against touchstone 0.14.5), `src/spec/` builder validation, pipeline stage predecessor checks and context resume path, `src/matrix/core.rs` invariants (Givens similarity is correct to ~1e-18 in every response-invariance test).

---

## B. P1 — silent wrong results

### P1-1 · `recover_channel_design` rescales the full-path D and destroys the fixed reflection zeros

- **File:** `src/synthesis/recovery.rs:42-50` (call sites: `src/duplexer.rs:700-705, 720-733`, `src/multiplexer.rs:178-188`)
- **Mechanism:** Both full-path callers pass `d_original = S_converged − F⁰` (S monic degree n, F⁰ monic degree n ⇒ D has degree n−1 with a **non-monic** leading coefficient). `recover_channel_design` unconditionally rescales `d_original` by `s_leading / d_leading`, which is a no-op only when `d_leading ≈ 1`. Measured on the §V golden spec (order 7, 8 iters, conv 3.5e-10): `|lead(S−F⁰)| = 1.2952`, rescale factor `0.7721`, and the recovered `F = S − (S−F⁰)·0.7721` deviates from F⁰ by max rel-err ~1.6e30 — the a priori reflection zeros are destroyed for **every** full-iteration run.
- **Why gates stay green:** junction observables come directly from `JunctionCommon`; per-channel `FilterDesign` (duplexer `tx`/`rx`, multiplexer `channels`) is built from the corrupted E/F but is not compared against the original F zeros by any gate.
- **History:** the WIP branch (commit `1779bf9`) passed `rescale_d = !full_iteration` with the doc "rescaling would corrupt the exact reconstruction"; commit `818d38c` dropped the flag justifying it as "a no-op for the monic S−F⁰ D" — that premise is false (S−F⁰ is *not* monic).
- **Correct behavior:** the simplified path is fine — there `d_original = (E−F)/2` has a cancelling (~0) leading coefficient and the existing `d_leading.norm() > 1e-15` guard already skips the rescale. The full path must use `D = S − F` **without rescaling** (normalization-independent; `F = S − D` then holds exactly).
- **Suggested fix:** add a `rescale_d: bool` parameter (full path passes `false`) or branch on degree; do not scale a D already expressed in S's monic normalization.

### P1-2 · `denormalize_resonator_frequency` mirrors the detuning sign vs `BandPassMapping`

- **File:** `src/matrix/scaling.rs:108-113`
- **Mechanism:** the inverse band-pass map is computed with the linear term **subtracted**: `center·(√(1+(n·fBW/2)²) − n·fBW/2)`. The authoritative `BandPassMapping::map_normalized_to_hz` (`src/freq.rs:178`) adds it: `(n·bw + √((n·bw)² + 4·center²))/2`. Numerically `scaling.denorm(n) == freq.denorm(−n)` exactly — a positive normalized detuning (resonator above center) is placed *below* center, so all physical resonator rest frequencies are mirrored.
- **Why gates stay green:** `bandpass_scaling_round_trips_internal_couplings` (`src/matrix/mod.rs:68-72`) passes because `normalize_resonator_frequency` (`scaling.rs:115-127`) is the self-consistent inverse of the wrong denorm (both sign-flipped vs `freq.rs`); the detuning is never compared against the canonical convention used everywhere else (`edges_to_center_bandwidth`, `band_edge_pair`, `FilterDesign::bandpass_edges`).
- **Reachable via:** `denormalize_bandpass` / `denormalize_bandpass_with_external_q` (`src/synthesis/engine.rs:90,101`, `src/design.rs:609,620`).
- **Suggested fix:** use `(normalized*bandwidth_hz + sqrt((normalized*bandwidth_hz)^2 + 4*center_hz^2))/2.0` and route `normalize_resonator_frequency` through `BandPassMapping::map_hz_to_normalized` so the two cannot drift; add a cross-check test against `freq.rs`.

### P1-3 · LU-backend group delay sums over port nodes

- **File:** `src/response/backend.rs:64-66`
- **Mechanism:** group delay is `−d arg(S21)/dω`. With `A(ω) = M + jω I_r − jR_S e₀e₀ᵀ − jR_L eₗeₗᵀ`, `dA/dω` is the identity on the **resonator block only**; by `A⁻¹′ = −A⁻¹ A′ A⁻¹` only resonator indices contribute. The numerator folds `(0..side)`, adding the spurious `Im(A⁻¹[0,0] + A⁻¹[L,L])`.
- **Evidence (reviewer, independent numpy central-difference):** order-2 test matrix at ω=0: code formula 2.67873 vs true 0.90498; the resonator-only sum `(1..side−1)` matches to 1e-10. Same O(1) error (1.0–1.9) across orders 2/4/6, lossless and Qu=300.
- **Why gates stay green:** group-delay tests assert `is_some()/is_none()/finiteness` only, never the value; the wrong number is exposed via Python `response()`, `src/output/format.rs`, and `src/verify/mod.rs`.
- **Suggested fix:** `let numerator = (1..side - 1).fold(...)` and add a value-level test against the numerical derivative.

### P1-4 · Overlapping-channel p0² clamped silently

- **File:** `src/synthesis/machiarella.rs:354-360 (N=2), 381 (N≥3)`
- **Mechanism:** the doc claims diagonal dominance, which holds only for well-separated channels. For adjacent/overlapping channels the p0² system has no all-positive solution; the LU solve yields a negative component that is clamped by `.max(1e-15)`. Measured on an overlapping N=3 spec (200 MHz bands, 100 MHz center spacing): `p0_sq = [1.112e-4, 1.000e-15, 9.039e-5]` — the middle value is exactly the floor (clamped from negative) with `converged = true`, i.e. that channel's junction S21 ≈ −160 dB is a silently corrupted design rather than a flagged non-physical configuration. The N=2 Cramer branch and its equal-split degenerate fallback share the pattern.
- **Why gates stay green:** all gated scenarios (well-separated N=2, non-overlapping N=3/4, Zhao N=4) return positive p0².
- **Suggested fix:** after the solve, reject any `p0² < −tol` with a typed `PreconditionViolation` (non-physical overlap must fail loudly); keep the floor only for the ~0 round-off case.

### P1-5 · `to_cascaded_quadruplets` "merge" is a no-op (CQ == CT byte-identical)

- **File:** `src/matrix/cascade.rs:104-122`
- **Mechanism:** the merge applies `diagonal_rotation_angle` (a 2×2 diagonal-equalizing Jacobi step) to the two trisection centers. For symmetric ±ω₀ pairs the two diagonals are equal ⇒ angle 0 ⇒ **no rotation**; the output is byte-identical to `to_cascaded_trisections` (reviewer ran the crate: N=6, zeros `[-1.5, 1.5, −2.0, 2.0]`, max|diff| = 0.0). For asymmetric pairs it applies the wrong operation (max|diff| ~0.3 vs CT) — still not a Cameron quadruplet. Contradicts the doc at `cascade.rs:87` and the comment at `:116-117`.
- **Public API:** `FilterDesign::to_cq` (`src/design.rs:355-357`).
- **Related (P2-4 below):** the CT "pull" loop in the same file also deviates from the verified extraction path, so *neither* cascade output is a realizable CT/CQ topology today (response is preserved — everything is a similarity transform — but the sparsity structure is wrong).
- **Suggested fix (sub-project):** implement the Cameron quadruplet merge (the tail-pivot rotation `extract_quadruplet` uses, applied across the two trisection centers) or route CQ through the verified `extract_quadruplet` path; add structural tests (CQ ≠ CT; trisection trident pattern present; mainline connected). See §D for the recommended wave split.

### P1-6 · Python `FilterDesign.tune` can panic the interpreter

- **File:** `src/python.rs:410-428`
- **Mechanism:** only `freqs.len() == target_s21_re.len()` is checked; the sample loop then indexes `target_s11_re[i]/target_s11_im[i]/target_s21_im[i]` over `0..freqs.len()`. A shorter s11 array → Rust index-out-of-bounds into Python; PyO3 does not catch Rust panics, so the interpreter crashes instead of raising `PyValueError`.
- **Crate convention proves the gap:** the sibling `extract()` (`src/python.rs:611-628`) validates all four companion arrays first.
- **Suggested fix:** extend the line-410 check to all five arrays with the same `PyValueError`.

---

## C. P2 / P3

### P2-1 · CLI `optimize` with a `.s2p` target cannot supply the matrix

- **File:** `src/bin/mfs_cli.rs:1041-1050`
- When `--target` is a `.s2p`, `json_input` is set to `String::new()`, then `parse_matrix_json(&json_input)` fails with "invalid matrix JSON" — the documented usage (matrix from `--input`/stdin + Touchstone target) is impossible. `run_extract` handles the same case correctly (reads JSON input unconditionally, line ~1109).
- **Fix:** always `read_input(cli)?` for the matrix; pass only the target path to `load_target`.

### P2-2 · CLI `design` silently mis-parses center-without-bandwidth

- **File:** `src/bin/mfs_cli.rs:432-458`
- `design -n 4 -r 20 --center 6.75e9` (no `--bandwidth`, no edges) falls into the wildcard arm returning `(Some(center), None)`; the `if let (Some(c), Some(bw))` guard then fails and a **normalized prototype** is synthesized — the user's center frequency is silently ignored.
- **Fix:** typed `InputParse` error for `(Some(_), None)` when start/stop were not both given.

### P2-3 · Pipeline schema/validator reject the API's topology case

- **File:** `src/pipeline/schema.rs:66-69, ~293`
- The generated JSON Schema and `validate_request` accept only lowercase `[transversal, folded, arrow]`, but the real request path deserializes through `MatrixTopology` (PascalCase, no serde rename), so `topology: "Folded"` — used by `tests/cli_integration.rs::full_request_json` and the `request.rs` doc example — is valid per the API and rejected by both the schema and the validator.
- **Fix:** make schema, validator, and serde agree on one case (or accept both) and fix the `request.rs` doc example.

### P2-4 · `to_cascaded_trisections` pull severs the mainline

- **File:** `src/matrix/cascade.rs:62-75`
- The pull rotates on `matrix.get(pos, side−1)` (fixed load column) instead of the shifting resonator diagonal. N=6 with 3 zeros: mainline `M[k+1][k+2] = [−0.88, 0.0, 0.40, 0.28, −0.82]` (bridge at index 1 severed), CT diagonals all 0 for symmetric zeros. The response-invariant **and** structurally-verified `extract_trisection` path on the same matrix yields a continuous mainline `[−0.88, −0.58, −0.56, −0.33, −0.46]`, matching the Cameron 2011 fixture convention (`tests/cameron2011_fixtures.json`). `ct_synthesis_produces_valid_matrix` only checks S21 magnitude, so the non-canonical topology is undetected.
- **Fix:** reuse the verified `extract_trisection` pull per zero; assert CT mainline connectivity.

### P2-5 · `DuplexerSynthesisResult` drops the shared loop's `converged` flag

- **File:** `src/duplexer.rs:227-245, 608-610`
- `run_full_iteration` returns `FullIterationState { converged, .. }` (`src/synthesis/machiarella.rs:128`); `MultiplexerSynthesis::synthesize` forwards it (`src/multiplexer.rs:195`) and the B3 gate asserts it, but `DuplexerSynthesisResult` has no `converged` field. A duplexer caller with a too-small `max_iterations` must heuristically compare `convergence_error` against the tolerance to distinguish convergence from budget exhaustion.
- **Fix:** add `pub converged: bool` and forward `state.converged`, matching the multiplexer.

### P3-1 · `FilterDesignBuilder::synthesize` masks missing bandwidth

- **File:** `src/design.rs:77-81` (pattern repeated at ~373, ~557, ~610)
- `bandwidth_hz.unwrap_or(0.0)` → `BandPassMapping::new` reports "bandwidth must be > 0, got 0" instead of saying bandwidth was never supplied.
- **Fix:** `ok_or_else` naming the missing parameter, matching the center check.

### P3-2 · Python `tune` rebuilds the grid as a uniform linspace

- **File:** `src/python.rs:432-437`
- `FrequencyGrid::linspace(freqs.first(), freqs.last(), freqs.len())` ignores the user's actual (possibly non-uniform) `freqs`; the optimizer then zips residuals against targets at the wrong frequencies — a silent mis-fit. `extract()` shares the construction but its docstring states the uniform-grid contract explicitly; `tune`'s does not.
- **Fix:** document the uniform-ascending-grid contract (or evaluate at the actual sample frequencies if the grid type allows).

### P3-3 · `FilterDesign::bandpass_edges` performs no validation

- **File:** `src/design.rs:191-196`
- `f_start_hz <= 0` ⇒ geometric-mean center is NaN; the input error surfaces one stage later as "center frequency must be > 0, got NaN". The CLI path guards it (`run_design` checks `f2 <= f1`), but the public/prelude/Python API does not.
- **Fix:** validate `0 < f_start < f_stop` at the constructor.

### P3-4 · `from_polynomials_and_matrix` double-unwrap / silent spec substitution

- **File:** `src/design.rs:262-263`
- `FilterSpec::new(order, rl).unwrap_or_else(|_| FilterSpec::new(order, 20.0).unwrap())` — panic if `order == 0`; and silently stores RL=20.0 when the supplied `return_loss_db` is invalid, so `order()/spec().return_loss_db()` can disagree with the polynomials/matrix. Sole caller is `recover_channel_design` (guarded today by `MultiplexerSynthesis` rejecting order 0, `src/multiplexer.rs:140`).
- **Fix:** drop the unwrap fallback; propagate the error.

### P3-5 · `PoleExpansionData::from_polynomials` lacks a y22 length check

- **File:** `src/response/pole_expansion.rs:81-105` (consumed at `src/synthesis/residues/transversal.rs:19-28`)
- Validates `y11` vs `y12` residue counts but not `y22`; `evaluate_at` indexes y22 over `poles.len()`, so a mismatch panics or drops poles. **Informational:** no production caller invokes `from_polynomials` (all use `from_matrix`); pre-dates the release window.
- **Fix:** add the `y22` length check (cheap, removes a future panic).

### P3-6 · Givens-similarity step-2 comment has flipped off-diagonal signs

- **File:** `src/matrix/core.rs:301-303`
- The code (304-308) correctly implements `M'' = M'·Rᵀ`; the comment states the signs flipped, describing `R·M·R`. Every response-invariance test passes to ~1e-18 because the *code* is right. A maintainer "correcting" the code to match the comment would silently break response invariance.
- **Fix:** comment-only rewrite to `col_a'' = c·col_a' − s·col_b'; col_b'' = s·col_a' + c·col_b'`.

---

## D. Adjudicated false positive (1)

**`src/design.rs:560` `response_lossy` "inverted Q formula" (S5, reported P1).**
The reported claim: `normalized_dissipation_q = unloaded_q * bw / center` is the reciprocal of the documented `δ = f0/(Qu·BW)`, understating loss by ~(Qu·fBW)². **Re-derivation:** the backend applies loss as `δ = 1/unloaded_q` on the resonator diagonals (`src/response/backend.rs:122-126`), so the intended normalized dissipation is `δ = center/(Qu·bw)`, i.e. `Qu_norm = Qu·bw/center` and `δ = f0/(Qu·BW)` — exactly the textbook normalized dissipation and exactly what the in-file comment says. The code is correct; the reviewer inverted the relation between "Q" and "dissipation". **Verdict: false positive, no change.**

---

## E. Slice S1 (approx/frequency) — manual re-review, clean

S1's review agent failed (model length-stop after 1h30m, no structured output). Re-reviewed manually:

- `multi_band.rs`: LM solve uses **successive-cost acceptance** (correct — the D2-opt `optimize_matrix` best-cost early-termination bug does not exist here), enforces in-band ordering/feasibility, `lambda` capped at 1e12, and explicitly rejects the overdetermined K≥3 case with a typed error. Oracle test `equiripple_solve_recovers_fixture_rz` recovers the Cameron 2011 Fig.-25 reflection zeros to 5e-4; `fixture_polynomials_are_equiripple_at_specified_rl` anchors the RL.
- `freq.rs`: band-pass / high-pass / band-stop mappings are algebraic inverses of each other (self-consistent round-trips); singular points (f=0, f=f0) clamp to `±SINGULAR_CLAMP` with documented rationale; `map_hz_to_normalized` / `map_normalized_to_hz` validation is complete.
- `complex_poly.rs`: three-stage root-finding fallback chain (faer companion eigen → Durand–Kerner → Aberth) with residual-check acceptance; nalgebra companion deliberately excluded (documented complex-branch defect).
- Literature fixtures (Cameron 1999/2011, Zhao) + `tests/synthesis_regression.rs` + `tests/synthesis_properties.rs` anchor the layer extensively.

**No new findings in S1.**

---

## F. Recommended remediation split (for a future wave)

1. **Wave "silent errors" (local, low ambiguity):** P1-1 (recovery flag), P1-2 (scaling sign + cross-check test), P1-3 (group-delay range + value test), P1-4 (p0² rejection), P1-6 (Python length checks), P2-1..P2-3, P2-5, all P3s. Each is a small change with a targeted regression test; the recovery and scaling fixes each need one new gate comparing the corrupted observable against its a priori / canonical value (that comparison is precisely what the green gates lack).
2. **Wave "CT/CQ cascade" (re-implementation):** P1-5 + P2-4 together — route both through the verified `extract_trisection`/`extract_quadruplet` machinery (the design is already specified in `docs/ct-cq-topology-design.md`) plus structural sparsity tests (trident pattern, connected mainline, CQ ≠ CT).
3. **Gate-gap note:** five of the six P1s are invisible to the current gates because gates anchor on junction observables and S21 magnitude. When fixing, each P1 fix should ship with a gate that evaluates the *previously unanchored* observable (per-channel F zeros vs a priori, physical detuning vs `freq.rs`, group-delay value vs numerical derivative, p0² sign, CQ/CT sparsity).
