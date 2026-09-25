# Work Plan — MFS Architecture Remediation

**Plan id:** `arch-remediation-20260925`
**Source of truth:** `docs/internal/architecture-review-2026-09-25.md` (every P-number below maps to that review).
**Repo:** `/mnt/c/Users/eynulai/dev/mfs` (Rust, single crate, edition 2021, MSRV 1.75).
**Goal:** Close the architectural gaps from the review in a safe, TDD-verified,
reversible sequence — fix the silent correctness defects first, remove the
velocity tax second, clean up last. No behavior regression at any gate.

This plan is **decision-complete**: a worker executes the waves in order with
no further interview. Each task names its category, skills, target file(s), the
atomic change, and a binary verification (RED→GREEN and/or a real-surface
artifact). Tasks in the same wave are independent and may run in parallel.

---

## Ground rules (bind every task)

- **TDD floor:** any behavior change follows RED → GREEN → (SURFACE). A
  production edit without a preceding failing test recorded in the notepad is a
  rejection; revert and redo. Exemptions (pure doc, version bump, comment,
  dead-code *deletion* with a pinned test already green) must be justified in
  `## Findings`.
- **Never delete/skip a failing test.** Fix code, not tests.
- **Gate before proceeding:** after each wave, `cargo test --all-targets` +
  `cargo clippy --all-targets -- -D warnings` + `cargo doc --no-deps` must be
  green (or note pre-existing failures, unchanged).
- **Commits:** one atomic commit per verified increment. Mimic repo style
  (subjects like `fix: ...`, `test: ...`, `docs: ...`). Load `git-master`.
- **Notepad:** durable at the path captured at start; append-only.

---

## Scenario contract (the binding acceptance set)

These are the observable outcomes that must ALL pass before the plan is "done."
Each has a binary check + the real-surface artifact.

| ID | Scenario (class) | Pass condition (binary) | Real-surface artifact |
|----|------------------|--------------------------|------------------------|
| **S1** | N-channel mutual loading now propagates (Happy+correctness) | `cargo test --test multiplexer` green AND a new `separation_limit` test asserts far-apart channels converge to independent synthesis (max\|Δg\|<1e-9) AND an `effect` test asserts adjacent channels differ from independent synthesis | `cargo test --test multiplexer -- --nocapture` transcript showing the two new tests PASS |
| **S2** | Engine sensitivity (regression guard for the P1 class) | New unit test: perturbing `e_s` roots changes `MatrixSynthesisEngine` output matrix (max\|Δ\|>tol) | `cargo test -p mfs --lib` transcript, test id + PASS |
| **S3** | Orchestration is single-sourced (P2) | `design.rs::synthesize`, `pipeline::run_full_pipeline`, `pipeline::run_stage`, and CLI all call one `orchestrate_*` core; `grep` shows the FilterSpec-build block exists in exactly 1 place | `grep -rn "FilterSpec::new" src` → 1 canonical definition + thin call sites; full `cargo test` green |
| **S4** | Band-center convention unified (P2/P9) | CLI `design` subcommand and `FilterDesign::bandpass_edges` produce identical center for the same edges | A small integration test asserting equal centers (or a doctest); `cargo test` green |
| **S5** | Pipeline `mapping` no longer a silent trap (P3) | Either (a) physical-domain evaluation implemented + a test asserting a physical-Hz grid yields the Hz-anchored response, or (b) an unrecognized/`mapping` key is rejected with a typed `MfsError` + a test asserting the error | `cargo test --test cli_integration` (or pipeline test) transcript |
| **S6** | Legacy free functions deprecated (F) | `lib.rs` root fns carry `#[deprecated(note=...)]`; prelude still re-exports them; `cargo build` emits deprecation warnings only for legacy callers (none in-tree) | `cargo build 2>&1` shows zero in-tree deprecation warnings; `cargo test` green |
| **S7** | Numeric stack de-duplicated (P4) | Exactly one companion-eigen root solver in the production `AdaptiveRootSolver` chain; the other demoted to a test-only cross-validation oracle; rustfft removed OR gated behind a feature (and the `use_full_iteration` path wired to it) | `grep -rn "faer\|CompanionMatrixRootSolver" src` shows the split; `cargo test` green; if bench added, `cargo bench` runs |
| **S8** | No silent degradation (cross-cutting policy) | `group_delay` on the pole-expansion fast path either computes or is explicitly flagged (e.g. `ResponseSample.group_delay` + a `degraded` marker), never a silent 0.0; pipeline `mapping` and multiperter full-iteration are explicit | Code inspection + a test asserting the flag/marker is set when the fast path zeroes group delay |
| **S9** | Docs + version aligned (P5/P6) | `docs/architecture.md` module list matches `src/` (all root modules present, incl. duplexer/multiplexer/optimize/spectral/touchstone/error/python); `Cargo.toml` version == `pyproject.toml` == README example | `cargo doc --no-deps` green; a CI/local check that the three versions match; diff of architecture.md module list vs `ls src` |
| **S10** | Full regression + CI parity (final gate) | `cargo test --all-targets` green; `cargo clippy --all-targets -- -D warnings` green; `cargo doc --no-deps` green (matches `.github/workflows/ci.yml`) | The three commands' exit-0 output captured |

**Stop condition:** I stop right away when S1–S10 are all PASS with both the
RED→GREEN proof and the real-surface artifact captured, the notepad is
updated, and `## Todo` is empty.

---

## Notepad (initialize at start, append-only)

```
# Ultrawork Notepad — MFS architecture remediation
Started: <ISO>
Plan: .omo/plans/arch-remediation-20260925.md
## Now
## Todo
## Findings   (incl. TDD-exemption justifications)
## Evidence   (S1..S10 artifacts, file:line)
## Learnings
```

---

## Wave 0 — Baseline & safety net  (serial, must finish first)

**T0.1 — Capture baseline.**  `category: quick`, skills: []
- Action: run `cargo test --all-targets`, `cargo clippy --all-targets -- -D warnings`, `cargo doc --no-deps`. Record pass/fail counts and any pre-existing failures verbatim into the notepad.
- Verify: notepad `## Findings` has a "BASELINE:" block with exact counts + the list of pre-existing failures (so later waves don't get blamed for them).
- Why first: every later "green" is judged against this baseline.

---

## Wave 1 — Independent, parallel (3 tasks)

> All three are independent of each other and of Wave 2. Run concurrently.

**T1.1 — P5/P6 stop-the-bleeding docs + version sync.**  `category: writing`, skills: []
- Action:
  1. In `docs/architecture.md`: add the 7 omitted root modules (duplexer, multiplexer, optimize, spectral, touchstone, error, python) + `bin/`; update the dependency table to match `Cargo.toml` (add faer, rustfft, touchstone, pyo3, criterion); remove the stale "Touchstone export not available" limitation (it exists at `src/touchstone.rs`).
  2. Pick the real version (git history shows v0.4.1) and set `Cargo.toml` == `pyproject.toml` == README `mfs = "x.y"` to the same value. (If the maintainer's intent is to stay 0.2.0, use 0.2.0 everywhere — but they MUST match. Default: align to the highest in git history, 0.4.1, and note the assumption in Findings.)
- Verify (S9 partial): `grep` architecture.md contains all 7 module names; `grep` the version in the 3 files is identical. No production code touched → TDD-exempt (justified: docs/version only).
- Commit: `docs: sync architecture.md module list + deps; align version across Cargo/pyproject/README`.

**T1.2 — P4a root-solver benchmark scaffold.**  `category: quick`, skills: []
- Action: add a `criterion` bench (new file under `benches/`, registered in `Cargo.toml`) that times `FaerRootSolver`, `CompanionMatrixRootSolver` (nalgebra), `DurandKernerRootSolver`, `AberthRootSolver` on polynomials of degree {10, 20, 30} with known roots, and asserts each recovers the roots within 1e-6 (so the bench also doubles as a correctness cross-check).
- Verify: `cargo bench --bench <name>` runs and prints per-solver times; the accuracy asserts pass. (This is data for the T2.x decision, not the decision itself.)
- Commit: `bench: add root-solver comparison bench (faer vs nalgebra companion vs DK vs Aberth)`.

**T1.3 — P1 RED test (lock the no-op contract).**  `category: deep`, skills: []
- Action (tests only, NO production change yet): in `tests/multiplexer.rs` add two tests:
  1. `separation_limit_converges_to_independent`: two channels far apart in frequency → assert the multiplexer's recovered channel matrices equal independently-synthesized matrices (max|Δg| < 1e-9).
  2. `adjacent_channels_show_loading`: two adjacent/overlapping channels → assert the recovered channel matrices DIFFER from independent synthesis beyond a tolerance (e.g. max|Δg| > 1e-4).
- Run them. **Record the RED output** (test 2 — the effect test — is expected to FAIL today because the perturbation is a no-op; test 1 may pass trivially). Paste RED into notepad with the exact assertion messages.
- Verify (RED proof): notepad shows `separation_limit...` result and `adjacent_channels_show_loading` FAIL with a "matrices identical / max|Δg| = 0" style message. **Do NOT touch production code in this task.**
- Commit: `test: add failing multiplexer mutual-loading propagation tests (RED)`.

---

## Wave 2 — Core fixes (2 tasks, parallel; both depend only on Wave 1)

**T2.1 — P1 GREEN: make mutual loading actually propagate.**  `category: deep`, skills: []
- Action (smallest change to flip T1.3 RED→GREEN):
  1. In `recover_filter_design` (`src/multiplexer.rs:285-325` and the `duplexer.rs` twin), after rebuilding E/F coefficients, **compute their roots** and write them to `new_polys.generalized.e_s` / `new_polys.generalized.f_s` (verify `P` invariance; re-root `p_s` only if the perturbation math changes it). This is option (i) from the review.
  2. Un-hardcode `use_full_iteration` (`src/duplexer.rs:391`): make it a field on the synthesis config (default `false` to preserve current behavior), so the full-complex path is reachable. Do NOT change its default behavior in this task (separate concern).
  3. Add the S2 boundary unit test: feed `MatrixSynthesisEngine` slightly perturbed `e_s` roots → assert the output matrix changes (max|Δ|>tol).
- Verify (S1, S2 GREEN): `cargo test --test multiplexer` green (both T1.3 tests now pass); `cargo test -p mfs --lib` shows the S2 boundary test PASS. Capture both transcripts as artifacts.
- Commit: `fix: propagate N-channel mutual loading into generalized.e_s/f_s (closes no-op)`.
- **If GREEN requires >~40 lines of production change, STOP** — the test was too coarse or the fix is bigger than expected; split it (per TDD rule) and note in Findings.

**T2.2 — P2: single-source the orchestration.**  `category: deep`, skills: ["refactor"]
- Action:
  1. Introduce one canonical orchestrator in `src/synthesis/orchestration.rs`, e.g. `fn synthesize_to_matrix(spec: &FilterSpec, opts: &SynthOptions) -> Result<SynthesisOutcome>` (and a `+response` variant), that performs the spec→approx→synthesize(→transform) sequence.
  2. Refactor the 4 call sites to delegate to it: `design.rs::FilterDesignBuilder::synthesize` (:84-119), `pipeline/request.rs::run_full_pipeline` (:146-248), `pipeline/execution.rs::run_stage` (:40-143), and `mfs_cli.rs` `build_context_from_request` (:655-678). Collapse the 3× copy-pasted FilterSpec-build into one shared builder helper.
  3. Unify the band-center convention (S4): one definition (geometric mean, matching `bandpass_edges`), used by both the CLI `design` subcommand and `FilterDesign::bandpass_edges`.
- Verify (S3, S4): `grep -rn "FilterSpec::new" src` shows the canonical build in 1 place + thin call sites; all four entry points call the new orchestrator; a test asserts CLI and `bandpass_edges` give identical centers; full `cargo test` green (no behavior change — this is a refactor, so existing tests are the characterization net; if any test was pinning the old arithmetic-mean CLI center, update THAT test to the unified value and justify in Findings).
- Commit: `refactor: single-source synthesis orchestration; unify band-center convention`.

---

## Wave 3 — Dependent on Wave 2 (3 tasks, parallel)

**T3.1 — P3: pipeline `mapping` (implement or reject-with-error).**  `category: deep`, skills: []
- Action: pick (a) implement physical-domain evaluation in `ResponseEvaluationStage` (consume `options.mapping`, map the grid, evaluate in Hz) OR (b) remove the dead `mapping` passthrough and reject an unrecognized/unsupported `mapping` key with a typed `MfsError`. Default to (a) if it is <~60 lines and reuses the canonical orchestrator's response path (from T2.2); otherwise (b). Add a test (S5) for the chosen behavior.
- Verify (S5): `cargo test --test cli_integration` (or the pipeline test) shows the mapping behavior is correct/explicit. Capture transcript.
- Commit: `fix: pipeline mapping — <physical-domain eval | reject-with-error>`.

**T3.2 — F: deprecate legacy root free functions.**  `category: quick`, skills: []
- Action: add `#[deprecated(note = "use mfs::design::FilterDesign / module APIs")]` to the crate-root legacy fns in `lib.rs` (:98-172). Keep them functional and keep the prelude re-exports. Fix any in-tree callers to the new API so the build is warning-clean.
- Verify (S6): `cargo build 2>&1` shows zero in-tree deprecation warnings; `cargo test` green. Capture the build output.
- Commit: `deprecate: legacy crate-root free functions (route to FilterDesign/module APIs)`.

**T3.3 — P4b: de-duplicate the numeric stack (uses T1.2 bench data).**  `category: deep`, skills: []
- Action: using the T1.2 bench results, keep **one** companion-eigen solver in the `AdaptiveRootSolver` production chain and **demote the other to a test-only cross-validation oracle** (a test that runs both on random polynomials and asserts they agree within tolerance). Remove `faer` from `Cargo.toml` if it is the demoted one (then delete `FaerRootSolver` at `complex_poly.rs:616-673` and its Strategy-1 slot at :684). For rustfft: either wire its cepstrum into the now-reachable full-iteration path (T2.1) behind a feature, or remove the dependency + the `spectral.rs` cepstrum half. Fix the 4-way doc inconsistency on fallback order (state it once, in `complex_poly.rs` + README, matching the final code).
- Verify (S7): `grep -rn "faer\|CompanionMatrixRootSolver" src` shows exactly one in the production chain + the other only under `#[cfg(test)]`; the cross-validation oracle test passes; `cargo test` green; README fallback-order text matches code. Capture grep + test transcript.
- Commit: `refactor: single root-solver in production chain; demote the other to test oracle`.

---

## Wave 4 — Cleanup (3 tasks, parallel; depend on Wave 2/3 being stable)

**T4.1 — P7: dead-code removal.**  `category: quick`, skills: []
- Action: delete `response/pole_expansion.rs:294-325` (`evaluate_with_pole_expansion_verified`, dead, different-tolerance duplicate), `spectral.rs:236-339` (dead balanced-companion duplicates + no-op `balance_matrix`), `duplexer.rs:663` (`wilson_spectral_factor`, zero callers). **Before deleting each**, confirm no callers (`grep`) and that existing tests stay green (they are the pinned behavior).
- Verify: `cargo test --all-targets` green after deletions; `cargo clippy --all-targets -- -D warnings` green (dead-code warnings gone). TDD-exempt (justified: pure dead-code deletion with green pinned tests).
- Commit: `chore: remove dead code (pole_expansion verified helper, spectral balanced companions, wilson_spectral_factor)`.

**T4.2 — P8: split oversized `residues.rs`.**  `category: deep`, skills: ["refactor"]
- Action: split `src/synthesis/residues.rs` (1279 LOC) along the engine seams into `residues/admittance.rs` (input parsing / `synthesize_admittance_polynomials`), `residues/expansion.rs` (residue computation), `residues/transversal.rs` (matrix assembly / `build_transversal_from_residues`), with `residues/mod.rs` re-exporting. Pure move — no logic change.
- Verify: `cargo test --all-targets` green (characterization net); `wc -l` shows no file >~600; module public surface unchanged (prelude/legacy still compile).
- Commit: `refactor: split synthesis/residues.rs into admittance/expansion/transversal`.

**T4.3 — S8: close the "silent degradation" policy.**  `category: deep`, skills: []
- Action: (1) `group_delay` on the pole-expansion fast path: either compute it or add an explicit marker (e.g. `ResponseSample.degraded: bool` or a `GroupDelay` enum with `NotComputed`) so it is never a silent 0.0. (2) Ensure the multiplexer full-iteration and pipeline mapping (T3.1) surface their state explicitly (no silent no-op). Add/adjust a test asserting the marker is set when the fast path does not compute group delay.
- Verify (S8): a test asserts the degradation flag/marker; `cargo test` green. Capture transcript.
- Commit: `fix: surface silent degradations (group_delay fast path, multiplexer iteration, pipeline mapping)`.

---

## Wave 5 — Final gate (serial)

**T5.1 — Full regression + CI parity.**  `category: quick`, skills: ["git-master"]
- Action: run the exact CI trio from `.github/workflows/ci.yml`: `cargo test --all-targets`, `cargo clippy --all-targets -- -D warnings`, `cargo doc --no-deps`. Fix anything this wave's changes broke (root cause, not suppression). Re-run the S1–S10 scenario list end-to-end and record PASS/FAIL for each with both artifacts.
- Verify (S10 + full S1–S10): all three commands exit 0; notepad `## Evidence` has every S-row marked PASS with its artifact path. 
- Commit: (only if fixes were needed) `fix: ...` for each; else no commit.

**T5.2 — Sync docs + CHANGELOG to the as-built state.**  `category: writing`, skills: []
- Action: update `docs/architecture.md` to reflect the new orchestration core, the final root-solver chain, the removed dead code, and the deprecation of legacy fns; add a `## [0.x.0]` entry to `CHANGELOG.md` summarizing the fixes (P1 no-op, P2 orchestration, P3 mapping, P4 stack, P7 dead code, deprecations) with file refs.
- Verify (S9): architecture.md matches the final `src/` tree; CHANGELOG entry present; `cargo doc --no-deps` green.
- Commit: `docs: sync architecture.md + CHANGELOG to as-built architecture`.

---

## Parallel task graph (waves & dependencies)

```
Wave 0  [T0.1 baseline]
Wave 1  [T1.1 docs+ver]  [T1.2 bench]  [T1.3 P1-RED]        (parallel)
Wave 2  [T2.1 P1-GREEN]  [T2.2 P2-orchestration]            (parallel; both need Wave 1)
Wave 3  [T3.1 P3-mapping]  [T3.2 deprecate]  [T3.3 P4-stack] (parallel; T3.1/T3.2 need T2.2; T3.3 needs T1.2)
Wave 4  [T4.1 dead-code]  [T4.2 split residues]  [T4.3 silent-degrade] (parallel)
Wave 5  [T5.1 regression] → [T5.2 docs/changelog]           (serial)
```

- **Critical path:** T0.1 → T1.3 → T2.1 (P1 correctness) → T5.1.
- **Biggest parallelism:** Wave 1 (3-way) and Wave 3/4 (3-way each).
- **Do not parallelize** the RED (T1.3) and GREEN (T2.1) of the same scenario — T2.1 must see T1.3's RED first.

## Reviewer gate (mandatory)

Because this touches 3+ files across the crate and is a refactor+correctness
effort, run a high-rigor reviewer after Wave 3 (post-P2/P3) and again at Wave 5
(final). Pass the reviewer: goal + S1–S10 contract + evidence + the diff + the
notepad path. Fix every criterion-cited blocker; re-run only the affected
scenarios; re-submit at most twice.

## Explicit non-goals (out of scope for this plan)

- No new filter families (Butterworth/Elliptic), no HighPass/BandStop mappings,
  no passband-TZ / repeated-pole support — these are feature roadmap, not
  remediation.
- No workspace split, no pyo3 API redesign.
- No change to the *math* of the perturbation beyond making it actually
  propagate (option i); a root-domain rewrite (option iii) is a future follow-up.
