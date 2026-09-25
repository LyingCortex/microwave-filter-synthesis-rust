# MFS Architecture Review

**Date:** 2026-09-25
**Scope:** Whole-crate architecture soundness + optimization opportunities.
**Method:** Parallel codebase exploration (module inventory, API-layering, dependency/numeric audit, duplexer/multiplexer correctness trace) cross-validated against `file:line` evidence, then an independent high-level architectural judgment pass.
**Companion docs:** `refactor-roadmap.md`, `current-code-to-target-architecture-migration-plan.md`, `multiplexer-code-review.md`.

This is an **assessment**, not a code change. Every finding below carries the
evidence it was verified against so it can be re-checked or acted on independently.

---

## A. Verdict

**The core architecture is sound. The failure mode is boundary discipline, not
layering.**

The layering (`spec → approx → matrix → transform → response`), a single
algorithmic implementation shared by every entry point, lean `Vec<f64>` data
types with algebra confined to the compute boundary, exemplary feature
isolation, and literature-fixture-anchored testing are exactly right for a
numerical synthesis library and more mature than typical 0.2-era code.

The risk is **not** "the architecture is wrong." It is two classes of
**silent failure** already in the tree:

1. Two advertised mechanisms — **N-channel mutual loading** and the pipeline
   **`mapping`** field — **silently do nothing**.
2. The **orchestration contract is copy-pasted 5×**, so the entry points have
   already begun to diverge in ways that ship wrong results without an error
   (CLI arithmetic-mean band center, `group_delay` silently zero on the fast
   path, a physical-Hz grid misread as normalized ω).

For a numerical tool, "plausible-looking output that is quietly wrong" is the
most dangerous defect class. The good news: all of this is fixable before 1.0
**without a redesign.**

---

## B. What is done well (keep)

| Dimension | Evidence |
|---|---|
| **No core-algorithm duplication** | All 4 entry surfaces (FilterDesign / module APIs / pipeline / legacy fns) converge on the same leaf functions: `approx::generalized_chebyshev_polynomials`, `MatrixSynthesisEngine::synthesize` (`synthesis/engine.rs:69`), `CouplingMatrix::transform_topology` (`matrix/core.rs:167`, `to_folded:332`, `to_arrow:397`), `ResponseSolver` (`response/mod.rs:91`). |
| **Lean data types** | `CouplingMatrix` stores `Vec<f64>` (`matrix/core.rs:36`); heavy deps (nalgebra/faer/rustfft) enter only at the compute boundary (`to_complex_dense`, `matrix/core.rs:182`), never in the types. |
| **Exemplary feature isolation** | Exactly one `#[cfg(feature = "python")]` in `src/lib.rs:64`; the CLI bin is gated at the *target* level via `required-features = ["cli"]` (`Cargo.toml:57-60`). No cfg scattering in core. |
| **Strong test culture** | 10 `proptest` properties, literature fixtures (Cameron 1999/2011, Macchiarella), regression suites, ~2,564 LOC of integration tests. |
| **Clean error model** | Typed `MfsError`, no panics in the synthesis path. |
| **Deliberate boundary refactors already done** | spec/freq separation, matrix/synthesis split, response backend isolation (per the maintainer's own roadmap). |

---

## C. Problems

### Architectural (compound as the project grows)

**P1 — N-channel mutual loading is a silent no-op. (Most severe.)**
The headline feature of `duplexer`/`multiplexer` — "iterative mutual-loading
compensation" — **never reaches the output.**
`recover_filter_design` (`multiplexer.rs:285-325`) rebuilds E/F coefficients
from the perturbed S-roots and writes `new_polys.e/f`, but it first clones the
original (`new_polys = original_polys.clone()`, `multiplexer.rs:311`) so
`new_polys.generalized` is unchanged. `MatrixSynthesisEngine::synthesize`
reads **only** `generalized.e_s/f_s/p_s`
(`synthesis/residues.rs:242-250`); it **never reads `.e`/`.f`**.
Consequence: output coupling matrices are **bit-identical** to independent
per-channel synthesis. `junction_coupling`, `alpha`, `max_iterations`, and
`convergence_tol` have no effect on the returned `channels`. **No regression
test guards this.** ~2,300 LOC of `duplexer`/`multiplexer` may rest on a false
premise.
- The full-complex junction path is additionally gated off by a hardcoded
  `let use_full_iteration = false; // TODO` (`duplexer.rs:391`).
- The `fully-reflective channel` note in `CHANGELOG.md:78-83` [Unreleased] is
  now **stale**: that symptom's root cause was a *core* odd-order residue bug,
  fixed by `residue_real_component` (`synthesis/residues.rs:72-78`); the guard
  test `multiplexer_recovered_channel_is_passband` is active and passing
  (`peak_db > -0.5`).

**P2 — Orchestration duplicated 4.5×. (Largest ongoing velocity tax.)**
The same `spec → approx → synthesize → transform` sequence is re-implemented in
5 places: `design.rs:84-119`, `synthesis/orchestration.rs:47-70`,
`pipeline/request.rs:146-248`, `pipeline/execution.rs:40-143`, and
`mfs_cli.rs:655-678` (`build_context_from_request`, a duplicate of
`request.rs:148-170`). The `FilterSpec` construction block is copy-pasted 3×
(`design.rs:101-107`, `request.rs:153-159`, `mfs_cli.rs:661-667`). Divergence
is already shipping: the CLI `design` subcommand uses the **arithmetic mean**
`(f1+f2)/2` for the band center (`mfs_cli.rs:308`) while
`FilterDesign::bandpass_edges` uses the **geometric mean** `√(f1·f2)`
(`design.rs:201`).

**P3 — Pipeline `mapping` field is dead + a silent correctness trap.**
`mapping` (`pipeline/context.rs:21`, `request.rs:90`) is accepted in JSON and
documented for "physical-frequency evaluation," but **nothing consumes
`options.mapping`** — the response stage always calls
`ResponseSolver::evaluate_normalized` on the raw grid
(`pipeline/stages.rs:144`, grid built raw at `request.rs:235`). So a physical-Hz
grid in the doc example is misread as normalized ω. The pipeline also silently
lacks capabilities the high-level API has: pole-expansion fast path, lossy Q,
band-pass scaling, touchstone export.

**P4 — Redundant numeric stack (faer + nalgebra + dead rustfft).**
faer and nalgebra **both** implement companion-matrix-eigenvalue root-finding
for the same polynomial. faer's `FaerRootSolver` is Strategy 1 of
`AdaptiveRootSolver` (`approx/complex_poly.rs:684`); the robust nalgebra
`CompanionMatrixRootSolver` (`complex_poly.rs:444-610`) sits **outside** the
adaptive chain. rustfft's only production use (the `spectral.rs` cepstrum) is
behind the hardcoded `use_full_iteration = false` flag (`duplexer.rs:391`) →
dead in the shipped path. The root-finding fallback order has **four
mutually-contradictory** descriptions: `complex_poly.rs:49-55` (nalgebra→DK→Aberth),
`complex_poly.rs:208-211` (DK→Aberth→companion), the actual code
`complex_poly.rs:679-699` (Faer→DK→Aberth), and the README
("Durand-Kerner → Aberth → Companion Matrix").

### Housekeeping (clean up, not architecture)

- **P5 — Documentation drift (severe).** `docs/architecture.md` omits 7 root
  modules (`duplexer`, `multiplexer`, `optimize`, `spectral`, `touchstone`,
  `error`, `python`); its dependency table is stale (missing faer/rustfft/
  touchstone/pyo3/criterion); it still lists "Touchstone export not available"
  although `src/touchstone.rs` (198 LOC) exists and is exported.
- **P6 — Version drift.** `Cargo.toml` 0.2.0 vs `pyproject.toml` 0.1.0 vs README
  `mfs = "0.1"`; git history has v0.3.0/v0.4.0/v0.4.1 commits but `Cargo.toml`
  was never bumped.
- **P7 — Dead code.** `response/pole_expansion.rs:294-325` (a second,
  different-tolerance fallback dispatcher, `#[allow(dead_code)]`);
  `spectral.rs:236-339` (duplicate companion construction + a no-op
  `balance_matrix`, so "balanced" is false advertising); `duplexer.rs:663`
  (`wilson_spectral_factor`, zero callers).
- **P8 — Oversized files.** `synthesis/residues.rs` 1279, `duplexer.rs` 1176,
  `design.rs` 1038, `approx/complex_poly.rs` 913, `bin/mfs_cli.rs` 1107.
- **P9 — Minor inconsistencies.** CLI arithmetic-mean vs geometric-mean band
  center (see P2); `group_delay` silently `0.0` on the pole-expansion fast path
  (`response/pole_expansion.rs:249,277`) while the LU path computes it.

---

## D. Prioritized optimization (by velocity-impact × risk-if-left)

| # | Problem | Effort | Minimal fix direction |
|---|---|---|---|
| **1** | P1 no-op mutual loading | Med (1–2 d) | See §E. Correctness defect; the sooner it is fixed the cheaper the retrofit. |
| **2** | P2 orchestration duplication | Med (1–2 d) | One canonical `run(spec, options) -> Result<DesignOutput>` in `synthesis/orchestration.rs`; make FilterDesign / pipeline / CLI thin adapters over it; collapse the 3× FilterSpec build into one builder. Prerequisite for deprecating the legacy fns. |
| **3** | P3 dead `mapping` + pipeline capability gap | Short (1–4 h, after #2) | Either implement physical-domain evaluation, or **remove the field and reject an unrecognized `mapping` key with a typed error** — never silently ignore. Capability parity comes free once #2 lands (the pipeline rides the canonical path). |
| **4** | P4 redundant numeric stack | Short (decision) / Med (if re-enabling cepstrum) | Benchmark both root solvers with criterion; keep one in the production chain, demote the other to a test-only oracle (cross-validate on random polynomials). Wire rustfft's cepstrum into the N-channel full iteration (ties to #1) or remove the dependency. Don't leave "adaptive" meaning "we haven't decided." |
| **5** | P5 + P6 docs & version drift | Short (<1 d, parallelizable) | Regenerate the `architecture.md` module list from code; state the fallback order exactly once; add a CI check `Cargo.toml version == pyproject == README`; a release checklist that bumps all three. |

**Folded in, not ranked:** while touching `residues.rs` for #1, split the
1,279-line file along the engine seams (input parsing / residue computation /
matrix assembly). The band-center convention (P2/P9) is unified under #2 (one
builder, one definition). `group_delay=0` is closed by the policy in §F.1.

---

## E. Fixing the no-op perturbation (P1)

**Choose (i): populate `new_polys.generalized.e_s/f_s` in `recover_filter_design`.**

- The engine's contract is "give me the **root sets** `e_s/f_s/p_s`; I return
  a matrix." The multiplexer is a preprocessor whose job is to produce
  *exactly that input*.
- **(ii)** — making the engine accept caller-supplied `.e`/`.f` coefficients —
  widens the engine boundary to a second representation, the opposite of the
  maintainer's demonstrated refactoring direction.
- **(iii)** — staying entirely in the root domain and skipping the coefficient
  round-trip — is the ideal end state (coefficient↔root conversion at high
  degree loses precision) but is a larger rewrite of the perturbation math; do
  it only after (i) proves the loop works end-to-end and tests pin the behavior.

**Concrete steps:**
1. In `recover_filter_design`, after rebuilding the E/F coefficients, compute
   their roots and write them to `new_polys.generalized.e_s/f_s`. Verify the
   perturbation math leaves `P` invariant; if it does not, re-root `p_s` too.
2. **Un-hardcode** `use_full_iteration` (`duplexer.rs:391`) — make it a spec
   option. The full path has evidently never run with effect, so expect
   debugging there.

**Regression tests (in priority order):**
1. **Separation-limit test** (strongest, no literature needed): two channels far
   apart in frequency → the iterative result must converge to the independent
   per-channel matrices (max |Δgᵢⱼ| < 1e-9). This verifies the perturbation
   *does something* and *vanishes where physics says it should* — it would have
   caught the no-op.
2. **Effect test:** adjacent/overlapping channels → output matrices differ from
   independent synthesis beyond tolerance, and both channels' reconstructed
   S-parameters show loading (passband edge shift).
3. **Fixture:** a published multi-channel example with tabulated coupling values
   (Cameron multiplexer / Macchiarella papers), an exact-value anchor in the
   style of the existing 1999/2011 fixtures.
4. **Boundary unit test:** feed the engine slightly perturbed `e_s` roots →
   assert the matrix changes. This tests engine sensitivity directly, so a
   future regression is caught at the seam, not 2,300 lines upstream.

---

## F. The four API surfaces — too many; consolidate to three (do not retire)

Target shape:
1. **Module API** = low-level contract, unchanged.
2. **`FilterDesign`** = the single high-level facade.
3. **Pipeline + CLI** = **one** shared execution path, a thin adapter over the
   canonical orchestration from D#2 (the CLI is just args→JSON→pipeline).
4. **Legacy root free functions** = `#[deprecated(note)]` pointing at the
   module APIs now; remove at 0.5/1.0.

Rationale: the pipeline earns its keep (JSON-driven automation/scripting; the
likely home for the Python story), so retiring it punishes users with
hand-rolled Rust loops — but only if it stops *re-implementing* and starts
*adapting*. The legacy functions earn nothing and are already drifting
(mean-center bug). **Litmus test:** each surface must add a capability, never
re-implement one.

---

## G. Cross-cutting items an experienced reviewer would flag

1. **The real defect is "silent degradation" as a pattern.** P1 (no-op
   loading), P3 (`mapping` ignored), P9 (`group_delay=0`, hardcoded flag) are
   the *same* bug: the library has a strong typed-error philosophy in the core
   but applies it inconsistently at exactly the new-subsystem boundaries. One
   policy closes all of them: *any capability the caller requested that cannot
   be honored must either error or be explicitly flagged in the output* (e.g. a
   `degraded: Vec<...>` field or an enum variant), never approximated away.
2. **Numeric robustness gaps.** No ill-conditioning detection —
   companion-matrix root-finding degrades at N ≳ 10–12 and users will hit it;
   add a root-residual check (‖p(r)‖) surfaced as a typed warning. Tolerances
   are internal; expose them. Pin eigenvalue-solver threading in CI so
   literature fixtures compare deterministically across machines.
3. **No benchmark suite.** For a numerical library, criterion benches on
   root-finding + synthesis are what settle D#4 with evidence and guard against
   silent perf regressions.
4. **Verify `duplexer` and `multiplexer` share one N-channel machine.** The
   zero-caller `wilson_spectral_factor` plus `spectral.rs`'s duplicate companion
   construction smells like the duplexer is a frozen early N=2 implementation
   and the multiplexer a later generalization. If so, the duplexer should be an
   N=2 call *into* the multiplexer, and the `spectral.rs` duplication disappears.
5. **`Vec<f64>` matrix.** The choice is right, but flat Vec + index arithmetic
   through three topology transforms is a bug factory. Ensure the transforms are
   pure permutation/reindexing (not numeric), add a dimension invariant at
   construction + `debug_assert` in the index helpers, and a property test that
   the S-response is invariant under transversal ↔ folded ↔ arrow round-trips.
6. **Keep the single crate.** 18.6k LOC with one `#[cfg(feature)]` in `lib.rs`
   is working; a workspace buys nothing now. Revisit only if the Python wheel
   and the library diverge in release cadence — and even then, a workspace with
   path deps, not a split.
7. (Minor) Check whether `optimize` (Levenberg–Marquardt) differentiates through
   `ResponseSolver` numerically — if so, that's a fragility/perf risk worth an
   analytic-gradient option later. Lowest priority here.

---

## H. Suggested sequence

1. Fix P1 (no-op loading) with the §E tests — restores a core promise.
2. Land the canonical orchestration (D#2) — removes the #1 velocity tax and is
   the prerequisite for everything below.
3. P3 (mapping: implement or reject-with-error) + deprecate legacy fns — comes
   nearly free after step 2.
4. P4 (numeric stack: benchmark, pick one, demote the other; wire or drop
   rustfft) — evidence-based, small.
5. P5 + P6 (docs regenerated from code; CI version check; release checklist).

After 1–5, re-run the full scenario/fixture suite and the `cargo clippy
--all-targets -- -D warnings` + `cargo doc` CI jobs before cutting the next
version.
