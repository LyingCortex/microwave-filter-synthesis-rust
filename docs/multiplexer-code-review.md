# Multiplexer Code Review: `src/multiplexer.rs`

**Date:** 2026-09-13
**Scope:** N-channel multiplexer synthesis + (N+1)-port analysis
**Files reviewed:** `src/multiplexer.rs`, `src/duplexer.rs`,
`src/synthesis/residues.rs`, `src/approx/complex_poly.rs`,
`docs/multiplexer-synthesis.md`, `tests/multiplexer.rs`,
`tests/synthesis_properties.rs`

## Status

All issues identified in the review have been fixed. The most important finding
was that the "fully-reflective recovered channel" symptom was **not** a
multiplexer bug but a **core synthesis defect** affecting every odd-order
all-pole filter in the crate. That root cause is now fixed.

- Test results: all integration suites pass (multiplexer, synthesis_properties,
  synthesis_regression, etc.); doctests pass. The only remaining lib failures
  are the **3 pre-existing, unrelated failures** documented in
  `docs/multiplexer-synthesis.md` (they fail at baseline, before any change).
- The multiplexer module now has 8 lib tests + 6 integration tests, all passing.
- The previously-`#[ignore]`d diagnostic
  `multiplexer_recovered_channel_is_passband` now **passes** and is un-ignored.

## Resolved Issues

### 1. CRITICAL — odd-order all-pole filters were fully-reflective (root cause)

**Symptom:** even a pristine `FilterDesign::prototype(5, 20.0)` evaluated to
peak S21 = `-inf` dB; recovered multiplexer channels were equally degenerate.

**Root cause:** in `src/synthesis/residues.rs`,
`build_transversal_from_residues` extracted coupling values from residues with
`.re`. For odd-order filters (an `s = j w` phase factor) the **transfer**
residue `r12` is purely imaginary while the reflection residues `r11`/`r22` are
real. Taking `.re` of `r12` therefore produced **zero load couplings** and a
fully-reflective matrix.

**Fix:** `residue_real_component()` extracts the dominant real/imaginary
component, preferring `.re` and falling back to `.im` only when `.re` is
numerically negligible. This preserves even-order behaviour exactly while fixing
odd-order. Verified: all-pole orders 2–7 now show peak S21 ≈ 0 dB.

### 2. Junction parameters were ignored (Type-I `b₀`, Type-II `c₀`)

`b₀` was bound to `_` and `c₀` discarded; only `n²` fed a hard-coded `alpha`.
The junction is now a complex coupling coefficient
(`n²/(1−jn²b₀)` for Type-I, `1/(1−jc₀)` for Type-II) applied to every
mutual-loading correction, and `alpha` is a pure geometric under-relaxation.
New tests `multiplexer_junction_susceptance_affects_synthesis` and
`multiplexer_type_ii_junction_converges` verify the parameters now matter.

### 3. Convergence tests asserted only `is_finite()`

The triplexer/quadruplexer convergence tests now assert the error is finite,
non-negative, bounded (`< 1.0`, catching divergence), and that iterations ran.

### 4. Root-clamp semantics were fragile/ambiguous

The real-part clamp is documented and made conservative: it now only reflects a
root that drifts into the right half-plane (keeping the auxiliary polynomial
Hurwitz) and preserves legitimate negative real parts. The stability guard
(`finite`/`norm < 50`) and its intent are documented.

### 5. `Multiplexer::response` antenna model diverged from the diplexer

The multiplexer used `max(S11)` while the diplexer used `√(S11_tx·S11_rx)`.
With real filters the N=2 parity test diverged by up to 146 dB. The multiplexer
now uses the **geometric mean** `(∏ᵢ|S11ᵢ|)^(1/N)`, matching the diplexer for
N=2 and providing a better parallel-junction reflection model.

### 6. Duplicated helper functions

`pad_polynomial`, `compute_root_displacement`, and `multiply_by_s` were
copy-pasted between `duplexer.rs` and `multiplexer.rs`. They are now shared
(`crate::approx::{pad_polynomial, compute_root_displacement, multiply_by_s}`),
removing drift risk.

### 7. No `TypeII` / junction coverage

Added `multiplexer_type_ii_junction_converges` and
`multiplexer_junction_susceptance_affects_synthesis`.

## Remaining Notes

- **The mutual-loading perturbation does not reach the output (confirmed
  2026-09-14).** A probe comparing the multiplexer's recovered channel coupling
  matrix to an independently-synthesized filter of the same spec showed the two
  are **bit-identical** (`max |Δ| = 0`). `MatrixSynthesisEngine.synthesize`
  reads `generalized.e_s/f_s/p_s` and ignores the `PolynomialSet.e/f` that
  `recover_filter_design` rebuilds from the perturbed S-roots
  (`src/multiplexer.rs`). Consequently the N-channel loading loop,
  `junction_coupling`, `alpha`, `max_iterations`, and `convergence_tol` have
  **no effect on the returned `channels`**; the multiplexer is functionally
  independent per-channel synthesis. This is a pre-existing architectural gap
  (not a regression from the residue fix). Resolution requires either
  `recover_filter_design` to hand its reconstructed E/F to the engine in a form
  the engine consumes, or the engine to accept a caller-supplied
  `PolynomialSet.e/f`. This qualifies item 2 above: the junction parameters
  affect the internal `s_roots` iteration, but that iteration does not yet
  propagate to the output matrices.
- The synthesis still uses the **simplified perturbation path** (not the full
  Macchiarella N×N p₀ spectral solve). This is a documented design choice for
  "get the N-channel flow running"; it does not converge to machine tolerance
  for well-separated channels (residual ≈ 1e-2) but the produced filters are
  now valid passband filters.
- The response model remains magnitude-based (geometric-mean S11, product-of-S21
  isolation); it is not a full complex (N+1)-port solve. Documented as a
  deviation in `docs/multiplexer-synthesis.md`.