# Folded-Topology Fidelity: flip_sign Fix and Real-Filter Verification

## Purpose

This note records a real correctness bug found in the folded-topology
reduction, its one-line root cause and fix, and the regression coverage added
to keep it from recurring. It also answers a practical question: does this
crate's synthesized + folded output reproduce the coupling matrices published
in Richard J. Cameron's *Microwave Filters for Communication Systems*?

## TL;DR

- **The bug:** `CouplingMatrix::flip_sign` (in `src/matrix/core.rs`) negated the
  diagonal element a *third* time, so any node whose main coupling was negative
  after the raw reduction got a wrong self-coupling sign. For even-order
  asymmetric filters this silently **dropped the transmission zeros** from the
  folded matrix's response.
- **The fix:** remove that redundant third negation. A proper similarity flip
  `R*M*Rᵀ` (with `R[i][i] = -1`) leaves the diagonal unchanged, which is exactly
  what the fuqtgui implementation and the Cameron book produce.
- **Verification:** after the fix, all tested orders/parities preserve the
  folded response to `≤ 1e-12`, and the folded matrices match the book's
  published matrices entry-for-entry.
- **Regression test:** `tests/folded_response_invariance.rs` synthesizes real
  filters and checks folded-vs-transversal response invariance. It fails on the
  four even-order asymmetric cases with the buggy code and passes with the fix.

## Background

The prior folding test (`tests/transform_invariance.rs`) only exercised a toy,
hand-built, all-positive transversal matrix with no self-couplings. That matrix
is so benign that folding it trivially preserves the response, so the test
could not detect the diagonal-sign bug. Real synthesized filters — especially
even-order ones with an asymmetric distribution of transmission zeros — expose
it immediately.

## The Bug

In `src/matrix/core.rs`, `flip_sign` implemented the sign flip that makes main
chain couplings positive after the Givens reduction:

```rust
// Negate row `diagonal_index`
for col in 0..side { result.data[diagonal_index * side + col] = -...; }
// Negate column `diagonal_index`
for row in 0..side { result.data[row * side + diagonal_index] = -...; }
// The diagonal element was negated twice, restore it
result.data[diagonal_index * side + diagonal_index] =
    -result.data[diagonal_index * side + diagonal_index];
```

After the row and column negations the diagonal element has been negated twice
and is back to its original value. The last line then negated it a *third* time,
turning it into the wrong sign. The comment claimed to "restore" it, but the
operation actually changed it.

The correct semantics of a similarity sign flip are: negate row `i`, negate
column `i`, and leave the diagonal unchanged (it is affected twice and cancels
out). The fuqtgui reference implementation does exactly this:

```python
R = np.eye(N); R[i+1, i+1] = -1
CM = R @ CM @ R.T   # diagonal unchanged
```

## Impact

The sign error only surfaces when the raw reduction leaves a main coupling
negative, which happens for **even-order asymmetric** filters (the cases where
the Cameron folding must be performed in the mirrored/load-side orientation).
Observed symptom: the folded matrix's response **lost the transmission zeros** —
the stopband flattened to ~ -27 dB with no deep nulls instead of showing the
expected zeros.

Measured response deviation (transversal vs folded, linear `|S21|`), with the
buggy code:

| Case | Order | max \|S21\| deviation |
|---|---|---|
| 4-2 (Fig 8.23b) | 4 | 0.139 |
| 6-3 (Fig 8.32) | 6 | 0.978 |
| 8-3 (Fig 9.22a) | 8 | 0.447 |
| 10-4 | 10 | 0.988 |
| 7-1-2 (Fig 8.13) | 7 | ~0 (odd, unaffected) |
| 5/7/3 odd, all-pole, symmetric-even | — | ~0 |

Same input matrix, correct implementation: the fuqtgui folded result matched the
book to `4.5e-5`; the buggy mfs result differed by `0.31` in the diagonal signs
(`M22` `+0.1439` vs `-0.1439`, `M44` `-0.1549` vs `+0.1549`).

## The Fix

Remove the redundant third negation so `flip_sign` performs a proper similarity
transformation:

```rust
// Negate row `diagonal_index`
for col in 0..side { result.data[diagonal_index * side + col] = -...; }
// Negate column `diagonal_index`
for row in 0..side { result.data[row * side + diagonal_index] = -...; }
// Diagonal is unchanged after the two negations; do NOT negate it again.
result
```

After the fix every tested case preserves the folded response to `≤ 1e-12`, and
the 4-2 / 6-3 / 8-3 folded matrices match the book entry-for-entry.

## Regression Coverage

`tests/folded_response_invariance.rs` synthesizes real generalized-Chebyshev
filters across orders and parities, folds them, and asserts the folded response
matches the transversal on a fine grid that reaches the transmission zeros.
Cases:

- 4-2 asymmetric (Fig 8.23b) — the exact case that used to fail
- 6-3 asymmetric (Fig 8.32)
- 8-3 asymmetric (Fig 9.22a)
- 6-2 symmetric (Fig 9.21)
- 7-1-2 asymmetric (Fig 8.13, odd order)
- 10-4 asymmetric (high-order stress)
- 4 all-pole
- 5 all-pole

Run with:

```
cargo test --test folded_response_invariance
```

With the buggy code this reports **4 failures** (the even-order asymmetric
cases); with the fix, **8 pass**.

## Related Findings

1. **Fully canonical synthesis was also rejected (now fixed).** Synthesizing a
   fully canonical filter (`nfz == order`, e.g. the 4-4 case of Fig 8.20/8.21)
   previously failed with `dimension mismatch: expected <order>, got <order+1>`.
   The guard in `src/approx/polynomial.rs` required `P(s)` to have *at most*
   `order` coefficients (degree `≤ order-1`), which excludes the fully canonical
   case where `P(s)` has degree `order` (nfz + 1 = `order + 1` coefficients).
   The bound was corrected to `≤ order + 1` (i.e. reject only `nfz > order`).
   The approximation layer already handled the fully canonical `eps_r`
   (`find_eps` in `src/approx/generalized_chebyshev_helpers.rs` uses
   `eps/sqrt(eps²-1)` when `finite_zero_count == order`). After the fix, the
   4-4 fully canonical case synthesizes and folds correctly, reproducing the
   book's Fig 8.21 matrix including the source-load coupling `MSL = 0.0151`.
   Regression coverage added in `tests/folded_response_invariance.rs`
   (`folded_preserves_response_order4_fully_canonical`,
   `folded_preserves_response_order6_fully_canonical`).

2. **Book "folded" layouts that are reconfigured topologies** (e.g. cul-de-sac,
   extended box) are not produced by the plain folded reduction. The reduction
   yields the standard folded (inline) arrangement; converting to those layouts
   requires additional similarity reconfigurations. The coupling *values* are
   correct; only the arrangement differs (a plane-rotation equivalence class).