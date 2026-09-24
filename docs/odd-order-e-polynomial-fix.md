# Odd-Order Generalized-Chebyshev Synthesis: E(s) and Admittance Parity Fix

## Purpose

This note records two correctness bugs that together broke **all odd-order**
generalized-Chebyshev synthesis (orders 3, 5, 7, 9, ...), including finite-zero
cases such as the 7p2z B42 winning configuration. For odd orders the
synthesized response carried spurious passband transmission zeros (S21 dipping
to −284…−313 dB in-band) and the return loss was destroyed (S11 ≈ 0 dB
instead of the specified ripple). Even orders appeared to work by accident.

## TL;DR

- **Bug 1** — `build_e_polynomial_stage`
  (`src/approx/generalized_chebyshev_helpers.rs`): the E(w) stage did not
  (a) normalize F(w) and P(w) to monic, (b) apply the `j` factor to the
  P(w)/eps term, and (c) normalize the reconstructed E(s) to monic. The
  `s = jw` transform rotates each coefficient by `j^k`, so after the transform
  the leading coefficient is `j^N`, not 1. For even N that is a real ±1 (roots
  unchanged), for odd N it is ±j (a 90° rotation that corrupts the F/P term
  balance). The missing `j` on the P term is wrong for every order with finite
  zeros.
- **Bug 2** — `synthesize_admittance_polynomials`
  (`src/synthesis/residues.rs`): the `E(-s)*` and `F(-s)*` terms built from
  `alternating_conjugate` were missing the global `(-1)^N` factor. For odd N
  the leading coefficient of the admittance denominator cancels to zero, so the
  poles land off the imaginary axis and the transversal builder rejects them.
- **Verification:** after both fixes, 11 filter configs (orders 3–8 all-pole
  plus 5 finite-zero cases) match the fuqtgui Python engine (`coupling_matrix.py`
  `find_Ps`/`find_Es`/`coupling_matrix` pipeline) to 0.000 dB on S11 and a
  max coupling-matrix deviation of 8.77e-15.
- **Regression test:** the order-4/2-TZ backward-compatibility reference values
  (`src/synthesis/residues.rs`) were re-captured from the corrected output and
  cross-checked against the Python engine.

## Background

The generalized-Chebyshev prototype is built from three polynomials in the
complex-frequency variable `s`, with the standard substitution `s = jw` used to
evaluate the literature recurrences in the helper variable `w`:

- `F(s)` — from the Cameron recurrence (all-pole part).
- `P(s)` — induced by the finite transmission zeros, with a parity phase
  applied when `(order − finite_count)` is even.
- `E(s)` — the denominator polynomial whose roots are the filter poles.

The reference implementation lives in
`fuqtgui/src/domain/coupling_matrix.py`. The Rust approximation layer must
reproduce `F(s)`, `P(s)`, and `E(s)` identically for the downstream admittance
and residue synthesis to be valid.

## The Bugs

### Bug 1 — `build_e_polynomial_stage`

Cameron's formula for the denominator in the helper variable `w` is

```
E(w) = F(w)/eps_r + j·P(w)/eps
```

where `F(w)` and `P(w)` are **monic**. The previous Rust code computed

```rust
let f_w = s_to_w_coefficients(f_s);          // leading coeff = j^N, NOT monic
let p_w = s_to_w_coefficients(p_s);          // leading coeff = j^NZ·phase, NOT monic
let e_w = f_w / eps_r + p_w / eps;           // missing the j on P(w)
...
let e_s = w_to_s_coefficients(e_w_from_roots); // leading coeff = (-j)^N, NOT monic
```

Three deviations from the reference:

1. `F(w)` and `P(w)` were not normalized to monic. After `s = jw`, the leading
   coefficient is `j^N`. For odd N this is `±j`, which rotates the F term by 90°
   relative to the P term and destroys the relative weighting that fixes the
   E(w) roots.
2. The `j` factor on `P(w)/eps` was absent. This is wrong for *every* order
   that has finite zeros (P is non-constant); for all-pole cases it is masked
   because P is a scalar.
3. `E(s)` was not re-normalized to monic after the `w → s` transform, which
   introduces a `(-j)^N` leading coefficient.

The Python reference does all three:

```python
Fw = Polynomial(new_coeffs / new_coeffs[-1])      # normalize to monic
Pw = Polynomial(new_coeffs / new_coeffs[-1])      # normalize to monic
Es3 = Fw / eps_r + 1j * Pw / eps                  # j factor on P
...
Es = new_coeffs / new_coeffs[-1]                  # normalize to monic
```

### Bug 2 — `synthesize_admittance_polynomials`

The lossless admittance denominator is

```
Y_d(s) = E(s) + F(s)/eps_r + (-1)^N·E(-s)* + (-1)^N·F(-s)* / eps_r
```

The Rust helper `alternating_conjugate` implements `P(-s)* = Σ conj(c_k)(-1)^k s^k`,
i.e. it omits the global `(-1)^N` factor. For even N that factor is +1 (no
effect). For odd N it is −1, so the leading coefficient of `Y_d` cancelled to
zero, pushing the poles out of the imaginary axis. The transversal builder then
rejected the result:

```
transversal synthesis expects poles close to the imaginary axis
(pole real parts 2.39e-01 exceed tolerance 1.00e-06 for order 3)
```

## The Fixes

```rust
// build_e_polynomial_stage
let f_w_raw = s_to_w_coefficients(f_s);
let p_w_raw = s_to_w_coefficients(p_s);
let f_w = normalize_to_monic(&f_w_raw)?;   // (1)
let p_w = normalize_to_monic(&p_w_raw)?;   // (1)
let j_over_eps = ComplexCoefficient::new(0.0, 1.0 / eps);
let e_w = f_w.scale(1.0 / eps_r)?
    .add(&p_w.scale(j_over_eps)?)?;        // (2) j factor
...
let e_s_raw = w_to_s_coefficients(&e_w_from_roots);
let e_s = normalize_to_monic(&e_s_raw)?;   // (3)
```

```rust
// synthesize_admittance_polynomials
let parity_sign = complex_from_real(if order % 2 == 0 { 1.0 } else { -1.0 });
let conjugated_e = e_normalized.alternating_conjugate()?.scale(parity_sign)?;
let conjugated_f = f_over_eps_r.alternating_conjugate()?.scale(parity_sign)?;
```

## Impact / Verification

End-to-end S-parameter and coupling-matrix comparison between the corrected
Rust engine and the fuqtgui Python engine (RL = 20 dB, 3.4–3.6 GHz):

| Config | S11 @ center (dB) | S21 @ center (dB) | CM max dev |
|---|---|---|---|
| order 3 all-pole | −47.318 | 0.000 | 0.000 dB |
| order 4 all-pole | −20.014 | −0.044 | 0.000 dB |
| order 5 all-pole | −42.886 | 0.000 | 0.000 dB |
| order 6 all-pole | −20.032 | −0.043 | 0.000 dB |
| order 7 all-pole | −39.971 | 0.000 | 0.000 dB |
| order 8 all-pole | −20.056 | −0.043 | 0.000 dB |
| order 6, zeros [−3,3] | −20.030 | −0.043 | 0.000 dB |
| order 7, zeros [−2.5,2.5] (B42 7p2z) | −40.180 | 0.000 | 0.000 dB |
| order 7, zeros [2.5] | −30.201 | −0.004 | 0.000 dB |
| order 5, zeros [−3,3] | −43.087 | 0.000 | 0.000 dB |
| order 3, zeros [2.5] | −28.801 | −0.006 | 0.000 dB |

Odd-order centers sit at a reflection zero (very low S11); even-order centers
sit at a ripple peak (S11 ≈ −RL). All match the Python reference. The order-4
two-TZ coupling matrix (M[0,1], M[0,2], …) now agrees with Python to 8.77e-15.

## Regression Coverage

- `src/synthesis/residues.rs::order_4_with_2_tzs_backward_compatibility`
  reference values re-captured from the corrected output and cross-validated
  against the fuqtgui Python engine.
- `src/fixtures/mod.rs::cameron_order3_generalized_pipeline_exact_case` — the
  order-3 E(s) anchor was a stale capture of the old un-normalized gauge
  (leading coefficient `j`); updated to the monic gauge, verified against the
  independent fuqtgui Python reference to ~1e-15.
- `tests/synthesis_properties.rs::property_10` — the old coefficient-parity
  assertion (even powers real, odd powers imaginary) was gauge-fragile: the
  admittance denominator is only defined up to a global ±1/±j rotation of
  E(s), so a fixed coefficient parity cannot hold across gauges (it had been
  failing pre-fix with a committed regression seed). Rewritten to the
  gauge-independent physical invariant it was meant to check: all denominator
  poles lie on the imaginary axis — the same invariant the transversal builder
  enforces via `pole_axis_tolerance`.
- `src/design.rs::bandpass_edges_equivalent` — the test compared
  `bandpass_edges(f1, f2)` (geometric center `√(f1·f2)`, which the bandpass
  transform requires for Ω(±1) at the edges) against `bandpass(6.75e9, …)`
  (arithmetic mean). Fixed the test to compare against the geometric center.
- `src/matrix/sections.rs::extract_trisection` — added the missing
  Arrow-topology precondition (the test
  `trisection_extraction_rejects_non_arrow_input` already expected it; all
  legitimate call paths convert to Arrow first).
- `cargo test`: **285 passed, 0 failed, 1 ignored** (full suite: lib unit,
  integration, and property tests).
