//! Shared channel-recovery helper for the Macchiarella-Tamiazzo iterative
//! duplexer/multiplexer synthesis.
//!
//! Both iterative syntheses converge a shared diplexer S polynomial and then
//! must rebuild each channel's `FilterDesign` from its updated S roots:
//! E = S + D, F = S − D (D is supplied in S's normalization by every caller),
//! then re-run the matrix synthesis engine.

use crate::approx::{pad_polynomial, ComplexPolynomial, PolynomialSet};
use crate::design::FilterDesign;
use crate::error::Result;
use crate::synthesis::MatrixSynthesisEngine;
use num_complex::Complex64;

/// Rebuilds a channel's `FilterDesign` after the Macchiarella-Tamiazzo
/// iteration has updated the channel's S-polynomial roots.
///
/// `s_roots` are the converged S-polynomial roots (left half-plane);
/// `d_original` is the channel's D polynomial, already expressed in S's
/// normalization; `original_polys` is the channel's initial `PolynomialSet`.
/// The perturbed E/F are written into both the plain `e`/`f` fields and
/// `generalized.e_s`/`f_s` — the synthesis engine reads `generalized.e_s/f_s`
/// (not `.e/.f`), so writing only the plain fields silently drops the loading
/// correction (see architecture-review P1).
pub(crate) fn recover_channel_design(
    s_roots: &[Complex64],
    d_original: &ComplexPolynomial,
    original_polys: &PolynomialSet,
    order: usize,
    return_loss_db: f64,
    center_hz: f64,
    bandwidth_hz: f64,
) -> Result<FilterDesign> {
    // Build S polynomial from converged roots (always monic, degree n).
    let s_poly = ComplexPolynomial::from_complex_roots(s_roots)?;

    // E = S + D,  F = S − D, with `d_original` used as-is.
    //
    // Every call site now supplies D already in S's normalization, so no
    // rescaling is applied:
    //   - full path: D = S − F⁰ with the channel's *fixed* reflection zeros
    //     F⁰; degree n − 1 with a non-monic leading coefficient.
    //   - simplified path: D = (Eᵢ − Fᵢ)/2 from the initial synthesis.
    //     `ComplexPolynomial::sub` trims the cancelled degree-n term, so the
    //     stored D is genuinely degree n − 1 with an O(1) leading coefficient.
    //     In the zero-perturbation limit `S + D = Eᵢ` and `S − D = Fᵢ` hold
    //     *exactly* with no rescale; the previous `rescale_d` branch rescaled
    //     D to monic (its `d_leading.norm() > 1e-15` guard always fired because
    //     the trimmed leading coefficient is O(1), not ~0), which changed F's
    //     roots and systematically destroyed the passband return loss
    //     (review P1-1). Using D as-is restores the exact reconstruction.
    let d_padded = pad_polynomial(d_original, s_poly.coefficients.len())?;

    let e_poly = s_poly.add(&d_padded)?;
    let f_poly = s_poly.sub(&d_padded)?;

    let mut new_polys = original_polys.clone();
    if let Some(g) = new_polys.generalized.as_mut() {
        g.e_s = Some(e_poly.clone());
        g.f_s = f_poly.clone();
    }
    new_polys.e = e_poly;
    new_polys.f = f_poly;

    let matrix = MatrixSynthesisEngine.synthesize(&new_polys)?;

    FilterDesign::from_polynomials_and_matrix(
        order,
        return_loss_db,
        center_hz,
        bandwidth_hz,
        new_polys,
        matrix,
    )
}
