//! Shared channel-recovery helper for the Macchiarella-Tamiazzo iterative
//! duplexer/multiplexer synthesis.
//!
//! Both iterative syntheses converge a shared diplexer S polynomial and then
//! must rebuild each channel's `FilterDesign` from its updated S roots:
//! E = S + D, F = S − D (D rescaled to S's normalization on the simplified
//! path, used as-is on the full path — see `rescale_d`), then re-run the
//! matrix synthesis engine.

use crate::approx::{pad_polynomial, ComplexPolynomial, PolynomialSet};
use crate::design::FilterDesign;
use crate::error::Result;
use crate::synthesis::MatrixSynthesisEngine;
use num_complex::Complex64;

/// Rebuilds a channel's `FilterDesign` after the Macchiarella-Tamiazzo
/// iteration has updated the channel's S-polynomial roots.
///
/// `s_roots` are the converged S-polynomial roots (left half-plane);
/// `d_original` is the D polynomial from the channel's initial synthesis;
/// `original_polys` is the channel's initial `PolynomialSet`. The perturbed
/// E/F are written into both the plain `e`/`f` fields and
/// `generalized.e_s`/`f_s` — the synthesis engine reads `generalized.e_s/f_s`
/// (not `.e/.f`), so writing only the plain fields silently drops the loading
/// correction (see architecture-review P1).
#[allow(clippy::too_many_arguments)]
pub(crate) fn recover_channel_design(
    s_roots: &[Complex64],
    d_original: &ComplexPolynomial,
    original_polys: &PolynomialSet,
    order: usize,
    return_loss_db: f64,
    center_hz: f64,
    bandwidth_hz: f64,
    rescale_d: bool,
) -> Result<FilterDesign> {
    // 8 args: the review-mandated `rescale_d` flag (P1-1) pushed an
    // already-at-limit function over the clippy 7-arg ceiling. The four
    // scalars (order/rl/center/bandwidth) pass straight through to
    // `FilterDesign::from_polynomials_and_matrix`, so grouping them would not
    // reduce the real arity at either end.
    // Build S polynomial from converged roots
    let s_poly = ComplexPolynomial::from_complex_roots(s_roots)?;

    // E = S + D,  F = S - D.
    // Simplified path: `d_original = (E − F)/2` from the initial synthesis,
    // whose leading coefficient cancels to ~0; when it is numerically
    // non-negligible it is rescaled to S's normalization (the historical
    // behavior). Full path: `d_original = S − F⁰` with the channel's *fixed*
    // reflection zeros F⁰; it is already expressed in S's monic
    // normalization (degree n − 1, non-monic leading coefficient), and
    // rescaling it corrupts the exact reconstruction F = S − D (the WIP
    // branch's rescale_d flag, re-added after commit 818d38c dropped it on
    // the false premise that S − F⁰ is monic).
    let d_scaled = if rescale_d {
        let s_leading = s_poly.leading_coefficient();
        let d_leading = d_original.leading_coefficient();
        if d_leading.norm() > 1e-15 {
            let scale_factor = s_leading / d_leading;
            d_original.scale(scale_factor)?
        } else {
            d_original.clone()
        }
    } else {
        d_original.clone()
    };

    let d_padded = pad_polynomial(&d_scaled, s_poly.coefficients.len())?;

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
