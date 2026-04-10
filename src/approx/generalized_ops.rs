//! Domain-specific helper operations for generalized-Chebyshev synthesis.
//!
//! These helpers stay narrower than `complex_poly`: they encode the `w <-> s`
//! transforms, recurrence-side real-polynomial arithmetic, and a few
//! generalized filter conventions used by the Cameron-style pipeline.

use crate::error::{MfsError, Result};

use super::{ComplexCoefficient, ComplexPolynomial};

const COMPLEX_ZERO: ComplexCoefficient = ComplexCoefficient::new(0.0, 0.0);

pub(crate) fn complex_from_real(value: f64) -> ComplexCoefficient {
    ComplexCoefficient::new(value, 0.0)
}

pub(crate) fn mul_i_pow(value: ComplexCoefficient, power: usize) -> ComplexCoefficient {
    match power % 4 {
        0 => value,
        1 => ComplexCoefficient::new(-value.im, value.re),
        2 => ComplexCoefficient::new(-value.re, -value.im),
        _ => ComplexCoefficient::new(value.im, -value.re),
    }
}

/// Convert a recurrence polynomial written in descending powers of `w`
/// into the crate's standard ascending-coefficient `s`-domain storage.
///
/// Cameron's recurrence is naturally written as `U(w) = sum a_k w^(n-k)`.
/// Public polynomial artifacts in this crate are stored as
/// `F(s) = sum c_k s^k`, so this helper applies `s = j w` term by term while
/// also flipping from descending to ascending coefficient order.
pub(crate) fn descending_w_to_ascending_s(descending: &[f64]) -> Result<ComplexPolynomial> {
    let order = descending.len().saturating_sub(1);
    let mut ascending = vec![COMPLEX_ZERO; descending.len()];

    for (index, coefficient) in descending.iter().copied().enumerate() {
        let power = order - index;
        ascending[power] = mul_i_pow(complex_from_real(coefficient), index);
    }

    ComplexPolynomial::new(ascending)
}

/// Rewrite an ascending `s`-domain polynomial into ascending `w`-domain
/// coefficients using `s = j w`.
///
/// If `P(s) = sum c_k s^k`, then `P(jw) = sum c_k j^k w^k`, so each
/// coefficient is rotated by `j^k` while keeping the ascending storage order.
pub(crate) fn s_to_w_coefficients(coefficients: &[ComplexCoefficient]) -> Vec<ComplexCoefficient> {
    coefficients
        .iter()
        .copied()
        .enumerate()
        .map(|(index, coefficient)| mul_i_pow(coefficient, index))
        .collect()
}

/// Rewrite an ascending `w`-domain polynomial back into ascending `s`-domain
/// coefficients.
///
/// This is the inverse of [`s_to_w_coefficients`]: if
/// `Q(w) = sum d_k w^k`, then `Q(-j s) = sum d_k (-j)^k s^k`.
pub(crate) fn w_to_s_coefficients(coefficients: &[ComplexCoefficient]) -> Vec<ComplexCoefficient> {
    coefficients
        .iter()
        .copied()
        .enumerate()
        .map(|(index, coefficient)| mul_i_pow(coefficient, (4 - index % 4) % 4))
        .collect()
}

/// Rotate a root from the helper `w` plane back into the public `s` plane.
///
/// Because `s = j w`, each `w`-plane root maps to `s_root = j * w_root`.
pub(crate) fn rotate_w_root_to_s(root: ComplexCoefficient) -> ComplexCoefficient {
    ComplexCoefficient::new(-root.im, root.re)
}

pub(crate) fn reflect_to_upper_half_plane(root: ComplexCoefficient) -> ComplexCoefficient {
    if root.im >= 0.0 {
        root
    } else {
        ComplexCoefficient::new(root.re, -root.im)
    }
}

pub(crate) fn convolve_descending(lhs: &[f64], rhs: &[f64]) -> Vec<f64> {
    let mut output = vec![0.0; lhs.len() + rhs.len() - 1];
    for (left_index, left) in lhs.iter().copied().enumerate() {
        for (right_index, right) in rhs.iter().copied().enumerate() {
            output[left_index + right_index] += left * right;
        }
    }
    output
}

pub(crate) fn add_descending(lhs: &[f64], rhs: &[f64]) -> Vec<f64> {
    let target_len = lhs.len().max(rhs.len());
    let mut output = vec![0.0; target_len];

    for (index, value) in lhs.iter().copied().enumerate() {
        output[target_len - lhs.len() + index] += value;
    }
    for (index, value) in rhs.iter().copied().enumerate() {
        output[target_len - rhs.len() + index] += value;
    }

    output
}

pub(crate) fn reciprocal_or_zero(value: f64) -> f64 {
    if value.is_infinite() {
        0.0
    } else {
        1.0 / value
    }
}

pub(crate) fn safe_sqrt_term(value: f64) -> Result<f64> {
    if value.is_infinite() {
        return Ok(1.0);
    }
    let term = 1.0 - 1.0 / value.powi(2);
    if term < 0.0 {
        return Err(MfsError::Unsupported(
            "current generalized Chebyshev helper only supports |zero| >= 1 real zeros".to_string(),
        ));
    }
    Ok(term.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(diff <= tol, "expected {lhs} ~= {rhs} within {tol}, diff={diff}");
    }

    #[test]
    fn s_to_w_and_back_preserve_coefficients_up_to_roundoff() {
        let coefficients = vec![
            ComplexCoefficient::new(1.0, 0.0),
            ComplexCoefficient::new(0.5, -0.25),
            ComplexCoefficient::new(-0.75, 1.0),
        ];

        let rotated = s_to_w_coefficients(&coefficients);
        let restored = w_to_s_coefficients(&rotated);

        for (expected, actual) in coefficients.iter().zip(restored.iter()) {
            approx_eq(expected.re, actual.re, 1e-12);
            approx_eq(expected.im, actual.im, 1e-12);
        }
    }
}
