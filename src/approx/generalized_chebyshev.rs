//! Cameron/generalized-Chebyshev helper pipeline.
//!
//! This module owns the domain-specific recurrence, ripple-parameter helpers,
//! and `w <-> s` polynomial transforms used by generalized filter synthesis.
//! Generic complex-polynomial storage and root solving live in `complex_poly`.

use crate::error::{MfsError, Result};

use super::generalized_ops::{
    add_descending, complex_from_real, convolve_descending, descending_w_to_ascending_s,
    reciprocal_or_zero, reflect_to_upper_half_plane, rotate_w_root_to_s, s_to_w_coefficients,
    safe_sqrt_term, w_to_s_coefficients,
};
use super::{ComplexCoefficient, ComplexPolynomial, DurandKernerRootSolver};
use super::complex_poly::multiply_by_monic_root;

const COMPLEX_ZERO: ComplexCoefficient = ComplexCoefficient::new(0.0, 0.0);
const COMPLEX_ONE: ComplexCoefficient = ComplexCoefficient::new(1.0, 0.0);

/// Transmission-zero list padded with infinities up to the target order.
#[derive(Debug, Clone, PartialEq)]
pub struct PaddedTransmissionZeros {
    /// Original zeros followed by as many infinities as required.
    pub padded: Vec<f64>,
    /// Number of explicit finite zeros before padding.
    pub finite_count: usize,
}

/// Intermediate polynomials produced by Cameron's recurrence.
#[derive(Debug, Clone, PartialEq)]
pub struct CameronRecurrence {
    /// `U` polynomial stored in descending powers of the helper variable.
    pub u_descending: Vec<f64>,
    /// `V` polynomial stored in descending powers of the helper variable.
    pub v_descending: Vec<f64>,
    /// Recurrence output converted into the `s` domain.
    pub f_s: ComplexPolynomial,
}

/// Extra data required for generalized Chebyshev synthesis paths.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedChebyshevData {
    /// Padded zero list used by the recurrence.
    pub padded_zeros: Vec<f64>,
    /// Number of finite zeros in the original specification.
    pub finite_zero_count: usize,
    /// `F(s)` polynomial from Cameron's recurrence.
    pub f_s: ComplexPolynomial,
    /// `P(s)` polynomial induced by the transmission zeros.
    pub p_s: ComplexPolynomial,
    /// Optional auxiliary `A(s)` polynomial.
    pub a_s: Option<ComplexPolynomial>,
    /// Detailed intermediate results for the auxiliary `A` stage when available.
    pub a_stage: Option<APolynomialStage>,
    /// Optional generalized denominator `E(s)` polynomial.
    pub e_s: Option<ComplexPolynomial>,
    /// Detailed intermediate results for the generalized `E` stage.
    pub e_stage: Option<EPolynomialStage>,
    /// Conventional ripple parameter.
    pub eps: f64,
    /// Ripple parameter adjusted for all-finite-zero cases.
    pub eps_r: f64,
}

/// Detailed intermediate results for the generalized `A`-polynomial stage.
#[derive(Debug, Clone, PartialEq)]
pub struct APolynomialStage {
    /// `P(w)` polynomial obtained from the `s -> w` coefficient transform.
    pub p_w: ComplexPolynomial,
    /// Summed `Psi(w)` contribution built from the remaining finite zeros.
    pub psi_sum: ComplexPolynomial,
    /// Final auxiliary `A(w)` polynomial before root rotation.
    pub a_w: ComplexPolynomial,
    /// Roots of `A(w)` before rotating them back into the `s` plane.
    pub a_w_roots: Vec<ComplexCoefficient>,
    /// Rotated roots in the `s` plane.
    pub a_s_roots: Vec<ComplexCoefficient>,
}

/// Detailed intermediate results for the generalized `E`-polynomial stage.
#[derive(Debug, Clone, PartialEq)]
pub struct EPolynomialStage {
    /// `F(w)` polynomial obtained from the `s -> w` coefficient transform.
    pub f_w: ComplexPolynomial,
    /// `P(w)` polynomial obtained from the `s -> w` coefficient transform.
    pub p_w: ComplexPolynomial,
    /// Raw `E(w)` polynomial before root reflection.
    pub e_w: ComplexPolynomial,
    /// Raw roots of `E(w)` before enforcing upper-half-plane stability.
    pub raw_roots: Vec<ComplexCoefficient>,
    /// Reflected roots used to reconstruct the stable polynomial.
    pub reflected_roots: Vec<ComplexCoefficient>,
    /// Stable `E(w)` reconstructed from reflected roots.
    pub e_w_from_roots: ComplexPolynomial,
    /// Final `E(s)` polynomial in the `s` plane.
    pub e_s: ComplexPolynomial,
    /// Rotated roots in the `s` plane.
    pub e_s_roots: Vec<ComplexCoefficient>,
}

/// Pads a finite-zero list with infinities so it matches the target order.
pub fn pad_transmission_zeros(
    order: usize,
    finite_zeros: &[f64],
) -> Result<PaddedTransmissionZeros> {
    if order == 0 {
        return Err(MfsError::InvalidOrder { order });
    }
    if finite_zeros.len() > order {
        return Err(MfsError::DimensionMismatch {
            expected: order,
            actual: finite_zeros.len(),
        });
    }
    if finite_zeros
        .iter()
        .any(|zero| !zero.is_finite() || *zero == 0.0)
    {
        return Err(MfsError::InvalidTransmissionZero(
            "finite transmission zeros must be non-zero and finite".to_string(),
        ));
    }

    let mut padded = finite_zeros.to_vec();
    padded.resize(order, f64::INFINITY);

    Ok(PaddedTransmissionZeros {
        padded,
        finite_count: finite_zeros.len(),
    })
}

/// Runs the Cameron recurrence used for generalized Chebyshev prototypes.
///
/// The literature derivation is expressed in the helper variable `w`, not in
/// the public `s` domain. The returned `u_descending`/`v_descending` vectors
/// therefore preserve the recurrence's native descending-`w` representation,
/// while `f_s` stores the same polynomial after applying `s = j w`.
pub fn cameron_recursive(padded_zeros: &[f64]) -> Result<CameronRecurrence> {
    if padded_zeros.is_empty() {
        return Err(MfsError::InvalidTransmissionZero(
            "cameron recursion requires at least one transmission zero".to_string(),
        ));
    }

    let first_zero = padded_zeros[0];
    if !first_zero.is_finite() && !first_zero.is_infinite() || first_zero == 0.0 {
        return Err(MfsError::InvalidTransmissionZero(
            "first transmission zero must be non-zero or infinity".to_string(),
        ));
    }

    let mut u = vec![1.0, -reciprocal_or_zero(first_zero)];
    let mut v = vec![safe_sqrt_term(first_zero)?];

    for &next_zero in padded_zeros.iter().skip(1) {
        if next_zero == 0.0 || !next_zero.is_finite() && !next_zero.is_infinite() {
            return Err(MfsError::InvalidTransmissionZero(
                "transmission zeros must be non-zero real values or infinity".to_string(),
            ));
        }

        let temp = safe_sqrt_term(next_zero)?;
        // Update the recurrence in descending-power form to match the reference derivation.
        let u2 = add_descending(
            &convolve_descending(&u, &[1.0, -reciprocal_or_zero(next_zero)]),
            &convolve_descending(
                &[1.0, 0.0, -1.0],
                &v.iter().map(|value| value * temp).collect::<Vec<_>>(),
            ),
        );
        let v2 = add_descending(
            &convolve_descending(&v, &[1.0, -reciprocal_or_zero(next_zero)]),
            &u.iter().map(|value| value * temp).collect::<Vec<_>>(),
        );
        u = u2;
        v = v2;
    }

    let f_s = descending_w_to_ascending_s(&u)?;
    Ok(CameronRecurrence {
        u_descending: u,
        v_descending: v,
        f_s,
    })
}

/// Builds the generalized `P(s)` polynomial from the padded zero list.
pub fn find_p_polynomial(
    order: usize,
    padded_zeros: &[f64],
    finite_count: usize,
) -> Result<ComplexPolynomial> {
    if order == 0 {
        return Err(MfsError::InvalidOrder { order });
    }
    let finite_zeros = padded_zeros
        .iter()
        .copied()
        .filter(|zero| zero.is_finite())
        .collect::<Vec<_>>();

        let mut coefficients = vec![COMPLEX_ONE];
    if order > 1 {
        for zero in finite_zeros {
            let root = ComplexCoefficient::new(0.0, zero);
            coefficients = multiply_by_monic_root(&coefficients, root);
        }
    }

    if (order - finite_count) % 2 == 0 {
        coefficients = coefficients
            .into_iter()
            .map(|coefficient| coefficient * ComplexCoefficient::new(0.0, 1.0))
            .collect();
    }

    ComplexPolynomial::new(coefficients)
}

/// Computes the generalized Chebyshev ripple parameters `eps` and `eps_r`.
pub fn find_eps(
    finite_zero_count: usize,
    p_s: &ComplexPolynomial,
    f_s: &ComplexPolynomial,
    return_loss_db: f64,
    order: usize,
) -> Result<(f64, f64)> {
    if order == 0 {
        return Err(MfsError::InvalidOrder { order });
    }
    if !return_loss_db.is_finite() || return_loss_db <= 0.0 {
        return Err(MfsError::InvalidReturnLoss { return_loss_db });
    }

    let rl_linear = 10_f64.powf(return_loss_db.abs() / 10.0) - 1.0;
    let numerator = p_s.evaluate(ComplexCoefficient::new(0.0, 1.0));
    let denominator = f_s.evaluate(ComplexCoefficient::new(0.0, 1.0));

    let denominator_norm = denominator.norm();
    if denominator_norm <= 1e-18 {
        return Err(MfsError::Unsupported(
            "F polynomial evaluates to zero at s=j".to_string(),
        ));
    }

    let eps = numerator.norm() / denominator_norm / rl_linear.sqrt();
    let eps_r = if finite_zero_count < order {
        1.0
    } else {
        eps / (eps * eps - 1.0).sqrt()
    };
    Ok((eps, eps_r))
}

/// Builds the auxiliary `A` polynomial and returns its transformed roots.
pub fn find_a_polynomial(
    padded_zeros: &[f64],
    order: usize,
    p_s: &ComplexPolynomial,
) -> Result<(Option<ComplexPolynomial>, Vec<ComplexCoefficient>)> {
    let stage = build_a_polynomial_stage(padded_zeros, order, p_s)?;
    Ok((
        stage.as_ref().map(|stage| stage.a_w.clone()),
        stage.map_or_else(
            || vec![ComplexCoefficient::new(0.0, f64::INFINITY); order],
            |stage| stage.a_s_roots,
        ),
    ))
}

/// Builds the auxiliary `A` stage and preserves intermediate helper artifacts.
///
/// `P` enters this function in the public `s` domain, but the Cameron-style
/// `Psi` summation is defined in `w`. We therefore transform `P(s)` into
/// `P(w)`, complete the helper-domain algebra there, solve for roots in `w`,
/// and only then rotate those roots back into the `s` plane.
pub fn build_a_polynomial_stage(
    padded_zeros: &[f64],
    order: usize,
    p_s: &ComplexPolynomial,
) -> Result<Option<APolynomialStage>> {
    let solver = DurandKernerRootSolver;
    let finite_zeros = padded_zeros
        .iter()
        .copied()
        .filter(|zero| zero.is_finite())
        .collect::<Vec<_>>();
    let finite_count = finite_zeros.len();

    if finite_count == 0 {
        return Ok(None);
    }

    let parity_factor = if (order - finite_count) % 2 == 0 {
        ComplexCoefficient::new(0.0, 1.0)
    } else {
        COMPLEX_ONE
    };

    let p_w = ComplexPolynomial::new(s_to_w_coefficients(&p_s.coefficients))?;

    let mut psi_sum = ComplexPolynomial::new(vec![COMPLEX_ZERO])?;
    for (index, zero) in finite_zeros.iter().copied().enumerate() {
        let rn = zero * safe_sqrt_term(zero)?;
        let remaining_zeros = finite_zeros
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(candidate_index, candidate_zero)| {
                (candidate_index != index).then_some(candidate_zero)
            })
            .collect::<Vec<_>>();
        let psi =
            ComplexPolynomial::from_real_roots(&remaining_zeros)?.scale(parity_factor * rn)?;
        psi_sum = psi_sum.add(&psi)?;
    }

    let a_w = p_w
        .scale(complex_from_real((order - finite_count) as f64))?
        .add(&psi_sum)?;
    let a_w_roots = a_w.roots_with(&solver)?;
    // Rotate the roots back into the `s` plane expected by downstream code.
    let a_s_roots = a_w_roots
        .iter()
        .copied()
        .map(rotate_w_root_to_s)
        .collect::<Vec<_>>();
    Ok(Some(APolynomialStage {
        p_w,
        psi_sum,
        a_w,
        a_w_roots,
        a_s_roots,
    }))
}

/// Builds the generalized `E(s)` polynomial and returns its transformed roots.
pub fn find_e_polynomial(
    f_s: &ComplexPolynomial,
    p_s: &ComplexPolynomial,
    eps: f64,
    eps_r: f64,
) -> Result<(ComplexPolynomial, Vec<ComplexCoefficient>)> {
    let stage = build_e_polynomial_stage(f_s, p_s, eps, eps_r)?;
    Ok((stage.e_s.clone(), stage.e_s_roots))
}

/// Builds the generalized `E` stage and preserves intermediate helper artifacts.
///
/// The validated literature formulas for this stage are written in the helper
/// variable `w`. Even though callers pass `F(s)` and `P(s)`, we must first
/// convert them into `F(w)` and `P(w)` using `s = j w`, form `E(w)` there,
/// enforce the stable half-plane in `w`, and then rotate the result back to
/// the public `E(s)` representation. Mixing domains inside this stage causes
/// the literature coefficients and pole set to drift.
pub fn build_e_polynomial_stage(
    f_s: &ComplexPolynomial,
    p_s: &ComplexPolynomial,
    eps: f64,
    eps_r: f64,
) -> Result<EPolynomialStage> {
    let solver = DurandKernerRootSolver;
    let f_w = ComplexPolynomial::new(s_to_w_coefficients(&f_s.coefficients))?;
    let p_w = ComplexPolynomial::new(s_to_w_coefficients(&p_s.coefficients))?;

    let e_w = f_w
        .scale(complex_from_real(1.0 / eps_r))?
        .add(&p_w.scale(complex_from_real(1.0 / eps))?)?;

    let raw_roots = e_w.roots_with(&solver)?;
    // Reflect roots into the stable half-plane before reconstructing the polynomial.
    let reflected_roots = raw_roots
        .iter()
        .copied()
        .map(reflect_to_upper_half_plane)
        .collect::<Vec<_>>();
    let e_w_from_roots = ComplexPolynomial::from_complex_roots(&reflected_roots)?;
    let e_s = ComplexPolynomial::new(w_to_s_coefficients(&e_w_from_roots.coefficients))?;
    let e_s_roots = reflected_roots
        .iter()
        .copied()
        .map(rotate_w_root_to_s)
        .collect::<Vec<_>>();

    Ok(EPolynomialStage {
        f_w,
        p_w,
        e_w,
        raw_roots,
        reflected_roots,
        e_w_from_roots,
        e_s,
        e_s_roots,
    })
}

/// Executes the current generalized Chebyshev helper pipeline end to end.
pub fn synthesize_generalized_chebyshev_data(
    order: usize,
    finite_zeros: &[f64],
    return_loss_db: f64,
) -> Result<GeneralizedChebyshevData> {
    let padded = pad_transmission_zeros(order, finite_zeros)?;
    let recurrence = cameron_recursive(&padded.padded)?;
    let f_s = normalize_to_monic(&recurrence.f_s)?;
    let p_s = find_p_polynomial(order, &padded.padded, padded.finite_count)?;
    let (eps, eps_r) = find_eps(
        padded.finite_count,
        &p_s,
        &f_s,
        return_loss_db,
        order,
    )?;
    let a_stage = build_a_polynomial_stage(&padded.padded, order, &p_s)?;
    let e_stage = build_e_polynomial_stage(&f_s, &p_s, eps, eps_r)?;

    Ok(GeneralizedChebyshevData {
        padded_zeros: padded.padded,
        finite_zero_count: padded.finite_count,
        f_s,
        p_s,
        a_s: a_stage.as_ref().map(|stage| stage.a_w.clone()),
        a_stage,
        e_s: Some(e_stage.e_s.clone()),
        e_stage: Some(e_stage),
        eps,
        eps_r,
    })
}

fn normalize_to_monic(polynomial: &ComplexPolynomial) -> Result<ComplexPolynomial> {
    let leading = polynomial.leading_coefficient();
    if leading.norm_sqr() <= 1e-24 {
        return Err(MfsError::Unsupported(
            "cannot normalize polynomial with zero leading coefficient".to_string(),
        ));
    }

    polynomial.scale(COMPLEX_ONE / leading)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn transmission_zeros_are_padded_to_order() -> Result<()> {
        let padded = pad_transmission_zeros(4, &[1.5, 2.0])?;

        assert_eq!(padded.finite_count, 2);
        assert_eq!(padded.padded.len(), 4);
        assert!(padded.padded[2].is_infinite());
        Ok(())
    }

    #[test]
    fn cameron_recursive_builds_expected_first_order_polynomial() -> Result<()> {
        let recurrence = cameron_recursive(&[2.0])?;

        assert_eq!(recurrence.u_descending, vec![1.0, -0.5]);
        approx_eq(recurrence.v_descending[0], (0.75_f64).sqrt(), 1e-12);
        approx_eq(recurrence.f_s.coefficients[0].re, 0.0, 1e-12);
        approx_eq(recurrence.f_s.coefficients[0].im, -0.5, 1e-12);
        approx_eq(recurrence.f_s.coefficients[1].re, 1.0, 1e-12);
        approx_eq(recurrence.f_s.coefficients[1].im, 0.0, 1e-12);
        Ok(())
    }

    #[test]
    fn p_polynomial_matches_python_style_construction() -> Result<()> {
        let padded = pad_transmission_zeros(3, &[2.0])?;
        let polynomial = find_p_polynomial(3, &padded.padded, padded.finite_count)?;

        approx_eq(polynomial.coefficients[0].re, 2.0, 1e-12);
        approx_eq(polynomial.coefficients[0].im, 0.0, 1e-12);
        approx_eq(polynomial.coefficients[1].re, 0.0, 1e-12);
        approx_eq(polynomial.coefficients[1].im, 1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn eps_matches_return_loss_expression() -> Result<()> {
        let f_s = ComplexPolynomial::new(vec![
            ComplexCoefficient::new(0.0, -0.5),
            ComplexCoefficient::new(1.0, 0.0),
        ])?;
        let p_s = ComplexPolynomial::new(vec![COMPLEX_ONE])?;

        let (eps, eps_r) = find_eps(0, &p_s, &f_s, 20.0, 2)?;
        approx_eq(eps, 0.20100756305184242, 1e-12);
        approx_eq(eps_r, 1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn generalized_chebyshev_data_pipeline_returns_core_polynomials() -> Result<()> {
        let data = synthesize_generalized_chebyshev_data(3, &[2.0], 20.0)?;

        assert_eq!(data.finite_zero_count, 1);
        assert_eq!(data.padded_zeros.len(), 3);
        assert!(data.a_s.is_some());
        assert!(data.a_stage.is_some());
        assert!(data.e_s.is_some());
        assert!(data.e_stage.is_some());
        assert!(data.eps > 0.0);
        assert!(data.eps_r > 0.0);
        Ok(())
    }

    #[test]
    fn detailed_stage_builders_expose_intermediate_generalized_artifacts() -> Result<()> {
        let padded = pad_transmission_zeros(3, &[2.0])?;
        let recurrence = cameron_recursive(&padded.padded)?;
        let p_s = find_p_polynomial(3, &padded.padded, padded.finite_count)?;
        let (eps, eps_r) = find_eps(
            padded.finite_count,
            &p_s,
            &recurrence.f_s,
            20.0,
            3,
        )?;

        let a_stage = build_a_polynomial_stage(&padded.padded, 3, &p_s)?
            .expect("finite-zero case should build an A stage");
        let e_stage = build_e_polynomial_stage(&recurrence.f_s, &p_s, eps, eps_r)?;

        assert_eq!(a_stage.a_w_roots.len(), 1);
        assert_eq!(a_stage.a_s_roots.len(), 1);
        assert_eq!(e_stage.raw_roots.len(), e_stage.reflected_roots.len());
        assert_eq!(e_stage.e_s_roots.len(), e_stage.reflected_roots.len());
        Ok(())
    }
}
