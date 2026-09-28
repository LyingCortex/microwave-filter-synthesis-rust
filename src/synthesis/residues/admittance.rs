use crate::approx::{ComplexCoefficient, ComplexPolynomial, PolynomialSet};
use crate::error::{MfsError, Result};

fn complex_from_real(value: f64) -> ComplexCoefficient {
    ComplexCoefficient::new(value, 0.0)
}

/// Polynomial-form Y parameters derived from the approximation output.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittancePolynomials {
    /// Common denominator polynomial used by all Y parameters.
    pub denominator: ComplexPolynomial,
    /// Numerator of normalized input admittance.
    pub y11: ComplexPolynomial,
    /// Numerator of transfer admittance.
    pub y12: ComplexPolynomial,
    /// Numerator of output admittance.
    pub y22: ComplexPolynomial,
}

/// Builds polynomial-form Y parameters from generalized Chebyshev helper data.
pub fn synthesize_admittance_polynomials(
    polynomials: &PolynomialSet,
) -> Result<AdmittancePolynomials> {
    let generalized = polynomials.generalized.as_ref().ok_or_else(|| {
        MfsError::PreconditionViolation(
            "admittance synthesis requires generalized Chebyshev helper data".to_string(),
        )
    })?;
    let e_s = generalized.e_s.as_ref().ok_or_else(|| {
        MfsError::PreconditionViolation("generalized helper data is missing E(s)".to_string())
    })?;
    let f_over_eps_r = generalized
        .f_s
        .scale(complex_from_real(1.0 / polynomials.eps_r))?;
    let p_transfer = generalized
        .p_s
        .scale(complex_from_real(-2.0 / polynomials.eps))?;

    // Normalize E(s) to prevent leading-coefficient cancellation with F(s)/eps_r.
    // The w_to_s transform introduces a (-j)^N factor on E's leading coefficient.
    // For N mod 4 == 2 (orders 2, 6, 10, ...), E_lead = -1 while F_lead = +1,
    // causing exact cancellation in E + F/eps_r. We detect this and negate E(s)
    // to make the leading coefficients additive rather than cancelling.
    let e_leading = e_s.leading_coefficient();
    let f_leading = f_over_eps_r.leading_coefficient();
    let sum_leading = e_leading + f_leading;
    let e_normalized = if sum_leading.norm_sqr() < 1e-20 && e_leading.re < -1e-15 {
        // Leading coefficients cancel because E_lead ≈ -F_lead.
        // Negate E(s) so that E_lead becomes +1, making the sum ≈ 2.
        e_s.scale(complex_from_real(-1.0))?
    } else {
        e_s.clone()
    };

    // Cameron's admittance formula uses E(-s)* and F(-s)*, which include a
    // global (-1)^N factor relative to the simple alternating_conjugate
    // (conj(c_k) * (-1)^k). Without this factor the leading coefficients of
    // the denominator cancel for odd N, producing poles off the imaginary
    // axis and a fully-reflective coupling matrix.
    let order = polynomials.order;
    let parity_sign = if order % 2 == 0 {
        complex_from_real(1.0)
    } else {
        complex_from_real(-1.0)
    };

    let conjugated_e = e_normalized.alternating_conjugate()?.scale(parity_sign)?;
    let conjugated_f = f_over_eps_r.alternating_conjugate()?.scale(parity_sign)?;

    let denominator = e_normalized
        .add(&f_over_eps_r)?
        .add(&conjugated_f)?
        .add(&conjugated_e)?;
    let y11 = e_normalized
        .sub(&f_over_eps_r)?
        .add(&conjugated_f)?
        .sub(&conjugated_e)?;
    let y22 = e_normalized
        .add(&f_over_eps_r)?
        .sub(&conjugated_f)?
        .sub(&conjugated_e)?;

    // Fallback: if the symmetric construction still causes degree loss
    // (e.g., due to numerical noise), use the simpler Cameron formulation.
    if denominator.degree() < order {
        let alt_denominator = e_normalized.add(&f_over_eps_r)?;
        let alt_y11 = e_normalized.sub(&f_over_eps_r)?;
        let alt_y22 = conjugated_e.sub(&conjugated_f)?;

        return Ok(AdmittancePolynomials {
            denominator: alt_denominator,
            y11: alt_y11,
            y12: p_transfer,
            y22: alt_y22,
        });
    }

    Ok(AdmittancePolynomials {
        denominator,
        y11,
        y12: p_transfer,
        y22,
    })
}
