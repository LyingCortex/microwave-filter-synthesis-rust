use crate::error::{MfsError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ComplexCoefficient {
    pub re: f64,
    pub im: f64,
}

impl ComplexCoefficient {
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };

    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn from_real(value: f64) -> Self {
        Self { re: value, im: 0.0 }
    }

    pub fn norm(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    fn mul_i_pow(self, power: usize) -> Self {
        match power % 4 {
            0 => self,
            1 => Self::new(-self.im, self.re),
            2 => Self::new(-self.re, -self.im),
            _ => Self::new(self.im, -self.re),
        }
    }
}

impl std::ops::Add for ComplexCoefficient {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::AddAssign for ComplexCoefficient {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl std::ops::Sub for ComplexCoefficient {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for ComplexCoefficient {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::Mul<f64> for ComplexCoefficient {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl std::ops::Mul<ComplexCoefficient> for f64 {
    type Output = ComplexCoefficient;

    fn mul(self, rhs: ComplexCoefficient) -> Self::Output {
        rhs * self
    }
}

impl std::ops::Div for ComplexCoefficient {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

impl std::ops::Neg for ComplexCoefficient {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ComplexPolynomial {
    pub coefficients: Vec<ComplexCoefficient>,
}

impl ComplexPolynomial {
    pub fn new(coefficients: Vec<ComplexCoefficient>) -> Result<Self> {
        if coefficients.is_empty() {
            return Err(MfsError::Unsupported(
                "complex polynomial must contain at least one coefficient".to_string(),
            ));
        }
        if coefficients
            .iter()
            .any(|coeff| !coeff.re.is_finite() || !coeff.im.is_finite())
        {
            return Err(MfsError::Unsupported(
                "complex polynomial coefficients must be finite".to_string(),
            ));
        }
        Ok(Self { coefficients })
    }

    pub fn evaluate(&self, x: ComplexCoefficient) -> ComplexCoefficient {
        self.coefficients
            .iter()
            .rev()
            .copied()
            .fold(ComplexCoefficient::ZERO, |acc, coeff| acc * x + coeff)
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len().saturating_sub(1)
    }

    pub fn scale(&self, scalar: ComplexCoefficient) -> Result<Self> {
        Self::new(
            self.coefficients
                .iter()
                .copied()
                .map(|coefficient| coefficient * scalar)
                .collect(),
        )
    }

    pub fn add(&self, rhs: &Self) -> Result<Self> {
        let target_len = self.coefficients.len().max(rhs.coefficients.len());
        let mut coefficients = vec![ComplexCoefficient::ZERO; target_len];

        for (index, coefficient) in self.coefficients.iter().copied().enumerate() {
            coefficients[index] += coefficient;
        }
        for (index, coefficient) in rhs.coefficients.iter().copied().enumerate() {
            coefficients[index] += coefficient;
        }

        Self::new(trim_trailing_complex_zeros(coefficients))
    }

    pub fn from_real_roots(roots: &[f64]) -> Result<Self> {
        let mut coefficients = vec![ComplexCoefficient::ONE];
        for &root in roots {
            coefficients =
                multiply_by_monic_root(&coefficients, ComplexCoefficient::from_real(root));
        }
        Self::new(coefficients)
    }

    pub fn from_complex_roots(roots: &[ComplexCoefficient]) -> Result<Self> {
        let mut coefficients = vec![ComplexCoefficient::ONE];
        for &root in roots {
            coefficients = multiply_by_monic_root(&coefficients, root);
        }
        Self::new(coefficients)
    }

    pub fn roots(&self) -> Result<Vec<ComplexCoefficient>> {
        let degree = self.degree();
        if degree == 0 {
            return Ok(Vec::new());
        }

        let leading = *self.coefficients.last().ok_or_else(|| {
            MfsError::Unsupported("polynomial is missing a leading coefficient".to_string())
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::Unsupported(
                "polynomial leading coefficient must be non-zero".to_string(),
            ));
        }

        let normalized = self
            .coefficients
            .iter()
            .copied()
            .map(|coefficient| coefficient / leading)
            .collect::<Vec<_>>();
        let radius = 1.0
            + normalized[..degree]
                .iter()
                .copied()
                .map(ComplexCoefficient::norm)
                .fold(0.0_f64, f64::max);

        let mut roots = (0..degree)
            .map(|index| {
                let angle = 2.0 * std::f64::consts::PI * index as f64 / degree as f64;
                ComplexCoefficient::new(radius * angle.cos(), radius * angle.sin())
            })
            .collect::<Vec<_>>();

        for _ in 0..128 {
            let mut max_delta = 0.0_f64;
            for index in 0..degree {
                let root = roots[index];
                let mut denominator = ComplexCoefficient::ONE;
                for (other_index, other_root) in roots.iter().copied().enumerate() {
                    if index != other_index {
                        denominator = denominator * (root - other_root);
                    }
                }

                if denominator.norm_sqr() <= 1e-24 {
                    continue;
                }

                let delta = evaluate_monic_polynomial(&normalized, root) / denominator;
                roots[index] = root - delta;
                max_delta = max_delta.max(delta.norm());
            }

            if max_delta <= 1e-12 {
                return Ok(roots);
            }
        }

        Err(MfsError::Unsupported(
            "complex polynomial root solver did not converge".to_string(),
        ))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PaddedTransmissionZeros {
    pub padded: Vec<f64>,
    pub finite_count: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CameronRecurrence {
    pub u_descending: Vec<f64>,
    pub v_descending: Vec<f64>,
    pub f_s: ComplexPolynomial,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedChebyshevData {
    pub padded_zeros: Vec<f64>,
    pub finite_zero_count: usize,
    pub f_s: ComplexPolynomial,
    pub p_s: ComplexPolynomial,
    pub a_s: Option<ComplexPolynomial>,
    pub e_s: Option<ComplexPolynomial>,
    pub eps: f64,
    pub eps_r: f64,
}

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

    let mut u = vec![1.0, reciprocal_or_zero(first_zero).neg()];
    let mut v = vec![safe_sqrt_term(first_zero)?];

    for &next_zero in padded_zeros.iter().skip(1) {
        if next_zero == 0.0 || !next_zero.is_finite() && !next_zero.is_infinite() {
            return Err(MfsError::InvalidTransmissionZero(
                "transmission zeros must be non-zero real values or infinity".to_string(),
            ));
        }

        let temp = safe_sqrt_term(next_zero)?;
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

    let mut coefficients = vec![ComplexCoefficient::ONE];
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

pub fn find_a_polynomial(
    padded_zeros: &[f64],
    order: usize,
    p_s: &ComplexPolynomial,
) -> Result<(Option<ComplexPolynomial>, Vec<ComplexCoefficient>)> {
    let finite_zeros = padded_zeros
        .iter()
        .copied()
        .filter(|zero| zero.is_finite())
        .collect::<Vec<_>>();
    let finite_count = finite_zeros.len();

    if finite_count == 0 {
        return Ok((
            None,
            vec![ComplexCoefficient::new(0.0, f64::INFINITY); order],
        ));
    }

    let parity_factor = if (order - finite_count) % 2 == 0 {
        ComplexCoefficient::new(0.0, 1.0)
    } else {
        ComplexCoefficient::ONE
    };

    let p_w = ComplexPolynomial::new(
        p_s.coefficients
            .iter()
            .copied()
            .enumerate()
            .map(|(index, coefficient)| coefficient.mul_i_pow(index))
            .collect(),
    )?;

    let mut psi_sum = ComplexPolynomial::new(vec![ComplexCoefficient::ZERO])?;
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
        .scale(ComplexCoefficient::from_real((order - finite_count) as f64))?
        .add(&psi_sum)?;
    let a_w_roots = a_w.roots()?;
    let a_s_roots = a_w_roots
        .iter()
        .copied()
        .map(|root| ComplexCoefficient::new(-root.im, root.re))
        .collect::<Vec<_>>();
    Ok((Some(a_w), a_s_roots))
}

pub fn find_e_polynomial(
    f_s: &ComplexPolynomial,
    p_s: &ComplexPolynomial,
    eps: f64,
    eps_r: f64,
) -> Result<(ComplexPolynomial, Vec<ComplexCoefficient>)> {
    let f_w = ComplexPolynomial::new(
        f_s.coefficients
            .iter()
            .copied()
            .enumerate()
            .map(|(index, coefficient)| coefficient.mul_i_pow(index))
            .collect(),
    )?;
    let p_w = ComplexPolynomial::new(
        p_s.coefficients
            .iter()
            .copied()
            .enumerate()
            .map(|(index, coefficient)| coefficient.mul_i_pow(index))
            .collect(),
    )?;

    let e_w = f_w
        .scale(ComplexCoefficient::from_real(1.0 / eps_r))?
        .add(&p_w.scale(ComplexCoefficient::new(0.0, 1.0 / eps))?)?;

    let raw_roots = e_w.roots()?;
    let reflected_roots = raw_roots
        .into_iter()
        .map(|root| {
            if root.im >= 0.0 {
                root
            } else {
                ComplexCoefficient::new(root.re, -root.im)
            }
        })
        .collect::<Vec<_>>();
    let e_w_from_roots = ComplexPolynomial::from_complex_roots(&reflected_roots)?;
    let e_s = ComplexPolynomial::new(
        e_w_from_roots
            .coefficients
            .iter()
            .copied()
            .enumerate()
            .map(|(index, coefficient)| coefficient.mul_i_pow((4 - index % 4) % 4))
            .collect(),
    )?;
    let e_s_roots = reflected_roots
        .iter()
        .copied()
        .map(|root| ComplexCoefficient::new(-root.im, root.re))
        .collect::<Vec<_>>();

    Ok((e_s, e_s_roots))
}

pub fn synthesize_generalized_chebyshev_data(
    order: usize,
    finite_zeros: &[f64],
    return_loss_db: f64,
) -> Result<GeneralizedChebyshevData> {
    let padded = pad_transmission_zeros(order, finite_zeros)?;
    let recurrence = cameron_recursive(&padded.padded)?;
    let p_s = find_p_polynomial(order, &padded.padded, padded.finite_count)?;
    let (eps, eps_r) = find_eps(
        padded.finite_count,
        &p_s,
        &recurrence.f_s,
        return_loss_db,
        order,
    )?;
    let (a_s, _) = find_a_polynomial(&padded.padded, order, &p_s)?;
    let (e_s, _) = find_e_polynomial(&recurrence.f_s, &p_s, eps, eps_r)?;

    Ok(GeneralizedChebyshevData {
        padded_zeros: padded.padded,
        finite_zero_count: padded.finite_count,
        f_s: recurrence.f_s,
        p_s,
        a_s,
        e_s: Some(e_s),
        eps,
        eps_r,
    })
}

fn descending_w_to_ascending_s(descending: &[f64]) -> Result<ComplexPolynomial> {
    let order = descending.len().saturating_sub(1);
    let mut ascending = vec![ComplexCoefficient::ZERO; descending.len()];

    for (index, coefficient) in descending.iter().copied().enumerate() {
        let power = order - index;
        ascending[power] = ComplexCoefficient::from_real(coefficient).mul_i_pow(index);
    }

    ComplexPolynomial::new(ascending)
}

fn multiply_by_monic_root(
    coefficients: &[ComplexCoefficient],
    root: ComplexCoefficient,
) -> Vec<ComplexCoefficient> {
    let mut next = vec![ComplexCoefficient::ZERO; coefficients.len() + 1];
    for (index, coefficient) in coefficients.iter().copied().enumerate() {
        next[index] = next[index] + coefficient * (-root);
        next[index + 1] = next[index + 1] + coefficient;
    }
    next
}

fn trim_trailing_complex_zeros(
    mut coefficients: Vec<ComplexCoefficient>,
) -> Vec<ComplexCoefficient> {
    while coefficients.len() > 1
        && coefficients
            .last()
            .is_some_and(|coefficient| coefficient.norm_sqr() <= 1e-24)
    {
        coefficients.pop();
    }
    coefficients
}

fn evaluate_monic_polynomial(
    coefficients: &[ComplexCoefficient],
    x: ComplexCoefficient,
) -> ComplexCoefficient {
    coefficients
        .iter()
        .rev()
        .copied()
        .fold(ComplexCoefficient::ZERO, |acc, coefficient| {
            acc * x + coefficient
        })
}

fn convolve_descending(lhs: &[f64], rhs: &[f64]) -> Vec<f64> {
    let mut output = vec![0.0; lhs.len() + rhs.len() - 1];
    for (left_index, left) in lhs.iter().copied().enumerate() {
        for (right_index, right) in rhs.iter().copied().enumerate() {
            output[left_index + right_index] += left * right;
        }
    }
    output
}

fn add_descending(lhs: &[f64], rhs: &[f64]) -> Vec<f64> {
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

fn reciprocal_or_zero(value: f64) -> f64 {
    if value.is_infinite() {
        0.0
    } else {
        1.0 / value
    }
}

fn safe_sqrt_term(value: f64) -> Result<f64> {
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

trait NegExt {
    fn neg(self) -> Self;
}

impl NegExt for f64 {
    fn neg(self) -> Self {
        -self
    }
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
        let p_s = ComplexPolynomial::new(vec![ComplexCoefficient::ONE])?;

        let (eps, eps_r) = find_eps(0, &p_s, &f_s, 20.0, 2)?;
        approx_eq(eps, 0.20100756305184242, 1e-12);
        approx_eq(eps_r, 1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn complex_polynomial_root_solver_recovers_known_roots() -> Result<()> {
        let polynomial = ComplexPolynomial::from_real_roots(&[1.0, 2.0])?;
        let mut roots = polynomial.roots()?;
        roots.sort_by(|lhs, rhs| {
            lhs.re
                .partial_cmp(&rhs.re)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        approx_eq(roots[0].re, 1.0, 1e-8);
        approx_eq(roots[0].im, 0.0, 1e-8);
        approx_eq(roots[1].re, 2.0, 1e-8);
        approx_eq(roots[1].im, 0.0, 1e-8);
        Ok(())
    }

    #[test]
    fn generalized_chebyshev_data_pipeline_returns_core_polynomials() -> Result<()> {
        let data = synthesize_generalized_chebyshev_data(3, &[2.0], 20.0)?;

        assert_eq!(data.finite_zero_count, 1);
        assert_eq!(data.padded_zeros.len(), 3);
        assert!(data.a_s.is_some());
        assert!(data.e_s.is_some());
        assert!(data.eps > 0.0);
        assert!(data.eps_r > 0.0);
        Ok(())
    }
}
