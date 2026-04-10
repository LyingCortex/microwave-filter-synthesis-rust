use super::{ComplexCoefficient, ComplexPolynomial, GeneralizedChebyshevData};
use crate::error::{MfsError, Result};

/// Polynomial data produced by the approximation stage.
#[derive(Debug, Clone, PartialEq)]
pub struct PolynomialSet {
    /// Prototype order carried through the rest of the pipeline.
    pub order: usize,
    /// Chebyshev ripple factor derived from return loss.
    pub ripple_factor: f64,
    /// Generalized ripple parameter used by coupling-matrix synthesis formulas.
    pub eps: f64,
    /// Adjusted ripple parameter for all-finite-zero cases.
    pub eps_r: f64,
    /// Finite transmission zeros expressed in normalized coordinates.
    pub transmission_zeros_normalized: Vec<f64>,
    /// Denominator-like polynomial coefficients in ascending-power order.
    pub e: ComplexPolynomial,
    /// Numerator-like polynomial coefficients in ascending-power order.
    pub f: ComplexPolynomial,
    /// Finite-zero polynomial coefficients in ascending-power order.
    pub p: ComplexPolynomial,
    /// Extra generalized Chebyshev data when finite zeros are present.
    pub generalized: Option<GeneralizedChebyshevData>,
}

impl PolynomialSet {
    /// Creates a validated polynomial bundle.
    pub fn new(
        order: usize,
        ripple_factor: f64,
        eps: f64,
        eps_r: f64,
        transmission_zeros_normalized: Vec<f64>,
        e: Vec<f64>,
        f: Vec<f64>,
        p: Vec<f64>,
    ) -> Result<Self> {
        let set = Self {
            order,
            ripple_factor,
            eps,
            eps_r,
            transmission_zeros_normalized,
            e: complex_polynomial_from_real_coefficients(e)?,
            f: complex_polynomial_from_real_coefficients(f)?,
            p: complex_polynomial_from_real_coefficients(p)?,
            generalized: None,
        };
        set.validate()?;
        Ok(set)
    }

    /// Attaches generalized Chebyshev helper data to the bundle.
    pub fn with_generalized(mut self, generalized: GeneralizedChebyshevData) -> Self {
        self.eps = generalized.eps;
        self.eps_r = generalized.eps_r;
        if let Some(e_s) = generalized.e_s.as_ref() {
            self.e = e_s.clone();
        }
        self.f = generalized.f_s.clone();
        self.p = generalized.p_s.clone();
        self.generalized = Some(generalized);
        self
    }

    /// Validates vector lengths and scalar ranges before downstream use.
    pub fn validate(&self) -> Result<()> {
        if self.order == 0 {
            return Err(MfsError::InvalidOrder { order: self.order });
        }
        if !self.ripple_factor.is_finite() || self.ripple_factor < 0.0 {
            return Err(MfsError::Unsupported(
                "ripple factor must be finite and non-negative".to_string(),
            ));
        }
        if !self.eps.is_finite() || self.eps <= 0.0 {
            return Err(MfsError::Unsupported(
                "eps must be finite and positive".to_string(),
            ));
        }
        if !self.eps_r.is_finite() || self.eps_r <= 0.0 {
            return Err(MfsError::Unsupported(
                "eps_r must be finite and positive".to_string(),
            ));
        }
        if self.e.coefficients.len() != self.order + 1 {
            return Err(MfsError::DimensionMismatch {
                expected: self.order + 1,
                actual: self.e.coefficients.len(),
            });
        }
        if self.f.coefficients.len() != self.order + 1 {
            return Err(MfsError::DimensionMismatch {
                expected: self.order + 1,
                actual: self.f.coefficients.len(),
            });
        }
        if self.p.coefficients.len() > self.order {
            return Err(MfsError::DimensionMismatch {
                expected: self.order,
                actual: self.p.coefficients.len(),
            });
        }
        if self.transmission_zeros_normalized.iter().any(|value| !value.is_finite())
            || self
                .e
                .coefficients
                .iter()
                .chain(self.f.coefficients.iter())
                .chain(self.p.coefficients.iter())
                .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(MfsError::Unsupported(
                "polynomial coefficients and transmission zeros must be finite".to_string(),
            ));
        }

        Ok(())
    }

    /// Returns the degree of the `E` polynomial.
    pub fn e_degree(&self) -> usize {
        self.e.coefficients.len().saturating_sub(1)
    }

    /// Returns the degree of the `F` polynomial.
    pub fn f_degree(&self) -> usize {
        self.f.coefficients.len().saturating_sub(1)
    }

    /// Returns the degree of the `P` polynomial.
    pub fn p_degree(&self) -> usize {
        self.p.coefficients.len().saturating_sub(1)
    }

    /// Returns the `E` polynomial projected onto the closest real-valued representation.
    pub fn e_real_projection(&self) -> Vec<f64> {
        normalize_and_project_complex_polynomial(&self.e, self.order + 1)
            .unwrap_or_else(|| vec![0.0; self.order + 1])
    }

    /// Returns the `F` polynomial projected onto the closest real-valued representation.
    pub fn f_real_projection(&self) -> Vec<f64> {
        normalize_and_project_complex_polynomial(&self.f, self.order + 1)
            .unwrap_or_else(|| vec![0.0; self.order + 1])
    }

    /// Returns the `P` polynomial projected onto the closest real-valued representation.
    pub fn p_real_projection(&self) -> Vec<f64> {
        normalize_and_project_complex_polynomial(&self.p, self.p.coefficients.len())
            .unwrap_or_else(|| vec![0.0; self.p.coefficients.len()])
    }
}

/// Converts passband return loss into the Chebyshev ripple factor.
pub fn chebyshev_ripple_factor(return_loss_db: f64) -> f64 {
    let power_ratio = 10_f64.powf(return_loss_db / 10.0) - 1.0;
    1.0 / power_ratio.sqrt()
}

/// Builds a monic real polynomial from the provided real roots.
pub fn monic_polynomial_from_real_roots(roots: &[f64]) -> Vec<f64> {
    let mut coefficients = vec![1.0];

    for &root in roots {
        let mut next = vec![0.0; coefficients.len() + 1];
        for (index, coeff) in coefficients.iter().copied().enumerate() {
            next[index] += coeff;
            next[index + 1] -= coeff * root;
        }
        coefficients = next;
    }

    coefficients
}

fn normalize_and_project_complex_polynomial(
    polynomial: &ComplexPolynomial,
    expected_len: usize,
) -> Option<Vec<f64>> {
    let coefficients = &polynomial.coefficients;
    if coefficients.is_empty() {
        return None;
    }

    let already_real = coefficients.iter().all(|coefficient| coefficient.im.abs() <= 1e-12);
    if already_real {
        let mut projected = coefficients.iter().map(|coefficient| coefficient.re).collect::<Vec<_>>();
        if projected.len() < expected_len {
            projected.resize(expected_len, 0.0);
        } else if projected.len() > expected_len {
            projected.truncate(expected_len);
        }
        return Some(projected);
    }

    let phase_candidates = [
        ComplexCoefficient::new(1.0, 0.0),
        ComplexCoefficient::new(0.0, -1.0),
        ComplexCoefficient::new(-1.0, 0.0),
        ComplexCoefficient::new(0.0, 1.0),
    ];

    let mut best_projection: Option<Vec<f64>> = None;
    let mut best_imag_residual = f64::INFINITY;

    for phase in phase_candidates {
        let rotated = coefficients
            .iter()
            .map(|coefficient| *coefficient * phase)
            .collect::<Vec<_>>();
        let real_weight = rotated.iter().map(|coefficient| coefficient.re.abs()).sum::<f64>();
        let imag_weight = rotated.iter().map(|coefficient| coefficient.im.abs()).sum::<f64>();
        if real_weight <= 1e-18 {
            continue;
        }
        if imag_weight < best_imag_residual {
            best_imag_residual = imag_weight;
            best_projection = Some(rotated.iter().map(|coefficient| coefficient.re).collect());
        }
    }

    let mut projected = best_projection?;

    if let Some(leading_nonzero) = projected.iter().rfind(|value| value.abs() > 1e-18) {
        if *leading_nonzero < 0.0 {
            for coefficient in &mut projected {
                *coefficient = -*coefficient;
            }
        }
    }

    if projected.len() < expected_len {
        projected.resize(expected_len, 0.0);
    } else if projected.len() > expected_len {
        projected.truncate(expected_len);
    }

    Some(projected)
}

fn complex_polynomial_from_real_coefficients(coefficients: Vec<f64>) -> Result<ComplexPolynomial> {
    ComplexPolynomial::new(
        coefficients
            .into_iter()
            .map(|coefficient| ComplexCoefficient::new(coefficient, 0.0))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::{
        ComplexCoefficient, ComplexPolynomial, GeneralizedChebyshevData,
    };

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn polynomial_set_reports_degrees() {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![-1.0, 1.0],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, -2.0],
        )
        .expect("valid polynomial set");

        assert_eq!(polynomials.e_degree(), 3);
        assert_eq!(polynomials.f_degree(), 3);
        assert_eq!(polynomials.p_degree(), 1);
        approx_eq(polynomials.eps, 0.1, 1e-12);
        approx_eq(polynomials.eps_r, 1.0, 1e-12);
        assert!(polynomials.generalized.is_none());
    }

    #[test]
    fn polynomial_set_rejects_mismatched_coefficient_lengths() {
        let error =
            PolynomialSet::new(3, 0.1, 0.1, 1.0, vec![], vec![1.0], vec![0.8, 0.6], vec![1.0])
                .expect_err("mismatched coefficient lengths must fail");

        assert!(matches!(error, MfsError::DimensionMismatch { .. }));
    }

    #[test]
    fn ripple_factor_matches_return_loss() {
        approx_eq(chebyshev_ripple_factor(20.0), 0.10050378152592121, 1e-12);
    }

    #[test]
    fn monic_polynomial_is_built_from_real_roots() {
        let coefficients = monic_polynomial_from_real_roots(&[-2.0, 1.5]);

        approx_eq(coefficients[0], 1.0, 1e-12);
        approx_eq(coefficients[1], 0.5, 1e-12);
        approx_eq(coefficients[2], -3.0, 1e-12);
    }

    #[test]
    fn generalized_attachment_overrides_eps_metadata() {
        let generalized = GeneralizedChebyshevData {
            padded_zeros: vec![2.0, f64::INFINITY, f64::INFINITY],
            finite_zero_count: 1,
            f_s: ComplexPolynomial::new(vec![ComplexCoefficient::new(1.0, 0.0)])
                .expect("valid polynomial"),
            p_s: ComplexPolynomial::new(vec![ComplexCoefficient::new(1.0, 0.0)])
                .expect("valid polynomial"),
            a_s: None,
            a_stage: None,
            e_s: None,
            e_stage: None,
            eps: 0.23,
            eps_r: 1.7,
        };
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![2.0],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, -2.0],
        )
        .expect("valid polynomial set")
        .with_generalized(generalized);

        approx_eq(polynomials.eps, 0.23, 1e-12);
        approx_eq(polynomials.eps_r, 1.7, 1e-12);
    }

    #[test]
    fn generalized_attachment_replaces_summary_polynomials_with_complex_forms() {
        let generalized = GeneralizedChebyshevData {
            padded_zeros: vec![2.0, f64::INFINITY, f64::INFINITY],
            finite_zero_count: 1,
            f_s: ComplexPolynomial::new(vec![
                ComplexCoefficient::new(0.0, 0.5),
                ComplexCoefficient::new(0.0, -0.25),
            ])
            .expect("valid polynomial"),
            p_s: ComplexPolynomial::new(vec![ComplexCoefficient::new(1.0, 0.0)])
                .expect("valid polynomial"),
            a_s: None,
            a_stage: None,
            e_s: Some(
                ComplexPolynomial::new(vec![
                    ComplexCoefficient::new(1.0, 0.0),
                    ComplexCoefficient::new(-0.5, 0.0),
                ])
                .expect("valid polynomial"),
            ),
            e_stage: None,
            eps: 0.23,
            eps_r: 1.7,
        };
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![2.0],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, -2.0],
        )
        .expect("valid polynomial set")
        .with_generalized(generalized);

        approx_eq(polynomials.e.coefficients[0].re, 1.0, 1e-12);
        approx_eq(polynomials.e.coefficients[1].re, -0.5, 1e-12);
        approx_eq(polynomials.f.coefficients[0].im, 0.5, 1e-12);
        approx_eq(polynomials.f.coefficients[1].im, -0.25, 1e-12);
        approx_eq(polynomials.p.coefficients[0].re, 1.0, 1e-12);
        let projected_f = polynomials.f_real_projection();
        approx_eq(projected_f[0], -0.5, 1e-12);
        approx_eq(projected_f[1], 0.25, 1e-12);
    }

    #[test]
    fn generalized_projection_normalizes_global_phase_before_projecting() {
        let polynomial = ComplexPolynomial::new(vec![
            ComplexCoefficient::new(0.0, 2.0),
            ComplexCoefficient::new(0.0, -1.0),
            ComplexCoefficient::new(0.0, 0.5),
        ])
        .expect("valid polynomial");

        let projected =
            normalize_and_project_complex_polynomial(&polynomial, 3).expect("projection exists");

        approx_eq(projected[0], 2.0, 1e-12);
        approx_eq(projected[1], -1.0, 1e-12);
        approx_eq(projected[2], 0.5, 1e-12);
    }
}
