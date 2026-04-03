use super::GeneralizedChebyshevData;
use crate::error::{MfsError, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct PolynomialSet {
    pub order: usize,
    pub ripple_factor: f64,
    pub transmission_zeros_normalized: Vec<f64>,
    pub e: Vec<f64>,
    pub f: Vec<f64>,
    pub p: Vec<f64>,
    pub generalized: Option<GeneralizedChebyshevData>,
}

impl PolynomialSet {
    pub fn new(
        order: usize,
        ripple_factor: f64,
        transmission_zeros_normalized: Vec<f64>,
        e: Vec<f64>,
        f: Vec<f64>,
        p: Vec<f64>,
    ) -> Result<Self> {
        let set = Self {
            order,
            ripple_factor,
            transmission_zeros_normalized,
            e,
            f,
            p,
            generalized: None,
        };
        set.validate()?;
        Ok(set)
    }

    pub fn with_generalized(mut self, generalized: GeneralizedChebyshevData) -> Self {
        self.generalized = Some(generalized);
        self
    }

    pub fn validate(&self) -> Result<()> {
        if self.order == 0 {
            return Err(MfsError::InvalidOrder { order: self.order });
        }
        if !self.ripple_factor.is_finite() || self.ripple_factor < 0.0 {
            return Err(MfsError::Unsupported(
                "ripple factor must be finite and non-negative".to_string(),
            ));
        }
        if self.e.len() != self.order + 1 {
            return Err(MfsError::DimensionMismatch {
                expected: self.order + 1,
                actual: self.e.len(),
            });
        }
        if self.f.len() != self.order + 1 {
            return Err(MfsError::DimensionMismatch {
                expected: self.order + 1,
                actual: self.f.len(),
            });
        }
        if self.p.len() > self.order {
            return Err(MfsError::DimensionMismatch {
                expected: self.order,
                actual: self.p.len(),
            });
        }
        if self
            .transmission_zeros_normalized
            .iter()
            .chain(self.e.iter())
            .chain(self.f.iter())
            .chain(self.p.iter())
            .any(|value| !value.is_finite())
        {
            return Err(MfsError::Unsupported(
                "polynomial coefficients and transmission zeros must be finite".to_string(),
            ));
        }

        Ok(())
    }

    pub fn e_degree(&self) -> usize {
        self.e.len().saturating_sub(1)
    }

    pub fn f_degree(&self) -> usize {
        self.f.len().saturating_sub(1)
    }

    pub fn p_degree(&self) -> usize {
        self.p.len().saturating_sub(1)
    }
}

pub fn chebyshev_ripple_factor(return_loss_db: f64) -> f64 {
    let power_ratio = 10_f64.powf(return_loss_db / 10.0) - 1.0;
    1.0 / power_ratio.sqrt()
}

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
    fn polynomial_set_reports_degrees() {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            vec![-1.0, 1.0],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, -2.0],
        )
        .expect("valid polynomial set");

        assert_eq!(polynomials.e_degree(), 3);
        assert_eq!(polynomials.f_degree(), 3);
        assert_eq!(polynomials.p_degree(), 1);
        assert!(polynomials.generalized.is_none());
    }

    #[test]
    fn polynomial_set_rejects_mismatched_coefficient_lengths() {
        let error = PolynomialSet::new(3, 0.1, vec![], vec![1.0], vec![0.8, 0.6], vec![1.0])
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
}
