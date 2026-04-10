//! Reusable complex-polynomial primitives for the approximation layer.
//!
//! This module deliberately stays domain-light: it owns coefficient storage,
//! basic polynomial algebra, and the current default root-solving strategy.
//! Filter-specific transforms remain in `generalized_chebyshev`.

use crate::error::{MfsError, Result};
use num_complex::Complex64;

/// Complex scalar used by approximation helpers.
pub type ComplexCoefficient = Complex64;

const COMPLEX_ZERO: ComplexCoefficient = ComplexCoefficient::new(0.0, 0.0);
const COMPLEX_ONE: ComplexCoefficient = ComplexCoefficient::new(1.0, 0.0);

/// Strategy interface for complex-polynomial root finding.
pub trait ComplexRootSolver {
    /// Estimates all roots of the given complex polynomial.
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>>;
}

/// Default root finder used by the current approximation helpers.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct DurandKernerRootSolver;

/// Dense polynomial with complex coefficients in ascending-power order.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexPolynomial {
    /// Coefficients ordered from constant term upward.
    pub coefficients: Vec<ComplexCoefficient>,
}

impl ComplexPolynomial {
    /// Creates a validated complex polynomial.
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

    /// Evaluates the polynomial at a complex point using Horner's rule.
    pub fn evaluate(&self, x: ComplexCoefficient) -> ComplexCoefficient {
        self.coefficients
            .iter()
            .rev()
            .copied()
            .fold(COMPLEX_ZERO, |acc, coeff| acc * x + coeff)
    }

    /// Returns the polynomial degree.
    pub fn degree(&self) -> usize {
        self.coefficients.len().saturating_sub(1)
    }

    /// Multiplies every coefficient by the given scalar.
    pub fn scale(&self, scalar: ComplexCoefficient) -> Result<Self> {
        Self::new(
            self.coefficients
                .iter()
                .copied()
                .map(|coefficient| coefficient * scalar)
                .collect(),
        )
    }

    /// Adds two polynomials, padding the shorter one with implicit zeros.
    pub fn add(&self, rhs: &Self) -> Result<Self> {
        let target_len = self.coefficients.len().max(rhs.coefficients.len());
        let mut coefficients = vec![COMPLEX_ZERO; target_len];

        for (index, coefficient) in self.coefficients.iter().copied().enumerate() {
            coefficients[index] += coefficient;
        }
        for (index, coefficient) in rhs.coefficients.iter().copied().enumerate() {
            coefficients[index] += coefficient;
        }

        Self::new(trim_trailing_complex_zeros(coefficients))
    }

    /// Subtracts another polynomial, padding the shorter one with implicit zeros.
    pub fn sub(&self, rhs: &Self) -> Result<Self> {
        let target_len = self.coefficients.len().max(rhs.coefficients.len());
        let mut coefficients = vec![COMPLEX_ZERO; target_len];

        for (index, coefficient) in self.coefficients.iter().copied().enumerate() {
            coefficients[index] += coefficient;
        }
        for (index, coefficient) in rhs.coefficients.iter().copied().enumerate() {
            coefficients[index] += -coefficient;
        }

        Self::new(trim_trailing_complex_zeros(coefficients))
    }

    /// Returns the formal derivative of the polynomial.
    pub fn derivative(&self) -> Result<Self> {
        if self.coefficients.len() == 1 {
            return Self::new(vec![COMPLEX_ZERO]);
        }

        Self::new(
            self.coefficients
                .iter()
                .copied()
                .enumerate()
                .skip(1)
                .map(|(power, coefficient)| coefficient * power as f64)
                .collect(),
        )
    }

    /// Applies coefficient conjugation with alternating signs, equivalent to `Q(-s)^*`.
    pub fn alternating_conjugate(&self) -> Result<Self> {
        Self::new(
            self.coefficients
                .iter()
                .copied()
                .enumerate()
                .map(|(power, coefficient)| {
                    let sign = if power % 2 == 0 { 1.0 } else { -1.0 };
                    coefficient.conj() * sign
                })
                .collect(),
        )
    }

    /// Returns the leading non-zero coefficient in ascending-power storage.
    pub fn leading_coefficient(&self) -> ComplexCoefficient {
        self.coefficients.last().copied().unwrap_or(COMPLEX_ZERO)
    }

    /// Builds a monic polynomial whose roots are all real.
    pub fn from_real_roots(roots: &[f64]) -> Result<Self> {
        let mut coefficients = vec![COMPLEX_ONE];
        for &root in roots {
            coefficients = multiply_by_monic_root(&coefficients, ComplexCoefficient::new(root, 0.0));
        }
        Self::new(coefficients)
    }

    /// Builds a monic polynomial whose roots may be complex.
    pub fn from_complex_roots(roots: &[ComplexCoefficient]) -> Result<Self> {
        let mut coefficients = vec![COMPLEX_ONE];
        for &root in roots {
            coefficients = multiply_by_monic_root(&coefficients, root);
        }
        Self::new(coefficients)
    }

    /// Estimates all roots with the default root solver.
    pub fn roots(&self) -> Result<Vec<ComplexCoefficient>> {
        DurandKernerRootSolver.roots_of(self)
    }

    /// Estimates all roots with an explicit complex root solver.
    pub fn roots_with<S: ComplexRootSolver>(
        &self,
        solver: &S,
    ) -> Result<Vec<ComplexCoefficient>> {
        solver.roots_of(self)
    }
}

impl ComplexRootSolver for DurandKernerRootSolver {
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>> {
        let degree = polynomial.degree();
        if degree == 0 {
            return Ok(Vec::new());
        }

        let leading = *polynomial.coefficients.last().ok_or_else(|| {
            MfsError::Unsupported("polynomial is missing a leading coefficient".to_string())
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::Unsupported(
                "polynomial leading coefficient must be non-zero".to_string(),
            ));
        }

        let normalized = polynomial
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
                let mut denominator = COMPLEX_ONE;
                for (other_index, other_root) in roots.iter().copied().enumerate() {
                    if index != other_index {
                        denominator *= root - other_root;
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

pub(crate) fn multiply_by_monic_root(
    coefficients: &[ComplexCoefficient],
    root: ComplexCoefficient,
) -> Vec<ComplexCoefficient> {
    let mut next = vec![COMPLEX_ZERO; coefficients.len() + 1];
    for (index, coefficient) in coefficients.iter().copied().enumerate() {
        next[index] += coefficient * (-root);
        next[index + 1] += coefficient;
    }
    next
}

fn trim_trailing_complex_zeros(mut coefficients: Vec<ComplexCoefficient>) -> Vec<ComplexCoefficient> {
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
        .fold(COMPLEX_ZERO, |acc, coefficient| acc * x + coefficient)
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
    fn explicit_root_solver_matches_default_roots_wrapper() -> Result<()> {
        let polynomial = ComplexPolynomial::from_real_roots(&[1.0, 2.0, 3.0])?;
        let solver = DurandKernerRootSolver;

        let mut via_wrapper = polynomial.roots()?;
        let mut via_solver = polynomial.roots_with(&solver)?;

        let by_real_then_imag = |lhs: &ComplexCoefficient, rhs: &ComplexCoefficient| {
            lhs.re
                .partial_cmp(&rhs.re)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    lhs.im
                        .partial_cmp(&rhs.im)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        };
        via_wrapper.sort_by(by_real_then_imag);
        via_solver.sort_by(by_real_then_imag);

        for (wrapper_root, solver_root) in via_wrapper.iter().zip(via_solver.iter()) {
            approx_eq(wrapper_root.re, solver_root.re, 1e-10);
            approx_eq(wrapper_root.im, solver_root.im, 1e-10);
        }
        Ok(())
    }
}
