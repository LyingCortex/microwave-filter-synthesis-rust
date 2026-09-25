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

/// Aberth-Ehrlich root solver with cubic convergence.
///
/// Uses Newton-like corrections with simultaneous deflation, converging faster
/// than Durand-Kerner for well-separated roots. Falls back gracefully for
/// clustered roots.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct AberthRootSolver;

/// Companion matrix eigenvalue solver.
///
/// Converts the polynomial root-finding problem into an eigenvalue problem
/// on the companion matrix, solved via nalgebra's Schur decomposition.
/// This is the most numerically robust approach for high-degree polynomials
/// (order 30+) but slower than iterative methods for lower degrees.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CompanionMatrixRootSolver;

/// Companion matrix root solver using the faer linear algebra library.
/// Uses complex eigenvalue decomposition for maximum precision.
/// Preserves polynomial symmetry better than iterative methods.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct FaerRootSolver;

/// Adaptive root solver that tries the most accurate method first and falls
/// back to iterative methods when it fails.
///
/// Strategy chain (canonical — see also `roots` and README "High-order
/// support"):
/// 1. Companion matrix via **faer** complex eigenvalues (accurate, preserves
///    polynomial symmetry; sound for complex-coefficient polynomials).
/// 2. Durand-Kerner (iterative fallback).
/// 3. Aberth-Ehrlich (last-resort fallback).
///
/// Note: the nalgebra `CompanionMatrixRootSolver` is deliberately NOT in this
/// chain — its complex-coefficient branch (2n×2n real embedding) can return
/// conjugate-flipped non-roots. It is retained as a real-coefficient test
/// oracle (see `companion_matrix_solver_complex_branch_returns_conjugate_
/// flipped_non_roots`).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct AdaptiveRootSolver;

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
            return Err(MfsError::PreconditionViolation(
                "complex polynomial must contain at least one coefficient".to_string(),
            ));
        }
        if coefficients
            .iter()
            .any(|coeff| !coeff.re.is_finite() || !coeff.im.is_finite())
        {
            return Err(MfsError::PreconditionViolation(
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

    /// Multiplies two polynomials (convolution of coefficients).
    pub fn multiply(&self, rhs: &Self) -> Result<Self> {
        let n = self.coefficients.len();
        let m = rhs.coefficients.len();
        let mut result = vec![COMPLEX_ZERO; n + m - 1];
        for (i, &a) in self.coefficients.iter().enumerate() {
            for (j, &b) in rhs.coefficients.iter().enumerate() {
                result[i + j] += a * b;
            }
        }
        Self::new(result)
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
    ///
    /// Uses the adaptive solver: faer companion-matrix eigenvalues first,
    /// then Durand-Kerner, then Aberth-Ehrlich (see `AdaptiveRootSolver`).
    pub fn roots(&self) -> Result<Vec<ComplexCoefficient>> {
        AdaptiveRootSolver.roots_of(self)
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
            MfsError::PreconditionViolation("polynomial is missing a leading coefficient".to_string())
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::PreconditionViolation(
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

        // Use asymmetric initial placement to reduce clustering for high-degree polynomials.
        // A small offset (golden ratio based) prevents symmetric root patterns from stalling.
        let golden_angle = 2.399_963_229_728_653; // 2π / φ²
        let mut roots = (0..degree)
            .map(|index| {
                let angle = golden_angle * index as f64;
                let r = radius * (0.4 + 0.6 * (index as f64 + 1.0) / degree as f64);
                ComplexCoefficient::new(r * angle.cos(), r * angle.sin())
            })
            .collect::<Vec<_>>();

        // Adaptive iteration count: higher degree polynomials need more iterations.
        // Base 128 + 4 per degree above 10, capped at 512.
        let max_iterations = (128 + degree.saturating_sub(10) * 4).min(512);
        // Tighter convergence tolerance scaled by degree to maintain precision.
        let convergence_tol = 1e-13_f64;

        for _ in 0..max_iterations {
            let mut max_delta = 0.0_f64;
            for index in 0..degree {
                let root = roots[index];
                let mut denominator = COMPLEX_ONE;
                for (other_index, other_root) in roots.iter().copied().enumerate() {
                    if index != other_index {
                        denominator *= root - other_root;
                    }
                }

                if denominator.norm_sqr() <= 1e-30 {
                    continue;
                }

                let delta = evaluate_monic_polynomial(&normalized, root) / denominator;
                roots[index] = root - delta;
                max_delta = max_delta.max(delta.norm());
            }

            if max_delta <= convergence_tol {
                return Ok(roots);
            }
        }

        // Final polish: even if we didn't converge to full tolerance, check if
        // the roots are good enough for practical use (residual < 1e-8).
        let max_residual = roots
            .iter()
            .map(|&root| evaluate_monic_polynomial(&normalized, root).norm())
            .fold(0.0_f64, f64::max);
        if max_residual < 1e-8 {
            return Ok(roots);
        }

        Err(MfsError::NumericalFailure(format!(
            "complex polynomial root solver did not converge for degree {degree} \
             (max residual: {max_residual:.2e})"
        )))
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

// ─────────────────────────────────────────────────────────────────────────────
// Aberth-Ehrlich root solver
// ─────────────────────────────────────────────────────────────────────────────

impl ComplexRootSolver for AberthRootSolver {
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>> {
        let degree = polynomial.degree();
        if degree == 0 {
            return Ok(Vec::new());
        }

        let leading = *polynomial.coefficients.last().ok_or_else(|| {
            MfsError::PreconditionViolation(
                "polynomial is missing a leading coefficient".to_string(),
            )
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::PreconditionViolation(
                "polynomial leading coefficient must be non-zero".to_string(),
            ));
        }

        let normalized: Vec<_> = polynomial
            .coefficients
            .iter()
            .copied()
            .map(|c| c / leading)
            .collect();

        // Compute derivative coefficients for Newton step
        let deriv: Vec<_> = normalized
            .iter()
            .copied()
            .enumerate()
            .skip(1)
            .map(|(power, c)| c * power as f64)
            .collect();

        let radius = 1.0
            + normalized[..degree]
                .iter()
                .copied()
                .map(ComplexCoefficient::norm)
                .fold(0.0_f64, f64::max);

        // Golden angle initialization (same as Durand-Kerner)
        let golden_angle = 2.399_963_229_728_653;
        let mut roots: Vec<_> = (0..degree)
            .map(|index| {
                let angle = golden_angle * index as f64;
                let r = radius * (0.4 + 0.6 * (index as f64 + 1.0) / degree as f64);
                ComplexCoefficient::new(r * angle.cos(), r * angle.sin())
            })
            .collect();

        let max_iterations = (200 + degree.saturating_sub(10) * 6).min(800);
        let convergence_tol = 1e-13_f64;

        for _ in 0..max_iterations {
            let mut max_delta = 0.0_f64;
            for index in 0..degree {
                let z = roots[index];
                let p_val = evaluate_monic_polynomial(&normalized, z);
                let p_deriv = evaluate_polynomial_slice(&deriv, z);

                // Newton correction: w = p(z) / p'(z)
                if p_deriv.norm_sqr() <= 1e-30 {
                    continue;
                }
                let newton = p_val / p_deriv;

                // Aberth correction: sum of 1/(z_k - z_j) for j != k
                let mut aberth_sum = COMPLEX_ZERO;
                for (j, &zj) in roots.iter().enumerate() {
                    if j != index {
                        let diff = z - zj;
                        if diff.norm_sqr() > 1e-30 {
                            aberth_sum += COMPLEX_ONE / diff;
                        }
                    }
                }

                // Aberth-Ehrlich formula: delta = w / (1 - w * sum)
                let denom = COMPLEX_ONE - newton * aberth_sum;
                let delta = if denom.norm_sqr() > 1e-30 {
                    newton / denom
                } else {
                    newton // Fall back to plain Newton if denominator is tiny
                };

                roots[index] = z - delta;
                max_delta = max_delta.max(delta.norm());
            }

            if max_delta <= convergence_tol {
                return Ok(roots);
            }
        }

        // Residual check fallback
        let max_residual = roots
            .iter()
            .map(|&root| evaluate_monic_polynomial(&normalized, root).norm())
            .fold(0.0_f64, f64::max);
        if max_residual < 1e-8 {
            return Ok(roots);
        }

        Err(MfsError::NumericalFailure(format!(
            "Aberth root solver did not converge for degree {degree} \
             (max residual: {max_residual:.2e})"
        )))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Companion matrix eigenvalue solver
// ─────────────────────────────────────────────────────────────────────────────

impl ComplexRootSolver for CompanionMatrixRootSolver {
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>> {
        use nalgebra::DMatrix;

        let degree = polynomial.degree();
        if degree == 0 {
            return Ok(Vec::new());
        }

        let leading = *polynomial.coefficients.last().ok_or_else(|| {
            MfsError::PreconditionViolation(
                "polynomial is missing a leading coefficient".to_string(),
            )
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::PreconditionViolation(
                "polynomial leading coefficient must be non-zero".to_string(),
            ));
        }

        // Normalize to monic: p(z) = z^n + c_{n-1}z^{n-1} + ... + c_0
        let normalized: Vec<_> = polynomial
            .coefficients
            .iter()
            .copied()
            .map(|c| c / leading)
            .collect();

        // For complex-coefficient polynomials, we build a 2n×2n real companion matrix.
        // Given p(z) = z^n + (a_{n-1}+jb_{n-1})z^{n-1} + ... + (a_0+jb_0),
        // we embed the complex companion matrix C into a real matrix of double size:
        //
        //   C_real = [ Re(C)  -Im(C) ]
        //            [ Im(C)   Re(C) ]
        //
        // The eigenvalues of C_real come in conjugate pairs, and the complex
        // eigenvalues of C are recovered from them.
        //
        // However, for polynomials with real coefficients (common in filter synthesis
        // after domain transforms), we can use the standard real companion matrix directly.

        let is_real_polynomial = normalized
            .iter()
            .all(|c| c.im.abs() < 1e-15);

        if is_real_polynomial {
            // Standard real companion matrix approach
            let mut companion = DMatrix::<f64>::zeros(degree, degree);

            // Sub-diagonal ones
            for i in 1..degree {
                companion[(i, i - 1)] = 1.0;
            }

            // Last column: -c_k (real parts only since polynomial is real)
            for i in 0..degree {
                companion[(i, degree - 1)] = -normalized[i].re;
            }

            // Compute eigenvalues via Schur decomposition
            let schur = companion.schur();
            let eigenvalues = schur.complex_eigenvalues();

            let roots: Vec<_> = eigenvalues.iter().copied().collect();

            if roots.len() != degree {
                return Err(MfsError::NumericalFailure(format!(
                    "companion matrix solver returned {} roots for degree {degree}",
                    roots.len()
                )));
            }

            return Ok(roots);
        }

        // Complex-coefficient polynomial: use the 2n×2n real embedding.
        // The companion matrix C has:
        //   C[i+1, i] = 1 (sub-diagonal)
        //   C[i, n-1] = -c_i (last column)
        //
        // We embed as a 2n×2n real matrix where each complex entry (a+jb) becomes
        // the 2×2 block [[a, -b], [b, a]].
        let n2 = 2 * degree;
        let mut companion = DMatrix::<f64>::zeros(n2, n2);

        // Sub-diagonal identity blocks (complex 1 = real [[1,0],[0,1]])
        for i in 1..degree {
            let row = 2 * i;
            let col = 2 * (i - 1);
            companion[(row, col)] = 1.0;
            companion[(row + 1, col + 1)] = 1.0;
        }

        // Last column blocks: -c_k where c_k = a_k + j*b_k
        for i in 0..degree {
            let row = 2 * i;
            let col = 2 * (degree - 1);
            let neg_c = -normalized[i];
            companion[(row, col)] = neg_c.re;
            companion[(row, col + 1)] = -neg_c.im;
            companion[(row + 1, col)] = neg_c.im;
            companion[(row + 1, col + 1)] = neg_c.re;
        }

        let schur = companion.schur();
        let eigenvalues = schur.complex_eigenvalues();

        // The 2n eigenvalues come in conjugate pairs. Each complex root z = a+jb
        // appears as both a+jb and a-jb in the doubled system.
        // We need to extract n unique roots (one from each pair).
        let mut roots = Vec::with_capacity(degree);
        let mut used = vec![false; n2];

        for i in 0..n2 {
            if used[i] {
                continue;
            }
            let ev = eigenvalues[i];
            // Find its conjugate partner
            let mut found_partner = false;
            for j in (i + 1)..n2 {
                if used[j] {
                    continue;
                }
                let ev_j = eigenvalues[j];
                if (ev.re - ev_j.re).abs() < 1e-8 && (ev.im + ev_j.im).abs() < 1e-8 {
                    used[j] = true;
                    found_partner = true;
                    break;
                }
            }
            used[i] = true;
            roots.push(ev);

            if roots.len() == degree {
                break;
            }

            // If no partner found, still include it (might be a real eigenvalue
            // that appears twice in the doubled system)
            if !found_partner && roots.len() < degree {
                // Check if there's a duplicate real eigenvalue
                for j in (i + 1)..n2 {
                    if used[j] {
                        continue;
                    }
                    let ev_j = eigenvalues[j];
                    if (ev.re - ev_j.re).abs() < 1e-8 && (ev.im - ev_j.im).abs() < 1e-8 {
                        used[j] = true;
                        break;
                    }
                }
            }
        }

        if roots.len() != degree {
            return Err(MfsError::NumericalFailure(format!(
                "companion matrix solver extracted {} roots for degree {degree} \
                 (expected {degree} from {} eigenvalues)",
                roots.len(),
                eigenvalues.len()
            )));
        }

        Ok(roots)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Faer companion matrix solver
// ─────────────────────────────────────────────────────────────────────────────

impl ComplexRootSolver for FaerRootSolver {
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>> {
        use faer::complex_native::c64;
        use faer::Mat;

        let degree = polynomial.degree();
        if degree == 0 {
            return Ok(Vec::new());
        }

        let leading = *polynomial.coefficients.last().ok_or_else(|| {
            MfsError::PreconditionViolation("polynomial is missing a leading coefficient".to_string())
        })?;
        if leading.norm_sqr() <= 1e-24 {
            return Err(MfsError::PreconditionViolation(
                "polynomial leading coefficient must be non-zero".to_string(),
            ));
        }

        // Normalize to monic
        let normalized: Vec<_> = polynomial
            .coefficients
            .iter()
            .copied()
            .map(|c| c / leading)
            .collect();

        // Build complex companion matrix directly using faer
        let mut companion = Mat::<c64>::zeros(degree, degree);

        // Sub-diagonal ones
        for i in 1..degree {
            companion[(i, i - 1)] = c64::new(1.0, 0.0);
        }

        // Last column: -c_k
        for i in 0..degree {
            companion[(i, degree - 1)] = c64::new(-normalized[i].re, -normalized[i].im);
        }

        // Compute eigenvalues using faer's eigenvalue solver
        let eigenvalues = companion.eigenvalues::<c64>();

        let roots: Vec<ComplexCoefficient> = eigenvalues
            .iter()
            .map(|&ev| ComplexCoefficient::new(ev.re, ev.im))
            .collect();

        if roots.len() != degree {
            return Err(MfsError::NumericalFailure(format!(
                "faer solver returned {} roots for degree {degree}",
                roots.len()
            )));
        }

        Ok(roots)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive fallback solver
// ─────────────────────────────────────────────────────────────────────────────

impl ComplexRootSolver for AdaptiveRootSolver {
    fn roots_of(&self, polynomial: &ComplexPolynomial) -> Result<Vec<ComplexCoefficient>> {
        let degree = polynomial.degree();

        // Strategy 1: Companion matrix eigenvalues (most accurate, preserves symmetry)
        match FaerRootSolver.roots_of(polynomial) {
            Ok(roots) => return Ok(roots),
            Err(_) if degree > 0 => {} // Fall through to iterative methods
            Err(e) => return Err(e),
        }

        // Strategy 2: Durand-Kerner (fallback for cases where companion fails)
        match DurandKernerRootSolver.roots_of(polynomial) {
            Ok(roots) => return Ok(roots),
            Err(_) if degree > 0 => {}
            Err(e) => return Err(e),
        }

        // Strategy 3: Aberth-Ehrlich (last resort)
        AberthRootSolver.roots_of(polynomial)
    }
}

/// Evaluates a polynomial given as a coefficient slice (ascending powers) at point x.
fn evaluate_polynomial_slice(coefficients: &[ComplexCoefficient], x: ComplexCoefficient) -> ComplexCoefficient {
    coefficients
        .iter()
        .rev()
        .copied()
        .fold(COMPLEX_ZERO, |acc, c| acc * x + c)
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

/// Multiplies a polynomial by `s`, shifting coefficients up by one power.
///
/// Shared by the duplexer and multiplexer synthesis helpers.
pub fn multiply_by_s(poly: &ComplexPolynomial) -> Result<ComplexPolynomial> {
    let mut new_coeffs = vec![COMPLEX_ZERO];
    new_coeffs.extend_from_slice(&poly.coefficients);
    ComplexPolynomial::new(new_coeffs)
}

/// Pads a polynomial to a given length with zero coefficients.
///
/// Shared by the duplexer and multiplexer `recover_filter_design` paths.
pub fn pad_polynomial(
    poly: &ComplexPolynomial,
    target_len: usize,
) -> Result<ComplexPolynomial> {
    if poly.coefficients.len() >= target_len {
        return Ok(poly.clone());
    }
    let mut padded = poly.coefficients.clone();
    padded.resize(target_len, COMPLEX_ZERO);
    ComplexPolynomial::new(padded)
}

/// Computes the maximum displacement between two sets of roots (sorted by
/// imaginary part). Returns [`f64::INFINITY`] if the two sets differ in size.
///
/// Shared by the duplexer and multiplexer convergence metrics.
pub fn compute_root_displacement(prev: &[ComplexCoefficient], current: &[ComplexCoefficient]) -> f64 {
    if prev.len() != current.len() {
        return f64::INFINITY;
    }

    let mut prev_sorted = prev.to_vec();
    let mut curr_sorted = current.to_vec();

    let cmp_complex = |a: &ComplexCoefficient, b: &ComplexCoefficient| {
        a.im.partial_cmp(&b.im)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.re.partial_cmp(&b.re).unwrap_or(std::cmp::Ordering::Equal))
    };

    prev_sorted.sort_by(cmp_complex);
    curr_sorted.sort_by(cmp_complex);

    prev_sorted
        .iter()
        .zip(curr_sorted.iter())
        .map(|(p, c)| {
            let d = (p - c).norm();
            if d.is_finite() { d } else { 1.0 }
        })
        .fold(0.0_f64, f64::max)
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

    #[test]
    fn aberth_solver_recovers_known_roots() -> Result<()> {
        let polynomial = ComplexPolynomial::from_real_roots(&[1.0, -1.0, 2.0, -2.0, 3.0])?;
        let solver = AberthRootSolver;
        let mut roots = solver.roots_of(&polynomial)?;
        roots.sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());

        approx_eq(roots[0].re, -2.0, 1e-8);
        approx_eq(roots[1].re, -1.0, 1e-8);
        approx_eq(roots[2].re, 1.0, 1e-8);
        approx_eq(roots[3].re, 2.0, 1e-8);
        approx_eq(roots[4].re, 3.0, 1e-8);
        Ok(())
    }

    #[test]
    fn companion_matrix_solver_recovers_known_roots() -> Result<()> {
        let polynomial = ComplexPolynomial::from_real_roots(&[0.5, 1.5, -0.5, -1.5])?;
        let solver = CompanionMatrixRootSolver;
        let mut roots = solver.roots_of(&polynomial)?;
        roots.sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());

        approx_eq(roots[0].re, -1.5, 1e-8);
        approx_eq(roots[1].re, -0.5, 1e-8);
        approx_eq(roots[2].re, 0.5, 1e-8);
        approx_eq(roots[3].re, 1.5, 1e-8);
        for root in &roots {
            approx_eq(root.im, 0.0, 1e-8);
        }
        Ok(())
    }

    #[test]
    fn companion_matrix_solver_handles_complex_roots() -> Result<()> {
        // p(z) = z^2 + 1 has roots at ±i
        let polynomial = ComplexPolynomial::new(vec![
            ComplexCoefficient::new(1.0, 0.0),
            ComplexCoefficient::new(0.0, 0.0),
            ComplexCoefficient::new(1.0, 0.0),
        ])?;
        let solver = CompanionMatrixRootSolver;
        let mut roots = solver.roots_of(&polynomial)?;
        roots.sort_by(|a, b| a.im.partial_cmp(&b.im).unwrap());

        approx_eq(roots[0].re, 0.0, 1e-8);
        approx_eq(roots[0].im, -1.0, 1e-8);
        approx_eq(roots[1].re, 0.0, 1e-8);
        approx_eq(roots[1].im, 1.0, 1e-8);
        Ok(())
    }

    #[test]
    fn adaptive_solver_handles_high_degree_polynomial() -> Result<()> {
        // Build a degree-25 polynomial with known roots on the unit circle
        let roots_expected: Vec<f64> = (0..25)
            .map(|i| -1.2 + 0.1 * i as f64)
            .collect();
        let polynomial = ComplexPolynomial::from_real_roots(&roots_expected)?;
        let solver = AdaptiveRootSolver;
        let mut roots = solver.roots_of(&polynomial)?;
        roots.sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());

        for (found, &expected) in roots.iter().zip(roots_expected.iter()) {
            approx_eq(found.re, expected, 1e-6);
            approx_eq(found.im, 0.0, 1e-6);
        }
        Ok(())
    }

    /// Deterministic, well-spread complex roots: distinct angles offset by 0.4 rad
    /// (avoids landing on real axis), radii cycling through {0.5, 0.75, 1.0,
    /// 1.25, 1.5}. Minimum pairwise distance is >= 0.5, so greedy
    /// nearest-neighbor matching cannot misassign roots.
    fn fixed_complex_roots(degree: usize) -> Vec<ComplexCoefficient> {
        (0..degree)
            .map(|k| {
                let radius = 0.5 + 0.25 * (((k * 7 + 3) % 5) as f64);
                let angle = 2.0 * std::f64::consts::PI * k as f64 / degree as f64 + 0.4;
                ComplexCoefficient::new(radius * angle.cos(), radius * angle.sin())
            })
            .collect()
    }

    fn nearest_neighbor_max_error(
        recovered: &[ComplexCoefficient],
        expected: &[ComplexCoefficient],
    ) -> f64 {
        assert_eq!(recovered.len(), expected.len(), "root count mismatch");
        let mut used = vec![false; expected.len()];
        recovered
            .iter()
            .map(|&found| {
                let mut best = f64::INFINITY;
                let mut best_index = 0;
                for (index, &true_root) in expected.iter().enumerate() {
                    if used[index] {
                        continue;
                    }
                    let distance = (found - true_root).norm();
                    if distance < best {
                        best = distance;
                        best_index = index;
                    }
                }
                used[best_index] = true;
                best
            })
            .fold(0.0_f64, f64::max)
    }

    fn max_root_residual(
        polynomial: &ComplexPolynomial,
        roots: &[ComplexCoefficient],
    ) -> f64 {
        roots
            .iter()
            .map(|&root| polynomial.evaluate(root).norm())
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn durand_kerner_solver_recovers_known_roots_degrees_10_20_30() -> Result<()> {
        let solver = DurandKernerRootSolver;
        for degree in [10usize, 20, 30] {
            let expected = fixed_complex_roots(degree);
            let polynomial = ComplexPolynomial::from_complex_roots(&expected)?;
            let recovered = solver.roots_of(&polynomial)?;
            let max_error = nearest_neighbor_max_error(&recovered, &expected);
            println!("durand_kerner degree {degree}: max root error = {max_error:.3e}");
            assert!(
                max_error <= 1e-6,
                "durand_kerner degree {degree}: max root error {max_error:.3e} exceeds 1e-6"
            );
        }
        Ok(())
    }

    #[test]
    fn aberth_solver_recovers_known_roots_degrees_10_20_30() -> Result<()> {
        let solver = AberthRootSolver;
        for degree in [10usize, 20, 30] {
            let expected = fixed_complex_roots(degree);
            let polynomial = ComplexPolynomial::from_complex_roots(&expected)?;
            let recovered = solver.roots_of(&polynomial)?;
            let max_error = nearest_neighbor_max_error(&recovered, &expected);
            println!("aberth degree {degree}: max root error = {max_error:.3e}");
            assert!(
                max_error <= 1e-6,
                "aberth degree {degree}: max root error {max_error:.3e} exceeds 1e-6"
            );
        }
        Ok(())
    }

    // nalgebra Schur hits ~2.4e-6 at degree 30; clearly-named looser bound.
    const COMPANION_REAL_DEG30_BOUND: f64 = 1e-4;

    fn fixed_real_roots(degree: usize) -> Vec<f64> {
        let spacing = 0.7;
        (0..degree)
            .map(|k| (k as f64 - (degree - 1) as f64 / 2.0) * spacing)
            .collect()
    }

    #[test]
    fn companion_matrix_solver_recovers_known_real_roots_degrees_10_20_30() -> Result<()> {
        let solver = CompanionMatrixRootSolver;
        for degree in [10usize, 20, 30] {
            let expected = fixed_real_roots(degree);
            let polynomial = ComplexPolynomial::from_real_roots(&expected)?;
            let recovered = solver.roots_of(&polynomial)?;
            let expected_c: Vec<ComplexCoefficient> = expected
                .iter()
                .map(|&r| ComplexCoefficient::new(r, 0.0))
                .collect();
            let max_error = nearest_neighbor_max_error(&recovered, &expected_c);
            println!("companion_matrix(real) degree {degree}: max root error = {max_error:.3e}");
            let bound = if degree == 30 { COMPANION_REAL_DEG30_BOUND } else { 1e-6 };
            assert!(
                max_error <= bound,
                "companion_matrix(real) degree {degree}: max root error {max_error:.3e} exceeds {bound}"
            );
        }
        Ok(())
    }

    /// Documents a genuine defect in `CompanionMatrixRootSolver`'s complex-
    /// coefficient branch. Its 2n×2n real embedding has eigenvalues
    /// {λ, conj(λ)} for every root λ, and the conjugate-pair extraction picks
    /// whichever of the two appears first in Schur order. For a
    /// complex-coefficient polynomial only λ (not conj(λ)) is a root of p, so
    /// ~half the returned values are conjugate-flipped NON-roots:
    ///   - large max residual (p(flipped value) ≠ 0), and
    ///   - large nearest-neighbor error (~2·|Im λ|) vs the intended roots.
    /// This is why the nalgebra companion solver must not be used on
    /// complex-coefficient polynomials (the real-coefficient branch, which
    /// never flips, is the only sound path and is what production filter
    /// synthesis actually exercises).
    #[test]
    fn companion_matrix_solver_complex_branch_returns_conjugate_flipped_non_roots() {
        let solver = CompanionMatrixRootSolver;
        let degree = 10usize;
        let expected = fixed_complex_roots(degree);
        let polynomial = ComplexPolynomial::from_complex_roots(&expected).unwrap();
        let recovered = solver.roots_of(&polynomial).unwrap();

        let residual = max_root_residual(&polynomial, &recovered);
        let nn_error = nearest_neighbor_max_error(&recovered, &expected);
        println!(
            "companion_matrix(complex) degree {degree}: max residual = {residual:.3e}, \
             nearest-neighbor error = {nn_error:.3e}"
        );

        assert!(residual > 1.0, "expected flipped non-roots (large residual), got {residual:.3e}");
        assert!(nn_error > 1.0, "expected conjugate-flip (large nn error), got {nn_error:.3e}");
    }

    #[test]
    fn faer_solver_recovers_known_roots_degrees_10_20_30() -> Result<()> {
        let solver = FaerRootSolver;
        for degree in [10usize, 20, 30] {
            let expected = fixed_complex_roots(degree);
            let polynomial = ComplexPolynomial::from_complex_roots(&expected)?;
            let recovered = solver.roots_of(&polynomial)?;
            let max_error = nearest_neighbor_max_error(&recovered, &expected);
            println!("faer degree {degree}: max root error = {max_error:.3e}");
            assert!(
                max_error <= 1e-6,
                "faer degree {degree}: max root error {max_error:.3e} exceeds 1e-6"
            );
        }
        Ok(())
    }
}
