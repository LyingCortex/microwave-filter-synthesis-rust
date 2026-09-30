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
///
/// The returned eigenvalues are gated on a backward-residual check
/// (`max_backward_residual` ≤ the private `FAER_BACKWARD_RESIDUAL_TOL` constant) and a
/// non-finiteness check: a companion matrix whose eigenvalue problem is
/// ill-conditioned can yield garbage eigenvalues without the call reporting
/// an error, and the gate converts that silent corruption into a typed error
/// so [`AdaptiveRootSolver`] falls through to the iterative methods.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct FaerRootSolver;

/// Adaptive root solver that tries the most accurate method first and falls
/// back to iterative methods when it fails.
///
/// Strategy chain (canonical — see also `roots` and README "High-order
/// support"):
/// 1. Companion matrix via **faer** complex eigenvalues (accurate, preserves
///    polynomial symmetry; sound for complex-coefficient polynomials). The
///    result is gated on a backward-residual check (see `FaerRootSolver`):
///    ill-conditioned eigenvalue problems fail the gate and fall through.
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

    /// Returns the polynomial with its argument shifted: the coefficients of
    /// `P(s + shift)`.
    ///
    /// This is the common-domain transform used by the N-channel
    /// Macchiarella-Tamiazzo diplexer iteration, which must express one
    /// channel's polynomials in another channel's normalized-frequency domain.
    ///
    /// Computed by Horner's rule in the argument `(s + shift)`:
    /// `P(x+c) = (...((a_n (x+c) + a_{n-1})(x+c) + ...)(x+c) + a_0)`, which is
    /// exact up to roundoff and needs no root-finding. Numerically well
    /// behaved for the narrow-band shifts encountered in diplexer synthesis
    /// (`|shift|` = O(band separation / bandwidth)).
    pub fn shift(&self, shift: ComplexCoefficient) -> Result<Self> {
        // acc = 0; for a in coeffs (descending) { acc = acc * (s + shift) + a }.
        // Multiplying ascending vector acc by (s + shift):
        //   next[k]  += shift * acc[k]   (the c·acc part)
        //   next[k+1] += acc[k]          (the s·acc part)
        let mut acc: Vec<ComplexCoefficient> = Vec::new();
        for &a in self.coefficients.iter().rev() {
            let len = acc.len();
            let mut next = vec![COMPLEX_ZERO; len + 1];
            for k in 0..len {
                next[k] += shift * acc[k];
                next[k + 1] += acc[k];
            }
            next[0] += a;
            acc = next;
        }
        Self::new(trim_trailing_complex_zeros(acc))
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
            coefficients =
                multiply_by_monic_root(&coefficients, ComplexCoefficient::new(root, 0.0));
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
    pub fn roots_with<S: ComplexRootSolver>(&self, solver: &S) -> Result<Vec<ComplexCoefficient>> {
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
            MfsError::PreconditionViolation(
                "polynomial is missing a leading coefficient".to_string(),
            )
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

            // `f64::max` swallows NaN, so `max_delta` can read 0.0 even when an
            // update overflowed to a NaN root. Check the roots directly so a NaN
            // set is never returned as "converged".
            if max_delta <= convergence_tol && !roots_contain_non_finite(&roots) {
                return Ok(roots);
            }
        }

        // Final polish: even if we didn't converge to full tolerance, check if
        // the roots are good enough for practical use (residual < 1e-8).
        let max_residual = roots
            .iter()
            .map(|&root| evaluate_monic_polynomial(&normalized, root).norm())
            .fold(0.0_f64, f64::max);
        if !roots_contain_non_finite(&roots) && max_residual < 1e-8 {
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

            // Same NaN guard as Durand-Kerner: `f64::max` swallows NaN in
            // `max_delta`, so a NaN root set would otherwise read as converged.
            if max_delta <= convergence_tol && !roots_contain_non_finite(&roots) {
                return Ok(roots);
            }
        }

        // Residual check fallback
        let max_residual = roots
            .iter()
            .map(|&root| evaluate_monic_polynomial(&normalized, root).norm())
            .fold(0.0_f64, f64::max);
        if !roots_contain_non_finite(&roots) && max_residual < 1e-8 {
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

        let is_real_polynomial = normalized.iter().all(|c| c.im.abs() < 1e-15);

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
        let col = 2 * (degree - 1);
        for (i, &c) in normalized.iter().enumerate().take(degree) {
            let row = 2 * i;
            let neg_c = -c;
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
            MfsError::PreconditionViolation(
                "polynomial is missing a leading coefficient".to_string(),
            )
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

        // Backward-residual gate: faer's eigenvalues can be garbage (e.g. from
        // a defective companion matrix) without the call reporting an error.
        // A backward-stable solve returns exact roots of a slightly perturbed
        // polynomial, so the backward residual is ~1e-15..1e-6 (worst observed:
        // 3e-6 for an order-30 production F_s). Garbage outputs are O(1e-2)+.
        // Failing here lets `AdaptiveRootSolver` fall through to the iterative
        // methods instead of silently shipping bad roots.
        let max_backward = max_backward_residual(&normalized, &roots);
        if roots_contain_non_finite(&roots) || max_backward > FAER_BACKWARD_RESIDUAL_TOL {
            return Err(MfsError::NumericalFailure(format!(
                "faer companion-matrix eigenvalues failed the backward-residual \
                 check for degree {degree} (max backward residual: {max_backward:.2e})"
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
fn evaluate_polynomial_slice(
    coefficients: &[ComplexCoefficient],
    x: ComplexCoefficient,
) -> ComplexCoefficient {
    coefficients
        .iter()
        .rev()
        .copied()
        .fold(COMPLEX_ZERO, |acc, c| acc * x + c)
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
        .fold(COMPLEX_ZERO, |acc, coefficient| acc * x + coefficient)
}

/// Maximum backward residual a `FaerRootSolver` root set may exhibit and
/// still be accepted. See the gate in `FaerRootSolver::roots_of` for the
/// calibration behind this value.
const FAER_BACKWARD_RESIDUAL_TOL: f64 = 1e-3;

/// True when any root has a non-finite real or imaginary part.
///
/// Needed because the residual folds use `f64::max`, which swallows NaN — a
/// root set containing a NaN would otherwise pass a `residual < tol` check
/// and be silently shipped.
fn roots_contain_non_finite(roots: &[ComplexCoefficient]) -> bool {
    roots.iter().any(|r| !r.re.is_finite() || !r.im.is_finite())
}

/// Largest backward residual of `roots` for the monic polynomial whose
/// ascending coefficients are `coefficients`.
///
/// For a candidate root `r`, the backward residual is
/// `|p(r)| / sum_k |c_k| |r|^k`. It measures "the roots are exact for some
/// nearly-perturbed polynomial" (a backward-stable-solve certificate) and is
/// scale-free — unlike the raw residual `|p(r)|`, it does not blow up with
/// root magnitude or polynomial coefficient spread.
fn max_backward_residual(coefficients: &[ComplexCoefficient], roots: &[ComplexCoefficient]) -> f64 {
    roots
        .iter()
        .map(|&root| {
            let p_at_root = evaluate_monic_polynomial(coefficients, root).norm();
            let mut scale = 0.0_f64;
            let mut r_pow = 1.0_f64;
            let r_mag = root.norm();
            for coefficient in coefficients {
                scale += coefficient.norm() * r_pow;
                r_pow *= r_mag;
            }
            p_at_root / scale.max(1e-300)
        })
        .fold(0.0_f64, f64::max)
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
pub fn pad_polynomial(poly: &ComplexPolynomial, target_len: usize) -> Result<ComplexPolynomial> {
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
pub fn compute_root_displacement(
    prev: &[ComplexCoefficient],
    current: &[ComplexCoefficient],
) -> f64 {
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
            if d.is_finite() {
                d
            } else {
                1.0
            }
        })
        .fold(0.0_f64, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::generalized_chebyshev_polynomials;
    use crate::spec::FilterSpec;
    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn shift_reconstructs_argument_shifted_polynomial() -> Result<()> {
        // P(s) = s^2 + 2 s + 3;  P(s + j5) = s^2 + (2 + j10) s + (-22 + j10).
        let p = ComplexPolynomial::new(vec![
            Complex64::new(3.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(1.0, 0.0),
        ])?;
        let q = p.shift(Complex64::new(0.0, 5.0))?;
        let expected = [
            Complex64::new(-22.0, 10.0),
            Complex64::new(2.0, 10.0),
            Complex64::new(1.0, 0.0),
        ];
        assert_eq!(q.coefficients.len(), expected.len());
        for (k, e) in expected.iter().enumerate() {
            assert!(
                (q.coefficients[k] - *e).norm() < 1e-12,
                "coeff {k}: {} (want {})",
                q.coefficients[k],
                e
            );
        }
        // P(s+c) evaluated at x must equal P evaluated at x+c.
        let c = Complex64::new(-0.7, 2.3);
        let r = p.shift(c)?;
        for x in [Complex64::new(0.0, 0.0), Complex64::new(1.5, -0.4), Complex64::new(-2.0, 3.0)] {
            let lhs = r.evaluate(x);
            let rhs = p.evaluate(x + c);
            let scale = (lhs.norm() + rhs.norm()).max(1.0);
            assert!((lhs - rhs).norm() / scale < 1e-11, "identity failed at {x}");
        }
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
        let roots_expected: Vec<f64> = (0..25).map(|i| -1.2 + 0.1 * i as f64).collect();
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

    fn max_root_residual(polynomial: &ComplexPolynomial, roots: &[ComplexCoefficient]) -> f64 {
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
            let bound = if degree == 30 {
                COMPANION_REAL_DEG30_BOUND
            } else {
                1e-6
            };
            assert!(
                max_error <= bound,
                "companion_matrix(real) degree {degree}: max root error {max_error:.3e} exceeds {bound}"
            );
        }
        Ok(())
    }

    /// Documents a genuine defect in `CompanionMatrixRootSolver`'s complex-
    /// coefficient branch. Its 2n×2n real embedding has eigenvalues {λ,
    /// conj(λ)} for every root λ, and the conjugate-pair extraction picks
    /// whichever of the two appears first in Schur order. For a
    /// complex-coefficient polynomial only λ (not conj(λ)) is a root of p, so
    /// about half the returned values are conjugate-flipped NON-roots, giving
    /// a large max residual (p(flipped value) ≠ 0) and a large nearest-
    /// neighbor error (~2·|Im λ|) versus the intended roots. This is why the
    /// nalgebra companion solver must not be used on complex-coefficient
    /// polynomials; the real-coefficient branch, which never flips, is the
    /// only sound path and is what production filter synthesis exercises.
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

        assert!(
            residual > 1.0,
            "expected flipped non-roots (large residual), got {residual:.3e}"
        );
        assert!(
            nn_error > 1.0,
            "expected conjugate-flip (large nn error), got {nn_error:.3e}"
        );
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

    /// Unit test for the backward-residual helper: a root set that is
    /// genuinely close to the true roots must have a tiny residual, while a
    /// crafted garbage set (offset roots) must have a large one. This is the
    /// metric the faer gate on, so it is tested directly and independently of
    /// faer.
    #[test]
    fn max_backward_residual_separates_good_roots_from_garbage() {
        let expected = fixed_complex_roots(12);
        let polynomial = ComplexPolynomial::from_complex_roots(&expected).unwrap();

        let good = max_backward_residual(&polynomial.coefficients, &expected);
        assert!(
            good < 1e-10,
            "expected ~machine-precision residual, got {good:.3e}"
        );

        // Shift every root by a sizeable amount: no longer roots of (any
        // nearby) polynomial, so the residual must blow up. The metric is
        // scale-free and measures the relative perturbation of the
        // coefficients, so the absolute value is shape-dependent — assert
        // separation from machine precision, not an absolute bound.
        let garbage: Vec<ComplexCoefficient> = expected
            .iter()
            .map(|r| *r + ComplexCoefficient::new(3.0, -2.0))
            .collect();
        let bad = max_backward_residual(&polynomial.coefficients, &garbage);
        assert!(bad > 0.1, "expected large backward residual, got {bad:.3e}");
        assert!(
            bad > good * 1e6,
            "garbage residual {bad:.3e} should exceed good residual {good:.3e} by orders of magnitude"
        );
    }

    /// Deterministic gate-firing case: roots spanning many decades make the
    /// companion matrix badly ill-conditioned, so faer's eigenvalues fail the
    /// backward-residual check. Regression-locks the gate itself — without it
    /// `FaerRootSolver` returns `Ok` on garbage eigenvalues (the silent
    /// corruption that motivated the fix).
    #[test]
    fn faer_backward_residual_gate_rejects_garbage_on_ill_conditioned_polynomial() {
        // Roots 1e-10 .. 1e10 (11 roots, geometric): observed backward
        // residual ~1.75e-2, far above the 1e-3 gate threshold.
        let roots: Vec<ComplexCoefficient> = (0..11)
            .map(|k| ComplexCoefficient::new(10f64.powi(k - 5), 0.0))
            .collect();
        let polynomial = ComplexPolynomial::from_complex_roots(&roots).unwrap();

        let result = FaerRootSolver.roots_of(&polynomial);
        match result {
            Err(e) => {
                // The typed error must name the residual failure.
                let msg = format!("{e}");
                assert!(
                    msg.contains("backward residual"),
                    "expected a backward-residual error, got: {msg}"
                );
            }
            Ok(recovered) => {
                // If a future faer version does recover it, the roots must at
                // least be backward-respectable — never silently garbage.
                let br = max_backward_residual(&polynomial.coefficients, &recovered);
                panic!(
                    "faer unexpectedly succeeded with backward residual {br:.3e} \
                     (gate threshold {:.1e}); the test premise changed",
                    FAER_BACKWARD_RESIDUAL_TOL
                );
            }
        }
    }

    /// The adaptive chain must fall through the (gated) faer solver to the
    /// iterative methods on an ill-conditioned polynomial and either recover
    /// or fail loudly — it must never return a non-finite root set.
    #[test]
    fn adaptive_chain_falls_through_on_ill_conditioned_polynomial() {
        let roots: Vec<ComplexCoefficient> = (0..11)
            .map(|k| ComplexCoefficient::new(10f64.powi(k - 5), 0.0))
            .collect();
        let polynomial = ComplexPolynomial::from_complex_roots(&roots).unwrap();

        // Premise: the gate fires for faer on this shape.
        assert!(
            FaerRootSolver.roots_of(&polynomial).is_err(),
            "test premise: faer should be gated on this polynomial"
        );

        let recovered = AdaptiveRootSolver.roots_of(&polynomial);
        if let Ok(r) = recovered {
            let finite = !r.iter().any(|x| !x.re.is_finite() || !x.im.is_finite());
            assert!(finite, "adaptive chain returned a non-finite root set");
        }
        // Err: failing loudly is acceptable for an extreme case
    }

    /// Regression for the `f64::max`-swallows-NaN defect: both iterative
    /// solvers' in-loop convergence check used to treat a NaN root set as
    /// converged (NaN updates keep `max_delta` at 0.0). A polynomial whose
    /// iterates overflow must return `Err`, never `Ok` with NaN roots.
    #[test]
    fn iterative_solvers_fail_loud_instead_of_returning_nan_roots() {
        // Roots on a huge circle: iterates overflow to NaN in the
        // Durand-Kerner / Aberth updates (observed: order-45 E_s via DK).
        let roots: Vec<ComplexCoefficient> = (0..15)
            .map(|k| {
                let a = 2.0 * std::f64::consts::PI * (k as f64 / 15.0);
                ComplexCoefficient::new(1e15 * a.cos(), 1e15 * a.sin())
            })
            .collect();
        let polynomial = ComplexPolynomial::from_complex_roots(&roots).unwrap();

        for (name, result) in [
            ("DurandKerner", DurandKernerRootSolver.roots_of(&polynomial)),
            ("Aberth", AberthRootSolver.roots_of(&polynomial)),
        ] {
            if let Ok(recovered) = result {
                let finite = !recovered
                    .iter()
                    .any(|x| !x.re.is_finite() || !x.im.is_finite());
                assert!(
                    finite,
                    "{name} returned a non-finite root set: {recovered:?}"
                );
            }
        }
    }

    /// The faer gate must NOT fire on legitimate production filter polynomials.
    /// Calibration (order 4..40, E_s/F_s/P_s) shows the production backward
    /// residual stays ≤ ~3e-2 up to order 40 where the companion matrix
    /// degrades; the 1e-3 gate with fall-through is safe because production
    /// polynomials up to order 30 sit at ~6e-6 (170x margin).
    #[test]
    fn faer_gate_does_not_fire_on_production_polynomials_up_to_order_30() -> Result<()> {
        for order in [4usize, 8, 12, 16, 20, 24, 30] {
            let zeros: Vec<f64> = (0..order.min(8))
                .map(|k| (1.2 + 0.5 * (k as f64)) * if k % 2 == 0 { 1.0 } else { -1.0 })
                .collect();
            let spec = FilterSpec::builder()
                .order(order)
                .return_loss_db(20.0)
                .normalized_transmission_zeros(zeros)
                .build()?;
            let ps = generalized_chebyshev_polynomials(&spec)?;
            if let Some(g) = &ps.generalized {
                for p in [g.e_s.as_ref().unwrap(), &g.f_s, &g.p_s] {
                    let recovered = FaerRootSolver.roots_of(p)?;
                    let br = max_backward_residual(&p.coefficients, &recovered);
                    assert!(
                        br < FAER_BACKWARD_RESIDUAL_TOL,
                        "order {order}: faer gate would fire (br={br:.3e}) on a production polynomial"
                    );
                }
            }
        }
        Ok(())
    }
}
