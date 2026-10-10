//! Spectral factorization methods for polynomial decomposition.
//!
//! Given a non-negative polynomial R(s) on the jω-axis such that R(jω) ≥ 0 ∀ω,
//! spectral factorization finds D(s) with all roots in the strict LHP such that:
//!   R(s) = D(s) · D*(-s)
//!
//! Two methods are provided:
//! 1. **Cepstrum method** (FFT-based): Numerically robust, O(N log N), avoids root-finding
//! 2. **Balanced companion matrix**: Improved eigenvalue-based root-finding with matrix balancing

use crate::approx::ComplexPolynomial;
use crate::error::{MfsError, Result};
use num_complex::Complex64;
use rustfft::{num_complex::Complex as FftComplex, FftPlanner};

// ═══════════════════════════════════════════════════════════════════════════════
// Cepstrum-based spectral factorization
// ═══════════════════════════════════════════════════════════════════════════════

/// Performs spectral factorization using the cepstrum (log-FFT) method.
///
/// Given R(jω) = |D(jω)|², finds D(s) with all roots in the LHP.
///
/// Algorithm:
/// 1. Evaluate R(jω) on a dense frequency grid via FFT
/// 2. Take log: log R(jω) = log|D(jω)|² = 2·Re{log D(jω)}
/// 3. Compute the real cepstrum (IFFT of log)
/// 4. Apply causal windowing (keep only n≥0 terms, multiply by 2 for n>0)
/// 5. FFT back → log D(jω) → D(jω) → fit polynomial coefficients
///
/// This method is numerically robust for polynomials up to degree ~100 and
/// avoids all root-finding operations.
pub fn cepstrum_spectral_factor(
    r_poly: &ComplexPolynomial,
    target_degree: usize,
) -> Result<ComplexPolynomial> {
    // Use a reasonable FFT size (balance between accuracy and speed)
    let n_fft = (target_degree * 8).max(128).next_power_of_two();

    // Step 1: Evaluate R(jω) on the unit circle in the z-domain
    // We use the bilinear transform: s = (z-1)/(z+1), so jω maps to the unit circle.
    // Alternatively, evaluate R(jω) directly for ω ∈ [-π, π] on a uniform grid.
    //
    // For polynomials in s, we evaluate at s = jω for ω uniformly spaced.
    // The evaluation range must cover the significant frequency content.
    let omega_max = estimate_bandwidth(r_poly);

    // Evaluate R(jω) on uniform grid
    let mut r_values: Vec<Complex64> = Vec::with_capacity(n_fft);
    for k in 0..n_fft {
        let omega = -omega_max + 2.0 * omega_max * (k as f64) / (n_fft as f64);
        let s = Complex64::new(0.0, omega);
        let val = r_poly.evaluate(s);
        r_values.push(val);
    }

    // Step 2: Take logarithm of |R(jω)|
    // R(jω) should be real and non-negative (it's |D(jω)|²).
    // Due to numerical errors it might be slightly complex/negative.
    let mut log_r: Vec<f64> = Vec::with_capacity(n_fft);
    for val in &r_values {
        let mag = val.norm(); // Use magnitude since R should be real positive
        if mag > 1e-300 {
            log_r.push(mag.ln());
        } else {
            log_r.push(-690.0); // log(1e-300)
        }
    }

    // Step 3: Compute real cepstrum via IFFT
    // Convert to complex for FFT
    let mut fft_input: Vec<FftComplex<f64>> =
        log_r.iter().map(|&x| FftComplex::new(x, 0.0)).collect();

    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(n_fft);
    ifft.process(&mut fft_input);

    // Normalize IFFT
    let scale = 1.0 / n_fft as f64;
    for v in fft_input.iter_mut() {
        v.re *= scale;
        v.im *= scale;
    }

    // Step 4: Apply minimum-phase (causal) windowing
    // cepstrum[0] stays as is (halved for the factor)
    // cepstrum[1..N/2] doubled (causal part)
    // cepstrum[N/2+1..] zeroed (anti-causal part)
    let mut min_phase_cepstrum: Vec<FftComplex<f64>> = vec![FftComplex::new(0.0, 0.0); n_fft];
    min_phase_cepstrum[0] = FftComplex::new(fft_input[0].re * 0.5, 0.0);
    // Full causal part
    min_phase_cepstrum[1..n_fft / 2].copy_from_slice(&fft_input[1..n_fft / 2]);
    if n_fft > 1 {
        min_phase_cepstrum[n_fft / 2] = FftComplex::new(fft_input[n_fft / 2].re * 0.5, 0.0);
    }

    // Step 5: FFT back to get log D(jω)
    let fft = planner.plan_fft_forward(n_fft);
    fft.process(&mut min_phase_cepstrum);

    // Step 6: Exponentiate to get D(jω)
    let d_values: Vec<Complex64> = min_phase_cepstrum
        .iter()
        .map(|&v| {
            let log_d = Complex64::new(v.re, v.im);
            log_d.exp()
        })
        .collect();

    // Step 7: Fit polynomial coefficients from D(jω) samples
    // Use DFT-based coefficient extraction (essentially polynomial interpolation
    // on the uniform frequency grid)
    fit_polynomial_from_freq_samples(&d_values, target_degree, omega_max)
}

/// Estimates the bandwidth of a polynomial (where its magnitude is significant).
fn estimate_bandwidth(poly: &ComplexPolynomial) -> f64 {
    // For a degree-n polynomial, the significant content is within ω ∈ [-n, n]
    // approximately. We use a slightly larger range for safety.
    let degree = poly.degree();
    (degree as f64 + 2.0) * 2.0
}

/// Fits a polynomial of given degree from frequency-domain samples.
///
/// Given D(jω_k) for uniformly spaced ω_k, extracts the polynomial coefficients
/// using least-squares fitting.
fn fit_polynomial_from_freq_samples(
    d_values: &[Complex64],
    degree: usize,
    omega_max: f64,
) -> Result<ComplexPolynomial> {
    use nalgebra::{DMatrix, DVector};

    let n_fft = d_values.len();
    let m = degree + 1; // number of coefficients to find

    // Use a subset of well-conditioned sample points
    // Pick points distributed across the frequency range
    let n_samples = (2 * m + 4).min(n_fft);
    let step = n_fft / n_samples;

    // Build the Vandermonde system: V * c = d
    // where V[i,j] = (jω_i)^j and c are the polynomial coefficients
    // Split into real/imaginary parts for real-valued LS solve
    let mut a_mat = DMatrix::<f64>::zeros(2 * n_samples, 2 * m);
    let mut b_vec = DVector::<f64>::zeros(2 * n_samples);

    for (sample_idx, k_idx) in (0..n_fft).step_by(step).take(n_samples).enumerate() {
        let omega = -omega_max + 2.0 * omega_max * (k_idx as f64) / (n_fft as f64);
        let s = Complex64::new(0.0, omega);
        let d_val = d_values[k_idx];

        // Build row of Vandermonde: [s^0, s^1, ..., s^degree]
        let mut s_power = Complex64::new(1.0, 0.0);
        for j in 0..m {
            // Real part equation
            a_mat[(2 * sample_idx, 2 * j)] = s_power.re;
            a_mat[(2 * sample_idx, 2 * j + 1)] = -s_power.im;
            // Imaginary part equation
            a_mat[(2 * sample_idx + 1, 2 * j)] = s_power.im;
            a_mat[(2 * sample_idx + 1, 2 * j + 1)] = s_power.re;

            s_power *= s;
        }

        b_vec[2 * sample_idx] = d_val.re;
        b_vec[2 * sample_idx + 1] = d_val.im;
    }

    // Solve via normal equations: (A^T A) x = A^T b
    let ata = a_mat.transpose() * &a_mat;
    let atb = a_mat.transpose() * &b_vec;

    let solution = ata.lu().solve(&atb).ok_or_else(|| {
        MfsError::NumericalFailure("Cepstrum polynomial fit: singular normal equations".into())
    })?;

    // Extract complex coefficients
    let coeffs: Vec<Complex64> = (0..m)
        .map(|j| Complex64::new(solution[2 * j], solution[2 * j + 1]))
        .collect();

    ComplexPolynomial::new(coeffs)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Balanced companion matrix eigenvalue solver
// ═══════════════════════════════════════════════════════════════════════════════

/// Spectral factorization via companion-matrix eigenvalue decomposition.
///
/// Finds all roots of `R(s)` with the production
/// [`AdaptiveRootSolver`](crate::approx::AdaptiveRootSolver) chain (faer complex companion
/// eigenvalues gated on a backward residual, falling through to Durand-Kerner and
/// Aberth-Ehrlich), then selects the LHP roots for the spectral factor.
///
/// Note: the nalgebra `CompanionMatrixRootSolver` is deliberately NOT used
/// here — its complex-coefficient branch (2n×2n real embedding) returns
/// conjugate-flipped non-roots (the 2n eigenvalues are
/// `{λ, conj(λ)} ∪ {eig(C)}`, so a root and its conjugate are indistinguishable
/// by conjugate-pair matching alone). See
/// `companion_matrix_solver_complex_branch_returns_conjugate_flipped_non_roots`.
pub fn companion_spectral_factor(
    r_poly: &ComplexPolynomial,
    target_degree: usize,
) -> Result<ComplexPolynomial> {
    use crate::approx::AdaptiveRootSolver;
    use crate::approx::ComplexRootSolver;

    // Find all roots of R(s) using the adaptive (faer-gated) chain
    let all_roots = AdaptiveRootSolver.roots_of(r_poly)?;

    // Select LHP roots for the spectral factor
    let stable_roots = select_lhp_roots(&all_roots, target_degree);

    if stable_roots.is_empty() {
        return ComplexPolynomial::new(vec![Complex64::new(1.0, 0.0)]);
    }

    // Build D(s) from stable roots
    let mut result = ComplexPolynomial::from_complex_roots(&stable_roots)?;

    // Scale so that D(s)*D*(-s) matches R(s) at s = 0
    let r_at_0 = r_poly.evaluate(Complex64::new(0.0, 0.0));
    let d_at_0 = result.evaluate(Complex64::new(0.0, 0.0));
    let d_conj = result.alternating_conjugate()?;
    let d_star_at_0 = d_conj.evaluate(Complex64::new(0.0, 0.0));
    let product_at_0 = d_at_0 * d_star_at_0;

    if product_at_0.norm() > 1e-15 && r_at_0.norm() > 1e-15 {
        let scale = (r_at_0 / product_at_0).sqrt();
        if scale.re.is_finite() && scale.im.is_finite() {
            result = result.scale(scale)?;
        }
    }

    Ok(result)
}
/// Selects LHP roots from a set of roots that come in (z, -z*) pairs.
fn select_lhp_roots(all_roots: &[Complex64], target_count: usize) -> Vec<Complex64> {
    let mut stable: Vec<Complex64> = Vec::new();
    let mut used = vec![false; all_roots.len()];

    // Sort by real part
    let mut indices: Vec<usize> = (0..all_roots.len()).collect();
    indices.sort_by(|&a, &b| {
        all_roots[a]
            .re
            .partial_cmp(&all_roots[b].re)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Take roots from LHP first
    for &idx in &indices {
        if used[idx] {
            continue;
        }
        if stable.len() >= target_count {
            break;
        }

        let root = all_roots[idx];
        if !root.re.is_finite() || !root.im.is_finite() {
            used[idx] = true;
            continue;
        }

        if root.re < -1e-12 {
            stable.push(root);
            used[idx] = true;

            // Mark mirror partner (-conj(z)) as used
            for &j in &indices {
                if used[j] {
                    continue;
                }
                let mirror = all_roots[j];
                if (mirror.re + root.re).abs() < 1e-6 * (1.0 + root.re.abs())
                    && (mirror.im + root.im).abs() < 1e-6 * (1.0 + root.im.abs())
                {
                    used[j] = true;
                    break;
                }
            }
        }
    }

    // Handle imaginary axis roots
    if stable.len() < target_count {
        for &idx in &indices {
            if used[idx] {
                continue;
            }
            if stable.len() >= target_count {
                break;
            }

            let root = all_roots[idx];
            if root.re.abs() <= 1e-12 && root.re.is_finite() && root.im.is_finite() {
                stable.push(Complex64::new(-1e-14, root.im));
                used[idx] = true;
                // Mark conjugate pair
                for &j in &indices {
                    if used[j] {
                        continue;
                    }
                    let m = all_roots[j];
                    if (m.re.abs()) <= 1e-10
                        && (m.im + root.im).abs() < 1e-6 * (1.0 + root.im.abs())
                    {
                        used[j] = true;
                        break;
                    }
                }
            }
        }
    }

    // Fill remaining with closest-to-LHP roots
    if stable.len() < target_count {
        let mut remaining: Vec<(usize, f64)> = indices
            .iter()
            .filter(|&&i| !used[i] && all_roots[i].re.is_finite())
            .map(|&i| (i, all_roots[i].re))
            .collect();
        remaining.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (idx, _) in remaining {
            if stable.len() >= target_count {
                break;
            }
            stable.push(all_roots[idx]);
        }
    }

    stable.truncate(target_count);
    stable
}

// ═══════════════════════════════════════════════════════════════════════════════
// Combined spectral factorization: tries cepstrum first, companion as fallback
// ═══════════════════════════════════════════════════════════════════════════════

/// Checks the para-product identity `D(jω)·D̃(jω) = R(jω)` (= `|D(jω)|² = R(jω)`
/// on the jω axis) at 8 axis points and returns `true` when the worst-case
/// relative residual is below `0.01`.
///
/// This is a **magnitude** check: it verifies `D` is a correct spectral factor of
/// `R` (the right S-parameter response), not a grossly mis-factored or mis-scaled
/// one. It does **not** by itself enforce minimum-phase — for a real `R`, the
/// minimum-phase factor and its maximum-phase anti-factor share `|D(jω)|²`, so
/// both pass (a single point `s = j` would pass even more wrong factors, so 8
/// points are sampled). Minimum-phase *selection* is `select_lhp_roots`' LHP
/// bias; this check is the backstop that rejects a factor whose magnitude
/// response is wrong, which would otherwise silently corrupt the downstream
/// Machiarella coupling matrices (review round-2 C2).
fn para_product_verified(
    r_poly: &ComplexPolynomial,
    d: &ComplexPolynomial,
) -> bool {
    let mut worst = 0.0_f64;
    for k in 1..=8usize {
        let omega = 0.375 * k as f64 * (r_poly.degree() as f64 + 2.0).sqrt();
        let s_j = Complex64::new(0.0, omega);
        let r_val = r_poly.evaluate(s_j);
        let product = d.evaluate(s_j).norm_sqr();
        let rel_err = if r_val.norm() > 1e-10 {
            (product - r_val).norm() / r_val.norm()
        } else {
            product.abs()
        };
        worst = worst.max(rel_err);
    }
    worst < 0.01
}

/// Robust spectral factorization using cepstrum method with companion matrix fallback.
///
/// Given R(s) = sum of polynomial products, finds D(s) such that D(s)D*(-s) = R(s)
/// with D(s) having all roots in the LHP.
pub fn spectral_factor(
    r_poly: &ComplexPolynomial,
    target_degree: usize,
) -> Result<ComplexPolynomial> {
    // Try cepstrum method first (most robust for well-conditioned problems)
    if let Ok(d) = cepstrum_spectral_factor(r_poly, target_degree) {
        if para_product_verified(r_poly, &d) {
            return monic(d);
        }
        // Fall through to companion matrix if the cepstrum result is poor,
        // and likewise when the cepstrum method errors out.
    }

    // Fallback: companion matrix (adaptive, faer-gated) approach. The result
    // must also pass the para-product verification. `select_lhp_roots` can admit
    // a right-half-plane root when the LHP count is short (its last-resort fill
    // has no sign check); a wrong-magnitude factor — the kind that corrupts the
    // downstream coupling matrices — only the `|D(jω)|² = R(jω)` identity check
    // distinguishes from a correct one. (Minimum-phase *selection* is the LHP
    // bias in `select_lhp_roots`; the verification is the magnitude backstop.)
    let d = companion_spectral_factor(r_poly, target_degree)?;
    if !para_product_verified(r_poly, &d) {
        return Err(MfsError::NumericalFailure(
            "spectral factorization: companion-matrix factor failed para-product \
             verification (|D(jω)|² ≠ R(jω)); factorization does not reproduce \
             R's magnitude response"
                .to_string(),
        ));
    }
    monic(d)
}

/// Rotates a spectral factor by a unit phase so its leading (highest-degree)
/// coefficient is real and positive — the paper's monic normalization. The LHP
/// factor is otherwise determined only up to a unit phase `k` (|k|² = 1), which
/// leaves R(s) = D·D̃ invariant but changes any linear combination like
/// `(a·N + b·D)/2` in the diplexer iteration.
fn monic(d: ComplexPolynomial) -> Result<ComplexPolynomial> {
    if d.degree() == 0 {
        return Ok(d);
    }
    let lead = *d.coefficients.last().unwrap();
    if lead.norm() < 1e-15 {
        return Ok(d);
    }
    // Rotation k = conj(lead)/|lead|: unit phase taking lead → |lead| > 0.
    let rotation = lead.conj() / lead.norm();
    d.scale(rotation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cepstrum_factors_simple_polynomial() -> Result<()> {
        // R(s) = (s+1)(s-1)(-s+1)(-s-1) = (1-s²)(1-s²) = (1-s²)²
        // But in s-domain: R(s) = D(s)*D*(-s) where D(s) = s+1 (degree 1)
        // Actually: (s+1)*(-s+1) = -(s²-1) = 1-s²...
        //
        // Simpler: R(s) = (s+1)(-s+1) = 1 - s² evaluated as poly:
        // coefficients [1, 0, -1] in ascending power
        //
        // D(s) = s + 1, D*(-s) = -s + 1 = 1 - s
        // D(s)*D*(-s) = (s+1)(1-s) = 1 - s²  ✓

        let r = ComplexPolynomial::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(-1.0, 0.0),
        ])?;

        let d = spectral_factor(&r, 1)?;
        assert_eq!(d.degree(), 1);

        // D should have root at s = -1 (LHP)
        let roots = d.roots()?;
        assert!(
            (roots[0].re + 1.0).abs() < 0.1 || (roots[0].re - (-1.0)).abs() < 0.1,
            "Root should be near -1, got {:.4}",
            roots[0]
        );

        Ok(())
    }

    #[test]
    fn companion_factors_degree4_complex_branch() -> Result<()> {
        // D(s) = s² + 2s + 2 (roots at -1±j)
        // D*(-s) = s² - 2s + 2 (roots at 1±j)
        // R(s) = D(s)*D*(-s) = (s²+2s+2)(s²-2s+2) = s⁴ + 4
        // coefficients: [4, 0, 0, 0, 1]
        //
        // This exercises the complex-coefficient root path end to end
        // (R(s) is real here, but the factor's roots -1±j are complex, so
        // LHP selection must pair them correctly). It previously hung /
        // mis-factored under the nalgebra 2n×2n real embedding; the
        // adaptive (faer-gated) chain resolves it.
        let r = ComplexPolynomial::new(vec![
            Complex64::new(4.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ])?;

        let d = companion_spectral_factor(&r, 2)?;
        assert_eq!(d.degree(), 2);

        // D roots should be at -1±j (LHP)
        let roots = d.roots()?;
        for root in &roots {
            assert!(
                root.re < 0.0 && (root.re + 1.0).abs() < 1e-6 && (root.im.abs() - 1.0).abs() < 1e-6,
                "Root should be -1±j: {root:.6}",
            );
        }

        // Para-product identity: D(v)·D̃(v) = R(v) (up to the s=0 scaling,
        // which is exact here since R is real and even).
        let d_conj = d.alternating_conjugate()?;
        for omega in [0.0_f64, 0.5, 1.0, 2.0, 5.0] {
            let v = Complex64::new(0.0, omega);
            let lhs = d.evaluate(v) * d_conj.evaluate(v);
            let rhs = r.evaluate(v);
            assert!(
                (lhs - rhs).norm() / rhs.norm().max(1e-300) < 1e-8,
                "para-product identity violated at ω={omega}: |{lhs} - {rhs}|",
            );
        }

        Ok(())
    }

    #[test]
    fn spectral_factor_complex_r_poly_matches_known_lhp_set() -> Result<()> {
        // Regression for the conjugate-flip defect: R(s) with complex
        // coefficients whose LHP set is known. Build R from a known LHP
        // factor D, then check the dispatcher recovers the same LHP roots.
        let d_roots: Vec<Complex64> = vec![
            Complex64::new(-1.6718464616, 1.3824014892),
            Complex64::new(-1.6718464616, 1.1176531574),
            Complex64::new(-4.1116410137, -3.3228232740),
            Complex64::new(-4.1116410137, 5.8228779206),
        ];
        let d = ComplexPolynomial::from_complex_roots(&d_roots)?;
        let d_conj = d.alternating_conjugate()?;
        let r = d.multiply(&d_conj)?;

        let d_rec = spectral_factor(&r, 4)?;
        assert_eq!(d_rec.degree(), 4);

        // Recovered roots (up to the s=0 scale, which is positive here) must
        // match the LHP set — specifically the -4.1116-3.3228j root, whose
        // conjugate-flipped twin -4.1116+3.3228j is NOT a root of R.
        //
        // Match setwise (each expected root to its nearest recovered root),
        // not elementwise: two roots of the fixture share a real part to
        // ~1e-15, so any (re, im) sort order of the recovered set can
        // legitimately differ from the literal's order while the sets are
        // identical.
        let rec = d_rec.roots()?;
        assert_eq!(rec.len(), d_roots.len());
        let mut used = vec![false; rec.len()];
        for want in &d_roots {
            let (best_i, best_d) = rec
                .iter()
                .enumerate()
                .filter(|(i, _)| !used[*i])
                .map(|(i, got)| (i, (*got - *want).norm()))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .expect("no recovered roots left");
            assert!(
                best_d < 1e-6,
                "root mismatch: want {want}, nearest recovered {best_d:.3e} away",
            );
            used[best_i] = true;
        }

        // Pin the factor's scale/phase: the para-product identity
        // D̃(v)·D(v) = R(v) must hold on the axis (D̃(v)=conj(D(v)) there),
        // up to the s=0 scaling. A D with the right LHP roots but a wrong
        // magnitude or phase would fail this.
        let d_rec_conj = d_rec.alternating_conjugate()?;
        for omega in [0.0_f64, 0.5, 1.0, 2.0, 5.0] {
            let v = Complex64::new(0.0, omega);
            let lhs = d_rec.evaluate(v) * d_rec_conj.evaluate(v);
            let rhs = r.evaluate(v);
            assert!(
                (lhs - rhs).norm() / rhs.norm().max(1e-300) < 1e-6,
                "para-product identity violated at ω={omega}: |{lhs} - {rhs}|",
            );
        }
        Ok(())
    }

    /// C2 regression: the shared para-product verification is a **magnitude**
    /// check. It accepts the correct spectral factor of R and rejects a
    /// same-degree factor with a *different* jω-axis magnitude — the kind of
    /// wrong-magnitude root set `select_lhp_roots` can admit when its
    /// last-resort fill has no `re < 0` check (it can pull in an RHP root whose
    /// jω-axis magnitude differs from the true LHP root). This is the backstop
    /// that keeps the companion fallback from returning a grossly wrong-magnitude
    /// factor.
    ///
    /// It is **not** a minimum-phase enforcer: for a real R, the min-phase and
    /// max-phase anti-factors share `|D(jω)|²`, so both pass. That limitation is
    /// pinned explicitly below (the test uses a complex R for the reject case,
    /// where a single reflected root *does* change the magnitude).
    #[test]
    fn para_product_verification_checks_magnitude_not_min_phase() -> Result<()> {
        // --- Accept + reject on a complex R (where root placement changes |.|²) ---
        // True LHP factor: D(s) = (s - (-1+0.5j))(s - (-1-0.5j)).
        let true_roots: Vec<Complex64> = vec![
            Complex64::new(-1.0, 0.5),
            Complex64::new(-1.0, -0.5),
        ];
        let d_true = ComplexPolynomial::from_complex_roots(&true_roots)?;
        let d_true_conj = d_true.alternating_conjugate()?;
        let r = d_true.multiply(&d_true_conj)?;
        assert!(
            para_product_verified(&r, &d_true),
            "the correct factor must pass verification"
        );

        // A wrong same-degree factor: reflect one root to the RHP
        // (-1+0.5j → 1-0.5j). Its jω-axis magnitude differs, so it must be
        // rejected.
        let wrong_roots: Vec<Complex64> = vec![
            Complex64::new(1.0, -0.5),
            Complex64::new(-1.0, -0.5),
        ];
        let d_wrong = ComplexPolynomial::from_complex_roots(&wrong_roots)?;
        assert!(
            !para_product_verified(&r, &d_wrong),
            "a wrong-magnitude factor must be rejected"
        );

        // --- Document the limitation: for a real R, the max-phase anti-factor
        //     passes too (it has the same |D(jω)|²). The verification is a
        //     magnitude check, not a min-phase enforcer.
        // R = s⁴ + 4 = (s²+2s+2)(s²-2s+2); min-phase factor s²+2s+2 (roots
        // -1±j), max-phase anti-factor s²-2s+2 (roots +1±j). Both pass.
        let r_real = ComplexPolynomial::new(vec![
            Complex64::new(4.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ])?;
        let d_min = ComplexPolynomial::new(vec![
            Complex64::new(2.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(1.0, 0.0),
        ])?;
        let d_max = ComplexPolynomial::new(vec![
            Complex64::new(2.0, 0.0),
            Complex64::new(-2.0, 0.0),
            Complex64::new(1.0, 0.0),
        ])?;
        assert!(
            para_product_verified(&r_real, &d_min),
            "min-phase factor of real R must pass"
        );
        assert!(
            para_product_verified(&r_real, &d_max),
            "max-phase anti-factor of real R ALSO passes (magnitude is identical)"
        );
        Ok(())
    }

}
