//! Pole-expansion (partial-fraction) frequency response evaluation.
//!
//! Instead of solving a linear system at each frequency point (O(N³) per point),
//! this module evaluates S-parameters directly from the partial-fraction expansion
//! of the admittance parameters:
//!
//!   Y12(s) = Σ r_k / (s - p_k) + constant
//!
//! This gives O(N) per frequency point — typically 100-500× faster than LU-based
//! evaluation for large frequency sweeps.

use num_complex::Complex64;

use crate::error::{MfsError, Result};
use crate::freq::{FrequencyGrid, FrequencyMapping};

use super::{ResponseSample, SParameterResponse};

/// Pre-computed pole-expansion data for fast frequency sweeps.
///
/// Once constructed from a transversal coupling matrix via
/// [`from_matrix`](Self::from_matrix), this can evaluate S-parameters at
/// arbitrary frequencies in O(N) per point.
#[derive(Debug, Clone)]
pub struct PoleExpansionData {
    /// Poles shared by all Y parameters (denominator roots).
    poles: Vec<Complex64>,
    /// Residues of Y11 at each pole (M_Sk²).
    residues_y11: Vec<Complex64>,
    /// Residues of Y12 at each pole (M_Sk * M_kL).
    residues_y12: Vec<Complex64>,
    /// Residues of Y22 at each pole (M_kL²).
    residues_y22: Vec<Complex64>,
    /// Laurent principal parts for repeated (multiplicity >= 2) poles, one
    /// `Option` per pole: `Some([b_0, ..., b_{m-1}])` where `b_k` multiplies
    /// `1/(jω - p)^(m - k)`, or `None` for a simple pole. Each Y parameter's
    /// principal part is kept with its own length: a shared pole's effective
    /// order can differ between the three after numerator-denominator
    /// cancellation, and `evaluate_at` sums each with its own length.
    laurent_y11: Vec<Option<Vec<Complex64>>>,
    /// Laurent principal parts for Y12 (see `laurent_y11`).
    laurent_y12: Vec<Option<Vec<Complex64>>>,
    /// Laurent principal parts for Y22 (see `laurent_y11`).
    laurent_y22: Vec<Option<Vec<Complex64>>>,
    /// Constant term of Y12 (direct source-load coupling).
    y12_constant: Complex64,
    /// Constant term of Y11 (reserved for future use).
    #[allow(dead_code)]
    y11_constant: Complex64,
    /// Source resistance (normalized).
    source_r: f64,
    /// Load resistance (normalized).
    load_r: f64,
}

impl PoleExpansionData {
    /// Builds pole-expansion data from the transversal coupling matrix.
    ///
    /// For a transversal coupling matrix, the response matrix at frequency ω is:
    ///   A(ω) = jωI_r + M - jR_S·e₀e₀ᵀ - jR_L·eₙeₙᵀ
    ///
    /// The transversal structure means each resonator k only couples to source and load:
    ///   `M[0,k]` = M_Sk (source coupling)
    ///   `M[k,N+1]` = M_kL (load coupling)  
    ///   `M[k,k]` = B_k (diagonal/detuning)
    ///
    /// The inverse elements needed for S-parameters can be written as partial fractions:
    ///   [A⁻¹]_{0,0} = Σ M_Sk² / (jω - jB_k) / D(ω)  (approximately)
    ///
    /// But the exact formula requires accounting for the port terminations.
    /// We use the direct partial-fraction form of the network function:
    ///   [A⁻¹]_{N+1,0} = Σ (M_Sk · M_kL) / (jω + B_k + port_correction)
    ///
    /// Actually, for the transversal matrix the exact closed-form is:
    ///   [A⁻¹]_{0,0} = (-jR_S + Σ M_Sk²/(jω+B_k)) / denominator_correction
    ///   [A⁻¹]_{N+1,0} = (Σ M_Sk·M_kL/(jω+B_k) + M_SL) / denominator_correction
    ///
    /// The simplest correct approach: extract poles and residues directly from
    /// the transversal matrix structure, matching what the LU solver computes.
    ///
    /// **Pole-domain convention (review P2-12):** the stored poles are the
    /// *real ω-domain* poles `−B_k` (resonator detuning), so that
    /// `evaluate_at` forms `denom = ω − pole` with a real
    /// `ω`. The transversal matrix is purely real, so `−B_k` is real; a
    /// non-transversal (folded/arrow) matrix has complex couplings and must be
    /// transformed to transversal form first.
    ///
    /// **Input contract (enforced):** the resonator block `M[1..=N][1..=N]`
    /// must be diagonal — each resonator couples only to source and load,
    /// never to another resonator — and the port diagonals `M[0,0]` /
    /// `M[N+1,N+1]` must be zero (pure series terminations, no port shunt).
    /// A non-transversal (folded/arrow/wheel) matrix carries inter-resonator
    /// couplings, and a port-shunted matrix carries a port diagonal, that this
    /// O(N) partial-fraction path cannot model; passing either would silently
    /// drop the unmodeled entries and return a wrong response (measured
    /// ~1e-1 in |S21| for a source shunt). `from_matrix` rejects both with
    /// `MfsError::PreconditionViolation` rather than computing a silent
    /// misanswer. Use the LU-based `ResponseSolver` for arbitrary symmetric
    /// matrices, or transform to transversal form with series terminations.
    pub fn from_matrix(matrix: &crate::matrix::CouplingMatrix) -> Result<Self> {
        let order = matrix.order();
        let side = matrix.side();
        Self::assert_transversal_structure(matrix)?;

        let mut poles = Vec::with_capacity(order);
        let mut residues_y11 = Vec::with_capacity(order);
        let mut residues_y12 = Vec::with_capacity(order);
        let mut residues_y22 = Vec::with_capacity(order);

        for k in 1..=order {
            let b_k = matrix.at(k, k).unwrap_or(0.0);
            let m_sk = matrix.at(0, k).unwrap_or(0.0);
            let m_kl = matrix.at(k, side - 1).unwrap_or(0.0);

            // The response matrix diagonal for resonator k is (ω + B_k).
            // The pole (where denominator = 0) is at ω = -B_k.
            poles.push(Complex64::new(-b_k, 0.0));
            // Residues for the partial-fraction sums
            residues_y11.push(Complex64::new(m_sk * m_sk, 0.0));
            residues_y12.push(Complex64::new(m_sk * m_kl, 0.0));
            residues_y22.push(Complex64::new(m_kl * m_kl, 0.0));
        }

        // Direct source-load coupling (M[0, N+1])
        let m_sl = matrix.at(0, side - 1).unwrap_or(0.0);
        let y12_constant = Complex64::new(m_sl, 0.0);

        Ok(Self {
            poles,
            residues_y11,
            residues_y12,
            residues_y22,
            laurent_y11: vec![None; order],
            laurent_y12: vec![None; order],
            laurent_y22: vec![None; order],
            y12_constant,
            y11_constant: Complex64::new(0.0, 0.0),
            source_r: 1.0,
            load_r: 1.0,
        })
    }

    /// Structural precondition for [`from_matrix`]: the resonator block
    /// `M[1..=N][1..=N]` must be diagonal (each resonator couples only to
    /// source and load, never to another resonator).
    ///
    /// This is a *structural* check on the actual entries, not a test of the
    /// topology label. Folded/arrow/wheel matrices carry inter-resonator
    /// couplings `M[k][k']` (k≠k') that the O(N) partial-fraction model
    /// cannot represent; dropping them would yield a physically wrong but
    /// silently computed response. The Schur-2×2 model further assumes the
    /// port diagonals `M[0,0]` / `M[N+1,N+1]` are zero (pure series
    /// terminations); a nonzero port diagonal is dropped silently by the
    /// model, so it is rejected as well. We reject any off-diagonal resonator
    /// entry or port diagonal above a scale-aware tolerance.
    ///
    /// The tolerance is `1e-9 × max(1, largest |M| entry)`. It sits far above
    /// f64 rounding on a normalized matrix (~1e-16) and far below any
    /// physically meaningful normalized coupling (O(1)), so it accepts a
    /// genuinely transversal (up to round-off) matrix while rejecting any real
    /// inter-resonator coupling. It is deliberately not the main path's 1e-3
    /// S-parameter threshold, which measures response error and is unrelated
    /// to this structural precondition.
    fn assert_transversal_structure(matrix: &crate::matrix::CouplingMatrix) -> Result<()> {
        let order = matrix.order();
        let side = matrix.side();
        let data = matrix.as_slice();
        let max_entry = data.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let tol = 1e-9 * max_entry.max(1.0);
        for (i, row) in data.chunks_exact(side).enumerate() {
            if i < 1 || i > order {
                continue; // source (row 0) and load (row side-1) are not resonators
            }
            for (j, &value) in row.iter().enumerate() {
                if j == i || j < 1 || j > order {
                    continue; // skip the diagonal and the port columns
                }
                if value.abs() > tol {
                    return Err(MfsError::PreconditionViolation(format!(
                        "PoleExpansionData::from_matrix requires a transversal matrix, \
                         but the resonator-resonator coupling M[{i}][{j}] = {value} \
                         (|M| > {tol}) is present. Folded/arrow and other non-transversal \
                         topologies carry inter-resonator couplings that the O(N) fast path \
                         cannot model; transform to transversal form first, or use the \
                         LU-based ResponseSolver, which handles arbitrary symmetric matrices."
                    )));
                }
            }
        }
        // Port diagonals: the Schur-2×2 model assumes M[0,0] = M[N+1,N+1] = 0
        // (pure series terminations; port shunts are not part of the
        // coupling-matrix convention this path models). A nonzero port
        // diagonal would be silently dropped by the model, so reject it with
        // the same scale-aware tolerance as the resonator block.
        for (i, j) in [(0, 0), (order + 1, order + 1)] {
            if data[i * side + j].abs() > tol {
                return Err(MfsError::PreconditionViolation(format!(
                    "PoleExpansionData::from_matrix requires a transversal \
                     matrix with pure series terminations, but the port \
                     diagonal M[{i}][{j}] = {} (|M| > {tol}) is nonzero. The \
                     O(N) fast path models only M[0,k], M[k,k], M[k,N+1] and \
                     M[0,N+1]; a port shunt is not modeled. Use the LU-based \
                     ResponseSolver, which handles arbitrary symmetric \
                     matrices, or zero the port diagonal if it is spurious.",
                    data[i * side + j]
                )));
            }
        }
        Ok(())
    }

    /// Evaluates S11, S21 and the group delay at a single normalized frequency
    /// point.
    ///
    /// For the normalized network (R_S = R_L = 1):
    ///   S11 = 1 + 2j · [A⁻¹]_{0,0}
    ///   S21 = -2j · [A⁻¹]_{N+1,0}
    ///
    /// Group delay is `−Im(d/dω ln S21)`, evaluated analytically from the
    /// same Schur 2×2 form: with S21 = `j·scale·s01/det` (j·scale a
    /// frequency-independent constant),
    ///   d/dω ln S21 = ds01/s01 − ddet/det,
    /// where the per-pole derivatives d/dω [r_k/(ω−p_k)] = −r_k/(ω−p_k)²
    /// (simple poles) and d/dω [b_k/(ω−p_k)^{m−k}] = −(m−k)·b_k/(ω−p_k)^{m−k+1}
    /// (Laurent parts) are accumulated in the
    /// same single pass, so it stays **O(N) per point**. This is the same
    /// quantity the per-point LU backend returns —
    /// `Im Σ_k A⁻¹[k,N+1]·A⁻¹[k,0] / A⁻¹[N+1,0]`, from
    /// `dA⁻¹/dω = −A⁻¹·(dA/dω)·A⁻¹` with `dA/dω = I_r` on the resonator block
    /// only — both reduce to `−Im(d ln S21/dω)`. At a transmission zero where
    /// S21 → 0 the expression is undefined, so it is returned as `None`
    /// (explicit, not a fake 0.0), mirroring the LU backend's degenerate guard.
    fn evaluate_at(&self, omega: f64) -> Result<(Complex64, Complex64, Option<f64>)> {
        // For a transversal coupling matrix, the response matrix A has the structure:
        //   A[0,0] = -jR_S,  A[0,k] = M_Sk,  A[0,N+1] = M_SL
        //   A[k,k] = jω + B_k (resonator diagonal)
        //   A[k,0] = M_Sk, A[k,N+1] = M_kL
        //   A[N+1,N+1] = -jR_L, A[N+1,k] = M_kL, A[N+1,0] = M_SL
        //
        // Using Schur complement on the resonator block D = diag(jω+B_k):
        //   S[0,0] = -jR_S - Σ M_Sk² / (jω + B_k)
        //   S[0,1] = M_SL - Σ M_Sk·M_kL / (jω + B_k)
        //   S[1,1] = -jR_L - Σ M_kL² / (jω + B_k)
        //   [A⁻¹]_{0,0} = S⁻¹[0,0],  [A⁻¹]_{N+1,0} = S⁻¹[1,0]
        let jw = Complex64::new(omega, 0.0); // omega is real frequency variable
        let neg_j_rs = Complex64::new(0.0, -self.source_r);
        let neg_j_rl = Complex64::new(0.0, -self.load_r);

        // Partial-fraction sums (resonator contributions) and their ω-derivatives.
        // A simple pole r/(jω−p) differentiates to −j·r/(jω−p)²; a repeated pole's
        // Laurent part b_k/(jω−p)^{m−k} differentiates to −j·(m−k)·b_k/(jω−p)^{m−k+1}.
        let mut sum_11 = Complex64::new(0.0, 0.0);
        let mut sum_12 = Complex64::new(0.0, 0.0);
        let mut sum_22 = Complex64::new(0.0, 0.0);
        let mut ds_11 = Complex64::new(0.0, 0.0);
        let mut ds_12 = Complex64::new(0.0, 0.0);
        let mut ds_22 = Complex64::new(0.0, 0.0);

        for i in 0..self.poles.len() {
            // poles[i] = -B_k (real), so jω − poles[i] = jω + B_k
            let denom = jw - self.poles[i];
            if denom.norm_sqr() > 1e-30 {
                // A repeated pole contributes its full Laurent principal part
                // to its own partial-fraction sum, summed with its *own* length
                // (a shared pole's effective order can differ between y11/y12/y22
                // after numerator-denominator cancellation).
                if let Some(l) = self.laurent_y11[i].as_ref() {
                    sum_11 += Self::laurent_sum(l, denom);
                    ds_11 += Self::laurent_derivative(l, denom);
                } else {
                    sum_11 += self.residues_y11[i] / denom;
                    ds_11 -= self.residues_y11[i] / (denom * denom);
                }
                if let Some(l) = self.laurent_y12[i].as_ref() {
                    sum_12 += Self::laurent_sum(l, denom);
                    ds_12 += Self::laurent_derivative(l, denom);
                } else {
                    sum_12 += self.residues_y12[i] / denom;
                    ds_12 -= self.residues_y12[i] / (denom * denom);
                }
                if let Some(l) = self.laurent_y22[i].as_ref() {
                    sum_22 += Self::laurent_sum(l, denom);
                    ds_22 += Self::laurent_derivative(l, denom);
                } else {
                    sum_22 += self.residues_y22[i] / denom;
                    ds_22 -= self.residues_y22[i] / (denom * denom);
                }
            }
        }

        // 2×2 Schur complement matrix
        let s00 = neg_j_rs - sum_11;
        let s01 = self.y12_constant - sum_12; // M_SL - Σ...
        let s10 = s01; // symmetric
        let s11 = neg_j_rl - sum_22;

        // ω-derivatives of the Schur elements (resistances and M_SL are
        // frequency-independent).
        let ds00 = -ds_11;
        let ds01 = -ds_12;
        let ds11 = -ds_22;

        // Invert the 2×2 matrix. A near-singular Schur determinant marks a
        // pole/resonance of the reduced network; the partial-fraction form is
        // not evaluable there, so surface a typed error (the caller falls back
        // to LU) rather than emitting a fake "(S11, S21) = (1, 0)" total-
        // reflection sample (review P3 #1).
        let det = s00 * s11 - s01 * s10;
        if det.norm_sqr() <= 1e-30 {
            return Err(MfsError::NumericalFailure(format!(
                "pole expansion: near-singular Schur determinant (|det|² = {}) at ω = {omega}; \
                 frequency is a network pole, use LU-based evaluation",
                det.norm_sqr()
            )));
        }

        let inv_00 = s11 / det;
        let inv_10 = -s10 / det;

        // S-parameters from inverse matrix elements
        let s_11 = Complex64::new(1.0, 0.0) + Complex64::new(0.0, 2.0 * self.source_r) * inv_00;
        let scale = 2.0 * (self.source_r * self.load_r).sqrt();
        let s_21 = Complex64::new(0.0, -scale) * inv_10;

        // Group delay = −Im(d/dω ln S21) = −Im(ds01/s01 − ddet/det), with
        // ddet = ds00·s11 + s00·ds11 − 2·s01·ds01. Undefined at a transmission
        // zero (|S21| → 0, i.e. |inv_10| → 0), so guard exactly like the LU
        // backend (`denominator.norm_sqr() <= 1e-18`).
        let ddet = ds00 * s11 + s00 * ds11 - 2.0 * s01 * ds01;
        let group_delay = if inv_10.norm_sqr() <= 1e-18 {
            None
        } else {
            Some(-(ds01 / s01 - ddet / det).im)
        };

        Ok((s_11, s_21, group_delay))
    }

    /// Derivative of a Laurent principal part `Σ_k b_k/denom^{m−k}` (where
    /// `m = l.len()`): `d/dω Σ_k b_k/(jω−p)^{m−k} = −j·Σ_k (m−k)·b_k/(jω−p)^{m−k+1}`.
    fn laurent_derivative(l: &[Complex64], denom: Complex64) -> Complex64 {
        let m = l.len();
        let mut sum = Complex64::new(0.0, 0.0);
        for (k, &b) in l.iter().enumerate() {
            let power = m - k; // exponent of b_k on denom
            sum -= (power as f64) * b / denom.powi(power as i32 + 1);
        }
        sum
    }

    /// Sums one pole's Laurent principal part:
    /// `b_0/denom^m + b_1/denom^{m-1} + ... + b_{m-1}/denom`, where
    /// `b_k` multiplies `1/denom^(m-k)` and `m = laurent.len()`.
    fn laurent_sum(laurent: &[Complex64], denom: Complex64) -> Complex64 {
        let m = laurent.len();
        // Accumulate 1/denom, 1/denom^2, ..., 1/denom^m, then reverse so
        // that index k holds 1/denom^{m-k} (b_0 gets the most singular).
        let mut powers = Vec::with_capacity(m);
        let mut inv = Complex64::new(1.0, 0.0) / denom;
        for _ in 0..m {
            powers.push(inv);
            inv /= denom;
        }
        powers.reverse();
        let mut sum = Complex64::new(0.0, 0.0);
        for (k, &b) in laurent.iter().enumerate() {
            sum += b * powers[k];
        }
        sum
    }

    /// Evaluates S-parameters on a normalized frequency grid.
    ///
    /// Returns an error if any grid point hits a near-singular Schur
    /// determinant (a network pole); the caller should fall back to LU-based
    /// evaluation (review P3 #1).
    pub fn evaluate_normalized(&self, grid: &FrequencyGrid) -> Result<SParameterResponse> {
        let samples = grid
            .as_slice()
            .iter()
            .map(|&omega| {
                let (s11, s21, group_delay) = self.evaluate_at(omega)?;
                Ok(ResponseSample {
                    frequency_hz: omega,
                    normalized_omega: omega,
                    group_delay,
                    s11_re: s11.re,
                    s11_im: s11.im,
                    s21_re: s21.re,
                    s21_im: s21.im,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(SParameterResponse { samples })
    }

    /// Evaluates S-parameters on a physical frequency grid with band-pass mapping.
    pub fn evaluate_bandpass(
        &self,
        grid: &FrequencyGrid,
        mapping: &impl FrequencyMapping,
    ) -> Result<SParameterResponse> {
        let samples = grid
            .as_slice()
            .iter()
            .map(|&freq_hz| {
                let normalized = mapping.map_hz_to_normalized(freq_hz)?;
                let omega = normalized.omega;
                let (s11, s21, group_delay) = self.evaluate_at(omega)?;
                Ok(ResponseSample {
                    frequency_hz: freq_hz,
                    normalized_omega: omega,
                    group_delay,
                    s11_re: s11.re,
                    s11_im: s11.im,
                    s21_re: s21.re,
                    s21_im: s21.im,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(SParameterResponse { samples })
    }
}

#[cfg(test)]
mod tests {
    use super::PoleExpansionData;
    use crate::design::FilterDesign;
    use crate::freq::FrequencyGrid;
    use crate::response::ResponseSolver;
    use crate::transform::{transform_matrix, TopologyKind};
    use num_complex::Complex64;

    /// Verify that a repeated pole contributes its full Laurent principal
    /// part to the partial-fraction sums, and that the simple-pole path is
    /// unchanged.
    ///
    /// Hand-computed case (omega = 0.5, all parameters real):
    /// - Double pole at p = 1, y12 Laurent part [b_0, b_1] = [2, 3]:
    ///   contribution to sum_12 = 2/(0.5-1)^2 + 3/(0.5-1) = 8 - 6 = 2
    ///   (y11/y22 laurent parts are [0,0] so sum_11 = sum_22 = 0)
    /// - Simple pole at p = 2, y12 residue = 5: contribution = 5/(0.5-2) = -10/3
    ///
    /// So sum_12 = 2 - 10/3 = -4/3, y12_constant = 0:
    ///   s00 = s11s = -j, s01 = 4/3
    ///   det = (-j)(-j) - (4/3)^2 = -25/9
    ///   inv_00 = (-j)/(-25/9) = 9j/25,  inv_10 = -(4/3)/(-25/9) = 12/25
    ///   s_11 = 1 + 2j * (9j/25) = 7/25
    ///   s_21 = -2j * (12/25) = -24j/25
    #[test]
    fn laurent_principal_part_contributes_to_sum() {
        let c = Complex64::new(0.0, 0.0);
        let data = PoleExpansionData {
            poles: vec![Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)],
            residues_y11: vec![c, c],
            residues_y12: vec![c, Complex64::new(5.0, 0.0)],
            residues_y22: vec![c, c],
            laurent_y11: vec![Some(vec![c, c]), None],
            laurent_y12: vec![
                Some(vec![Complex64::new(2.0, 0.0), Complex64::new(3.0, 0.0)]),
                None,
            ],
            laurent_y22: vec![Some(vec![c, c]), None],
            y12_constant: c,
            y11_constant: c,
            source_r: 1.0,
            load_r: 1.0,
        };

        let (s11, s21, _group_delay) = data.evaluate_at(0.5).expect("evaluable point");

        let expected_s11 = Complex64::new(7.0 / 25.0, 0.0);
        assert!((s11 - expected_s11).norm() < 1e-12, "s11 = {s11}");
        let expected_s21 = Complex64::new(0.0, -24.0 / 25.0);
        assert!((s21 - expected_s21).norm() < 1e-12, "s21 = {s21}");
    }

    /// A shared repeated pole can have a *different effective Laurent order*
    /// in each Y parameter after numerator-denominator cancellation: here y11
    /// and y12 keep their double-pole principal parts, but y22's numerator
    /// cancels the pole so it collapses to a simple residue (`laurent: None`).
    /// Each parameter must be summed with its *own* principal part — a single
    /// shared `m` (the old code used `l11.len()` for all three) drops or
    /// over-runs the shorter parameters.
    ///
    /// omega = 2, double pole at p = 1 -> denom = 1:
    ///   sum_11 = 2 + 3 = 5      (y11 [2,3])
    ///   sum_12 = 4 + 5 = 9      (y12 [4,5])
    ///   sum_22 = 6 / 1 = 6      (y22 simple, residue 6)
    ///   s00 = -5 - j, s01 = -9, s11 = -6 - j
    ///   det = (-5-j)(-6-j) - 81 = -52 + 11j
    ///   s_11 = 1 + 2j*s11/det  = 0.91646018 + 0.21309734j
    ///   s_21 = -2j*(-s10)/det  = -0.07008850 + 0.33132744j
    #[test]
    fn repeated_pole_with_cancellation_uses_each_y_own_order() {
        let z = Complex64::new(0.0, 0.0);
        let data = PoleExpansionData {
            poles: vec![Complex64::new(1.0, 0.0)],
            residues_y11: vec![Complex64::new(2.0, 0.0)],
            residues_y12: vec![Complex64::new(4.0, 0.0)],
            residues_y22: vec![Complex64::new(6.0, 0.0)],
            laurent_y11: vec![Some(vec![
                Complex64::new(2.0, 0.0),
                Complex64::new(3.0, 0.0),
            ])],
            laurent_y12: vec![Some(vec![
                Complex64::new(4.0, 0.0),
                Complex64::new(5.0, 0.0),
            ])],
            laurent_y22: vec![None],
            y12_constant: z,
            y11_constant: z,
            source_r: 1.0,
            load_r: 1.0,
        };

        let (s11, s21, _group_delay) = data.evaluate_at(2.0).expect("evaluable point");

        let expected_s11 = Complex64::new(0.91646018, 0.21309734);
        assert!((s11 - expected_s11).norm() < 1e-6, "s11 = {s11}");
        let expected_s21 = Complex64::new(-0.07008850, 0.33132744);
        assert!((s21 - expected_s21).norm() < 1e-6, "s21 = {s21}");
    }
    /// P3 #1: a near-singular Schur determinant (the reduced 2×2 network has
    /// no finite solution) must surface a typed error so the caller falls back
    /// to LU — not a fake `(1, 0)` total-reflection sample. Construct a data
    /// set whose `det = s00*s11 - s01²` is exactly zero at the chosen ω.
    #[test]
    fn near_singular_schur_determinant_is_an_error() {
        let z = Complex64::new(0.0, 0.0);
        // One real pole at -0.5, evaluated at ω = 0.5 (denom = 1.0). A
        // pure-imaginary source residue makes the partial-fraction sum exactly
        // `sum_11 = -j`, cancelling the `-jR_S` diagonal so the (0,0) Schur
        // element is exactly zero: s00 = 0. With no source-load path (s01 = 0)
        // the determinant s00·s11 - s01² is exactly 0 — a genuine network pole.
        // This is the case the old code answered with a fake `(1, 0)` sample;
        // it must now surface a typed error (review P3 #1).
        let data = PoleExpansionData {
            poles: vec![Complex64::new(-0.5, 0.0)],
            residues_y11: vec![Complex64::new(0.0, -1.0)],
            residues_y12: vec![z],
            residues_y22: vec![z],
            laurent_y11: vec![None],
            laurent_y12: vec![None],
            laurent_y22: vec![None],
            y12_constant: z,
            y11_constant: z,
            source_r: 1.0,
            load_r: 1.0,
        };
        let err = data.evaluate_at(0.5).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::NumericalFailure(_)),
            "expected NumericalFailure, got: {err:?}"
        );
        let msg = format!("{err}");
        assert!(msg.contains("near-singular"), "message: {msg}");
    }

    /// A genuinely transversal 2-resonator matrix: the resonator block
    /// `M[1..=2][1..=2]` is diagonal and the ports are series-terminated
    /// (`M[0,0] = M[3,3] = 0`), so `from_matrix` must succeed.
    #[test]
    fn transversal_matrix_is_accepted() -> Result<(), Box<dyn std::error::Error>> {
        use crate::matrix::CouplingMatrix;
        // side = 4 (2 resonators + source + load). Row-major.
        // Resonator block (rows/cols 1,2) is diagonal: M[1][2] = M[2][1] = 0.
        // Port diagonals are zero (pure series terminations): the O(N) fast
        // path's Schur form is exact only when M[0,0] = M[N+1,N+1] = 0, and a
        // nonzero port diagonal is outside its domain (rejected → LU).
        let data = vec![
            // col0(src) col1(res1) col2(res2) col3(load)
            0.0, 0.40, 0.30, 0.0, // row0 (source)
            0.40, 0.10, 0.0, 0.20, // row1 (res 1)
            0.30, 0.0, -0.10, 0.15, // row2 (res 2)
            0.0, 0.20, 0.15, 0.0, // row3 (load) — M[3][3] = 0 (series)
        ];
        let m = CouplingMatrix::new(2, data)?;
        assert_eq!(m.topology(), crate::matrix::MatrixTopology::Transversal);
        PoleExpansionData::from_matrix(&m)?; // must not error
        Ok(())
    }

    /// A transversal matrix whose resonator block carries only f64 round-off
    /// level off-diagonal entries (~1e-12, the residue of a similarity
    /// rotation that happened to land transversal) must be accepted — the
    /// check must not treat `!= 0.0` as "non-transversal".
    #[test]
    fn tiny_roundoff_residual_is_not_rejected() -> Result<(), Box<dyn std::error::Error>> {
        use crate::matrix::CouplingMatrix;
        // Identical to the accepted matrix above except M[1][2] = M[2][1] = 1e-12
        // — far above exact zero, far below any real coupling, well inside the
        // 1e-9 * max(1, max|M|) structural tolerance.
        let data = vec![
            0.0, 0.40, 0.30, 0.0, // row0 (source)
            0.40, 0.10, 1e-12, 0.20, // row1 (res 1)
            0.30, 1e-12, -0.10, 0.15, // row2 (res 2)
            0.0, 0.20, 0.15, 0.0, // row3 (load) — M[3][3] = 0 (series)
        ];
        let m = CouplingMatrix::new(2, data)?;
        PoleExpansionData::from_matrix(&m)?; // round-off residual is accepted

        // The same matrix with a real O(1) resonator-resonator coupling must be
        // rejected — proving the check is not a no-op. (Port diagonals stay
        // zero so the rejection is attributable to M[1][2], not a port shunt.)
        let data_bad = vec![
            0.0, 0.40, 0.30, 0.0, // row0 (source)
            0.40, 0.10, 0.5, 0.20, // row1 (res 1)
            0.30, 0.5, -0.10, 0.15, // row2 (res 2)
            0.0, 0.20, 0.15, 0.0, // row3 (load) — M[3][3] = 0 (series)
        ];
        let m_bad = CouplingMatrix::new(2, data_bad)?;
        let err = PoleExpansionData::from_matrix(&m_bad).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "a real O(1) inter-resonator coupling must be rejected, got: {err:?}"
        );
        Ok(())
    }

    /// A folded topology matrix (carries nearest-neighbor resonator
    /// couplings the O(N) fast path cannot model) must be rejected with a
    /// typed `PreconditionViolation`, not silently answered.
    #[test]
    fn folded_matrix_is_rejected_with_typed_error() {
        let d = FilterDesign::prototype(4, 20.0)
            .zeros([1.5, -1.5])
            .synthesize()
            .expect("synthesis");
        let folded = transform_matrix(d.matrix(), TopologyKind::Folded)
            .expect("folded transform")
            .matrix;
        assert_eq!(folded.topology(), crate::matrix::MatrixTopology::Folded);
        let err = PoleExpansionData::from_matrix(&folded).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "expected PreconditionViolation for a folded matrix, got: {err:?}"
        );
        let msg = format!("{err}");
        assert!(
            msg.contains("transversal"),
            "error message should name the transversal precondition: {msg}"
        );
    }

    /// An arrow topology matrix (dense upper-triangular resonator couplings)
    /// must also be rejected with a typed `PreconditionViolation`.
    #[test]
    fn arrow_matrix_is_rejected_with_typed_error() {
        let d = FilterDesign::prototype(4, 20.0)
            .zeros([1.5, -1.5])
            .synthesize()
            .expect("synthesis");
        let arrow = transform_matrix(d.matrix(), TopologyKind::Arrow)
            .expect("arrow transform")
            .matrix;
        assert_eq!(arrow.topology(), crate::matrix::MatrixTopology::Arrow);
        let err = PoleExpansionData::from_matrix(&arrow).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "expected PreconditionViolation for an arrow matrix, got: {err:?}"
        );
    }

    /// A transversal resonator block with a nonzero *source* port diagonal
    /// `M[0,0]` (a source shunt) is not modeled by the Schur-2×2 fast path:
    /// `from_matrix` must reject it with a typed `PreconditionViolation`
    /// rather than silently dropping the shunt (round-3 review P2-2).
    #[test]
    fn source_port_diagonal_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
        use crate::matrix::CouplingMatrix;
        // Same base as `transversal_matrix_is_accepted`, but M[0,0] = 4.0
        // (a source shunt). The resonator block stays diagonal, so only the
        // port-diagonal check can catch this.
        let data = vec![
            4.0, 0.40, 0.30, 0.0, // row0 (source) — M[0,0] nonzero
            0.40, 0.10, 0.0, 0.20, // row1 (res 1)
            0.30, 0.0, -0.10, 0.15, // row2 (res 2)
            0.0, 0.20, 0.15, 0.05, // row3 (load)
        ];
        let m = CouplingMatrix::new(2, data)?;
        let err = PoleExpansionData::from_matrix(&m).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "a nonzero source port diagonal must be rejected, got: {err:?}"
        );
        let msg = format!("{err}");
        assert!(
            msg.contains("M[0][0]"),
            "error message should name the offending entry: {msg}"
        );
        Ok(())
    }

    /// The same contract applies to the *load* port diagonal `M[N+1,N+1]` —
    /// not just the source side (round-3 review P2-2, requirement: both ends).
    #[test]
    fn load_port_diagonal_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
        use crate::matrix::CouplingMatrix;
        let data = vec![
            0.0, 0.40, 0.30, 0.0, // row0 (source)
            0.40, 0.10, 0.0, 0.20, // row1 (res 1)
            0.30, 0.0, -0.10, 0.15, // row2 (res 2)
            0.0, 0.20, 0.15, 2.5, // row3 (load) — M[3][3] nonzero
        ];
        let m = CouplingMatrix::new(2, data)?;
        let err = PoleExpansionData::from_matrix(&m).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "a nonzero load port diagonal must be rejected, got: {err:?}"
        );
        let msg = format!("{err}");
        assert!(
            msg.contains("M[3][3]"),
            "error message should name the offending entry: {msg}"
        );
        Ok(())
    }

    /// Tolerance boundary for the port diagonal: a round-off level residual
    /// (1e-12, inside `1e-9 · max(1, max|M|)`) is accepted, a physical shunt
    /// above the tolerance is rejected. Same policy as the resonator block.
    #[test]
    fn port_diagonal_tolerance_boundary() -> Result<(), Box<dyn std::error::Error>> {
        use crate::matrix::CouplingMatrix;
        // Only the SOURCE port diagonal M[0,0] varies here; the load diagonal
        // is series-terminated (M[3][3] = 0) so it is never the cause of a
        // reject. tol = 1e-9 * max(1, max|M|) = 1e-9 (max|M| = 0.40).
        let tiny = vec![
            1e-12, 0.40, 0.30, 0.0, // M[0,0] = 1e-12 < tol: accepted
            0.40, 0.10, 0.0, 0.20,
            0.30, 0.0, -0.10, 0.15,
            0.0, 0.20, 0.15, 0.0, // M[3][3] = 0 (series)
        ];
        let m_tiny = CouplingMatrix::new(2, tiny)?;
        PoleExpansionData::from_matrix(&m_tiny)?; // round-off accepted

        // 1e-8 > 1e-9: a real (if small) source shunt.
        let real = vec![
            1e-8, 0.40, 0.30, 0.0, // M[0,0] = 1e-8 > tol: rejected
            0.40, 0.10, 0.0, 0.20,
            0.30, 0.0, -0.10, 0.15,
            0.0, 0.20, 0.15, 0.0, // M[3][3] = 0 (series)
        ];
        let m_real = CouplingMatrix::new(2, real)?;
        let err = PoleExpansionData::from_matrix(&m_real).unwrap_err();
        assert!(
            matches!(err, crate::error::MfsError::PreconditionViolation(_)),
            "a 1e-8 port shunt exceeds the 4e-9 tolerance and must be rejected"
        );
        Ok(())
    }

    /// The main-flow `response_normalized` on a design whose stored matrix is
    /// *folded* must still return a correct response: the structural check now
    /// rejects the O(N) fast path up front, so the call falls through to the
    /// LU backend. This pins the fallback-protection behaviour end to end.
    #[test]
    fn main_path_on_folded_design_matches_lu() -> Result<(), Box<dyn std::error::Error>> {
        let d = FilterDesign::prototype(4, 20.0)
            .zeros([1.5, -1.5])
            .synthesize()?;
        let folded = transform_matrix(d.matrix(), TopologyKind::Folded)?
            .matrix
            .clone();
        // Build a design whose stored matrix is the folded one.
        let folded_design = FilterDesign::from_polynomials_and_matrix(
            d.order(),
            20.0,
            6.75e9,
            300e6,
            d.polynomials().clone(),
            folded.clone(),
        )?;

        // Main-flow response on the folded-matrix design.
        let main = folded_design.response_normalized(-3.0, 3.0, 201)?;
        // Independent LU reference on the identical folded matrix.
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 201)?;
        let lu = ResponseSolver.evaluate_normalized(&folded, &grid)?;

        assert_eq!(main.samples.len(), lu.samples.len());
        for (a, b) in main.samples.iter().zip(lu.samples.iter()) {
            if b.s21_mag() > 1e-10 {
                let rel = ((a.s21_mag() - b.s21_mag()) / b.s21_mag()).abs();
                assert!(
                    rel < 1e-6,
                    "main path on folded design diverges from LU at ω={}: main={}, lu={}",
                    b.normalized_omega,
                    a.s21_mag(),
                    b.s21_mag()
                );
            }
        }
        Ok(())
    }
}
