use crate::error::{MfsError, Result};
use crate::matrix::CouplingMatrix;

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

use super::{ResponseSample, ResponseSettings, SParameterResponse};

/// Evaluates a response when the supplied grid is already normalized.
pub(super) fn evaluate_normalized_response(
    matrix: &CouplingMatrix,
    normalized_omegas: &[f64],
    settings: ResponseSettings,
) -> Result<SParameterResponse> {
    evaluate_response(matrix, normalized_omegas, normalized_omegas, settings)
}

/// Evaluates the response on paired physical and normalized frequency axes.
pub(super) fn evaluate_response(
    matrix: &CouplingMatrix,
    frequencies_hz: &[f64],
    normalized_omegas: &[f64],
    settings: ResponseSettings,
) -> Result<SParameterResponse> {
    validate_settings(settings)?;
    if frequencies_hz.len() != normalized_omegas.len() {
        return Err(MfsError::DimensionMismatch {
            expected: frequencies_hz.len(),
            actual: normalized_omegas.len(),
        });
    }

    let side = matrix.side();
    let source = settings.source_resistance;
    let load = settings.load_resistance;
    let transmission_scale = 2.0 * (source * load).sqrt();

    // Pre-compute the base complex matrix (frequency-independent part).
    // Only the resonator diagonals change with omega, so we build the base once
    // and add omega*I to the resonator block at each frequency point.
    let base_matrix = build_base_matrix(matrix, settings);

    let samples = frequencies_hz
        .iter()
        .copied()
        .zip(normalized_omegas.iter().copied())
        .map(|(frequency_hz, omega)| {
            // For S-parameter extraction we only need:
            //   - inverse[(0, 0)]       for S11
            //   - inverse[(side-1, 0)]  for S21
            //   - sum_k inverse[(side-1, k)] * inverse[(k, 0)] for group delay
            //
            // We solve two linear systems instead of computing the full inverse:
            //   A * x = e_0  (first column of inverse)
            //   A^T * y = e_{N-1}  (last row of inverse, transposed)
            // Then group_delay uses dot(y, x).
            let response = shifted_matrix(&base_matrix, omega, side);
            let (col_first, row_last) = solve_s_parameter_columns(&response, side)?;

            let s11 = Complex64::new(1.0, 0.0)
                + Complex64::new(0.0, 2.0 * source) * col_first[0];
            let s21 = Complex64::new(0.0, -transmission_scale) * col_first[side - 1];

            // Group delay: Im(sum_k y[k] * x[k] / x[side-1])
            let numerator = (0..side).fold(Complex64::new(0.0, 0.0), |acc, index| {
                acc + row_last[index] * col_first[index]
            });
            let denominator = col_first[side - 1];
            let group_delay = if denominator.norm_sqr() <= 1e-18 {
                0.0
            } else {
                (numerator / denominator).im
            };

            Ok(ResponseSample {
                frequency_hz,
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

/// Validates normalized source and load terminations before solving.
fn validate_settings(settings: ResponseSettings) -> Result<()> {
    if !settings.source_resistance.is_finite() || settings.source_resistance <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "source resistance must be > 0, got {}",
            settings.source_resistance
        )));
    }
    if !settings.load_resistance.is_finite() || settings.load_resistance <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "load resistance must be > 0, got {}",
            settings.load_resistance
        )));
    }
    Ok(())
}

/// Builds the frequency-independent base matrix (coupling matrix + port terminations + loss).
/// The resonator diagonal shift (omega) is applied separately per frequency point.
fn build_base_matrix(
    matrix: &CouplingMatrix,
    settings: ResponseSettings,
) -> DMatrix<Complex64> {
    let side = matrix.side();
    let mut response = matrix.to_complex_dense();

    // Add port terminations (frequency-independent imaginary parts)
    response[(0, 0)] += Complex64::new(0.0, -settings.source_resistance);
    response[(side - 1, side - 1)] += Complex64::new(0.0, -settings.load_resistance);

    // Add dissipation loss to resonator diagonals: delta = j/Qu (imaginary part)
    // This models the finite unloaded Q of each resonator.
    // The response matrix becomes: A(ω) = (jω + j/Qu)I_r + M - jR_S·e₀e₀ᵀ - jR_L·eₙeₙᵀ
    // The j/Qu term shifts the poles off the imaginary axis, causing dissipation.
    if settings.unloaded_q.is_finite() && settings.unloaded_q > 0.0 {
        let dissipation = 1.0 / settings.unloaded_q;
        for index in 1..(side - 1) {
            response[(index, index)] += Complex64::new(0.0, -dissipation);
        }
    }

    response
}

/// Returns the response matrix shifted by omega on the resonator diagonals.
/// This avoids cloning the full matrix — we modify in place and could restore,
/// but since DMatrix is cheap to clone for small sizes (N+2 ≤ 32), we clone.
#[inline]
fn shifted_matrix(
    base: &DMatrix<Complex64>,
    omega: f64,
    side: usize,
) -> DMatrix<Complex64> {
    let mut response = base.clone();
    let omega_c = Complex64::new(omega, 0.0);
    for index in 1..(side - 1) {
        response[(index, index)] += omega_c;
    }
    response
}

/// Solves for the first column and last row of the inverse using LU decomposition.
///
/// For a symmetric coupling matrix M, the response matrix A(ω) = M + jωI_r - jR terms
/// is also symmetric (A = A^T), because:
///   - M is symmetric (coupling matrix property)
///   - jωI_r is diagonal (symmetric)
///   - port termination terms are diagonal (symmetric)
///
/// Therefore we can solve both systems with a single LU factorization:
///   A * x = e_0       → x is the first column of A^{-1}
///   A * y = e_{N-1}   → y is the last column of A^{-1}
///
/// Since A = A^T, the last row of A^{-1} equals the last column of A^{-1},
/// so y directly gives us the needed row.
///
/// This halves the LU decomposition cost compared to factoring both A and A^T.
fn solve_s_parameter_columns(
    response: &DMatrix<Complex64>,
    side: usize,
) -> Result<(DVector<Complex64>, DVector<Complex64>)> {
    let lu = response.clone().lu();

    // Solve A * x = e_0 (first column of inverse)
    let mut rhs_first = DVector::zeros(side);
    rhs_first[0] = Complex64::new(1.0, 0.0);
    let col_first = lu.solve(&rhs_first).ok_or_else(|| {
        MfsError::NumericalFailure("response matrix became singular during solve".to_string())
    })?;

    // Solve A * y = e_{N-1} (last column of inverse = last row, since A is symmetric)
    let mut rhs_last = DVector::zeros(side);
    rhs_last[side - 1] = Complex64::new(1.0, 0.0);
    let row_last = lu.solve(&rhs_last).ok_or_else(|| {
        MfsError::NumericalFailure("response matrix became singular during solve (last column)".to_string())
    })?;

    Ok((col_first, row_last))
}
