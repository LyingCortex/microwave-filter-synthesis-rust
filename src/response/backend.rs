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

    let samples = frequencies_hz
        .iter()
        .copied()
        .zip(normalized_omegas.iter().copied())
        .map(|(frequency_hz, omega)| {
            let inverse = solve_inverse(matrix, omega, settings)?;

            let s11 = Complex64::new(1.0, 0.0) + Complex64::new(0.0, 2.0 * source) * inverse[(0, 0)];
            let s21 = Complex64::new(0.0, -transmission_scale) * inverse[(side - 1, 0)];

            let numerator = (0..side).fold(Complex64::new(0.0, 0.0), |acc, index| {
                acc + inverse[(side - 1, index)] * inverse[(index, 0)]
            });
            let denominator = inverse[(side - 1, 0)];
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

/// Builds the complex response matrix for one frequency sample.
fn build_response_matrix(
    matrix: &CouplingMatrix,
    omega: f64,
    settings: ResponseSettings,
) -> DMatrix<Complex64> {
    let side = matrix.side();
    let mut response = matrix.to_complex_dense();

    for index in 0..side {
        if index != 0 && index != side - 1 {
            // Resonator diagonals are shifted by the normalized frequency sample.
            response[(index, index)] += Complex64::new(omega, 0.0);
        }
    }

    response[(0, 0)] += Complex64::new(0.0, -settings.source_resistance);
    response[(side - 1, side - 1)] += Complex64::new(0.0, -settings.load_resistance);
    response
}

/// Solves the matrix inverse through the backend linear solver.
fn solve_inverse(
    matrix: &CouplingMatrix,
    omega: f64,
    settings: ResponseSettings,
) -> Result<DMatrix<Complex64>> {
    let side = matrix.side();
    let response = build_response_matrix(matrix, omega, settings);
    let lu = response.lu();
    let identity = DMatrix::from_diagonal(&DVector::from_element(side, Complex64::new(1.0, 0.0)));

    lu.solve(&identity).ok_or_else(|| {
        MfsError::Unsupported("response matrix became singular during solve".to_string())
    })
}
