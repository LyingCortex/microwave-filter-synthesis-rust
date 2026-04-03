use crate::error::{MfsError, Result};
use crate::matrix::CouplingMatrix;

use super::{ResponseSample, ResponseSettings, SParameterResponse};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Complex64 {
    re: f64,
    im: f64,
}

impl Complex64 {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };
    const ONE: Self = Self { re: 1.0, im: 0.0 };

    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn from_real(value: f64) -> Self {
        Self { re: value, im: 0.0 }
    }

    fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl std::ops::Add for Complex64 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::Sub for Complex64 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for Complex64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::Mul<f64> for Complex64 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl std::ops::Div for Complex64 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.norm_sqr();
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

impl std::ops::Div<f64> for Complex64 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

impl std::ops::Neg for Complex64 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

pub(super) fn evaluate_normalized_response(
    matrix: &CouplingMatrix,
    normalized_omegas: &[f64],
    settings: ResponseSettings,
) -> Result<SParameterResponse> {
    evaluate_response(matrix, normalized_omegas, normalized_omegas, settings)
}

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

    let side = matrix.shape().rows;
    let source = settings.source_resistance;
    let load = settings.load_resistance;
    let transmission_scale = 2.0 * (source * load).sqrt();

    let samples = frequencies_hz
        .iter()
        .copied()
        .zip(normalized_omegas.iter().copied())
        .map(|(frequency_hz, omega)| {
            let inverse = invert_response_matrix(matrix, omega, settings)?;

            let s11 = Complex64::ONE + Complex64::new(0.0, 2.0 * source) * inverse[0][0];
            let s21 = Complex64::new(0.0, -transmission_scale) * inverse[side - 1][0];

            let last_row = &inverse[side - 1];
            let numerator = last_row
                .iter()
                .copied()
                .zip(inverse.iter().map(|row| row[0]))
                .fold(Complex64::ZERO, |acc, (lhs, rhs)| acc + lhs * rhs);
            let denominator = inverse[side - 1][0];
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

fn invert_response_matrix(
    matrix: &CouplingMatrix,
    omega: f64,
    settings: ResponseSettings,
) -> Result<Vec<Vec<Complex64>>> {
    let side = matrix.shape().rows;
    let mut a = vec![vec![Complex64::ZERO; side]; side];
    let mut inv = vec![vec![Complex64::ZERO; side]; side];

    for row in 0..side {
        for col in 0..side {
            let mut value = Complex64::from_real(matrix.at(row, col).unwrap_or_default());
            if row == col {
                if row != 0 && row != side - 1 {
                    value = value + Complex64::from_real(omega);
                }
                if row == 0 {
                    value = value + Complex64::new(0.0, -settings.source_resistance);
                } else if row == side - 1 {
                    value = value + Complex64::new(0.0, -settings.load_resistance);
                }
            }
            a[row][col] = value;
        }
        inv[row][row] = Complex64::ONE;
    }

    gauss_jordan_inverse(&mut a, &mut inv)?;
    Ok(inv)
}

fn gauss_jordan_inverse(a: &mut [Vec<Complex64>], inv: &mut [Vec<Complex64>]) -> Result<()> {
    let n = a.len();
    for pivot_index in 0..n {
        let pivot_row = (pivot_index..n)
            .max_by(|&lhs, &rhs| {
                a[lhs][pivot_index]
                    .norm_sqr()
                    .partial_cmp(&a[rhs][pivot_index].norm_sqr())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .ok_or_else(|| MfsError::Unsupported("empty matrix inversion".to_string()))?;

        if a[pivot_row][pivot_index].norm_sqr() <= 1e-18 {
            return Err(MfsError::Unsupported(
                "response matrix became singular during inversion".to_string(),
            ));
        }

        if pivot_row != pivot_index {
            a.swap(pivot_row, pivot_index);
            inv.swap(pivot_row, pivot_index);
        }

        let pivot = a[pivot_index][pivot_index];
        for column in 0..n {
            a[pivot_index][column] = a[pivot_index][column] / pivot;
            inv[pivot_index][column] = inv[pivot_index][column] / pivot;
        }

        for row in 0..n {
            if row == pivot_index {
                continue;
            }
            let factor = a[row][pivot_index];
            if factor.norm_sqr() <= 1e-24 {
                continue;
            }
            for column in 0..n {
                a[row][column] = a[row][column] - factor * a[pivot_index][column];
                inv[row][column] = inv[row][column] - factor * inv[pivot_index][column];
            }
        }
    }

    Ok(())
}
