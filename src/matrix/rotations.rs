use crate::error::Result;

use super::core::CouplingMatrix;

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub(crate) enum RotationAxis {
    Row,
    Column,
}

#[allow(dead_code)]
pub(crate) fn rotation_matrix_basic(order: usize, theta: f64, row: usize, col: usize) -> Result<CouplingMatrix> {
    let mut rotation = CouplingMatrix::identity(order)?;
    let cosine = theta.cos();
    let sine = theta.sin();
    rotation.set_entry(row, row, cosine);
    rotation.set_entry(col, col, cosine);
    rotation.set_entry(row, col, -sine);
    rotation.set_entry(col, row, sine);
    Ok(rotation)
}

pub(crate) fn safe_angle(y: f64, x: f64) -> f64 {
    if x.abs() < 1e-10 {
        if y.is_sign_positive() {
            -std::f64::consts::FRAC_PI_2
        } else {
            std::f64::consts::FRAC_PI_2
        }
    } else {
        (y / x).atan()
    }
}

pub(crate) fn diagonal_rotation_angle(matrix: &CouplingMatrix, index_a: usize, index_b: usize) -> f64 {
    let diagonal_delta = matrix.get(index_b, index_b)
        - matrix.get(index_a, index_a);
    if diagonal_delta.abs() < 1e-10 {
        0.0
    } else {
        let ratio = (2.0 * matrix.get(index_a, index_b)) / diagonal_delta;
        0.5 * safe_angle(ratio, 1.0)
    }
}
