use crate::error::{MfsError, Result};

use super::core::{CouplingMatrix, MatrixTopology};
use super::rotations::{safe_angle, RotationAxis};

impl CouplingMatrix {
    /// Extracts one trisection from the tail of the matrix and moves it to the requested center.
    ///
    /// `center_resonator` uses 1-based resonator numbering, excluding source and load.
    /// For example, `center_resonator = 2` targets the resonator triplet `(1, 2, 3)`.
    pub(crate) fn extract_triplet(
        &self,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<Self> {
        validate_triplet_center(self.order(), center_resonator)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "triplet transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.get(tail, tail);
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.get(tail - 1, tail), denominator)
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 1, tail, cosine, sine);

        let move_steps = order - center_resonator - 1;
        for step in 0..move_steps {
            let pivot_a = order - step - 2;
            let pivot_b = order - step - 1;
            matrix = matrix.rotate_matrix_with_indices(
                order - step,
                pivot_a,
                pivot_b,
                1.0,
                RotationAxis::Row,
            );
        }

        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Extracts two neighboring trisections and merges them into a quadruplet.
    ///
    /// `position` matches the first triplet center in 1-based resonator numbering.
    /// `common_resonator` must be `1` or `4`, matching the two elimination formulas
    /// used by the Python prototype.
    pub(crate) fn extract_quadruplet(
        &self,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<Self> {
        validate_quadruplet_position(self.order(), position)?;
        if common_resonator != 1 && common_resonator != 4 {
            return Err(MfsError::PreconditionViolation(
                "common resonator for quadruplet extraction must be 1 or 4".to_string(),
            ));
        }

        let mut matrix = if swap_zero_order {
            self.extract_triplet(second_zero, position)?
                .extract_triplet(first_zero, position + 1)?
        } else {
            self.extract_triplet(first_zero, position)?
                .extract_triplet(second_zero, position + 1)?
        };

        let tail = position + 2;
        let theta = if common_resonator == 4 {
            -safe_angle(
                matrix.get(tail - 1, tail - 3),
                matrix.get(tail - 3, tail - 2),
            )
        } else {
            safe_angle(
                matrix.get(tail, tail - 2),
                matrix.get(tail - 1, tail),
            )
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 2, tail - 1, cosine, sine);
        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Converts an arrow-style matrix into a trisection-centered topology.
    ///
    /// `zero_positions` uses 1-based resonator numbering and must span exactly
    /// one center resonator, for example `(2, 4)` to target a trisection centered
    /// on resonator `3`.
    pub(crate) fn extract_trisection(
        &self,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<Self> {
        if self.topology() != MatrixTopology::Arrow {
            return Err(MfsError::PreconditionViolation(format!(
                "trisection extraction requires Arrow input, got {:?}",
                self.topology()
            )));
        }
        validate_trisection_positions(self.order(), zero_positions)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "trisection transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.get(tail, tail);
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.get(tail - 1, tail), denominator)
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 1, tail, cosine, sine);

        let center_resonator = (zero_positions.0 + zero_positions.1) / 2;
        let pull_steps = order - 1 - center_resonator;
        for step in 0..pull_steps {
            matrix = matrix.rotate_matrix_with_indices(
                order - step,
                order - step - 2,
                order - step - 1,
                1.0,
                RotationAxis::Row,
            );
        }

        matrix.clean_small_values();
        Ok(matrix)
    }
}

fn validate_triplet_center(order: usize, center_resonator: usize) -> Result<()> {
    if center_resonator < 2 || center_resonator >= order {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "triplet center must be in [2, {}), got {center_resonator}",
            order
        )));
    }

    Ok(())
}

fn validate_quadruplet_position(order: usize, position: usize) -> Result<()> {
    if position < 2 || position + 1 >= order {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "quadruplet position must leave room for two adjacent triplets, got {position} for order {order}"
        )));
    }

    Ok(())
}

fn validate_trisection_positions(order: usize, zero_positions: (usize, usize)) -> Result<()> {
    let (start, end) = zero_positions;
    if start < 1 || end > order || start >= end {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "trisection zero positions must be ordered resonator indices within 1..={order}, got ({start}, {end})"
        )));
    }
    if end - start != 2 {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "trisection zero positions must differ by exactly 2, got ({start}, {end})"
        )));
    }

    Ok(())
}
