//! Cascaded Trisection (CT) and Cascaded Quadruplet (CQ) synthesis.
//!
//! Implements the Cameron method (2011) for synthesizing coupling matrices
//! as cascades of trisections and quadruplets from the arrow canonical form.

use crate::error::{MfsError, Result};
use super::core::CouplingMatrix;
use super::rotations::{safe_angle, diagonal_rotation_angle};

impl CouplingMatrix {
    /// Synthesizes a Cascaded Trisection (CT) configuration from an arrow matrix.
    ///
    /// Each transmission zero is realized by one trisection. The trisections are
    /// extracted one by one from the arrow matrix and pulled to their target positions.
    ///
    /// `zeros` are the transmission zero frequencies in the order they should be
    /// assigned to trisections (from load end towards source end).
    ///
    /// Returns the CT coupling matrix with trisection topology.
    pub fn to_cascaded_trisections(&self, zeros: &[f64]) -> Result<Self> {
        if zeros.is_empty() {
            return Err(MfsError::PreconditionViolation(
                "at least one transmission zero is required for CT synthesis".into(),
            ));
        }

        let order = self.order();
        if zeros.len() > order - 2 {
            return Err(MfsError::PreconditionViolation(format!(
                "too many zeros ({}) for order {} (max {})",
                zeros.len(), order, order - 2
            )));
        }

        // Start from arrow form
        let mut matrix = self.transform_topology(super::core::MatrixTopology::Arrow)?;
        let side = matrix.side();

        // For each transmission zero, create and position a trisection
        for (tz_idx, &tz) in zeros.iter().enumerate() {
            let n = order; // current effective size (resonator indices 1..=N)

            // Conditioning rotation at pivot [N-1, N] with angle from (13) in Cameron 2011:
            //   θ₀₁ = atan(M_{N-1,N} / (ω₀₁ + M_{N,N}))
            let tail = n; // last resonator index
            let denominator = tz + matrix.get(tail, tail);
            let theta = if denominator.abs() < 1e-10 {
                std::f64::consts::FRAC_PI_2
            } else {
                safe_angle(matrix.get(tail - 1, tail), denominator)
            };
            matrix.apply_givens_similarity_inplace(tail - 1, tail, theta.cos(), theta.sin());

            // Pull the trisection up the diagonal to its target position
            // Target: trisection should end up at position (2*tz_idx+1, 2*tz_idx+2)
            // from the source end. We pull from bottom towards the target.
            let target_pos = order - 1 - tz_idx; // target middle resonator of trisection
            let current_pos = n - 1; // currently at N-1

            for pos in (target_pos..current_pos).rev() {
                // Rotation to pull trisection from pos+1 to pos
                let pivot_a = pos;
                let pivot_b = pos + 1;
                let numerator = matrix.get(pivot_a, side - 1); // element in load column
                let denom_val = matrix.get(pivot_b, side - 1);

                if denom_val.abs() < 1e-15 && numerator.abs() < 1e-15 {
                    continue;
                }

                let angle = safe_angle(numerator, denom_val);
                matrix.apply_givens_similarity_inplace(pivot_a, pivot_b, angle.cos(), angle.sin());
            }
        }

        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Synthesizes a Cascaded Quadruplet (CQ) configuration by merging trisection pairs.
    ///
    /// `zero_pairs` contains pairs of transmission zeros that should be merged into
    /// quartets. Each pair creates one quadruplet section.
    ///
    /// For symmetric zero pairs (±ω₀), the diagonal coupling in the quartet will be zero.
    pub fn to_cascaded_quadruplets(&self, zero_pairs: &[(f64, f64)]) -> Result<Self> {
        if zero_pairs.is_empty() {
            return Err(MfsError::PreconditionViolation(
                "at least one zero pair is required for CQ synthesis".into(),
            ));
        }

        // First create all trisections
        let all_zeros: Vec<f64> = zero_pairs.iter()
            .flat_map(|&(z1, z2)| vec![z1, z2])
            .collect();

        let mut matrix = self.to_cascaded_trisections(&all_zeros)?;
        let order = matrix.order();

        // Merge adjacent trisection pairs into quartets
        // Each pair of adjacent trisections (2k-1, 2k) and (2k+1, 2k+2) is merged
        // by a cross-pivot rotation that eliminates one main-line coupling
        for (pair_idx, &(z1, z2)) in zero_pairs.iter().enumerate() {
            // The two trisections for this pair are at positions starting from the
            // load end. Merge them using the quadruplet extraction we already have.
            let position = order - 2 * (pair_idx + 1);

            if position < 2 || position + 1 >= order {
                continue; // Skip if position is out of range
            }

            // Cross-pivot rotation to merge: annihilate the coupling between
            // the two trisections, creating the quartet diagonal coupling
            let pivot_a = position;
            let pivot_b = position + 1;
            let angle = diagonal_rotation_angle(&matrix, pivot_a, pivot_b);
            matrix.apply_givens_similarity_inplace(pivot_a, pivot_b, angle.cos(), angle.sin());
        }

        matrix.clean_small_values();
        Ok(matrix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::FilterDesign;
    use crate::freq::FrequencyGrid;
    use crate::response::ResponseSolver;

    #[test]
    fn ct_synthesis_produces_valid_matrix() -> Result<()> {
        let design = FilterDesign::prototype(6, 23.0)
            .zeros([-2.0, -1.2, 1.5])
            .synthesize()?;

        let ct = design.matrix().to_cascaded_trisections(&[-2.0, -1.2, 1.5])?;
        assert_eq!(ct.order(), 6);

        // Verify response is preserved (similarity transform)
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig_resp = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let ct_resp = ResponseSolver.evaluate_normalized(&ct, &grid)?;

        for (o, c) in orig_resp.samples.iter().zip(ct_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(diff < 1e-6, "CT response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega);
        }
        Ok(())
    }

    #[test]
    fn cq_synthesis_with_symmetric_zeros() -> Result<()> {
        let design = FilterDesign::prototype(6, 23.0)
            .zeros([-1.5, 1.5, -2.0, 2.0])
            .synthesize()?;

        let cq = design.matrix().to_cascaded_quadruplets(&[(-1.5, 1.5), (-2.0, 2.0)])?;
        assert_eq!(cq.order(), 6);

        // Verify response preservation
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig_resp = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let cq_resp = ResponseSolver.evaluate_normalized(&cq, &grid)?;

        for (o, c) in orig_resp.samples.iter().zip(cq_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(diff < 1e-4, "CQ response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega);
        }
        Ok(())
    }
}
