//! Cascaded Trisection (CT) and Cascaded Quadruplet (CQ) synthesis.
//!
//! Implements the Cameron method (2011) for synthesizing coupling matrices
//! as cascades of trisections and quadruplets from the arrow canonical form.

use super::core::CouplingMatrix;
use crate::error::{MfsError, Result};

impl CouplingMatrix {
    /// Synthesizes a Cascaded Trisection (CT) configuration from a filter matrix.
    ///
    /// The matrix is first reduced to arrow form, which already encodes every
    /// transmission zero. A single condition-at-load-tail + peel-to-source (the
    /// verified `extract_triplet` primitive) then
    /// reveals the full trisection cascade: `k` zeros become `k` source/load-
    /// connected trisections (3-resonator cliques) linked by a connected,
    /// single-bridge mainline.
    ///
    /// `zeros[0]` is the primary zero: it conditions the load-end trisection
    /// during the peel and becomes the source-side trisection's zero (for the
    /// Cameron 2011 Fig. 21 anchor, `zeros[0]` is `omega_01 = 1.8082`). The
    /// remaining zeros are already encoded in the arrow matrix and are used
    /// only for the count validation.
    ///
    /// The former per-zero loop re-conditioned the tail on every pass and
    /// rotated on the fixed load column, which severed the mainline and left a
    /// dense folded/wheel matrix (review P2-4).
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
                zeros.len(),
                order,
                order - 2
            )));
        }

        let arrow = self.transform_topology(super::core::MatrixTopology::Arrow)?;
        arrow.extract_triplet(zeros[0], 1)
    }

    /// Synthesizes a Cascaded Quadruplet (CQ) configuration.
    ///
    /// `zero_pairs` contains pairs of transmission zeros, each realized as one
    /// quadruplet (a 4-resonator section). The quadruplets are cascaded from
    /// the load end toward the source, linked by a single-bridge mainline.
    ///
    /// The merge is the verified `extract_quadruplet` tail-pivot Givens
    /// rotation. The former `diagonal_rotation_angle` merge was a no-op for
    /// symmetric pairs (the rotation vanishes), which left the CQ
    /// byte-identical to the CT (review P1-5).
    pub fn to_cascaded_quadruplets(&self, zero_pairs: &[(f64, f64)]) -> Result<Self> {
        if zero_pairs.is_empty() {
            return Err(MfsError::PreconditionViolation(
                "at least one zero pair is required for CQ synthesis".into(),
            ));
        }

        let order = self.order();
        let mut matrix = self.transform_topology(super::core::MatrixTopology::Arrow)?;

        // Peel quadruplets from the load end toward the source so adjacent
        // quartets stay connected by the single-bridge mainline. The verified
        // `extract_cascade_quadruplets` extracts the highest-position quartet
        // first; `common_resonator = 4` selects the load-end elimination
        // formula (the order-4 case is common-invariant).
        for (idx, &(z1, z2)) in zero_pairs.iter().enumerate() {
            let position = order - 2 * (idx + 1);
            if position < 2 {
                break;
            }
            matrix = matrix.extract_quadruplet(z1, z2, position, 4, false)?;
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

        let ct = design
            .matrix()
            .to_cascaded_trisections(&[-2.0, -1.2, 1.5])?;
        assert_eq!(ct.order(), 6);

        // Verify response is preserved (similarity transform)
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig_resp = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let ct_resp = ResponseSolver.evaluate_normalized(&ct, &grid)?;

        for (o, c) in orig_resp.samples.iter().zip(ct_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(
                diff < 1e-6,
                "CT response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega
            );
        }
        // P2-4 gate: the resonator mainline must stay connected. The former
        // per-zero loop re-conditioned the tail on each pass and rotated on
        // the fixed load column, severing the mainline (S-2 = 0, 3-L = 0).
        for i in 1..6 {
            assert!(
                ct.get(i, i + 1).abs() > 1e-3,
                "CT mainline broken at resonator pair ({i}, {}): {}",
                i + 1,
                ct.get(i, i + 1)
            );
        }
        Ok(())
    }

    #[test]
    fn cq_synthesis_with_symmetric_zeros() -> Result<()> {
        let design = FilterDesign::prototype(6, 23.0)
            .zeros([-1.5, 1.5, -2.0, 2.0])
            .synthesize()?;

        let cq = design
            .matrix()
            .to_cascaded_quadruplets(&[(-1.5, 1.5), (-2.0, 2.0)])?;
        assert_eq!(cq.order(), 6);

        // Verify response preservation
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig_resp = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let cq_resp = ResponseSolver.evaluate_normalized(&cq, &grid)?;

        for (o, c) in orig_resp.samples.iter().zip(cq_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(
                diff < 1e-4,
                "CQ response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega
            );
        }
        // P1-5 gate: the CQ must be structurally distinct from the CT. The
        // former `diagonal_rotation_angle` merge was a no-op for symmetric
        // pairs, which left the CQ byte-identical to the CT (dM = 0). The
        // `extract_quadruplet` merge produces a different realization.
        let ct = design
            .matrix()
            .to_cascaded_trisections(&[-1.5, 1.5, -2.0, 2.0])?;
        let dmc: f64 = (0..8usize)
            .flat_map(|i| (0..8).map(move |j| (i, j)))
            .map(|(i, j)| (cq.get(i, j) - ct.get(i, j)).abs())
            .fold(0.0, f64::max);
        assert!(
            dmc > 0.01,
            "P1-5: CQ is not distinct from CT (dM = {dmc:.2e}); the quadruplet merge is a no-op"
        );

        // Mainline must stay connected across the two quartets.
        for i in 1..6 {
            assert!(
                cq.get(i, i + 1).abs() > 1e-3,
                "CQ mainline broken at ({i}, {}): {}",
                i + 1,
                cq.get(i, i + 1)
            );
        }
        Ok(())
    }

    #[test]
    fn ct_synthesis_reproduces_fig21_cameron2011() -> Result<()> {
        // Cameron 2011 Fig. 21: order 4, 22 dB, TZ at 1.8082 (source side)
        // and 1.3217 (load side). A single condition-at-load-tail +
        // peel-to-source (`extract_triplet`) must reproduce the paper's
        // printed matrix: source S connects to resonators 1 and 2, load L
        // to 3 and 4, with the (2,3) bridge linking the two trisections.
        //
        // This is the P2-4 acceptance: the peel produces the edge-connected
        // trisections, not the dense/wheel matrix the old per-zero loop made.
        let design = FilterDesign::prototype(4, 22.0)
            .zeros([1.8082, 1.3217])
            .synthesize()?;
        let ct = design
            .matrix()
            .to_cascaded_trisections(&[1.8082, 1.3217])?;
        assert_eq!(ct.order(), 4);

        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let ct_resp = ResponseSolver.evaluate_normalized(&ct, &grid)?;
        for (o, c) in orig.samples.iter().zip(ct_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(
                diff < 1e-6,
                "Fig21 CT response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega
            );
        }

        // The trisection "trident" sparsity pattern.
        let want_nonzero = [(0usize, 1usize), (0, 2), (1, 2), (2, 3), (3, 4), (3, 5), (4, 5)];
        for (i, j) in want_nonzero {
            assert!(
                ct.get(i, j).abs() > 0.05,
                "Fig21 CT: expected nonzero coupling ({i},{j}), got {}",
                ct.get(i, j)
            );
        }
        let want_zero = [(0usize, 3usize), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 4), (2, 5)];
        for (i, j) in want_zero {
            assert!(
                ct.get(i, j).abs() < 1e-3,
                "Fig21 CT: expected zero coupling ({i},{j}), got {}",
                ct.get(i, j)
            );
        }

        // Value anchor: the printed Fig 21 matrix. Off-diagonals are
        // sign-insensitive (Givens rotations fix magnitudes, not signs); the
        // 4-decimal print quantum plus 22 dB synthesis error bounds the
        // tolerance at 5e-3 (measured deviation is ~1e-4).
        let fig21: [[f64; 6]; 6] = [
            [0.0, 0.9221, 0.5921, 0.0, 0.0, 0.0],
            [0.9221, -0.8333, 0.6259, 0.0, 0.0, 0.0],
            [0.5921, 0.6259, 0.6022, 0.8634, 0.0, 0.0],
            [0.0, 0.0, 0.8634, 0.5812, 0.2590, 0.8471],
            [0.0, 0.0, 0.0, 0.2590, -1.1092, 0.6952],
            [0.0, 0.0, 0.0, 0.8471, 0.6952, 0.0],
        ];
        for (i, row) in fig21.iter().enumerate() {
            for (j, expected) in row.iter().enumerate() {
                let d = (ct.get(i, j).abs() - expected.abs()).abs();
                assert!(
                    d < 5e-3,
                    "Fig21 CT value mismatch at ({i},{j}): got {}, expected {}",
                    ct.get(i, j),
                    expected
                );
            }
        }
        Ok(())
    }

    #[test]
    fn cq_synthesis_single_quadruplet_structure() -> Result<()> {
        // Order 4 with one symmetric zero pair: a single quadruplet. The
        // P1-5 acceptance: source connects to resonator 1 only, load to
        // resonator 4 only, the mainline 1-2-3-4 is connected, and the
        // quadruplet's signature 2-skip coupling M[1,4] is present while the
        // resonator diagonals are zeroed. This is structurally distinct from
        // the CT (which connects source to both 1 and 2).
        let design = FilterDesign::prototype(4, 23.0)
            .zeros([-1.5, 1.5])
            .synthesize()?;
        let cq = design
            .matrix()
            .to_cascaded_quadruplets(&[(-1.5, 1.5)])?;
        assert_eq!(cq.order(), 4);

        let grid = FrequencyGrid::linspace(-3.0, 3.0, 21)?;
        let orig = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;
        let cq_resp = ResponseSolver.evaluate_normalized(&cq, &grid)?;
        for (o, c) in orig.samples.iter().zip(cq_resp.samples.iter()) {
            let diff = (o.s21_mag() - c.s21_mag()).abs();
            assert!(
                diff < 1e-6,
                "CQ response mismatch at ω={:.2}: diff={diff:.2e}",
                o.normalized_omega
            );
        }

        // Source to resonator 1 only; load to resonator 4 only.
        assert!(cq.get(0, 1).abs() > 0.05, "CQ: source-res1 missing: {}", cq.get(0, 1));
        assert!(cq.get(0, 2).abs() < 1e-3, "CQ: source-res2 should be 0: {}", cq.get(0, 2));
        assert!(cq.get(3, 5).abs() < 1e-3, "CQ: load-res3 should be 0: {}", cq.get(3, 5));
        assert!(cq.get(4, 5).abs() > 0.05, "CQ: load-res4 missing: {}", cq.get(4, 5));

        // Mainline connected.
        for i in 1..4 {
            assert!(
                cq.get(i, i + 1).abs() > 0.05,
                "CQ mainline broken at ({i}, {}): {}",
                i + 1,
                cq.get(i, i + 1)
            );
        }

        // The quadruplet 2-skip coupling is present; no other 2-skip.
        assert!(cq.get(1, 4).abs() > 0.05, "CQ: quadruplet 2-skip M[1,4] missing: {}", cq.get(1, 4));
        assert!(cq.get(1, 3).abs() < 1e-3, "CQ: unexpected 2-skip M[1,3]: {}", cq.get(1, 3));
        assert!(cq.get(2, 4).abs() < 1e-3, "CQ: unexpected 2-skip M[2,4]: {}", cq.get(2, 4));

        // Resonator diagonals zeroed by the peel.
        for i in 1..5 {
            assert!(
                cq.get(i, i).abs() < 1e-3,
                "CQ: resonator {} diagonal not zeroed: {}",
                i,
                cq.get(i, i)
            );
        }
        Ok(())
    }
}
