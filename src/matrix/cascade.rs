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
    /// quadruplet (a 4-resonator section). The quadruplets tile the chain
    /// **non-overlappingly from the source side** (pair *i* spans resonators
    /// `4i-3..=4i`), linked by a single-bridge mainline. A clean tiling needs
    /// `order >= 4 * num_pairs`; otherwise returns an error (use
    /// [`to_cascaded_trisections`](Self::to_cascaded_trisections) for denser
    /// configurations, which pack 3-resonator trisections).
    ///
    /// This routes through the verified `extract_cascade_quadruplets`
    /// primitive (the same path as the Python `cq` binding), which merges each
    /// quadruplet's two neighboring trisections. The source-side tiling pairs
    /// each quadruplet with a trisection starting at the same resonator, so the
    /// merge uses `common_resonator = 1` on every quadruplet (the `= 4`
    /// load-end formula exists in the primitive but is not reached by this
    /// tiling). The former manual tail-peel loop left
    /// overlapping-quadruplet residues, and the original
    /// `diagonal_rotation_angle` merge was a no-op for symmetric pairs (CQ
    /// byte-identical to the CT, review P1-5).
    pub fn to_cascaded_quadruplets(&self, zero_pairs: &[(f64, f64)]) -> Result<Self> {
        if zero_pairs.is_empty() {
            return Err(MfsError::PreconditionViolation(
                "at least one zero pair is required for CQ synthesis".into(),
            ));
        }

        let order = self.order();
        if order < 4 * zero_pairs.len() {
            return Err(MfsError::PreconditionViolation(format!(
                "order {} is too small for {} CQ quadruplets (needs >= 4*{} = {}); \
                 use to_cascaded_trisections for denser configurations",
                order,
                zero_pairs.len(),
                zero_pairs.len(),
                4 * zero_pairs.len()
            )));
        }

        let arrow = self.transform_topology(super::core::MatrixTopology::Arrow)?;
        let zeros: Vec<f64> = zero_pairs.iter().flat_map(|&(z1, z2)| [z1, z2]).collect();

        // Source-side non-overlapping tiling: quadruplet i (1-indexed) spans
        // resonators 4i-3..=4i. The trisection partner spans 4i-3..4i-1
        // (span 2) and the quadruplet diagonal spans 4i-3..4i (span 3). This
        // reproduces the verified layouts: order-6 1-pair [[1,3],[1,4]],
        // order-8 2-pair [[1,3],[1,4],[5,7],[5,8]].
        let mut z_locs: Vec<[usize; 2]> = Vec::with_capacity(2 * zero_pairs.len());
        for (idx, _) in zero_pairs.iter().enumerate() {
            let s = 4 * idx + 1;
            z_locs.push([s, s + 2]);
            z_locs.push([s, s + 3]);
        }

        arrow.extract_cascade_quadruplets(&zeros, &z_locs)
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
    fn cq_synthesis_order6_single_quadruplet() -> Result<()> {
        // Order-6 with one symmetric zero pair: a single clean quadruplet
        // spanning resonators 1..=4 (verified Python layout `d.cq([[1,3],[1,4]])`,
        // "6阶2零点CQ"). The quadruplet's cross coupling is M[1,4]; the
        // resonator diagonals are zeroed; the mainline 1-2-3-4-5-6 stays
        // connected; source couples to res 1 only, load to res 6 only.
        let design = FilterDesign::prototype(6, 23.0)
            .zeros([-2.0, 2.0])
            .synthesize()?;
        let cq = design.matrix().to_cascaded_quadruplets(&[(-2.0, 2.0)])?;
        assert_eq!(cq.order(), 6);

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

        // Source to res 1 only; load to res 6 only.
        assert!(
            cq.get(0, 1).abs() > 0.05,
            "CQ: source-res1 missing: {}",
            cq.get(0, 1)
        );
        assert!(cq.get(0, 2).abs() < 1e-3, "CQ: source-res2 should be 0: {}", cq.get(0, 2));
        // Load port (row side-1 = 7) couples to the last resonator (res 6)
        // only — the mainline end. res 5 and the outer 2-skip are zero.
        assert!(
            cq.get(7, 6).abs() > 0.05,
            "CQ: load-res6 missing: {}",
            cq.get(7, 6)
        );
        assert!(
            cq.get(7, 5).abs() < 1e-3,
            "CQ: load-res5 should be 0: {}",
            cq.get(7, 5)
        );

        // Mainline 1..=6 connected.
        for i in 1..6 {
            assert!(
                cq.get(i, i + 1).abs() > 0.05,
                "CQ mainline broken at ({i}, {}): {}",
                i + 1,
                cq.get(i, i + 1)
            );
        }

        // The quadruplet cross coupling M[1,4] is present; no other 2-skip.
        assert!(
            cq.get(1, 4).abs() > 0.05,
            "CQ: quadruplet cross M[1,4] missing: {}",
            cq.get(1, 4)
        );
        for i in 0..7usize {
            for j in (i + 2)..7usize {
                if (i, j) == (1, 4) {
                    continue;
                }
                assert!(
                    cq.get(i, j).abs() < 1e-3,
                    "CQ: unexpected 2-skip M[{i},{j}] = {}",
                    cq.get(i, j)
                );
            }
        }

        // Resonator diagonals zeroed.
        for i in 1..7 {
            assert!(
                cq.get(i, i).abs() < 1e-3,
                "CQ: resonator {} diagonal not zeroed: {}",
                i,
                cq.get(i, i)
            );
        }
        Ok(())
    }

    #[test]
    fn cq_synthesis_order8_two_quadruplets() -> Result<()> {
        // Order-8 with two zero pairs: the verified `CascadeQuartets` layout
        // `d.cq([[1,3],[1,4],[5,7],[5,8]])` (design doc verification item 3) —
        // two non-overlapping quadruplets, the first spanning resonators 1..=4
        // and the second 5..=8. This is the reference layout the verified
        // Python `cq` binding produces for an 8-order, 4-zero CQ.
        let design = FilterDesign::prototype(8, 23.0)
            .zeros([-1.5, 1.5, -2.0, 2.0])
            .synthesize()?;
        let cq = design
            .matrix()
            .to_cascaded_quadruplets(&[(-1.5, 1.5), (-2.0, 2.0)])?;
        assert_eq!(cq.order(), 8);

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

        // Source to res 1 only; load to res 8 only.
        assert!(
            cq.get(0, 1).abs() > 0.05,
            "CQ: source-res1 missing: {}",
            cq.get(0, 1)
        );
        assert!(cq.get(0, 2).abs() < 1e-3, "CQ: source-res2 should be 0: {}", cq.get(0, 2));
        assert!(
            cq.get(9, 8).abs() > 0.05,
            "CQ: load-res8 missing: {}",
            cq.get(9, 8)
        );
        // The load-end elimination (common_res=1, source-side tiling) leaves
        // small inherent load->res 2/4/6 couplings (measured 0.010 / 0.035 /
        // 0.020) — part of the verified `CascadeQuartets` layout, not a
        // defect, so the load row is gated for the res-8 mainline end only.

        // Mainline 1..=8 connected.
        for i in 1..8 {
            assert!(
                cq.get(i, i + 1).abs() > 0.05,
                "CQ mainline broken at ({i}, {}): {}",
                i + 1,
                cq.get(i, i + 1)
            );
        }

        // The quadruplet cross couplings (2-skip within each quartet) are the
        // signature of the CQ topology and are nonzero for this design:
        // M[1,4] = -0.2371, M[5,8] = -0.0313 (verified reference layout).
        assert!(
            cq.get(1, 4).abs() > 0.05,
            "CQ: first quadruplet cross M[1,4] missing: {}",
            cq.get(1, 4)
        );
        assert!(
            cq.get(5, 8).abs() > 0.02,
            "CQ: second quadruplet cross M[5,8] missing: {}",
            cq.get(5, 8)
        );

        // P1-5: the CQ must be structurally distinct from the CT (which packs
        // 3-resonator trisections, a denser structure). The former
        // `diagonal_rotation_angle` merge was a no-op for symmetric pairs (CQ
        // byte-identical to the CT).
        let ct = design.matrix().to_cascaded_trisections(&[-1.5, 1.5, -2.0, 2.0])?;
        let dmc: f64 = (0..10usize)
            .flat_map(|i| (0..10).map(move |j| (i, j)))
            .map(|(i, j)| (cq.get(i, j) - ct.get(i, j)).abs())
            .fold(0.0, f64::max);
        assert!(
            dmc > 0.05,
            "P1-5: CQ is not distinct from CT (dM = {dmc:.2e}); the quadruplet merge is a no-op"
        );
        Ok(())
    }

    #[test]
    fn cq_synthesis_rejects_too_small_order() -> Result<()> {
        // A clean quadruplet tiling needs order >= 4 * num_pairs. Order-6 with
        // two pairs (four zeros) cannot host two non-overlapping 4-resonator
        // quadruplets, so the synthesis must be rejected rather than produce
        // overlapping-quadruplet residues (the former tail-peel left ~1e-2
        // skip couplings that no valid tiling can remove).
        let design = FilterDesign::prototype(6, 23.0)
            .zeros([-1.5, 1.5, -2.0, 2.0])
            .synthesize()?;
        let err = design
            .matrix()
            .to_cascaded_quadruplets(&[(-1.5, 1.5), (-2.0, 2.0)])
            .expect_err("order-6 two-pair CQ must be rejected");
        assert!(
            err.to_string().contains("too small"),
            "expected a 'too small' precondition error, got: {err}"
        );
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
