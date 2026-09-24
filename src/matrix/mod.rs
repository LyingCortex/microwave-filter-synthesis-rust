//! Dense coupling-matrix data structures and low-level helpers.
//!
//! This module remains the home of the core [`CouplingMatrix`] object and
//! matrix-oriented utilities such as builders, indexing, and low-level
//! operations.
//!
//! Advanced section extraction is now routed through [`crate::transform`] and
//! [`crate::synthesis`] facades instead of being part of the intended public
//! `matrix` surface.
//!
//! For new high-level workflows, prefer:
//!
//! - [`crate::synthesis::MatrixSynthesisEngine`]
//! - [`crate::synthesis::SectionSynthesis`]
//! - [`crate::transform::transform_matrix`]
//! - [`crate::verify`] helpers for structural and response checks

mod builder;
mod cascade;
mod core;
pub(crate) mod rotations;
mod scaling;
mod sections;

pub use builder::CouplingMatrixBuilder;
pub use core::{BandPassScaledCouplingMatrix, CouplingMatrix, MatrixShape, MatrixTopology};

#[cfg(test)]
mod tests {
    use crate::error::Result;
    use crate::freq::BandPassMapping;

    use super::{CouplingMatrixBuilder, MatrixTopology};

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn builder_can_set_symmetric_entries() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.75)?
            .build()?;

        approx_eq(matrix.at(0, 1).unwrap_or_default(), 0.75, 1e-12);
        approx_eq(matrix.at(1, 0).unwrap_or_default(), 0.75, 1e-12);
        assert_eq!(matrix.topology(), MatrixTopology::Transversal);
        Ok(())
    }

    #[test]
    fn bandpass_scaling_round_trips_internal_couplings() -> Result<()> {
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let matrix = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 0.85)?
            .set(1, 1, -0.3)?
            .set_symmetric(1, 2, 0.42)?
            .set(2, 2, 0.15)?
            .set_symmetric(2, 3, 0.37)?
            .set(3, 3, -0.2)?
            .set_symmetric(3, 4, 0.9)?
            .build()?;

        let denormalized = matrix.denormalize_bandpass(&mapping)?;
        let renormalized = denormalized.normalize_bandpass(&mapping)?;

        for (left, right) in matrix.as_slice().iter().zip(renormalized.as_slice()) {
            approx_eq(*left, *right, 1e-9);
        }
        Ok(())
    }

    #[test]
    fn bandpass_external_q_conversion_round_trips_ports() -> Result<()> {
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.9)?
            .set_symmetric(1, 2, 0.45)?
            .set_symmetric(2, 3, 1.1)?
            .build()?;

        let scaled = matrix.denormalize_bandpass_with_external_q(&mapping)?;
        assert!(scaled.source_external_q() > 0.0);
        assert!(scaled.load_external_q() > 0.0);

        let restored = scaled.matrix_hz().normalize_bandpass_with_external_q(&mapping)?;
        approx_eq(restored.at(0, 1).unwrap_or_default(), 0.9, 1e-9);
        approx_eq(restored.at(2, 3).unwrap_or_default(), 1.1, 1e-9);
        Ok(())
    }

    #[test]
    fn triplet_extraction_moves_zero_to_requested_center() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?;

        let extracted = matrix.extract_triplet(-1.3, 2)?;
        assert!(extracted.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(extracted.at(1, 3).unwrap_or_default().abs() > 1e-6);
        Ok(())
    }

    #[test]
    fn quadruplet_extraction_eliminates_one_internal_cross_coupling() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(6)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.84)?
            .set_symmetric(2, 3, 0.78)?
            .set_symmetric(3, 4, 0.72)?
            .set_symmetric(4, 5, 0.66)?
            .set_symmetric(5, 6, 0.61)?
            .set_symmetric(6, 7, 1.0)?
            .set(5, 5, 0.25)?
            .set(6, 6, -0.15)?
            .build()?;

        let extracted = matrix.extract_quadruplet(-1.1, 1.35, 2, 1, false)?;
        assert!(extracted.at(4, 2).unwrap_or_default().abs() <= 1e-6);
        assert!(extracted.at(3, 1).unwrap_or_default().abs() > 1e-6);
        Ok(())
    }

    #[test]
    fn trisection_extraction_pulls_tail_triplet_to_requested_window() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .topology(MatrixTopology::Arrow)
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.86)?
            .set_symmetric(1, 5, 0.25)?
            .set_symmetric(2, 3, 0.78)?
            .set_symmetric(3, 4, 0.69)?
            .set_symmetric(4, 5, 0.58)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.18)?
            .build()?;

        let extracted = matrix.extract_trisection(-1.25, (2, 4))?;
        assert!(extracted.at(2, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(extracted.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(extracted.at(1, 3).unwrap_or_default().abs() > 1e-6);
        assert!(extracted.at(1, 5).unwrap_or_default().abs() > 1e-6);
        Ok(())
    }

    #[test]
    fn trisection_extraction_rejects_non_arrow_input() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.86)?
            .set_symmetric(1, 5, 0.25)?
            .set_symmetric(2, 3, 0.78)?
            .set_symmetric(3, 4, 0.69)?
            .set_symmetric(4, 5, 0.58)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.18)?
            .build()?;

        let error = matrix
            .extract_trisection(-1.25, (2, 4))
            .expect_err("non-arrow matrix should be rejected");
        assert!(matches!(error, crate::error::MfsError::PreconditionViolation(_)));
        Ok(())
    }

    /// Reference triple-loop multiply for equivalence testing against nalgebra.
    fn reference_multiply(
        lhs: &super::CouplingMatrix,
        rhs: &super::CouplingMatrix,
    ) -> Vec<f64> {
        let side = lhs.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for col in 0..side {
                let mut acc = 0.0;
                for inner in 0..side {
                    acc += lhs.get(row, inner) * rhs.get(inner, col);
                }
                data[row * side + col] = acc;
            }
        }
        data
    }

    /// Reference transpose for equivalence testing against nalgebra.
    fn reference_transpose(matrix: &super::CouplingMatrix) -> Vec<f64> {
        let side = matrix.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for col in 0..side {
                data[col * side + row] = matrix.get(row, col);
            }
        }
        data
    }

    #[test]
    fn nalgebra_multiply_matches_reference_implementation() -> Result<()> {
        // Test case 1: Identity multiplication (order 3, side 5)
        let identity = super::CouplingMatrix::identity(3)?;
        let matrix_a = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 0.85)?
            .set(1, 1, -0.3)?
            .set_symmetric(1, 2, 0.42)?
            .set(2, 2, 0.15)?
            .set_symmetric(2, 3, 0.37)?
            .set(3, 3, -0.2)?
            .set_symmetric(3, 4, 0.9)?
            .build()?;

        let nalgebra_result = identity.multiply(&matrix_a);
        let reference_result = reference_multiply(&identity, &matrix_a);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "multiply deviation {diff} exceeds 1e-12");
        }

        // Test case 2: Non-trivial multiplication (order 4, side 6)
        let matrix_b = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 1.0)?
            .set(2, 2, 0.12)?
            .set(3, 3, -0.08)?
            .build()?;

        let matrix_c = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 0.95)?
            .set_symmetric(1, 2, 0.77)?
            .set_symmetric(2, 3, 0.63)?
            .set_symmetric(3, 4, 0.55)?
            .set_symmetric(4, 5, 0.88)?
            .set(1, 1, 0.05)?
            .set(4, 4, -0.11)?
            .build()?;

        let nalgebra_result = matrix_b.multiply(&matrix_c);
        let reference_result = reference_multiply(&matrix_b, &matrix_c);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "multiply deviation {diff} exceeds 1e-12");
        }

        // Test case 3: Small matrix (order 2, side 4)
        let matrix_d = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.9)?
            .set_symmetric(1, 2, 0.45)?
            .set_symmetric(2, 3, 1.1)?
            .set(1, 1, 0.3)?
            .set(2, 2, -0.2)?
            .build()?;

        let nalgebra_result = matrix_d.multiply(&matrix_d);
        let reference_result = reference_multiply(&matrix_d, &matrix_d);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "multiply deviation {diff} exceeds 1e-12");
        }

        Ok(())
    }

    #[test]
    fn nalgebra_transpose_matches_reference_implementation() -> Result<()> {
        // Test case 1: Symmetric matrix (transpose should be identity operation)
        let matrix_a = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 0.85)?
            .set(1, 1, -0.3)?
            .set_symmetric(1, 2, 0.42)?
            .set(2, 2, 0.15)?
            .set_symmetric(2, 3, 0.37)?
            .set(3, 3, -0.2)?
            .set_symmetric(3, 4, 0.9)?
            .build()?;

        let nalgebra_result = matrix_a.transpose();
        let reference_result = reference_transpose(&matrix_a);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "transpose deviation {diff} exceeds 1e-12");
        }

        // Test case 2: Non-symmetric matrix (order 4)
        let data: Vec<f64> = (0..36).map(|i| (i as f64) * 0.1 - 1.5).collect();
        let matrix_b = super::CouplingMatrix::new_with_topology(
            4,
            MatrixTopology::Transversal,
            data,
        )?;

        let nalgebra_result = matrix_b.transpose();
        let reference_result = reference_transpose(&matrix_b);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "transpose deviation {diff} exceeds 1e-12");
        }

        // Test case 3: Larger matrix (order 6, side 8)
        let data: Vec<f64> = (0..64).map(|i| ((i as f64) * 0.37).sin()).collect();
        let matrix_c = super::CouplingMatrix::new_with_topology(
            6,
            MatrixTopology::Transversal,
            data,
        )?;

        let nalgebra_result = matrix_c.transpose();
        let reference_result = reference_transpose(&matrix_c);
        for (na_val, ref_val) in nalgebra_result.as_slice().iter().zip(reference_result.iter()) {
            let diff = (na_val - ref_val).abs();
            assert!(diff < 1e-12, "transpose deviation {diff} exceeds 1e-12");
        }

        Ok(())
    }
}
