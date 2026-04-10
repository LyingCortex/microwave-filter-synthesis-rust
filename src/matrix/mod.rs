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
//! - [`crate::synthesis::CanonicalMatrixSynthesis`]
//! - [`crate::synthesis::SectionSynthesis`]
//! - [`crate::transform::TransformEngine`]
//! - [`crate::verify`] helpers for structural and response checks

mod builder;
mod coupling_matrix;

pub use builder::CouplingMatrixBuilder;
pub use coupling_matrix::{BandPassScaledCouplingMatrix, CouplingMatrix, MatrixShape, MatrixTopology};

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
        assert!(matches!(error, crate::error::MfsError::Unsupported(_)));
        Ok(())
    }
}
