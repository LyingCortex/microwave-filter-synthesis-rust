mod builder;
mod coupling_matrix;
mod synthesis;

pub use builder::CouplingMatrixBuilder;
pub use coupling_matrix::{BandPassScaledCouplingMatrix, CouplingMatrix, MatrixShape, MatrixTopology};
pub use synthesis::{
    AdmittancePolynomials, CouplingMatrixSynthesizer, MatrixSynthesisMethod,
    MatrixSynthesisOutcome, ResidueExpansion, ResiduePole,
};

#[cfg(test)]
mod tests {
    use crate::approx::{ApproximationEngine, PolynomialSet};
    use crate::error::Result;
    use crate::freq::{BandPassMapping, LowPassMapping};
    use crate::spec::{FilterParameter, TransmissionZero};

    use super::{CouplingMatrixBuilder, CouplingMatrixSynthesizer, MatrixShape, MatrixTopology};

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
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_builds_matrix_from_polynomials() -> Result<()> {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![-1.5, 1.5],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, 0.5, -2.25],
        )?;

        let matrix = CouplingMatrixSynthesizer.synthesize(&polynomials)?;
        assert_eq!(matrix.order(), 3);
        assert_eq!(matrix.side(), 5);
        assert_eq!(matrix.shape(), MatrixShape { rows: 5, cols: 5 });
        approx_eq(matrix.at(0, 1).unwrap_or_default(), 0.8, 1e-12);
        approx_eq(matrix.at(1, 2).unwrap_or_default(), 0.4, 1e-12);
        approx_eq(matrix.at(2, 2).unwrap_or_default(), 0.5, 1e-12);
        approx_eq(matrix.at(3, 3).unwrap_or_default(), -2.25, 1e-12);
        approx_eq(matrix.at(3, 4).unwrap_or_default(), 1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_apply_topology_during_synthesis() -> Result<()> {
        let polynomials = PolynomialSet::new(
            4,
            0.1,
            0.1,
            1.0,
            vec![-1.5, 1.5],
            vec![1.0, 0.9, 0.7, 0.5, 0.3],
            vec![0.8, 0.6, 0.4, 0.2, 0.1],
            vec![0.2, -0.1, 0.05],
        )?;

        let arrow = CouplingMatrixSynthesizer.synthesize_with_topology(&polynomials, MatrixTopology::Arrow)?;
        assert!(arrow.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        assert!(arrow.at(1, 3).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_emit_bandpass_scaled_matrix() -> Result<()> {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            0.1,
            1.0,
            vec![-1.5, 1.5],
            vec![1.0, 0.5, 0.3, 0.2],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![0.15, -0.08, 0.02],
        )?;
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let scaled =
            CouplingMatrixSynthesizer.synthesize_bandpass(&polynomials, MatrixTopology::Folded, &mapping)?;
        assert_eq!(scaled.order(), 3);
        assert!(scaled.at(1, 1).unwrap_or_default().abs() > 1.0e9);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_emit_external_q_view() -> Result<()> {
        let polynomials = PolynomialSet::new(
            2,
            0.1,
            0.1,
            1.0,
            vec![],
            vec![1.0, 0.4, 0.2],
            vec![0.9, 0.5, 0.3],
            vec![0.12, -0.03],
        )?;
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let scaled = CouplingMatrixSynthesizer
            .synthesize_bandpass_with_external_q(&polynomials, MatrixTopology::Transversal, &mapping)?;
        assert!(scaled.source_external_q() > 0.0);
        assert!(scaled.load_external_q() > 0.0);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_triplet() -> Result<()> {
        let polynomials = PolynomialSet::new(
            5,
            0.1,
            0.1,
            1.0,
            vec![-1.3],
            vec![1.0, 0.9, 0.8, 0.7, 0.6, 0.5],
            vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4],
            vec![0.2, -0.1, 0.05],
        )?;

        let matrix = CouplingMatrixSynthesizer.synthesize_triplet(&polynomials, -1.3, 2)?;
        assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_quadruplet() -> Result<()> {
        let polynomials = PolynomialSet::new(
            6,
            0.1,
            0.1,
            1.0,
            vec![-1.1, 1.35],
            vec![1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4],
            vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
            vec![0.25, -0.15, 0.08],
        )?;

        let matrix =
            CouplingMatrixSynthesizer.synthesize_quadruplet(&polynomials, -1.1, 1.35, 2, 1, false)?;
        assert!(matrix.at(4, 2).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_trisection() -> Result<()> {
        let polynomials = PolynomialSet::new(
            5,
            0.1,
            0.1,
            1.0,
            vec![-1.25],
            vec![1.0, 0.92, 0.84, 0.76, 0.68, 0.6],
            vec![0.95, 0.87, 0.79, 0.71, 0.63, 0.55],
            vec![0.18, -0.07, 0.03],
        )?;

        let matrix = CouplingMatrixSynthesizer.synthesize_trisection(&polynomials, -1.25, (2, 4))?;
        assert!(matrix.at(2, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_build_admittance_polynomials_from_generalized_data() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = crate::approx::ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let admittance = CouplingMatrixSynthesizer.synthesize_admittance_polynomials(&polynomials)?;
        assert!(admittance.denominator.degree() >= 3);
        assert!(admittance.denominator.degree() <= 4);
        assert!(admittance.y11.degree() <= admittance.denominator.degree());
        assert_eq!(admittance.y12.degree(), 2);
        assert!(admittance.y22.degree() <= admittance.denominator.degree());
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_expand_admittance_residues() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = crate::approx::ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let (y11, y12, y22) = CouplingMatrixSynthesizer.synthesize_residue_expansions(&polynomials)?;
        assert_eq!(y11.residues.len(), 4);
        assert_eq!(y12.residues.len(), 4);
        assert_eq!(y22.residues.len(), 4);
        assert!(y12.constant_term.is_none());
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_reports_residue_path_when_supported() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = crate::approx::ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let outcome = CouplingMatrixSynthesizer.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.method, super::MatrixSynthesisMethod::ResidueExpansion);
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_accepts_all_pole_chebyshev_with_generalized_metadata() -> Result<()> {
        let spec = FilterParameter::chebyshev(3, 20.0)?;
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = crate::approx::ChebyshevApproximation.synthesize(&spec, &mapping)?;

        assert!(polynomials.generalized.is_some());
        let outcome = CouplingMatrixSynthesizer.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.method, super::MatrixSynthesisMethod::PlaceholderFallback);
        assert_eq!(outcome.matrix.order(), 3);
        assert!(outcome.matrix.at(0, 1).unwrap_or_default().abs() > 1e-6);
        assert!(outcome.matrix.at(3, 4).unwrap_or_default().abs() > 1e-6);
        Ok(())
    }

    #[test]
    fn folded_transform_zeroes_far_source_couplings() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.4)?
            .set_symmetric(0, 3, 0.3)?
            .set_symmetric(0, 4, 0.2)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        let folded = matrix.transform_topology(MatrixTopology::Folded)?;
        assert_eq!(folded.shape(), matrix.shape());
        assert!(folded.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        assert!(folded.at(0, 3).unwrap_or_default().abs() <= 1e-6);
        assert!(folded.at(0, 4).unwrap_or_default().abs() <= 1e-6);
        approx_eq(
            folded.at(1, 0).unwrap_or_default(),
            folded.at(0, 1).unwrap_or_default(),
            1e-12,
        );
        Ok(())
    }

    #[test]
    fn arrow_transform_zeroes_non_arrow_source_couplings() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.4)?
            .set_symmetric(0, 3, 0.3)?
            .set_symmetric(0, 4, 0.2)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        let arrow = matrix.transform_topology(MatrixTopology::Arrow)?;
        assert_eq!(arrow.shape(), matrix.shape());
        for row in 0..(matrix.order().saturating_sub(1)) {
            for col in (row + 2)..=matrix.order() {
                assert!(
                    arrow.at(row, col).unwrap_or_default().abs() <= 1e-6,
                    "expected arrow reduction to clear ({row}, {col}), got {}",
                    arrow.at(row, col).unwrap_or_default()
                );
            }
        }
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
}
