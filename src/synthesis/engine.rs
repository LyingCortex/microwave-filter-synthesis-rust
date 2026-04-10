use crate::approx::PolynomialSet;
use crate::error::{MfsError, Result};
use crate::freq::BandPassMapping;
use crate::matrix::{BandPassScaledCouplingMatrix, CouplingMatrix, MatrixTopology};

use super::{
    build_transversal_from_residues, synthesize_admittance_polynomials,
    synthesize_placeholder_matrix, synthesize_residue_expansions, AdmittancePolynomials,
    ResidueExpansion, SectionSynthesis,
};

/// Indicates which matrix-construction path produced the final matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixSynthesisMethod {
    /// Matrix was reconstructed from Y-parameter poles and residues.
    ResidueExpansion,
    /// Matrix fell back to the placeholder chain builder.
    PlaceholderFallback,
}

/// Detailed result of one coupling-matrix synthesis attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixSynthesisOutcome {
    /// Final synthesized matrix.
    pub matrix: CouplingMatrix,
    /// Construction path that produced the matrix.
    pub method: MatrixSynthesisMethod,
}

/// Main synthesis engine for canonical and section-oriented matrix workflows.
#[derive(Debug, Default, Clone, Copy)]
pub struct MatrixSynthesisEngine;

impl MatrixSynthesisEngine {
    /// Builds polynomial-form Y parameters from generalized Chebyshev helper data.
    pub fn synthesize_admittance_polynomials(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<AdmittancePolynomials> {
        synthesize_admittance_polynomials(polynomials)
    }

    /// Splits the Y-parameter numerators into simple residues over the shared denominator.
    pub fn synthesize_residue_expansions(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<(ResidueExpansion, ResidueExpansion, ResidueExpansion)> {
        synthesize_residue_expansions(polynomials)
    }

    /// Produces a coupling matrix and reports which synthesis path succeeded.
    pub fn synthesize_with_details(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<MatrixSynthesisOutcome> {
        if polynomials.generalized.is_some() {
            match self
                .synthesize_residue_expansions(polynomials)
                .and_then(|(y11, y12, y22)| {
                    build_transversal_from_residues(polynomials, &y11, &y12, &y22)
                }) {
                Ok(matrix) => {
                    return Ok(MatrixSynthesisOutcome {
                        matrix,
                        method: MatrixSynthesisMethod::ResidueExpansion,
                    })
                }
                Err(MfsError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }

        Ok(MatrixSynthesisOutcome {
            matrix: synthesize_placeholder_matrix(polynomials)?,
            method: MatrixSynthesisMethod::PlaceholderFallback,
        })
    }

    /// Produces the canonical matrix implied by the approximation output.
    pub fn synthesize(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        Ok(self.synthesize_with_details(polynomials)?.matrix)
    }

    /// Synthesizes a matrix and immediately applies one topology transformation.
    pub fn synthesize_with_topology(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
    ) -> Result<CouplingMatrix> {
        self.synthesize(polynomials)?.transform_topology(topology)
    }

    /// Synthesizes a matrix, applies the requested topology, and scales it into band-pass units.
    pub fn synthesize_bandpass(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
        mapping: &BandPassMapping,
    ) -> Result<CouplingMatrix> {
        self.synthesize_with_topology(polynomials, topology)?
            .denormalize_bandpass(mapping)
    }

    /// Synthesizes a matrix, applies topology, and returns the band-pass/Qe representation.
    pub fn synthesize_bandpass_with_external_q(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
        mapping: &BandPassMapping,
    ) -> Result<BandPassScaledCouplingMatrix> {
        self.synthesize_with_topology(polynomials, topology)?
            .denormalize_bandpass_with_external_q(mapping)
    }

    /// Synthesizes a matrix and extracts one triplet section at the requested center.
    pub fn synthesize_triplet(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<CouplingMatrix> {
        SectionSynthesis::default().synthesize_triplet(polynomials, transmission_zero, center_resonator)
    }

    /// Synthesizes a matrix and extracts a quadruplet section from two adjacent triplets.
    pub fn synthesize_quadruplet(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<CouplingMatrix> {
        SectionSynthesis::default().synthesize_quadruplet(
            polynomials,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )
    }

    /// Synthesizes a matrix and pulls one trisection into the requested resonator window.
    pub fn synthesize_trisection(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<CouplingMatrix> {
        SectionSynthesis::default().synthesize_trisection(polynomials, transmission_zero, zero_positions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet};
    use crate::fixtures::load_filter_database_end_to_end_fixture;
    use crate::freq::{BandPassMapping, LowPassMapping};
    use crate::spec::{FilterSpec, TransmissionZero};

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
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

        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        assert_eq!(matrix.order(), 3);
        assert_eq!(matrix.side(), 5);
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

        let arrow = MatrixSynthesisEngine.synthesize_with_topology(&polynomials, MatrixTopology::Arrow)?;
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
            MatrixSynthesisEngine.synthesize_bandpass(&polynomials, MatrixTopology::Folded, &mapping)?;
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

        let scaled = MatrixSynthesisEngine
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

        let matrix = MatrixSynthesisEngine.synthesize_triplet(&polynomials, -1.3, 2)?;
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
            MatrixSynthesisEngine.synthesize_quadruplet(&polynomials, -1.1, 1.35, 2, 1, false)?;
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

        let matrix = MatrixSynthesisEngine.synthesize_trisection(&polynomials, -1.25, (2, 4))?;
        assert!(matrix.at(2, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_build_admittance_polynomials_from_generalized_data() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let admittance = MatrixSynthesisEngine.synthesize_admittance_polynomials(&polynomials)?;
        assert!(admittance.denominator.degree() >= 3);
        assert!(admittance.denominator.degree() <= 4);
        assert!(admittance.y11.degree() <= admittance.denominator.degree());
        assert_eq!(admittance.y12.degree(), 2);
        assert!(admittance.y22.degree() <= admittance.denominator.degree());
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_expand_admittance_residues() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

        let (y11, y12, y22) = MatrixSynthesisEngine.synthesize_residue_expansions(&polynomials)?;
        assert_eq!(y11.residues.len(), 4);
        assert_eq!(y12.residues.len(), 4);
        assert_eq!(y22.residues.len(), 4);
        assert!(y12.constant_term.is_none());
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_reports_residue_path_when_supported() -> Result<()> {
        let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
        let polynomials = ChebyshevApproximation.synthesize(&fixture.spec, &fixture.mapping)?;

        let outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.method, MatrixSynthesisMethod::ResidueExpansion);
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_accepts_all_pole_chebyshev_without_generalized_metadata() -> Result<()> {
        let spec = FilterSpec::chebyshev(3, 20.0)?;
        let mapping = LowPassMapping::new(1.0)?;
        let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

        assert!(polynomials.generalized.is_none());
        let outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.method, MatrixSynthesisMethod::PlaceholderFallback);
        assert_eq!(outcome.matrix.order(), 3);
        assert!(outcome.matrix.at(0, 1).unwrap_or_default().abs() > 1e-6);
        assert!(outcome.matrix.at(3, 4).unwrap_or_default().abs() > 1e-6);
        Ok(())
    }
}
