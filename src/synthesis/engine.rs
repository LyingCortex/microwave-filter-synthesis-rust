use crate::approx::PolynomialSet;
use crate::error::{MfsError, Result};
use crate::freq::BandPassMapping;
use crate::matrix::{BandPassScaledCouplingMatrix, CouplingMatrix, MatrixTopology};

use super::{
    build_transversal_from_residues, synthesize_admittance_polynomials,
    synthesize_residue_expansions, AdmittancePolynomials, ResidueExpansion, SectionSynthesis,
};

/// Indicates which matrix-construction path produced the final matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixSynthesisMethod {
    /// Matrix was reconstructed from Y-parameter poles and residues.
    ResidueExpansion,
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
        let _generalized = polynomials.generalized.as_ref().ok_or_else(|| {
            MfsError::PreconditionViolation(
                "residue expansion requires generalized Chebyshev data".to_string(),
            )
        })?;

        let (y11, y12, y22) = synthesize_residue_expansions(polynomials)?;
        let matrix = build_transversal_from_residues(polynomials, &y11, &y12, &y22)?;

        Ok(MatrixSynthesisOutcome {
            matrix,
            method: MatrixSynthesisMethod::ResidueExpansion,
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
        SectionSynthesis {}.synthesize_triplet_matrix(
            polynomials,
            transmission_zero,
            center_resonator,
        )
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
        SectionSynthesis {}.synthesize_quadruplet_matrix(
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
        SectionSynthesis {}.synthesize_trisection_matrix(
            polynomials,
            transmission_zero,
            zero_positions,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::{generalized_chebyshev_polynomials, PolynomialSet, PolynomialSetParams};
    use crate::fixtures::load_filter_database_end_to_end_fixture;
    use crate::freq::BandPassMapping;
    use crate::spec::FilterSpec;

    #[test]
    fn matrix_synthesizer_builds_matrix_from_generalized_polynomials() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        assert_eq!(matrix.order(), 3);
        assert_eq!(matrix.side(), 5);
        // Source and load couplings must be non-zero
        assert!(matrix.at(0, 1).unwrap_or_default().abs() > 1e-12);
        assert!(matrix.at(3, 4).unwrap_or_default().abs() > 1e-12);
        Ok(())
    }

    /// S2 boundary guard: the engine must respond to `generalized.e_s`, the
    /// field it actually consumes. Perturbing `e_s` must change the output
    /// matrix. This is the seam-level test for the P1 no-op class — callers
    /// that rebuild a polynomial set must update `generalized.e_s/f_s` (not
    /// just `.e`/`.f`), and this test fails if the engine ever stops reading
    /// the generalized data.
    #[test]
    fn matrix_synthesizer_is_sensitive_to_generalized_e_s() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;
        let base = generalized_chebyshev_polynomials(&spec)?;
        let base_matrix = MatrixSynthesisEngine.synthesize(&base)?;

        let mut perturbed = base.clone();
        if let Some(g) = perturbed.generalized.as_mut() {
            if let Some(e_s) = g.e_s.as_mut() {
                let c = e_s.coefficients.len();
                // Perturb the second-lowest coefficient (avoid the leading
                // term, which the synthesis renormalizes).
                let idx = c.min(2).saturating_sub(1);
                let k = e_s.coefficients[idx];
                e_s.coefficients[idx] = k * 1.001;
            }
        }
        let perturbed_matrix = MatrixSynthesisEngine.synthesize(&perturbed)?;

        let mut max_delta = 0.0_f64;
        let n = base_matrix.side();
        for i in 0..n {
            for j in 0..n {
                let d = (base_matrix.at(i, j).unwrap_or(0.0)
                    - perturbed_matrix.at(i, j).unwrap_or(0.0))
                .abs();
                max_delta = max_delta.max(d);
            }
        }
        assert!(
            max_delta > 1.0e-9,
            "engine must respond to a perturbed generalized.e_s (max|Δ| = {max_delta})"
        );
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_returns_precondition_violation_without_generalized_data() -> Result<()> {
        let polynomials = PolynomialSet::new(&PolynomialSetParams {
            order: 3,
            ripple_factor: 0.1,
            eps: 0.1,
            eps_r: 1.0,
            transmission_zeros_normalized: vec![-1.5, 1.5],
            e: vec![1.0, 0.2, 0.3, 0.4],
            f: vec![0.8, 0.6, 0.4, 0.2],
            p: vec![1.0, 0.5, -2.25],
        })?;

        let error = MatrixSynthesisEngine.synthesize(&polynomials).unwrap_err();
        assert!(matches!(error, MfsError::PreconditionViolation(_)));
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_apply_topology_during_synthesis() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let arrow =
            MatrixSynthesisEngine.synthesize_with_topology(&polynomials, MatrixTopology::Arrow)?;
        assert!(arrow.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        assert!(arrow.at(1, 3).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_emit_bandpass_scaled_matrix() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?.with_normalized_transmission_zeros(vec![-1.5, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let scaled = MatrixSynthesisEngine.synthesize_bandpass(
            &polynomials,
            MatrixTopology::Folded,
            &mapping,
        )?;
        assert_eq!(scaled.order(), 3);
        assert!(scaled.at(1, 1).unwrap_or_default().abs() > 1.0e9);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_emit_external_q_view() -> Result<()> {
        let spec = FilterSpec::new(2, 20.0)?;
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let scaled = MatrixSynthesisEngine.synthesize_bandpass_with_external_q(
            &polynomials,
            MatrixTopology::Transversal,
            &mapping,
        )?;
        assert!(scaled.source_external_q() > 0.0);
        assert!(scaled.load_external_q() > 0.0);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_triplet() -> Result<()> {
        let spec = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![-1.3]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let matrix = MatrixSynthesisEngine.synthesize_triplet(&polynomials, -1.3, 2)?;
        assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_quadruplet() -> Result<()> {
        let spec = FilterSpec::new(6, 20.0)?.with_normalized_transmission_zeros(vec![-1.1, 1.35]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let matrix =
            MatrixSynthesisEngine.synthesize_quadruplet(&polynomials, -1.1, 1.35, 2, 1, false)?;
        assert!(matrix.at(4, 2).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_extract_trisection() -> Result<()> {
        let spec = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![-1.25]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        let matrix = MatrixSynthesisEngine.synthesize_trisection(&polynomials, -1.25, (2, 4))?;
        assert!(matrix.at(2, 5).unwrap_or_default().abs() <= 1e-6);
        assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_can_build_admittance_polynomials_from_generalized_data() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

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
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

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
        let polynomials = generalized_chebyshev_polynomials(&fixture.spec)?;

        let outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.method, MatrixSynthesisMethod::ResidueExpansion);
        assert_eq!(outcome.matrix.order(), 4);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_accepts_all_pole_generalized_output() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;
        let polynomials = generalized_chebyshev_polynomials(&spec)?;

        assert!(polynomials.generalized.is_some());
        let outcome = MatrixSynthesisEngine.synthesize_with_details(&polynomials)?;
        assert_eq!(outcome.matrix.order(), 3);
        assert_eq!(outcome.matrix.side(), 5);
        Ok(())
    }
}
