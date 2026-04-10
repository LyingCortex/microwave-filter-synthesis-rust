use crate::approx::PolynomialSet;
use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::transform::{
    extract_quadruplet_section, extract_quadruplet_section_with_report,
    extract_quadruplet_section_with_response_check, extract_triplet_section,
    extract_triplet_section_with_report, extract_triplet_section_with_response_check,
    extract_trisection_section, extract_trisection_section_with_report,
    extract_trisection_section_with_response_check,
};
use crate::verify::{ResponseCheckReport, ResponseTolerance, SectionVerificationReport};

use super::CanonicalMatrixSynthesis;

/// Synthesis facade for section-oriented workflows built on top of a canonical matrix.
#[derive(Debug, Default, Clone, Copy)]
pub struct SectionSynthesis {
    canonical: CanonicalMatrixSynthesis,
}

/// Result of a section-oriented synthesis call paired with structural verification.
#[derive(Debug, Clone, PartialEq)]
pub struct VerifiedSectionSynthesis {
    /// Extracted matrix after the requested section workflow.
    pub matrix: CouplingMatrix,
    /// Structural verification summary for the extracted section.
    pub verification: SectionVerificationReport,
    /// Optional response-invariance summary for the extracted section workflow.
    pub response: ResponseCheckReport,
    /// Short notes about the transform path taken underneath the synthesis facade.
    pub notes: Vec<String>,
}

impl VerifiedSectionSynthesis {
    /// Returns whether the attached structure and electrical checks passed.
    pub fn passes(&self) -> bool {
        self.verification.passes() && self.response.passes()
    }
}

impl SectionSynthesis {
    /// Synthesizes a matrix and extracts one triplet section at the requested center.
    pub fn synthesize_triplet(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<CouplingMatrix> {
        let canonical = self.canonical.synthesize(polynomials)?;
        extract_triplet_section(&canonical, transmission_zero, center_resonator)
    }

    /// Synthesizes and verifies one triplet extraction workflow.
    pub fn synthesize_triplet_with_report(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<VerifiedSectionSynthesis> {
        let canonical = self.canonical.synthesize(polynomials)?;
        let outcome =
            extract_triplet_section_with_report(&canonical, transmission_zero, center_resonator)?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
    }

    /// Synthesizes, verifies, and checks response invariance for one triplet workflow.
    pub fn synthesize_triplet_with_response_check(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<VerifiedSectionSynthesis> {
        let canonical = self.canonical.synthesize(polynomials)?;
        let outcome = extract_triplet_section_with_response_check(
            &canonical,
            transmission_zero,
            center_resonator,
            grid,
            tolerance,
        )?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
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
        let canonical = self.canonical.synthesize(polynomials)?;
        extract_quadruplet_section(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )
    }

    /// Synthesizes and verifies a quadruplet extraction workflow.
    pub fn synthesize_quadruplet_with_report(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<VerifiedSectionSynthesis> {
        let canonical = self.canonical.synthesize(polynomials)?;
        let outcome = extract_quadruplet_section_with_report(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
    }

    /// Synthesizes, verifies, and checks response invariance for one quadruplet workflow.
    pub fn synthesize_quadruplet_with_response_check(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<VerifiedSectionSynthesis> {
        let canonical = self.canonical.synthesize(polynomials)?;
        let outcome = extract_quadruplet_section_with_response_check(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
            grid,
            tolerance,
        )?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
    }

    /// Synthesizes a matrix and pulls one trisection into the requested resonator window.
    pub fn synthesize_trisection(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<CouplingMatrix> {
        let arrow = self
            .canonical
            .synthesize(polynomials)?
            .transform_topology(MatrixTopology::Arrow)?;
        extract_trisection_section(&arrow, transmission_zero, zero_positions)
    }

    /// Synthesizes and verifies one trisection extraction workflow.
    pub fn synthesize_trisection_with_report(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<VerifiedSectionSynthesis> {
        let arrow = self
            .canonical
            .synthesize(polynomials)?
            .transform_topology(MatrixTopology::Arrow)?;
        let outcome =
            extract_trisection_section_with_report(&arrow, transmission_zero, zero_positions)?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
    }

    /// Synthesizes, verifies, and checks response invariance for one trisection workflow.
    pub fn synthesize_trisection_with_response_check(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<VerifiedSectionSynthesis> {
        let arrow = self
            .canonical
            .synthesize(polynomials)?
            .transform_topology(MatrixTopology::Arrow)?;
        let outcome = extract_trisection_section_with_response_check(
            &arrow,
            transmission_zero,
            zero_positions,
            grid,
            tolerance,
        )?;
        Ok(VerifiedSectionSynthesis {
            matrix: outcome.matrix,
            verification: outcome.verification,
            response: outcome.response,
            notes: outcome.notes,
        })
    }
}
