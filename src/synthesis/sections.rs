use crate::approx::PolynomialSet;
use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::transform::{
    extract_quadruplet_section, extract_quadruplet_section_matrix,
    extract_quadruplet_section_with_response_check, extract_triplet_section,
    extract_triplet_section_matrix, extract_triplet_section_with_response_check,
    extract_trisection_section, extract_trisection_section_matrix,
    extract_trisection_section_with_response_check, SectionTransformOutcome,
};
use crate::verify::ResponseTolerance;

use super::MatrixSynthesisEngine;

/// Synthesis facade for section-oriented workflows built on top of a canonical matrix.
#[derive(Debug, Default, Clone, Copy)]
pub struct SectionSynthesis;

impl SectionSynthesis {
    fn synthesize_canonical_matrix(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize(polynomials)
    }

    fn synthesize_arrow_matrix(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        self.synthesize_canonical_matrix(polynomials)?
            .transform_topology(MatrixTopology::Arrow)
    }

    /// Synthesizes a matrix and extracts one triplet section at the requested center.
    pub fn synthesize_triplet_matrix(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<CouplingMatrix> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_triplet_section_matrix(&canonical, transmission_zero, center_resonator)
    }

    /// Synthesizes and verifies one triplet extraction workflow.
    pub fn synthesize_triplet(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<SectionTransformOutcome> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_triplet_section(&canonical, transmission_zero, center_resonator)
    }

    /// Synthesizes, verifies, and checks response invariance for one triplet workflow.
    pub fn synthesize_triplet_with_response_check(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<SectionTransformOutcome> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_triplet_section_with_response_check(
            &canonical,
            transmission_zero,
            center_resonator,
            grid,
            tolerance,
        )
    }

    /// Synthesizes a matrix and extracts a quadruplet section from two adjacent triplets.
    pub fn synthesize_quadruplet_matrix(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<CouplingMatrix> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_quadruplet_section_matrix(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )
    }

    /// Synthesizes and verifies a quadruplet extraction workflow.
    pub fn synthesize_quadruplet(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<SectionTransformOutcome> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_quadruplet_section(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )
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
    ) -> Result<SectionTransformOutcome> {
        let canonical = self.synthesize_canonical_matrix(polynomials)?;
        extract_quadruplet_section_with_response_check(
            &canonical,
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
            grid,
            tolerance,
        )
    }

    /// Synthesizes a matrix and pulls one trisection into the requested resonator window.
    pub fn synthesize_trisection_matrix(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<CouplingMatrix> {
        let arrow = self.synthesize_arrow_matrix(polynomials)?;
        extract_trisection_section_matrix(&arrow, transmission_zero, zero_positions)
    }

    /// Synthesizes and verifies one trisection extraction workflow.
    pub fn synthesize_trisection(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<SectionTransformOutcome> {
        let arrow = self.synthesize_arrow_matrix(polynomials)?;
        extract_trisection_section(&arrow, transmission_zero, zero_positions)
    }

    /// Synthesizes, verifies, and checks response invariance for one trisection workflow.
    pub fn synthesize_trisection_with_response_check(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<SectionTransformOutcome> {
        let arrow = self.synthesize_arrow_matrix(polynomials)?;
        extract_trisection_section_with_response_check(
            &arrow,
            transmission_zero,
            zero_positions,
            grid,
            tolerance,
        )
    }
}
