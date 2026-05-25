use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::CouplingMatrix;
use crate::verify::{
    verify_quadruplet_extraction, verify_triplet_extraction,
    verify_trisection_extraction, MatrixPatternTolerance, ResponseCheckReport,
    ResponseTolerance, SectionVerificationReport,
};

use super::attach_response_invariance_check;

/// Result of a section-extraction transform plus its structural verification.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionTransformOutcome {
    /// Transformed matrix after the requested section extraction.
    pub matrix: CouplingMatrix,
    /// Structural verification summary for the extracted section.
    pub verification: SectionVerificationReport,
    /// Shared electrical-check summary for this section transform.
    pub response: ResponseCheckReport,
    /// Short notes about the transform path that was applied.
    pub notes: Vec<String>,
}

impl SectionTransformOutcome {
    /// Returns whether all attached checks passed.
    pub fn passes(&self) -> bool {
        self.verification.passes() && self.response.passes()
    }
}

fn extract_triplet_impl(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
) -> Result<SectionTransformOutcome> {
    let matrix = extract_triplet_section_matrix(matrix, transmission_zero, center_resonator)?;
    let verification =
        verify_triplet_extraction(&matrix, center_resonator, MatrixPatternTolerance::default())?;
    Ok(SectionTransformOutcome {
        matrix,
        verification,
        response: ResponseCheckReport::skipped(),
        notes: vec!["triplet extraction used the current matrix backend".to_string()],
    })
}

fn extract_quadruplet_impl(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
) -> Result<SectionTransformOutcome> {
    let matrix = extract_quadruplet_section_matrix(
        matrix,
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
    )?;
    let verification =
        verify_quadruplet_extraction(&matrix, position, MatrixPatternTolerance::default())?;
    Ok(SectionTransformOutcome {
        matrix,
        verification,
        response: ResponseCheckReport::skipped(),
        notes: vec!["quadruplet extraction used the current matrix backend".to_string()],
    })
}

fn extract_trisection_impl(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
) -> Result<SectionTransformOutcome> {
    let matrix = extract_trisection_section_matrix(matrix, transmission_zero, zero_positions)?;
    let verification =
        verify_trisection_extraction(&matrix, zero_positions, MatrixPatternTolerance::default())?;
    Ok(SectionTransformOutcome {
        matrix,
        verification,
        response: ResponseCheckReport::skipped(),
        notes: vec!["trisection extraction used the current matrix backend".to_string()],
    })
}

/// Convenience wrapper for triplet extraction with structural verification.
pub fn extract_triplet_section(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
) -> Result<SectionTransformOutcome> {
    extract_triplet_impl(matrix, transmission_zero, center_resonator)
}

/// Convenience wrapper for triplet extraction returning only the matrix.
pub fn extract_triplet_section_matrix(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
) -> Result<CouplingMatrix> {
    matrix.extract_triplet(transmission_zero, center_resonator)
}

/// Convenience wrapper for reported triplet extraction plus response checking.
pub fn extract_triplet_section_with_response_check(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<SectionTransformOutcome> {
    let mut outcome = extract_triplet_impl(matrix, transmission_zero, center_resonator)?;
    attach_response_check(&mut outcome, matrix, grid, tolerance)?;
    Ok(outcome)
}

/// Convenience wrapper for quadruplet extraction with structural verification.
pub fn extract_quadruplet_section(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
) -> Result<SectionTransformOutcome> {
    extract_quadruplet_impl(
        matrix,
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
    )
}

/// Convenience wrapper for quadruplet extraction returning only the matrix.
pub fn extract_quadruplet_section_matrix(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
) -> Result<CouplingMatrix> {
    matrix.extract_quadruplet(
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
    )
}

/// Convenience wrapper for reported quadruplet extraction plus response checking.
pub fn extract_quadruplet_section_with_response_check(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<SectionTransformOutcome> {
    let mut outcome = extract_quadruplet_impl(
        matrix,
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
    )?;
    attach_response_check(&mut outcome, matrix, grid, tolerance)?;
    Ok(outcome)
}

/// Convenience wrapper for trisection extraction with structural verification.
pub fn extract_trisection_section(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
) -> Result<SectionTransformOutcome> {
    extract_trisection_impl(matrix, transmission_zero, zero_positions)
}

/// Convenience wrapper for trisection extraction returning only the matrix.
pub fn extract_trisection_section_matrix(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
) -> Result<CouplingMatrix> {
    matrix.extract_trisection(transmission_zero, zero_positions)
}

/// Convenience wrapper for reported trisection extraction plus response checking.
pub fn extract_trisection_section_with_response_check(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<SectionTransformOutcome> {
    let mut outcome = extract_trisection_impl(matrix, transmission_zero, zero_positions)?;
    attach_response_check(&mut outcome, matrix, grid, tolerance)?;
    Ok(outcome)
}

fn attach_response_check(
    outcome: &mut SectionTransformOutcome,
    baseline_matrix: &CouplingMatrix,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<()> {
    let (report, note) =
        attach_response_invariance_check(&outcome.matrix, baseline_matrix, grid, tolerance)?;
    outcome.response = report;
    outcome.notes.push(note);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::FrequencyGrid;
    use crate::matrix::{CouplingMatrixBuilder, MatrixTopology};

    #[test]
    fn triplet_section_extraction_reports_structure() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?;

        let outcome = extract_triplet_section(&matrix, -1.3, 2)?;
        assert!(outcome.verification.passes());
        Ok(())
    }

    #[test]
    fn trisection_section_extraction_reports_structure() -> Result<()> {
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

        let outcome = extract_trisection_section(&matrix, -1.25, (2, 4))?;
        assert!(outcome.verification.passes());
        Ok(())
    }

    #[test]
    fn triplet_section_extraction_can_attach_response_summary() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;

        let outcome = extract_triplet_section_with_response_check(
            &matrix,
            -1.3,
            2,
            &grid,
            ResponseTolerance::default(),
        )?;
        assert!(outcome.passes());
        assert_eq!(outcome.response.invariant, Some(true));
        assert!(outcome.response.comparison.is_some());
        Ok(())
    }

    #[test]
    fn quadruplet_section_extraction_can_attach_response_summary() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(6)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.88)?
            .set_symmetric(2, 3, 0.81)?
            .set_symmetric(3, 4, 0.74)?
            .set_symmetric(4, 5, 0.67)?
            .set_symmetric(5, 6, 0.6)?
            .set_symmetric(6, 7, 1.0)?
            .set(6, 6, 0.16)?
            .build()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;

        let outcome = extract_quadruplet_section_with_response_check(
            &matrix,
            -1.1,
            1.35,
            2,
            1,
            false,
            &grid,
            ResponseTolerance::default(),
        )?;
        assert!(outcome.passes());
        assert_eq!(outcome.response.invariant, Some(true));
        assert!(outcome.response.comparison.is_some());
        Ok(())
    }

    #[test]
    fn trisection_section_extraction_can_attach_response_summary() -> Result<()> {
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
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;

        let outcome = extract_trisection_section_with_response_check(
            &matrix,
            -1.25,
            (2, 4),
            &grid,
            ResponseTolerance::default(),
        )?;
        assert!(outcome.passes());
        assert_eq!(outcome.response.invariant, Some(true));
        assert!(outcome.response.comparison.is_some());
        Ok(())
    }
}
