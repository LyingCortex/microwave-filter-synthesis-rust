use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::CouplingMatrix;
use crate::response::ResponseSolver;
use crate::verify::{
    compare_responses, verify_quadruplet_extraction, verify_triplet_extraction,
    verify_trisection_extraction, MatrixPatternTolerance, ResponseCheckReport,
    ResponseTolerance, SectionVerificationReport,
};

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

/// Transform facade for section-oriented extraction workflows.
#[derive(Debug, Default, Clone, Copy)]
pub struct SectionTransformEngine;

impl SectionTransformEngine {
    /// Extracts one triplet section at the requested center.
    pub fn extract_triplet(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<CouplingMatrix> {
        matrix.extract_triplet(transmission_zero, center_resonator)
    }

    /// Extracts one triplet section and returns a structural verification report.
    pub fn extract_triplet_with_report(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<SectionTransformOutcome> {
        let matrix = self.extract_triplet(matrix, transmission_zero, center_resonator)?;
        let verification =
            verify_triplet_extraction(&matrix, center_resonator, MatrixPatternTolerance::default())?;
        Ok(SectionTransformOutcome {
            matrix,
            verification,
            response: ResponseCheckReport::skipped(),
            notes: vec!["triplet extraction used the current matrix backend".to_string()],
        })
    }

    /// Extracts one triplet section and checks response invariance on a normalized sweep.
    pub fn extract_triplet_with_response_check(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        center_resonator: usize,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<SectionTransformOutcome> {
        let mut outcome = self.extract_triplet_with_report(matrix, transmission_zero, center_resonator)?;
        attach_response_check(&mut outcome, matrix, grid, tolerance)?;
        Ok(outcome)
    }

    /// Extracts a quadruplet from two adjacent triplets.
    pub fn extract_quadruplet(
        &self,
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

    /// Extracts a quadruplet and returns a structural verification report.
    pub fn extract_quadruplet_with_report(
        &self,
        matrix: &CouplingMatrix,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<SectionTransformOutcome> {
        let matrix = self.extract_quadruplet(
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

    /// Extracts a quadruplet and checks response invariance on a normalized sweep.
    pub fn extract_quadruplet_with_response_check(
        &self,
        matrix: &CouplingMatrix,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<SectionTransformOutcome> {
        let mut outcome = self.extract_quadruplet_with_report(
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

    /// Extracts one trisection from an arrow-style matrix.
    pub fn extract_trisection(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<CouplingMatrix> {
        matrix.extract_trisection(transmission_zero, zero_positions)
    }

    /// Extracts one trisection and returns a structural verification report.
    pub fn extract_trisection_with_report(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<SectionTransformOutcome> {
        let matrix = self.extract_trisection(matrix, transmission_zero, zero_positions)?;
        let verification = verify_trisection_extraction(
            &matrix,
            zero_positions,
            MatrixPatternTolerance::default(),
        )?;
        Ok(SectionTransformOutcome {
            matrix,
            verification,
            response: ResponseCheckReport::skipped(),
            notes: vec!["trisection extraction used the current matrix backend".to_string()],
        })
    }

    /// Extracts one trisection and checks response invariance on a normalized sweep.
    pub fn extract_trisection_with_response_check(
        &self,
        matrix: &CouplingMatrix,
        transmission_zero: f64,
        zero_positions: (usize, usize),
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<SectionTransformOutcome> {
        let mut outcome = self.extract_trisection_with_report(matrix, transmission_zero, zero_positions)?;
        attach_response_check(&mut outcome, matrix, grid, tolerance)?;
        Ok(outcome)
    }
}

/// Convenience wrapper for triplet extraction.
pub fn extract_triplet_section(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
) -> Result<CouplingMatrix> {
    SectionTransformEngine.extract_triplet(matrix, transmission_zero, center_resonator)
}

/// Convenience wrapper for reported triplet extraction.
pub fn extract_triplet_section_with_report(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
) -> Result<SectionTransformOutcome> {
    SectionTransformEngine.extract_triplet_with_report(matrix, transmission_zero, center_resonator)
}

/// Convenience wrapper for reported triplet extraction plus response checking.
pub fn extract_triplet_section_with_response_check(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    center_resonator: usize,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<SectionTransformOutcome> {
    SectionTransformEngine.extract_triplet_with_response_check(
        matrix,
        transmission_zero,
        center_resonator,
        grid,
        tolerance,
    )
}

/// Convenience wrapper for quadruplet extraction.
pub fn extract_quadruplet_section(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
) -> Result<CouplingMatrix> {
    SectionTransformEngine.extract_quadruplet(
        matrix,
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
    )
}

/// Convenience wrapper for reported quadruplet extraction.
pub fn extract_quadruplet_section_with_report(
    matrix: &CouplingMatrix,
    first_zero: f64,
    second_zero: f64,
    position: usize,
    common_resonator: usize,
    swap_zero_order: bool,
) -> Result<SectionTransformOutcome> {
    SectionTransformEngine.extract_quadruplet_with_report(
        matrix,
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
    SectionTransformEngine.extract_quadruplet_with_response_check(
        matrix,
        first_zero,
        second_zero,
        position,
        common_resonator,
        swap_zero_order,
        grid,
        tolerance,
    )
}

/// Convenience wrapper for trisection extraction.
pub fn extract_trisection_section(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
) -> Result<CouplingMatrix> {
    SectionTransformEngine.extract_trisection(matrix, transmission_zero, zero_positions)
}

/// Convenience wrapper for reported trisection extraction.
pub fn extract_trisection_section_with_report(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
) -> Result<SectionTransformOutcome> {
    SectionTransformEngine.extract_trisection_with_report(matrix, transmission_zero, zero_positions)
}

/// Convenience wrapper for reported trisection extraction plus response checking.
pub fn extract_trisection_section_with_response_check(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    zero_positions: (usize, usize),
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<SectionTransformOutcome> {
    SectionTransformEngine.extract_trisection_with_response_check(
        matrix,
        transmission_zero,
        zero_positions,
        grid,
        tolerance,
    )
}

fn attach_response_check(
    outcome: &mut SectionTransformOutcome,
    baseline_matrix: &CouplingMatrix,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<()> {
    let solver = ResponseSolver::default();
    let baseline = solver.evaluate_normalized(baseline_matrix, grid)?;
    let transformed = solver.evaluate_normalized(&outcome.matrix, grid)?;
    let comparison = compare_responses(&baseline, &transformed)?;
    outcome.response = ResponseCheckReport::from_comparison(comparison, tolerance);
    if outcome.response.passes() {
        outcome
            .notes
            .push("response invariance check passed on the supplied normalized grid".to_string());
    } else {
        outcome
            .notes
            .push("response invariance check failed on the supplied normalized grid".to_string());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::FrequencyGrid;
    use crate::matrix::{CouplingMatrixBuilder, MatrixTopology};

    #[test]
    fn section_transform_engine_reports_triplet_structure() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?;

        let outcome = extract_triplet_section_with_report(&matrix, -1.3, 2)?;
        assert!(outcome.verification.passes());
        Ok(())
    }

    #[test]
    fn section_transform_engine_reports_trisection_structure() -> Result<()> {
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

        let outcome = extract_trisection_section_with_report(&matrix, -1.25, (2, 4))?;
        assert!(outcome.verification.passes());
        Ok(())
    }

    #[test]
    fn section_transform_engine_can_attach_triplet_response_summary() -> Result<()> {
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
    fn section_transform_engine_can_attach_quadruplet_response_summary() -> Result<()> {
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
    fn section_transform_engine_can_attach_trisection_response_summary() -> Result<()> {
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
