mod arrow;
mod folded;
mod sections;
mod wheel;

use crate::error::Result;
use crate::freq::FrequencyGrid;
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::response::ResponseSolver;
use crate::verify::{matches_topology_pattern, MatrixPatternTolerance};
use crate::verify::{compare_responses, ResponseCheckReport, ResponseTolerance};

pub use arrow::to_arrow;
pub use folded::to_folded;
pub use sections::{
    extract_quadruplet_section, extract_quadruplet_section_with_report,
    extract_quadruplet_section_with_response_check, extract_triplet_section,
    extract_triplet_section_with_report, extract_triplet_section_with_response_check,
    extract_trisection_section, extract_trisection_section_with_report,
    extract_trisection_section_with_response_check, SectionTransformEngine,
    SectionTransformOutcome,
};
pub use wheel::to_wheel;

/// Public topology name used by the transform facade.
pub type TopologyKind = MatrixTopology;

/// Lightweight outcome for topology conversion requests.
#[derive(Debug, Clone, PartialEq)]
pub struct TransformOutcome {
    /// Transformed coupling matrix.
    pub matrix: CouplingMatrix,
    /// Requested output topology.
    pub topology: TopologyKind,
    /// Minimal report describing the transform result.
    pub report: TransformReport,
}

/// Minimal transform report aligned with the public architecture notes.
#[derive(Debug, Clone, PartialEq)]
pub struct TransformReport {
    /// Topology declared by the input matrix metadata.
    pub source_topology: TopologyKind,
    /// Topology requested by the caller.
    pub requested_topology: TopologyKind,
    /// Topology attached to the output matrix metadata.
    pub result_topology: TopologyKind,
    /// Whether the output matrix matches the expected pattern for the requested topology.
    pub pattern_verified: bool,
    /// Shared electrical-check summary for this transform.
    pub response: ResponseCheckReport,
    /// Short implementation notes about the transform path that was used.
    pub notes: Vec<String>,
}

impl TransformReport {
    /// Returns whether the output matrix metadata and structure agree with the requested topology.
    pub fn passes(&self) -> bool {
        self.requested_topology == self.result_topology
            && self.pattern_verified
            && self.response.passes()
    }
}

/// Minimal facade that separates topology-conversion intent from matrix storage.
#[derive(Debug, Default, Clone, Copy)]
pub struct TransformEngine;

impl TransformEngine {
    /// Converts a matrix into the requested topology when supported.
    pub fn transform(
        &self,
        matrix: &CouplingMatrix,
        topology: TopologyKind,
    ) -> Result<TransformOutcome> {
        match topology {
            TopologyKind::Transversal => Ok(build_transform_outcome(
                matrix,
                topology,
                matrix.clone(),
                vec!["no topology conversion was applied".to_string()],
            )),
            TopologyKind::Folded => to_folded(matrix),
            TopologyKind::Arrow => to_arrow(matrix),
            TopologyKind::Wheel => to_wheel(matrix),
        }
    }

    /// Converts a matrix and records response deviation on a normalized sweep.
    pub fn transform_with_response_check(
        &self,
        matrix: &CouplingMatrix,
        topology: TopologyKind,
        grid: &FrequencyGrid,
        tolerance: ResponseTolerance,
    ) -> Result<TransformOutcome> {
        let mut outcome = self.transform(matrix, topology)?;
        let solver = ResponseSolver::default();
        let baseline = solver.evaluate_normalized(matrix, grid)?;
        let transformed = solver.evaluate_normalized(&outcome.matrix, grid)?;
        let comparison = compare_responses(&baseline, &transformed)?;
        outcome.report.response = ResponseCheckReport::from_comparison(comparison, tolerance);
        if outcome.report.response.passes() {
            outcome
                .report
                .notes
                .push("response invariance check passed on the supplied normalized grid".to_string());
        } else {
            outcome
                .report
                .notes
                .push("response invariance check failed on the supplied normalized grid".to_string());
        }
        Ok(outcome)
    }
}

/// Convenience wrapper for topology conversion through the transform facade.
pub fn transform_matrix(
    matrix: &CouplingMatrix,
    topology: TopologyKind,
) -> Result<TransformOutcome> {
    TransformEngine.transform(matrix, topology)
}

/// Convenience wrapper for topology conversion plus normalized response checking.
pub fn transform_matrix_with_response_check(
    matrix: &CouplingMatrix,
    topology: TopologyKind,
    grid: &FrequencyGrid,
    tolerance: ResponseTolerance,
) -> Result<TransformOutcome> {
    TransformEngine.transform_with_response_check(matrix, topology, grid, tolerance)
}

pub(crate) fn build_transform_outcome(
    input: &CouplingMatrix,
    topology: TopologyKind,
    matrix: CouplingMatrix,
    notes: Vec<String>,
) -> TransformOutcome {
    let pattern_verified = match topology {
        TopologyKind::Transversal => matrix.topology() == TopologyKind::Transversal,
        _ => matches_topology_pattern(&matrix, topology, MatrixPatternTolerance::default()),
    };

    TransformOutcome {
        topology,
        report: TransformReport {
            source_topology: input.topology(),
            requested_topology: topology,
            result_topology: matrix.topology(),
            pattern_verified,
            response: ResponseCheckReport::skipped(),
            notes,
        },
        matrix,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::FrequencyGrid;
    use crate::matrix::CouplingMatrixBuilder;
    use crate::verify::{matches_arrow_pattern, matches_folded_pattern, MatrixPatternTolerance};

    #[test]
    fn transform_engine_wraps_folded_conversion() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.3)?
            .set_symmetric(0, 3, 0.2)?
            .set_symmetric(1, 2, 0.6)?
            .set_symmetric(2, 3, 0.5)?
            .set_symmetric(3, 4, 1.0)?
            .build()?;

        let outcome = transform_matrix(&matrix, TopologyKind::Folded)?;
        assert_eq!(outcome.topology, TopologyKind::Folded);
        assert_eq!(outcome.report.source_topology, TopologyKind::Transversal);
        assert_eq!(outcome.report.result_topology, TopologyKind::Folded);
        assert!(outcome.report.passes());
        assert!(outcome.matrix.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn transform_engine_wraps_arrow_conversion() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.3)?
            .set_symmetric(0, 3, 0.2)?
            .set_symmetric(1, 2, 0.6)?
            .set_symmetric(2, 3, 0.5)?
            .set_symmetric(3, 4, 1.0)?
            .build()?;

        let outcome = transform_matrix(&matrix, TopologyKind::Arrow)?;
        assert_eq!(outcome.topology, TopologyKind::Arrow);
        assert_eq!(outcome.report.source_topology, TopologyKind::Transversal);
        assert_eq!(outcome.report.result_topology, TopologyKind::Arrow);
        assert!(outcome.report.passes());
        assert!(outcome.matrix.at(0, 2).unwrap_or_default().abs() <= 1e-6);
        Ok(())
    }

    #[test]
    fn transform_engine_wraps_wheel_conversion() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.3)?
            .set_symmetric(0, 3, 0.2)?
            .set_symmetric(1, 2, 0.6)?
            .set_symmetric(2, 3, 0.5)?
            .set_symmetric(3, 4, 1.0)?
            .build()?;

        let outcome = transform_matrix(&matrix, TopologyKind::Wheel)?;
        assert_eq!(outcome.topology, TopologyKind::Wheel);
        assert_eq!(outcome.report.result_topology, TopologyKind::Wheel);
        assert!(outcome.report.passes());
        Ok(())
    }

    #[test]
    fn transform_engine_can_attach_response_invariance_summary() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(3)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.3)?
            .set_symmetric(0, 3, 0.2)?
            .set_symmetric(1, 2, 0.6)?
            .set_symmetric(2, 3, 0.5)?
            .set_symmetric(3, 4, 1.0)?
            .build()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;

        let outcome = TransformEngine.transform_with_response_check(
            &matrix,
            TopologyKind::Folded,
            &grid,
            ResponseTolerance::default(),
        )?;

        assert!(outcome.report.response.comparison.is_some());
        assert_eq!(outcome.report.response.invariant, Some(true));
        assert!(outcome.report.passes());
        Ok(())
    }

    #[test]
    fn folded_transform_matches_folded_shape_rules() -> Result<()> {
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

        let folded = to_folded(&matrix)?.matrix;
        assert_eq!(folded.shape(), matrix.shape());
        assert!(matches_folded_pattern(
            &folded,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    fn arrow_transform_matches_arrow_shape_rules() -> Result<()> {
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

        let arrow = to_arrow(&matrix)?.matrix;
        assert_eq!(arrow.shape(), matrix.shape());
        assert!(matches_arrow_pattern(
            &arrow,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    fn wheel_transform_matches_current_wheel_shape_rules() -> Result<()> {
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

        let wheel = to_wheel(&matrix)?.matrix;
        assert_eq!(wheel.shape(), matrix.shape());
        assert!(matches_topology_pattern(
            &wheel,
            TopologyKind::Wheel,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }
}
