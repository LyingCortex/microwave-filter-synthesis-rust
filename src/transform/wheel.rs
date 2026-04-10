use crate::error::Result;

use super::{build_transform_outcome, TopologyKind, TransformOutcome};
use crate::matrix::CouplingMatrix;

/// Converts a matrix into the current wheel-reduction target.
///
/// The present implementation shares the same reduction path as arrow form but
/// keeps `Wheel` explicit at the API layer so later work can specialize it
/// without renaming user workflows.
pub fn to_wheel(matrix: &CouplingMatrix) -> Result<TransformOutcome> {
    let transformed = matrix.transform_topology(TopologyKind::Wheel)?;
    Ok(build_transform_outcome(
        matrix,
        TopologyKind::Wheel,
        transformed,
        vec!["currently reuses the arrow-style reduction path".to_string()],
    ))
}
