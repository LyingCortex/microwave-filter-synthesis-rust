use crate::error::Result;
use crate::matrix::CouplingMatrix;

use super::{build_transform_outcome, TopologyKind, TransformOutcome};

/// Converts a matrix into folded topology through the current backend.
pub fn to_folded(matrix: &CouplingMatrix) -> Result<TransformOutcome> {
    let transformed = matrix.transform_topology(TopologyKind::Folded)?;
    Ok(build_transform_outcome(
        matrix,
        TopologyKind::Folded,
        transformed,
        vec!["used the current folded reduction backend".to_string()],
    ))
}
