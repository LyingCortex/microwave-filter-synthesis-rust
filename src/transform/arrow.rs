use crate::error::Result;
use crate::matrix::CouplingMatrix;

use super::{build_transform_outcome, TopologyKind, TransformOutcome};

/// Converts a matrix into arrow topology through the current backend.
pub fn to_arrow(matrix: &CouplingMatrix) -> Result<TransformOutcome> {
    let transformed = matrix.transform_topology(TopologyKind::Arrow)?;
    Ok(build_transform_outcome(
        matrix,
        TopologyKind::Arrow,
        transformed,
        vec!["used the current arrow reduction backend".to_string()],
    ))
}
