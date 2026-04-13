use std::fs;
use std::path::PathBuf;

use mfs::prelude::*;

fn main() -> Result<()> {
    // `filter_spec(...)` expects normalized prototype zeros.
    let spec = filter_spec(4, 20.0, [-1.5], None)?;
    let synthesis = generalized_chebyshev(&spec)?;
    let matrix = synthesis.matrix.clone();
    let grid = FrequencyGrid::linspace(0.5, 1.5, 21)?;
    let transform = transform_matrix_with_response_check(
        &matrix,
        TopologyKind::Folded,
        &grid,
        ResponseTolerance::default(),
    )?;
    let response = ResponseSolver::default().evaluate_normalized(&transform.matrix, &grid)?;
    let output_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("docs")
        .join("quickstart_report_output.md");
    let report =
        render_markdown_synthesis_report(&spec, &synthesis, &transform, &response, false)?;
    fs::write(&output_path, report).map_err(|error| {
        MfsError::Unsupported(format!(
            "failed to write quickstart report {}: {error}",
            output_path.display()
        ))
    })?;

    println!("quickstart markdown report written to {}", output_path.display());
    Ok(())
}
