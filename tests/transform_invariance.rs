use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrixBuilder;
use mfs::response::ResponseSolver;
use mfs::transform::{to_arrow, to_folded, to_wheel, TopologyKind};
use mfs::verify::{compare_responses, ResponseTolerance};

fn sample_transversal_matrix() -> Result<mfs::CouplingMatrix> {
    CouplingMatrixBuilder::new(4)?
        .set_symmetric(0, 1, 1.0)?
        .set_symmetric(0, 2, 0.4)?
        .set_symmetric(0, 3, 0.3)?
        .set_symmetric(0, 4, 0.2)?
        .set_symmetric(1, 2, 0.7)?
        .set_symmetric(2, 3, 0.6)?
        .set_symmetric(3, 4, 0.5)?
        .set_symmetric(4, 5, 1.0)?
        .build()
}

fn assert_invariant_after_transform(
    topology: TopologyKind,
    transformed: &mfs::CouplingMatrix,
) -> Result<()> {
    let matrix = sample_transversal_matrix()?;
    let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;

    let solver = ResponseSolver::default();
    let baseline = solver.evaluate_normalized(&matrix, &grid)?;
    let transformed = solver.evaluate_normalized(transformed, &grid)?;
    let comparison = compare_responses(&baseline, &transformed)?;

    assert!(
        comparison.passes(ResponseTolerance::default()),
        "response changed after {topology:?} transform: {comparison:?}"
    );
    Ok(())
}

#[test]
fn folded_transform_preserves_sampled_response_for_manual_case() -> Result<()> {
    let matrix = sample_transversal_matrix()?;
    let folded = to_folded(&matrix)?.matrix;
    assert_invariant_after_transform(TopologyKind::Folded, &folded)
}

#[test]
fn arrow_transform_preserves_sampled_response_for_manual_case() -> Result<()> {
    let matrix = sample_transversal_matrix()?;
    let arrow = to_arrow(&matrix)?.matrix;
    assert_invariant_after_transform(TopologyKind::Arrow, &arrow)
}

#[test]
fn wheel_transform_preserves_sampled_response_for_manual_case() -> Result<()> {
    let matrix = sample_transversal_matrix()?;
    let wheel = to_wheel(&matrix)?.matrix;
    assert_invariant_after_transform(TopologyKind::Wheel, &wheel)
}
