use mfs::approx::generalized_chebyshev_polynomials;
use mfs::error::Result;
use mfs::prelude::{FilterSpec, TopologyKind};
use mfs::verify::{matches_folded_pattern, MatrixPatternTolerance};

#[test]
fn top_level_helpers_support_canonical_and_topology_workflows() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .normalized_transmission_zeros(vec![-2.0, 1.5])
        .build()?;
    let polynomials = generalized_chebyshev_polynomials(&spec)?;

    let matrix = mfs::synthesize_canonical_matrix(&polynomials)?;
    let folded = mfs::synthesize_matrix_with_topology(&polynomials, TopologyKind::Folded)?;

    assert_eq!(matrix.order(), 4);
    assert!(matches_folded_pattern(
        &folded,
        MatrixPatternTolerance::default()
    ));
    Ok(())
}
