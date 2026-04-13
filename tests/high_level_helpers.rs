use mfs::prelude::*;

#[test]
fn top_level_helpers_support_canonical_and_topology_workflows() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .normalized_transmission_zeros(vec![-2.0, 1.5])
        .build()?;
    let polynomials = GeneralizedChebyshevApproximation.synthesize(&spec)?;

    let matrix = mfs::synthesize_canonical_matrix(&polynomials)?;
    let folded = mfs::synthesize_matrix_with_topology(&polynomials, TopologyKind::Folded)?;

    assert_eq!(matrix.order(), 4);
    assert!(matches_folded_pattern(
        &folded,
        MatrixPatternTolerance::default()
    ));
    Ok(())
}
