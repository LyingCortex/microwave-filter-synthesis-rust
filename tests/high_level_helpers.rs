use mfs::prelude::*;

#[test]
fn top_level_helpers_support_canonical_and_topology_workflows() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ])
        .build()?;
    let mapping = LowPassMapping::new(1.0)?;
    let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

    let matrix = mfs::synthesize_canonical_matrix(&polynomials)?;
    let folded = mfs::synthesize_matrix_with_topology(&polynomials, TopologyKind::Folded)?;

    assert_eq!(matrix.order(), 4);
    assert!(matches_folded_pattern(
        &folded,
        MatrixPatternTolerance::default()
    ));
    Ok(())
}
