use mfs::error::Result;
use mfs::fixtures::load_filter_database_end_to_end_fixture;
use mfs::synthesis::ApproximationStageKind;
use mfs::synthesize_and_evaluate_chebyshev_with_details;

#[test]
fn filter_database_case_can_drive_end_to_end_synthesis_flow() -> Result<()> {
    let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
    let outcome = synthesize_and_evaluate_chebyshev_with_details(
        &fixture.spec,
        &fixture.mapping,
        &fixture.grid,
    )?;

    assert_eq!(fixture.case_id, "Cameron_passband_symmetry_4_2");
    assert_eq!(fixture.spec.order, 4);
    assert_eq!(fixture.spec.transmission_zeros.len(), 2);
    assert_eq!(outcome.polynomials.order, fixture.spec.order);
    assert_eq!(
        outcome.approximation_kind,
        ApproximationStageKind::GeneralizedChebyshev
    );
    assert_eq!(outcome.matrix.order(), fixture.spec.order);
    assert_eq!(outcome.response.samples.len(), 21);
    Ok(())
}
