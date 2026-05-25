use mfs::error::Result;
use mfs::fixtures::load_filter_database_end_to_end_fixture;
use mfs::generalized_chebyshev_with_response;

#[test]
fn filter_database_case_can_drive_end_to_end_synthesis_flow() -> Result<()> {
    let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
    let outcome = generalized_chebyshev_with_response(&fixture.spec, &fixture.mapping, &fixture.grid)?;

    assert_eq!(fixture.case_id, "Cameron_passband_symmetry_4_2");
    assert_eq!(fixture.spec.order, 4);
    assert_eq!(fixture.spec.transmission_zeros.len(), 2);
    assert_eq!(outcome.synthesis.polynomials.order, fixture.spec.order);
    assert_eq!(outcome.synthesis.approximation, mfs::synthesis::ApproximationKind::GeneralizedChebyshev);
    assert_eq!(outcome.synthesis.matrix.order(), fixture.spec.order);
    assert_eq!(outcome.response.samples.len(), 21);
    Ok(())
}
