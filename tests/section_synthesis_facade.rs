use mfs::approx::generalized_chebyshev_polynomials;
use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::spec::FilterSpec;
use mfs::synthesis::SectionSynthesis;
use mfs::verify::ResponseTolerance;

#[test]
fn section_synthesis_supports_triplet_workflow() -> Result<()> {
    let spec = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![-1.3]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;

    let triplet = SectionSynthesis::default().synthesize_triplet(&polynomials, -1.3, 2)?;
    assert!(triplet.verification.passes());
    assert!(triplet.matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
    Ok(())
}

#[test]
fn section_synthesis_report_supports_triplet_and_trisection_workflows() -> Result<()> {
    let spec = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![-1.25]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let sections = SectionSynthesis::default();

    let triplet = sections.synthesize_triplet(&polynomials, -1.25, 2)?;
    assert!(triplet.verification.passes());

    let trisection = sections.synthesize_trisection(&polynomials, -1.25, (2, 4))?;
    assert!(trisection.verification.passes());
    Ok(())
}

#[test]
fn section_synthesis_can_attach_response_summaries() -> Result<()> {
    let spec = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![-1.25]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let grid = FrequencyGrid::linspace(-2.0, 2.0, 41)?;
    let sections = SectionSynthesis::default();

    let triplet = sections.synthesize_triplet_with_response_check(
        &polynomials,
        -1.25,
        2,
        &grid,
        ResponseTolerance::default(),
    )?;
    assert!(triplet.passes());
    assert_eq!(triplet.response.invariant, Some(true));

    let trisection = sections.synthesize_trisection_with_response_check(
        &polynomials,
        -1.25,
        (2, 4),
        &grid,
        ResponseTolerance::default(),
    )?;
    assert!(trisection.passes());
    assert_eq!(trisection.response.invariant, Some(true));
    Ok(())
}
