use mfs::approx::{ApproximationEngine, ChebyshevApproximation, PolynomialSet};
use mfs::error::Result;
use mfs::freq::{FrequencyGrid, LowPassMapping};
use mfs::spec::{FilterSpec, TransmissionZero};
use mfs::synthesis::SectionSynthesis;
use mfs::verify::ResponseTolerance;

#[test]
fn section_synthesis_supports_triplet_workflow() -> Result<()> {
    let spec = FilterSpec::chebyshev(5, 20.0)?
        .with_transmission_zeros(vec![TransmissionZero::normalized(-1.3)]);
    let mapping = LowPassMapping::new(1.0)?;
    let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

    let matrix = SectionSynthesis::default().synthesize_triplet(&polynomials, -1.3, 2)?;
    assert!(matrix.at(3, 5).unwrap_or_default().abs() <= 1e-6);
    Ok(())
}

#[test]
fn section_synthesis_report_supports_triplet_and_trisection_workflows() -> Result<()> {
    let triplet_polynomials = PolynomialSet::new(
        5,
        0.1,
        0.1,
        1.0,
        vec![-1.25],
        vec![1.0, 0.92, 0.84, 0.76, 0.68, 0.6],
        vec![0.95, 0.87, 0.79, 0.71, 0.63, 0.55],
        vec![0.18, -0.07, 0.03],
    )?;
    let sections = SectionSynthesis::default();

    let triplet = sections.synthesize_triplet_with_report(&triplet_polynomials, -1.25, 2)?;
    assert!(triplet.verification.passes());

    let trisection =
        sections.synthesize_trisection_with_report(&triplet_polynomials, -1.25, (2, 4))?;
    assert!(trisection.verification.passes());
    Ok(())
}

#[test]
fn section_synthesis_can_attach_response_summaries() -> Result<()> {
    let polynomials = PolynomialSet::new(
        5,
        0.1,
        0.1,
        1.0,
        vec![-1.25],
        vec![1.0, 0.92, 0.84, 0.76, 0.68, 0.6],
        vec![0.95, 0.87, 0.79, 0.71, 0.63, 0.55],
        vec![0.18, -0.07, 0.03],
    )?;
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
