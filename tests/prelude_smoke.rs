use mfs::prelude::*;

#[test]
fn prelude_supports_basic_builder_and_synthesis_workflow() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .transmission_zeros(vec![TransmissionZero::normalized(-1.5)])
        .build()?;
    let mapping = LowPassMapping::new(1.0)?;
    let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;
    let matrix = CanonicalMatrixSynthesis::default().synthesize(&polynomials)?;

    assert_eq!(matrix.order(), 4);
    Ok(())
}
