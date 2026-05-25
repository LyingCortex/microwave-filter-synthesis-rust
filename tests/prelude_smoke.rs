use mfs::prelude::*;

#[test]
fn prelude_supports_basic_builder_and_synthesis_workflow() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .normalized_transmission_zeros(vec![-1.5])
        .build()?;
    let synthesis = generalized_chebyshev(&spec)?;

    assert_eq!(synthesis.matrix.order(), 4);
    Ok(())
}

#[test]
fn prelude_supports_explicit_infinite_zero_workflow() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .transmission_zeros(vec![
            TransmissionZero::infinite(),
            TransmissionZero::infinite(),
            TransmissionZero::infinite(),
            TransmissionZero::infinite(),
        ])
        .build()?;
    let synthesis = generalized_chebyshev(&spec)?;
    let polynomials = synthesis.polynomials;

    assert_eq!(synthesis.matrix.order(), 4);
    assert!(polynomials
        .transmission_zeros_normalized
        .iter()
        .all(|zero| zero.is_infinite()));
    Ok(())
}
