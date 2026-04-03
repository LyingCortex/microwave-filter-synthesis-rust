use mfs::{BandPassMapping, ChebyshevSynthesis, FilterParameter, FrequencyGrid, TransmissionZero};

fn main() -> mfs::Result<()> {
    // Define a small generalized Chebyshev-style example in normalized zero space.
    let spec = FilterParameter::chebyshev(6, 23.0)?.with_transmission_zeros(vec![
        TransmissionZero::normalized(-2.0),
        TransmissionZero::normalized(-1.2),
        TransmissionZero::normalized(1.5),
    ]);

    // Map the physical band into the low-pass prototype domain used by synthesis.
    let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 11)?;

    // Run the end-to-end flow and print a few representative outputs.
    let outcome = ChebyshevSynthesis::default()
        .synthesize_and_evaluate_with_mapping(&spec, &mapping, &grid)?;

    println!("order: {}", outcome.polynomials.order);
    println!(
        "normalized transmission zeros: {:?}",
        outcome.polynomials.transmission_zeros_normalized
    );
    println!("matrix shape: {:?}", outcome.matrix.shape());
    println!("sample count: {}", outcome.response.samples.len());
    println!(
        "center sample s21 ~= {} + j{}",
        outcome.response.samples[outcome.response.samples.len() / 2].s21_re,
        outcome.response.samples[outcome.response.samples.len() / 2].s21_im
    );

    Ok(())
}
