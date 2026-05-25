use mfs::{
    bandpass, filter_spec, generalized_chebyshev_with_response, normalize_transmission_zeros_hz,
    FrequencyGrid,
};

fn main() -> mfs::Result<()> {
    // `FilterSpec` stores normalized prototype zeros, so physical Hz zeros
    // must be mapped into the prototype domain before building the spec.
    let mapping = bandpass(6.75e9, 300.0e6)?;
    let zeros = normalize_transmission_zeros_hz([6.4e9, 6.5e9, 7.0e9], &mapping)?;
    let spec = filter_spec(6, 23.0, zeros, None)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 11)?;

    // Run the end-to-end flow and print a few representative outputs.
    let outcome = generalized_chebyshev_with_response(&spec, &mapping, &grid)?;
    let synthesis = &outcome.synthesis;

    println!("order: {}", synthesis.polynomials.order);
    println!(
        "normalized transmission zeros: {:?}",
        synthesis.polynomials.transmission_zeros_normalized
    );
    println!("matrix shape: {:?}", synthesis.matrix.shape());
    println!("sample count: {}", outcome.response.samples.len());
    println!(
        "center sample s21 ~= {} + j{}",
        outcome.response.samples[outcome.response.samples.len() / 2].s21_re,
        outcome.response.samples[outcome.response.samples.len() / 2].s21_im
    );

    Ok(())
}
