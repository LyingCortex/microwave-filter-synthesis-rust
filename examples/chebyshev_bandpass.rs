use mfs::{BandPassPlan, ChebyshevSynthesis, FilterSpec, FrequencyGrid, TransmissionZero};

fn main() -> mfs::Result<()> {
    let spec = FilterSpec::chebyshev(6, 23.0)?.with_transmission_zeros(vec![
        TransmissionZero::normalized(-2.0),
        TransmissionZero::normalized(-1.2),
        TransmissionZero::normalized(1.5),
    ]);

    let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
    let grid = FrequencyGrid::linspace(6.0e9, 7.5e9, 11)?;

    let outcome = ChebyshevSynthesis::default().synthesize_and_evaluate(&spec, &plan, &grid)?;

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
