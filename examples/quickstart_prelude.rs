use mfs::prelude::*;

fn main() -> Result<()> {
    let spec = FilterSpec::builder()
        .order(4)
        .return_loss_db(20.0)
        .transmission_zeros(vec![TransmissionZero::normalized(-1.5)])
        .build()?;

    let mapping = LowPassMapping::new(1.0)?;
    let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;
    let matrix = CanonicalMatrixSynthesis::default().synthesize(&polynomials)?;
    let grid = FrequencyGrid::linspace(0.5, 1.5, 21)?;
    let transform = transform_matrix_with_response_check(
        &matrix,
        TopologyKind::Folded,
        &grid,
        ResponseTolerance::default(),
    )?;
    let transformed = transform.matrix;
    let response = ResponseSolver::default().evaluate_normalized(&transformed, &grid)?;

    println!("order: {}", transformed.order());
    println!("samples: {}", response.samples.len());
    println!("first normalized omega: {}", response.samples[0].normalized_omega);
    println!("transform passed: {}", transform.report.passes());
    Ok(())
}
