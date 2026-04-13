use mfs::fixtures::{
    cameron_generalized_order4_info, cameron_generalized_order4_spec,
    cameron_style_section_info, cameron_style_section_polynomials,
    cameron_style_trisection_request, literature_reference_grid,
};
use mfs::prelude::*;

fn main() -> Result<()> {
    let generalized_info = cameron_generalized_order4_info();
    let (spec, _) = cameron_generalized_order4_spec()?;
    let synthesis = generalized_chebyshev(&spec)?;

    println!("fixture: {}", generalized_info.key);
    println!("reference: {}", generalized_info.source);
    println!("approximation stage: {}", synthesis.approximation_kind());
    println!("matrix method: {:?}", synthesis.matrix_method);
    println!("matrix order: {}", synthesis.matrix.order());

    let section_info = cameron_style_section_info();
    let polynomials = cameron_style_section_polynomials()?;
    let (zero, positions) = cameron_style_trisection_request();
    let grid = literature_reference_grid()?;
    let verified = SectionSynthesis::default().synthesize_trisection_with_response_check(
        &polynomials,
        zero,
        positions,
        &grid,
        ResponseTolerance::default(),
    )?;

    println!("fixture: {}", section_info.key);
    println!("reference: {}", section_info.source);
    println!("section verification passed: {}", verified.verification.passes());
    println!("response check passed: {}", verified.response.passes());
    println!("overall section workflow passed: {}", verified.passes());

    Ok(())
}
