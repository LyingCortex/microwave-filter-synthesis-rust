use mfs::approx::{
    ComplexCoefficient, ComplexPolynomial, cameron_recursive, find_eps, find_p_polynomial,
    pad_transmission_zeros, synthesize_generalized_chebyshev_data,
};
use mfs::error::Result;
use mfs::fixtures::{
    cameron_order3_generalized_pipeline_exact_case,
    cameron_order3_generalized_pipeline_exact_info, cameron_single_zero_exact_case,
    cameron_single_zero_exact_info, load_filter_database_case_from_repo,
    load_filter_database_end_to_end_fixture,
    cameron_style_section_info, cameron_style_section_polynomials,
    cameron_style_trisection_request, cameron_style_triplet_request,
    literature_reference_grid, literature_reference_grid_info,
};
use mfs::synthesis::{MatrixSynthesisMethod, SectionSynthesis};
use mfs::generalized_chebyshev;
use mfs::verify::ResponseTolerance;

fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
    let diff = (lhs - rhs).abs();
    assert!(
        diff <= tol,
        "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
    );
}

#[test]
fn literature_generalized_fixture_drives_generalized_main_flow() -> Result<()> {
    let case = load_filter_database_case_from_repo("Cameron_passband_symmetry_4_2")?;
    let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
    let outcome = generalized_chebyshev(&fixture.spec)?;

    assert!(case.case_id.contains("Cameron"));
    assert_eq!(outcome.approximation_kind(), "GeneralizedChebyshev");
    assert_eq!(outcome.matrix_method, MatrixSynthesisMethod::ResidueExpansion);
    Ok(())
}

#[test]
fn literature_triplet_fixture_supports_reported_response_checked_synthesis() -> Result<()> {
    let info = cameron_style_section_info();
    let grid_info = literature_reference_grid_info();
    let polynomials = cameron_style_section_polynomials()?;
    let (zero, center) = cameron_style_triplet_request();
    let grid = literature_reference_grid()?;

    let outcome = SectionSynthesis::default().synthesize_triplet_with_response_check(
        &polynomials,
        zero,
        center,
        &grid,
        ResponseTolerance::default(),
    )?;

    assert!(info.source.contains("ZTE"));
    assert!(grid_info
        .expected_behavior
        .contains(&"provides repeatable normalized-sweep regression coverage"));
    assert!(outcome.passes());
    assert_eq!(outcome.response.invariant, Some(true));
    Ok(())
}

#[test]
fn literature_trisection_fixture_supports_reported_response_checked_synthesis() -> Result<()> {
    let polynomials = cameron_style_section_polynomials()?;
    let (zero, positions) = cameron_style_trisection_request();
    let grid = literature_reference_grid()?;

    let outcome = SectionSynthesis::default().synthesize_trisection_with_response_check(
        &polynomials,
        zero,
        positions,
        &grid,
        ResponseTolerance::default(),
    )?;

    assert!(outcome.passes());
    assert_eq!(outcome.response.invariant, Some(true));
    Ok(())
}

#[test]
fn exact_cameron_single_zero_fixture_matches_recurrence_and_eps_values() -> Result<()> {
    let info = cameron_single_zero_exact_info();
    let fixture = cameron_single_zero_exact_case();

    assert!(info.source.contains("Cameron"));
    assert!(info
        .expected_behavior
        .contains(&"anchors the helper epsilon computation to a fixed numeric case"));

    let recurrence = cameron_recursive(&[fixture.finite_zero])?;
    assert_eq!(recurrence.u_descending, fixture.expected_u_descending);
    approx_eq(recurrence.v_descending[0], fixture.expected_v0, 1e-12);
    approx_eq(recurrence.f_s.coefficients[0].re, fixture.expected_f_constant.0, 1e-12);
    approx_eq(recurrence.f_s.coefficients[0].im, fixture.expected_f_constant.1, 1e-12);
    approx_eq(recurrence.f_s.coefficients[1].re, fixture.expected_f_linear.0, 1e-12);
    approx_eq(recurrence.f_s.coefficients[1].im, fixture.expected_f_linear.1, 1e-12);

    let padded = pad_transmission_zeros(3, &[fixture.finite_zero])?;
    let p_s = find_p_polynomial(3, &padded.padded, padded.finite_count)?;
    approx_eq(p_s.coefficients[0].re, fixture.expected_p_constant.0, 1e-12);
    approx_eq(p_s.coefficients[0].im, fixture.expected_p_constant.1, 1e-12);
    approx_eq(p_s.coefficients[1].re, fixture.expected_p_linear.0, 1e-12);
    approx_eq(p_s.coefficients[1].im, fixture.expected_p_linear.1, 1e-12);

    let f_s = ComplexPolynomial::new(vec![
        ComplexCoefficient::new(
            fixture.expected_f_constant.0,
            fixture.expected_f_constant.1,
        ),
        ComplexCoefficient::new(fixture.expected_f_linear.0, fixture.expected_f_linear.1),
    ])?;
    let p_unit = ComplexPolynomial::new(vec![ComplexCoefficient::new(1.0, 0.0)])?;
    let (eps, eps_r) = find_eps(0, &p_unit, &f_s, fixture.return_loss_db, 2)?;
    approx_eq(eps, fixture.expected_eps, 1e-12);
    approx_eq(eps_r, fixture.expected_eps_r, 1e-12);
    Ok(())
}

#[test]
fn exact_order3_generalized_pipeline_fixture_matches_helper_outputs() -> Result<()> {
    let info = cameron_order3_generalized_pipeline_exact_info();
    let fixture = cameron_order3_generalized_pipeline_exact_case();

    assert!(info.source.contains("Cameron"));
    assert!(info
        .expected_behavior
        .contains(&"stabilizes F(s), P(s), A(s), and E(s) against accidental drift"));

    let data = synthesize_generalized_chebyshev_data(
        fixture.order,
        &[fixture.finite_zero],
        fixture.return_loss_db,
    )?;

    assert_eq!(data.padded_zeros, fixture.expected_padded_zeros);
    approx_eq(data.eps, fixture.expected_eps, 1e-12);
    approx_eq(data.eps_r, fixture.expected_eps_r, 1e-12);

    for (actual, expected) in data.f_s.coefficients.iter().zip(fixture.expected_f_s.iter()) {
        approx_eq(actual.re, expected.0, 1e-12);
        approx_eq(actual.im, expected.1, 1e-12);
    }
    for (actual, expected) in data.p_s.coefficients.iter().zip(fixture.expected_p_s.iter()) {
        approx_eq(actual.re, expected.0, 1e-12);
        approx_eq(actual.im, expected.1, 1e-12);
    }

    let a_s = data.a_s.as_ref().expect("fixture should provide A(s)");
    let a_stage = data.a_stage.as_ref().expect("fixture should provide A-stage details");
    for (actual, expected) in a_s.coefficients.iter().zip(fixture.expected_a_s.iter()) {
        approx_eq(actual.re, expected.0, 1e-12);
        approx_eq(actual.im, expected.1, 1e-12);
    }
    assert_eq!(a_stage.a_s_roots.len(), 1);
    assert_eq!(a_stage.a_w_roots.len(), 1);

    let e_s = data.e_s.as_ref().expect("fixture should provide E(s)");
    let e_stage = data.e_stage.as_ref().expect("fixture should provide E-stage details");
    for (actual, expected) in e_s.coefficients.iter().zip(fixture.expected_e_s.iter()) {
        approx_eq(actual.re, expected.0, 1e-12);
        approx_eq(actual.im, expected.1, 1e-12);
    }
    assert_eq!(e_stage.raw_roots.len(), e_stage.reflected_roots.len());
    assert_eq!(e_stage.e_s_roots.len(), e_stage.reflected_roots.len());

    Ok(())
}
