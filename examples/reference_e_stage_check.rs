use mfs::approx::{ComplexCoefficient, ComplexPolynomial, build_e_polynomial_stage};
use mfs::fixtures::load_filter_database_case_from_repo;
use num_complex::Complex64;

fn main() -> mfs::Result<()> {
    let case = load_filter_database_case_from_repo("Cameron_passband_symmetry_4_2")?;
    let model = case
        .mathematical_model
        .as_ref()
        .ok_or_else(|| mfs::MfsError::Unsupported("fixture missing mathematical_model".to_string()))?;
    let singularities = model
        .singularities
        .as_ref()
        .ok_or_else(|| mfs::MfsError::Unsupported("fixture missing singularities".to_string()))?;

    let f_s = ComplexPolynomial::new(
        model.polynomial_coefficients.f.iter().map(to_complex).collect(),
    )?;
    let p_s = ComplexPolynomial::new(
        model.polynomial_coefficients.p.iter().map(to_complex).collect(),
    )?;
    let stage = build_e_polynomial_stage(
        &f_s,
        &p_s,
        singularities.epsilon.unwrap_or(0.0),
        singularities.epsilon_r.unwrap_or(1.0),
    )?;

    println!("reference F = {:?}", f_s.coefficients);
    println!("reference P = {:?}", p_s.coefficients);
    println!("computed E = {:?}", stage.e_s.coefficients);
    println!("raw roots = {:?}", stage.raw_roots);
    println!("reflected roots = {:?}", stage.reflected_roots);
    println!("e_s roots = {:?}", stage.e_s_roots);

    let f_w = ComplexPolynomial::new(rotate_s_to_w(&f_s.coefficients))?;
    let p_w = ComplexPolynomial::new(rotate_s_to_w(&p_s.coefficients))?;
    let eps = singularities.epsilon.unwrap_or(0.0);
    let plus_j = f_w.add(&p_w.scale(ComplexCoefficient::new(0.0, 1.0 / eps))?)?;
    let minus_j = f_w.add(&p_w.scale(ComplexCoefficient::new(0.0, -1.0 / eps))?)?;
    let plus_real = f_w.add(&p_w.scale(ComplexCoefficient::new(1.0 / eps, 0.0))?)?;

    println!("plus_j roots = {:?}", plus_j.roots()?);
    println!("minus_j roots = {:?}", minus_j.roots()?);
    println!("plus_real roots = {:?}", plus_real.roots()?);

    Ok(())
}

fn to_complex(value: &mfs::fixtures::ComplexValue) -> ComplexCoefficient {
    ComplexCoefficient::new(value.re, value.im)
}

fn rotate_s_to_w(coefficients: &[Complex64]) -> Vec<Complex64> {
    coefficients
        .iter()
        .copied()
        .enumerate()
        .map(|(index, coefficient)| mul_i_pow(coefficient, index))
        .collect()
}

fn mul_i_pow(value: Complex64, power: usize) -> Complex64 {
    match power % 4 {
        0 => value,
        1 => Complex64::new(-value.im, value.re),
        2 => Complex64::new(-value.re, -value.im),
        _ => Complex64::new(value.im, -value.re),
    }
}
