use std::fs;
use std::path::PathBuf;

use mfs::approx::{ComplexCoefficient, ComplexPolynomial};
use mfs::fixtures::{
    ComplexValue, load_filter_database_case_from_repo,
    load_filter_database_end_to_end_fixture,
};
use mfs::prelude::render_terminal_filter_database_report;
use mfs::generalized_chebyshev_with_response;

fn main() -> mfs::Result<()> {
    let reference_case = load_filter_database_case_from_repo("Cameron_passband_symmetry_4_2")?;
    let fixture = load_filter_database_end_to_end_fixture("Cameron_passband_symmetry_4_2")?;
    let outcome =
        generalized_chebyshev_with_response(&fixture.spec, &fixture.mapping, &fixture.grid)?;
    let synthesis = &outcome.synthesis;
    let generalized = outcome
        .synthesis
        .polynomials
        .generalized
        .as_ref()
        .ok_or_else(|| mfs::MfsError::Unsupported("expected generalized helper data".to_string()))?;
    let reference_model = reference_case
        .mathematical_model
        .as_ref()
        .ok_or_else(|| mfs::MfsError::Unsupported("fixture is missing mathematical_model".to_string()))?;

    let center_index = outcome.response.samples.len() / 2;
    let center_sample = &outcome.response.samples[center_index];

    let output_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("docs")
        .join("filter_database_case_output.md");
    let actual_e_values = generalized
        .e_s
        .as_ref()
        .map(|poly| complex_coefficients_to_values(&poly.coefficients))
        .unwrap_or_default();
    let actual_e = format_complex_values(&actual_e_values);
    let actual_f_values = complex_coefficients_to_values(&generalized.f_s.coefficients);
    let actual_f = format_complex_values(&actual_f_values);
    let actual_p_values = complex_coefficients_to_values(&generalized.p_s.coefficients);
    let actual_p = format_complex_values(&actual_p_values);
    let reference_e = format_complex_values(&reference_model.polynomial_coefficients.e);
    let reference_f = format_complex_values(&reference_model.polynomial_coefficients.f);
    let reference_p = format_complex_values(&reference_model.polynomial_coefficients.p);
    let actual_e_roots_values = generalized
        .e_stage
        .as_ref()
        .map(|stage| complex_coefficients_to_values(&stage.e_s_roots))
        .unwrap_or_default();
    let actual_e_roots = format_complex_values(&actual_e_roots_values);
    let actual_f_roots_values = generalized
        .f_s
        .roots()?
        .into_iter()
        .map(|root| ComplexValue {
            re: root.re,
            im: root.im,
        })
        .collect::<Vec<_>>();
    let actual_f_roots_values = sort_complex_values(&actual_f_roots_values);
    let actual_f_roots = format_complex_values(&actual_f_roots_values);
    let actual_a_roots_values = generalized
        .a_stage
        .as_ref()
        .map(|stage| complex_coefficients_to_values(&stage.a_s_roots))
        .unwrap_or_default();
    let actual_a_roots = format_complex_values(&actual_a_roots_values);
    let actual_raw_w_roots_values = generalized
        .e_stage
        .as_ref()
        .map(|stage| complex_coefficients_to_values(&stage.raw_roots))
        .unwrap_or_default();
    let actual_raw_w_roots = format_complex_values(&actual_raw_w_roots_values);
    let actual_reflected_w_roots_values = generalized
        .e_stage
        .as_ref()
        .map(|stage| complex_coefficients_to_values(&stage.reflected_roots))
        .unwrap_or_default();
    let actual_reflected_w_roots = format_complex_values(&actual_reflected_w_roots_values);
    let reconstructed_f_from_zeros: Option<Vec<ComplexValue>> = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_zeros.as_ref())
        .map(|zeros| reconstruct_polynomial_from_roots(zeros.as_slice()))
        .transpose()?;
    let reconstructed_e_from_poles: Option<Vec<ComplexValue>> = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_poles.as_ref())
        .map(|poles| reconstruct_polynomial_from_roots(poles.as_slice()))
        .transpose()?;
    let reference_poles = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_poles.as_ref())
        .map(|poles| format_complex_values(&sort_complex_values(poles)))
        .unwrap_or_else(|| "None".to_string());
    let reference_poles_values = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_poles.as_ref())
        .map(|values| sort_complex_values(values))
        .unwrap_or_default();
    let reference_zeros = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_zeros.as_ref())
        .map(|zeros| format_complex_values(&sort_complex_values(zeros)))
        .unwrap_or_else(|| "None".to_string());
    let reference_zeros_values = reference_model
        .singularities
        .as_ref()
        .and_then(|singularities| singularities.reflection_zeros.as_ref())
        .map(|values| sort_complex_values(values))
        .unwrap_or_default();
    let reconstructed_f = reconstructed_f_from_zeros
        .as_ref()
        .map(|values| format_complex_values(values))
        .unwrap_or_else(|| "None".to_string());
    let reconstructed_e = reconstructed_e_from_poles
        .as_ref()
        .map(|values| format_complex_values(values))
        .unwrap_or_else(|| "None".to_string());
    let e_comparison = format_comparison_table(
        "E",
        &reference_model.polynomial_coefficients.e,
        actual_e_values.clone(),
    );
    let reconstructed_e_comparison = reconstructed_e_from_poles
        .as_ref()
        .map(|values| {
            format_comparison_table(
                "E(s) from reflection_poles vs reference E(s)",
                &reference_model.polynomial_coefficients.e,
                values.clone(),
            )
        })
        .unwrap_or_else(|| "### E(s) from reflection_poles vs reference E(s)\nreflection_poles: None".to_string());
    let reconstructed_f_comparison = reconstructed_f_from_zeros
        .as_ref()
        .map(|values| {
            format_comparison_table(
                "F(s) from reflection_zeros vs reference F(s)",
                &reference_model.polynomial_coefficients.f,
                values.clone(),
            )
        })
        .unwrap_or_else(|| "### F(s) from reflection_zeros vs reference F(s)\nreflection_zeros: None".to_string());
    let f_comparison = format_comparison_table(
        "F",
        &reference_model.polynomial_coefficients.f,
        actual_f_values.clone(),
    );
    let p_comparison = format_comparison_table(
        "P",
        &reference_model.polynomial_coefficients.p,
        actual_p_values.clone(),
    );
    let e_summary = summarize_alignment(
        "E(s)",
        &reference_model.polynomial_coefficients.e,
        &actual_e_values,
    );
    let e_from_poles_summary = reconstructed_e_from_poles
        .as_ref()
        .map(|values| summarize_alignment("E(s) reconstructed from reflection_poles", &reference_model.polynomial_coefficients.e, values))
        .unwrap_or_else(|| "E(s) reconstructed from reflection_poles: no reflection_poles in reference".to_string());
    let f_summary = summarize_alignment(
        "F(s)",
        &reference_model.polynomial_coefficients.f,
        &actual_f_values,
    );
    let f_from_zeros_summary = reconstructed_f_from_zeros
        .as_ref()
        .map(|values| summarize_alignment("F(s) reconstructed from reflection_zeros", &reference_model.polynomial_coefficients.f, values))
        .unwrap_or_else(|| "F(s) reconstructed from reflection_zeros: no reflection_zeros in reference".to_string());
    let p_summary = summarize_alignment(
        "P(s)",
        &reference_model.polynomial_coefficients.p,
        &actual_p_values,
    );
    let report = format!(
        "# Filter Database Case Output\n\n\
case_id: `{}`\n\n\
## Fixture\n\n\
- order: `{}`\n\
- return_loss_db: `{}`\n\
- center_hz: `{}`\n\
- bandwidth_hz: `{}`\n\
- normalized_transmission_zeros: `{:?}`\n\n\
## Transmission-Zero Convention\n\n\
- fixture specs store normalized prototype transmission zeros\n\
- the band-pass mapping is used for the physical response grid, not for zero normalization inside `FilterSpec`\n\n\
## Synthesis\n\n\
- approximation_kind: `{}`\n\
- matrix_method: `{:?}`\n\
- matrix_order: `{}`\n\
- response_samples: `{}`\n\
- eps: `{}`\n\
- eps_r: `{}`\n\n\
## Coefficient Ordering\n\n\
- all polynomial coefficient lists in this report use **ascending powers**\n\
- index `0` is the constant term\n\
- the last entry is the highest-order coefficient\n\
- for this Cameron case, monic polynomials therefore end with `+1`\n\n\
## PolynomialSet\n\n\
```text\n\
e = {}\n\
f = {}\n\
p = {}\n\
```\n\n\
## Reference Mathematical Model\n\n\
```text\n\
reference E = {}\n\
reference F = {}\n\
reference P = {}\n\
reconstructed F from reflection_zeros = {}\n\
reconstructed E from reflection_poles = {}\n\
```\n\n\
## Actual Generalized Helper Polynomials\n\n\
```text\n\
actual e_s = {}\n\
actual f_s = {}\n\
actual p_s = {}\n\
```\n\n\
## Coefficient Comparison\n\n\
{}\n\
{}\n\
{}\n\
{}\n\
{}\n\
\n\
## Root Comparison\n\n\
```text\n\
reference reflection zeros (roots of F(s)) = {}\n\
actual f_s roots = {}\n\
\n\
actual raw e_w roots = {}\n\
actual reflected e_w roots = {}\n\
actual a_s roots = {}\n\
reference reflection poles (roots of E(s)) = {}\n\
actual e_s roots = {}\n\
```\n\n\
## Comparison Summary\n\n\
- {}\n\
- {}\n\
- {}\n\
- {}\n\
- {}\n\
- reference E length: `{}`\n\
- actual e_s length: `{}`\n\
- reference F length: `{}`\n\
- actual f_s length: `{}`\n\
- reference P length: `{}`\n\
- actual p_s length: `{}`\n\n\
## Center Response Sample\n\n\
```text\n\
frequency_hz = {}\n\
normalized_omega = {}\n\
s11 = {}\n\
s21 = {}\n\
```\n",
        fixture.case_id,
        fixture.spec.order,
        fixture.spec.return_loss_db(),
        fixture.mapping.center_hz(),
        fixture.mapping.bandwidth_hz(),
        synthesis.polynomials.transmission_zeros_normalized,
        synthesis.approximation_kind(),
        synthesis.matrix_method,
        synthesis.matrix.order(),
        outcome.response.samples.len(),
        synthesis.polynomials.eps,
        synthesis.polynomials.eps_r,
        format_complex_values(&complex_coefficients_to_values(&synthesis.polynomials.e.coefficients)),
        format_complex_values(&complex_coefficients_to_values(&synthesis.polynomials.f.coefficients)),
        format_complex_values(&complex_coefficients_to_values(&synthesis.polynomials.p.coefficients)),
        reference_e,
        reference_f,
        reference_p,
        reconstructed_f,
        reconstructed_e,
        actual_e,
        actual_f,
        actual_p,
        e_comparison,
        reconstructed_e_comparison,
        reconstructed_f_comparison,
        f_comparison,
        p_comparison,
        reference_zeros,
        actual_f_roots,
        actual_raw_w_roots,
        actual_reflected_w_roots,
        actual_a_roots,
        reference_poles,
        actual_e_roots,
        e_summary,
        e_from_poles_summary,
        f_summary,
        f_from_zeros_summary,
        p_summary,
        reference_model.polynomial_coefficients.e.len(),
        generalized.e_s.as_ref().map_or(0, |poly| poly.coefficients.len()),
        reference_model.polynomial_coefficients.f.len(),
        generalized.f_s.coefficients.len(),
        reference_model.polynomial_coefficients.p.len(),
        generalized.p_s.coefficients.len(),
        center_sample.frequency_hz,
        center_sample.normalized_omega,
        format_complex_scalar(ComplexValue {
            re: center_sample.s11_re,
            im: center_sample.s11_im,
        }),
        format_complex_scalar(ComplexValue {
            re: center_sample.s21_re,
            im: center_sample.s21_im,
        }),
    );
    fs::write(&output_path, report).map_err(|error| {
        mfs::MfsError::Unsupported(format!(
            "failed to write example report {}: {error}",
            output_path.display()
        ))
    })?;

    let terminal_report = render_terminal_filter_database_report(
        &fixture.case_id,
        &fixture,
        &outcome,
        reference_model,
        &reference_model.polynomial_coefficients.e,
        &actual_e_values,
        &reference_model.polynomial_coefficients.f,
        &actual_f_values,
        &reference_model.polynomial_coefficients.p,
        &actual_p_values,
        &reference_zeros_values,
        &actual_f_roots_values,
        &reference_poles_values,
        &actual_e_roots_values,
        &actual_raw_w_roots_values,
        &actual_reflected_w_roots_values,
        &actual_a_roots_values,
        &e_summary,
        &e_from_poles_summary,
        &f_summary,
        &f_from_zeros_summary,
        &p_summary,
    );
    print!("{terminal_report}");
    println!("report written to {}", output_path.display());

    Ok(())
}

fn complex_coefficients_to_values(coefficients: &[num_complex::Complex64]) -> Vec<ComplexValue> {
    coefficients
        .iter()
        .map(|coefficient| ComplexValue {
            re: coefficient.re,
            im: coefficient.im,
        })
        .collect()
}

fn format_complex_values(values: &[ComplexValue]) -> String {
    let entries = values
        .iter()
        .map(|value| format_complex_scalar(value.clone()))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{}]", entries)
}

fn format_complex_scalar(value: ComplexValue) -> String {
    let tolerance = 1e-12;
    let real = round_to_4(if value.re.abs() < tolerance { 0.0 } else { value.re });
    let imag = round_to_4(if value.im.abs() < tolerance { 0.0 } else { value.im });

    if imag == 0.0 {
        return format_signed_number(real);
    }
    if real == 0.0 {
        return format!(
            "{sign}j{}",
            format_unsigned_number(imag.abs()),
            sign = if imag.is_sign_negative() { "-" } else { "+" }
        );
    }

    format!(
        "{} {sign} j{}",
        format_signed_number(real),
        format_unsigned_number(imag.abs()),
        sign = if imag.is_sign_negative() { "-" } else { "+" }
    )
}

fn round_to_4(value: f64) -> f64 {
    let rounded = (value * 10_000.0).round() / 10_000.0;
    if rounded.abs() < 1e-12 { 0.0 } else { rounded }
}

fn format_signed_number(value: f64) -> String {
    if value.is_sign_negative() {
        format!("-{}", format_unsigned_number(value.abs()))
    } else {
        format!("+{}", format_unsigned_number(value))
    }
}

fn format_unsigned_number(value: f64) -> String {
    let rounded = round_to_4(value);
    if (rounded.fract()).abs() < 1e-12 {
        return format!("{}", rounded as i64);
    }

    let mut text = format!("{rounded:.4}");
    while text.contains('.') && text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    text
}

fn format_comparison_table(
    name: &str,
    reference: &[ComplexValue],
    actual: Vec<ComplexValue>,
) -> String {
    let row_count = reference.len().max(actual.len());
    let mut rows = Vec::with_capacity(row_count + 3);
    rows.push(format!("### {}(s)", name));
    rows.push(String::from("| idx | reference | actual | delta |"));
    rows.push(String::from("| --- | --- | --- | --- |"));

    for index in 0..row_count {
        let reference_value = reference.get(index).cloned();
        let actual_value = actual.get(index).cloned();
        let delta = difference(reference_value.clone(), actual_value.clone());
        rows.push(format!(
            "| {} | {} | {} | {} |",
            index,
            reference_value
                .map(format_complex_scalar)
                .unwrap_or_else(|| "None".to_string()),
            actual_value
                .map(format_complex_scalar)
                .unwrap_or_else(|| "None".to_string()),
            delta
                .map(format_complex_scalar)
                .unwrap_or_else(|| "None".to_string()),
        ));
    }

    rows.join("\n")
}

fn difference(reference: Option<ComplexValue>, actual: Option<ComplexValue>) -> Option<ComplexValue> {
    match (reference, actual) {
        (Some(reference), Some(actual)) => Some(ComplexValue {
            re: actual.re - reference.re,
            im: actual.im - reference.im,
        }),
        (None, Some(actual)) => Some(actual),
        (Some(reference), None) => Some(ComplexValue {
            re: -reference.re,
            im: -reference.im,
        }),
        (None, None) => None,
    }
}

fn summarize_alignment(name: &str, reference: &[ComplexValue], actual: &[ComplexValue]) -> String {
    if reference.len() != actual.len() {
        return format!(
            "{}: length mismatch (reference {}, actual {}), needs review",
            name,
            reference.len(),
            actual.len()
        );
    }

    let max_delta = reference
        .iter()
        .zip(actual.iter())
        .map(|(reference, actual)| {
            let delta_re = actual.re - reference.re;
            let delta_im = actual.im - reference.im;
            (delta_re * delta_re + delta_im * delta_im).sqrt()
        })
        .fold(0.0_f64, f64::max);

    let verdict = if max_delta <= 1e-3 {
        "basic match"
    } else if max_delta <= 1e-1 {
        "close but worth checking"
    } else {
        "clear mismatch"
    };

    format!(
        "{}: {}, max |delta| ~= {}",
        name,
        verdict,
        format_unsigned_number(max_delta)
    )
}

fn reconstruct_polynomial_from_roots(roots: &[ComplexValue]) -> mfs::Result<Vec<ComplexValue>> {
    let roots = roots
        .iter()
        .map(|root| ComplexCoefficient::new(root.re, root.im))
        .collect::<Vec<_>>();
    let polynomial = ComplexPolynomial::from_complex_roots(&roots)?;
    Ok(complex_coefficients_to_values(&polynomial.coefficients))
}

fn sort_complex_values(values: &[ComplexValue]) -> Vec<ComplexValue> {
    let mut sorted = values.to_vec();
    sorted.sort_by(|lhs, rhs| {
        lhs.re
            .partial_cmp(&rhs.re)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                lhs.im
                    .partial_cmp(&rhs.im)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    sorted
}
