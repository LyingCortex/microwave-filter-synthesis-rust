
use crate::error::Result;
use crate::fixtures::{ComplexValue, EndToEndFixture, MathematicalModel};
use crate::output::{
    format_aligned_summary_with_width, format_box_table, format_complex_scalar_parts,
    format_decimal_scalar, format_markdown_table,
    format_matrix_table_data, format_out_of_band_attenuation_table_data,
    format_polynomial_table_data, format_response_samples_table_data,
    format_reference_actual_polynomial_table_data, format_root_comparison_table_data,
    format_singularity_table_data,
};
use crate::response::SParameterResponse;
use crate::spec::FilterSpec;
use crate::synthesis::{EvaluationOutcome, SynthesisOutcome};
use crate::transform::TransformOutcome;

pub fn render_terminal_synthesis_report(
    spec: &FilterSpec,
    synthesis: &SynthesisOutcome,
    transform: &TransformOutcome,
    response: &SParameterResponse,
    debug: bool,
) -> Result<String> {
    let polynomials = &synthesis.polynomials;
    let matrix = &synthesis.matrix;
    let transformed = &transform.matrix;
    let mut singularity_table =
        format_singularity_table_data(spec, polynomials.generalized.as_ref())?;
    let mut response_table = format_response_samples_table_data(&response.samples);
    let out_of_band_attenuation_table = spec
        .out_of_band_attenuation
        .as_ref()
        .map(format_out_of_band_attenuation_table_data);
    let mut spec_summary_rows = vec![
        ("order", spec.order.to_string()),
        ("return_loss_db", spec.return_loss_db().to_string()),
        (
            "normalized_transmission_zeros",
            format!("{:?}", polynomials.transmission_zeros_normalized),
        ),
    ];
    if let Some(unloaded_q) = spec.unloaded_q {
        spec_summary_rows.push(("unloaded_q", format_decimal_scalar(unloaded_q)));
    }
    if let Some(out_of_band_attenuation) = spec.out_of_band_attenuation.as_ref() {
        spec_summary_rows.push((
            "out_of_band_windows",
            out_of_band_attenuation.windows.len().to_string(),
        ));
    }
    let shared_summary_width = [
        "order",
        "return_loss_db",
        "normalized_transmission_zeros",
        "unloaded_q",
        "out_of_band_windows",
        "ripple_factor",
        "eps",
        "eps_r",
        "approximation_kind",
        "matrix_method",
        "topology",
        "shape",
        "requested_topology",
        "source_topology",
        "result_topology",
        "pattern_verified",
        "response_check_passed",
        "overall_passed",
        "sample_count",
    ]
    .iter()
    .map(|item| item.len())
    .max()
    .unwrap_or(0);
    singularity_table.headers = vec![
        "i".to_string(),
        "Reflection Zeros\n(Roots of F(s))".to_string(),
        "Transmission Zeros\n(Prescribed)".to_string(),
        "Transmission/Reflection Poles\n(Roots of E(s))".to_string(),
    ];
    response_table.headers = vec![
        "sample".to_string(),
        "frequency\n(Hz)".to_string(),
        "normalized\nomega".to_string(),
        "s11".to_string(),
        "s21".to_string(),
        "group\ndelay".to_string(),
    ];
    let polynomial_block = if debug {
        let polynomial_table = format_polynomial_table_data(polynomials);
        format!(
            "{polynomial_title}\n\
{polynomial_view}\n\n\
",
            polynomial_title = indent_line("Transfer and Reflection Function Polynomials", 3),
            polynomial_view = indent_block(&format_box_table(&polynomial_table), 3),
        )
    } else {
        String::new()
    };

    Ok(format!(
        "{spec_title}\n\
{spec_summary}\n\n\
{out_of_band_block}\
{prototype_title}\n\
{prototype_summary}\n\n\
{polynomial_block}\
{singularity_title}\n\
{singularity_view}\n\n\
{synthesis_title}\n\
{synthesis_summary}\n\n\
{matrix_title}\n\
{matrix_view}\n\n\
{transform_title}\n\
{transform_summary}\n\
{notes_block}\n\
{transformed_title}\n\
{transformed_view}\n\n\
{response_title}\n\
{response_summary}\n\n\
{response_samples_title}\n\
{response_table}\n",
        spec_title = format_numbered_heading(1, "Spec"),
        spec_summary = format_aligned_summary_with_width(&spec_summary_rows, shared_summary_width),
        out_of_band_block = format_optional_out_of_band_terminal_block(
            out_of_band_attenuation_table.as_ref()
        ),
        prototype_title = format_numbered_heading(2, "Prototype"),
        prototype_summary = format_aligned_summary_with_width(&[
            ("order", polynomials.order.to_string()),
            ("ripple_factor", format_decimal_scalar(polynomials.ripple_factor)),
            ("eps", format_decimal_scalar(polynomials.eps)),
            ("eps_r", format_decimal_scalar(polynomials.eps_r)),
        ], shared_summary_width),
        polynomial_block = polynomial_block,
        singularity_title = indent_line("Corresponding Singularities", 3),
        singularity_view = indent_block(&format_box_table(&singularity_table), 3),
        synthesis_title = format_numbered_heading(3, "Canonical Synthesis"),
        synthesis_summary = format_aligned_summary_with_width(&[
            ("approximation_kind", synthesis.approximation_kind().to_string()),
            ("matrix_method", format!("{:?}", synthesis.matrix_method)),
            ("topology", format!("{:?}", matrix.topology())),
            ("shape", format!("{:?}", matrix.shape())),
        ], shared_summary_width),
        matrix_title = indent_line("Canonical Coupling Matrix", 3),
        matrix_view = indent_block(&format_box_table(&format_matrix_table_data(matrix)), 3),
        transform_title = format_numbered_heading(4, "Transform"),
        transform_summary = format_aligned_summary_with_width(&[
            ("requested_topology", format!("{:?}", transform.topology)),
            ("source_topology", format!("{:?}", transform.report.source_topology)),
            ("result_topology", format!("{:?}", transform.report.result_topology)),
            ("pattern_verified", transform.report.pattern_verified.to_string()),
            (
                "response_check_passed",
                transform.report.response.passes().to_string(),
            ),
            ("overall_passed", transform.report.passes().to_string()),
        ], shared_summary_width),
        notes_block = format_notes_block(&transform.report.notes),
        transformed_title = indent_line("Transformed Coupling Matrix", 3),
        transformed_view = indent_block(&format_box_table(&format_matrix_table_data(transformed)), 3),
        response_title = format_numbered_heading(5, "Response Summary"),
        response_summary = format_aligned_summary_with_width(&[(
            "sample_count",
            response.samples.len().to_string(),
        )], shared_summary_width),
        response_samples_title = indent_line("Response Samples", 3),
        response_table = indent_block(&format_box_table(&response_table), 3),
    ))
}

pub fn print_terminal_synthesis_report(
    spec: &FilterSpec,
    synthesis: &SynthesisOutcome,
    transform: &TransformOutcome,
    response: &SParameterResponse,
    debug: bool,
) -> Result<()> {
    print!(
        "{}",
        render_terminal_synthesis_report(spec, synthesis, transform, response, debug)?
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn render_terminal_filter_database_report(
    case_id: &str,
    fixture: &EndToEndFixture,
    outcome: &EvaluationOutcome,
    reference_model: &MathematicalModel,
    reference_e: &[ComplexValue],
    actual_e: &[ComplexValue],
    reference_f: &[ComplexValue],
    actual_f: &[ComplexValue],
    reference_p: &[ComplexValue],
    actual_p: &[ComplexValue],
    reference_reflection_zeros: &[ComplexValue],
    actual_f_roots: &[ComplexValue],
    reference_reflection_poles: &[ComplexValue],
    actual_e_roots: &[ComplexValue],
    actual_raw_w_roots: &[ComplexValue],
    actual_reflected_w_roots: &[ComplexValue],
    actual_a_roots: &[ComplexValue],
    e_summary: &str,
    e_from_poles_summary: &str,
    f_summary: &str,
    f_from_zeros_summary: &str,
    p_summary: &str,
) -> String {
    let synthesis = &outcome.synthesis;
    let shared_summary_width = [
        "case_id",
        "order",
        "return_loss_db",
        "center_hz",
        "bandwidth_hz",
        "normalized_transmission_zeros",
        "approximation_kind",
        "matrix_method",
        "matrix_order",
        "response_samples",
        "eps",
        "eps_r",
        "reference_e_length",
        "actual_e_length",
        "reference_f_length",
        "actual_f_length",
        "reference_p_length",
        "actual_p_length",
        "frequency_hz",
        "normalized_omega",
        "s11",
        "s21",
    ]
    .iter()
    .map(|item| item.len())
    .max()
    .unwrap_or(0);
    let center_index = outcome.response.samples.len() / 2;
    let center_sample = &outcome.response.samples[center_index];
    let mut coefficient_table = format_reference_actual_polynomial_table_data(
        reference_e,
        actual_e,
        reference_f,
        actual_f,
        reference_p,
        actual_p,
    );
    let mut root_table = format_root_comparison_table_data(
        reference_reflection_zeros,
        actual_f_roots,
        reference_reflection_poles,
        actual_e_roots,
    );
    coefficient_table.headers = vec![
        "i".to_string(),
        "reference\nE(s)".to_string(),
        "actual\nE(s)".to_string(),
        "reference\nF(s)".to_string(),
        "actual\nF(s)".to_string(),
        "reference\nP(s)".to_string(),
        "actual\nP(s)".to_string(),
    ];
    root_table.headers = vec![
        "i".to_string(),
        "reference\nreflection zeros".to_string(),
        "actual\nF(s) roots".to_string(),
        "reference\nreflection poles".to_string(),
        "actual\nE(s) roots".to_string(),
    ];
    let root_details = indent_block(
        &format!(
            "actual raw e_w roots        {}\nactual reflected e_w roots  {}\nactual a_s roots            {}",
            format_complex_list(actual_raw_w_roots),
            format_complex_list(actual_reflected_w_roots),
            format_complex_list(actual_a_roots),
        ),
        3,
    );

    format!(
        "{fixture_title}\n\
{fixture_summary}\n\n\
{synthesis_title}\n\
{synthesis_summary}\n\n\
{coefficient_title}\n\
{coefficient_table}\n\n\
{root_title}\n\
{root_table}\n\
{root_details}\n\
\n\
{comparison_title}\n\
{comparison_summary}\n\
{length_summary}\n\n\
{center_title}\n\
{center_summary}\n",
        fixture_title = format_numbered_heading(1, "Fixture"),
        fixture_summary = format_aligned_summary_with_width(
            &[
                ("case_id", case_id.to_string()),
                ("order", fixture.spec.order.to_string()),
                ("return_loss_db", fixture.spec.return_loss_db().to_string()),
                ("center_hz", fixture.mapping.center_hz().to_string()),
                ("bandwidth_hz", fixture.mapping.bandwidth_hz().to_string()),
                (
                    "normalized_transmission_zeros",
                    format!("{:?}", synthesis.polynomials.transmission_zeros_normalized),
                ),
            ],
            shared_summary_width,
        ),
        synthesis_title = format_numbered_heading(2, "Synthesis"),
        synthesis_summary = format_aligned_summary_with_width(
            &[
                ("approximation_kind", synthesis.approximation_kind().to_string()),
                ("matrix_method", format!("{:?}", synthesis.matrix_method)),
                ("matrix_order", synthesis.matrix.order().to_string()),
                ("response_samples", outcome.response.samples.len().to_string()),
                ("eps", format_decimal_scalar(synthesis.polynomials.eps)),
                ("eps_r", format_decimal_scalar(synthesis.polynomials.eps_r)),
            ],
            shared_summary_width,
        ),
        coefficient_title = indent_line("Polynomial Comparison", 3),
        coefficient_table = indent_block(&format_box_table(&coefficient_table), 3),
        root_title = indent_line("Root Comparison", 3),
        root_table = indent_block(&format_box_table(&root_table), 3),
        root_details = root_details,
        comparison_title = format_numbered_heading(3, "Comparison Summary"),
        comparison_summary = indent_block(
            &format!(
                "- {e_summary}\n- {e_from_poles_summary}\n- {f_summary}\n- {f_from_zeros_summary}\n- {p_summary}"
            ),
            3,
        ),
        length_summary = format_aligned_summary_with_width(
            &[
                (
                    "reference_e_length",
                    reference_model.polynomial_coefficients.e.len().to_string(),
                ),
                ("actual_e_length", actual_e.len().to_string()),
                (
                    "reference_f_length",
                    reference_model.polynomial_coefficients.f.len().to_string(),
                ),
                ("actual_f_length", actual_f.len().to_string()),
                (
                    "reference_p_length",
                    reference_model.polynomial_coefficients.p.len().to_string(),
                ),
                ("actual_p_length", actual_p.len().to_string()),
            ],
            shared_summary_width,
        ),
        center_title = format_numbered_heading(4, "Center Response Sample"),
        center_summary = format_aligned_summary_with_width(
            &[
                ("frequency_hz", format_decimal_scalar(center_sample.frequency_hz)),
                (
                    "normalized_omega",
                    format_decimal_scalar(center_sample.normalized_omega),
                ),
                (
                    "s11",
                    format_complex_scalar_parts(center_sample.s11_re, center_sample.s11_im),
                ),
                (
                    "s21",
                    format_complex_scalar_parts(center_sample.s21_re, center_sample.s21_im),
                ),
            ],
            shared_summary_width,
        ),
    )
}

fn format_complex_list(values: &[ComplexValue]) -> String {
    let entries = values
        .iter()
        .map(|value| format_complex_scalar_parts(value.re, value.im))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{entries}]")
}

pub fn render_markdown_synthesis_report(
    spec: &FilterSpec,
    synthesis: &SynthesisOutcome,
    transform: &TransformOutcome,
    response: &SParameterResponse,
    debug: bool,
) -> Result<String> {
    let polynomials = &synthesis.polynomials;
    let matrix = &synthesis.matrix;
    let transformed = &transform.matrix;
    let singularity_table = format_singularity_table_data(spec, polynomials.generalized.as_ref())?;
    let out_of_band_attenuation_table = spec
        .out_of_band_attenuation
        .as_ref()
        .map(format_out_of_band_attenuation_table_data);
    let mut spec_markdown_lines = vec![
        format!("- order: `{}`", spec.order),
        format!("- return_loss_db: `{}`", spec.return_loss_db()),
        format!(
            "- normalized_transmission_zeros: `{:?}`",
            polynomials.transmission_zeros_normalized
        ),
    ];
    if let Some(unloaded_q) = spec.unloaded_q {
        spec_markdown_lines.push(format!(
            "- unloaded_q: `{}`",
            format_decimal_scalar(unloaded_q)
        ));
    }
    if let Some(out_of_band_attenuation) = spec.out_of_band_attenuation.as_ref() {
        spec_markdown_lines.push(format!(
            "- out_of_band_windows: `{}`",
            out_of_band_attenuation.windows.len()
        ));
    }
    let spec_markdown_block = spec_markdown_lines.join("\n");
    let polynomial_block = if debug {
        let polynomial_table = format_polynomial_table_data(polynomials);
        format!(
            "### Transfer and Reflection Function Polynomials\n\n{}\n\n",
            format_markdown_table(&polynomial_table)
        )
    } else {
        String::new()
    };
    Ok(format!(
        "# Quickstart Report Output\n\n\
## Spec\n\n\
{}\n\n\
{}\
## Default Prototype Convention\n\n\
- `generalized_chebyshev(&spec)` runs normalized prototype synthesis without a frequency mapping\n\
- `FilterSpec` transmission zeros must already be normalized prototype values\n\
- physical Hz zeros should be converted ahead of time with `normalize_transmission_zeros_hz(...)`\n\
- coefficient lists below use ascending powers\n\
- index `0` is the constant term\n\n\
## Prototype\n\n\
- order: `{}`\n\
- ripple_factor: `{}`\n\
- eps: `{}`\n\
- eps_r: `{}`\n\n\
{}\\
### Corresponding Singularities\n\n\
{}\n\n\
## Canonical Synthesis\n\n\
- approximation_kind: `{}`\n\
- matrix_method: `{:?}`\n\
- topology: `{:?}`\n\
- shape: `{:?}`\n\n\
```text\n\
{}\n\
```\n\n\
## Transform\n\n\
- requested_topology: `{:?}`\n\
- source_topology: `{:?}`\n\
- result_topology: `{:?}`\n\
- pattern_verified: `{}`\n\
- response_check_passed: `{}`\n\
- overall_passed: `{}`\n\
- notes: `{:?}`\n\n\
```text\n\
{}\n\
```\n\n\
## Response Summary\n\n\
- sample_count: `{}`\n\n\
{}\n",
        spec_markdown_block,
        format_optional_out_of_band_markdown_block(out_of_band_attenuation_table.as_ref()),
        polynomials.order,
        polynomials.ripple_factor,
        polynomials.eps,
        polynomials.eps_r,
        polynomial_block,
        format_markdown_table(&singularity_table),
        synthesis.approximation_kind(),
        synthesis.matrix_method,
        matrix.topology(),
        matrix.shape(),
        format_box_table(&format_matrix_table_data(matrix)),
        transform.topology,
        transform.report.source_topology,
        transform.report.result_topology,
        transform.report.pattern_verified,
        transform.report.response.passes(),
        transform.report.passes(),
        transform.report.notes,
        format_box_table(&format_matrix_table_data(transformed)),
        response.samples.len(),
        format_markdown_table(&format_response_samples_table_data(&response.samples)),
    ))
}

fn format_notes_block(notes: &[String]) -> String {
    if notes.is_empty() {
        return String::new();
    }

    let mut block = String::from("   Notes\n");
    for note in notes {
        block.push_str("     - ");
        block.push_str(note);
        block.push('\n');
    }
    block.push('\n');
    block
}

fn format_numbered_heading(index: usize, title: &str) -> String {
    format!("{index}. {title}")
}

fn format_optional_out_of_band_terminal_block(
    out_of_band_attenuation_table: Option<&crate::output::Table>,
) -> String {
    match out_of_band_attenuation_table {
        Some(table) => format!(
            "{}\n{}\n\n",
            indent_line("Out-of-Band Attenuation", 3),
            indent_block(&format_box_table(table), 3)
        ),
        None => String::new(),
    }
}

fn format_optional_out_of_band_markdown_block(
    out_of_band_attenuation_table: Option<&crate::output::Table>,
) -> String {
    match out_of_band_attenuation_table {
        Some(table) => format!(
            "### Out-of-Band Attenuation\n\n{}\n\n",
            format_markdown_table(table)
        ),
        None => String::new(),
    }
}

fn indent_block(text: &str, spaces: usize) -> String {
    let prefix = " ".repeat(spaces);
    text.lines()
        .map(|line| format!("{prefix}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn indent_line(text: &str, spaces: usize) -> String {
    format!("{}{}", " ".repeat(spaces), text)
}
