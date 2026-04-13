use crate::approx::{GeneralizedChebyshevData, PolynomialSet};
use crate::fixtures::ComplexValue;
use crate::freq::validated_transmission_zero;
use crate::matrix::CouplingMatrix;
use crate::response::ResponseSample;
use crate::spec::{FilterSpec, OutOfBandAttenuationSpec};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

pub fn format_complex_scalar_parts(re: f64, im: f64) -> String {
    let tolerance = 1e-12;
    let real = round_to_4(if re.abs() < tolerance { 0.0 } else { re });
    let imag = round_to_4(if im.abs() < tolerance { 0.0 } else { im });

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

pub fn format_decimal_scalar(value: f64) -> String {
    format_unsigned_number(value)
}

pub fn format_key_value_table_data(rows: &[(&str, String)]) -> Table {
    Table {
        headers: vec!["field".to_string(), "value".to_string()],
        rows: rows
            .iter()
            .map(|(key, value)| vec![(*key).to_string(), value.clone()])
            .collect(),
    }
}

pub fn format_aligned_summary(rows: &[(&str, String)]) -> String {
    if rows.is_empty() {
        return String::new();
    }

    let width = rows.iter().map(|(key, _)| key.len()).max().unwrap_or(0);
    format_aligned_summary_with_width(rows, width)
}

pub fn format_aligned_summary_with_width(rows: &[(&str, String)], width: usize) -> String {
    if rows.is_empty() {
        return String::new();
    }

    rows.iter()
        .map(|(key, value)| format!("   {key:<width$}  {value}", width = width))
        .collect::<Vec<_>>()
        .join("\n")
}

pub fn format_matrix_table_data(matrix: &CouplingMatrix) -> Table {
    let side = matrix.side();
    let headers = std::iter::once(String::new())
        .chain(matrix_headers(side))
        .collect::<Vec<_>>();
    let rows = (0..side)
        .map(|row| {
            let mut values = Vec::with_capacity(side + 1);
            values.push(header_label(row, side));
            for col in 0..side {
                values.push(format_signed_number(matrix.at(row, col).unwrap_or_default()));
            }
            values
        })
        .collect::<Vec<_>>();

    Table { headers, rows }
}

pub fn format_response_samples_table_data(samples: &[ResponseSample]) -> Table {
    let indices = sample_indices(samples.len());
    let rows = indices
        .into_iter()
        .map(|index| {
            let sample = &samples[index];
            vec![
                index.to_string(),
                format_unsigned_number(sample.frequency_hz),
                format_signed_number(sample.normalized_omega),
                format_complex_scalar_parts(sample.s11_re, sample.s11_im),
                format_complex_scalar_parts(sample.s21_re, sample.s21_im),
                format_signed_number(sample.group_delay),
            ]
        })
        .collect::<Vec<_>>();

    Table {
        headers: vec![
            "sample".to_string(),
            "frequency_hz".to_string(),
            "normalized_omega".to_string(),
            "s11".to_string(),
            "s21".to_string(),
            "group_delay".to_string(),
        ],
        rows,
    }
}

pub fn format_polynomial_table_data(polynomials: &PolynomialSet) -> Table {
    let e = polynomials
        .e
        .coefficients
        .iter()
        .map(|value| format_complex_scalar_parts(value.re, value.im))
        .collect::<Vec<_>>();
    let f = polynomials
        .f
        .coefficients
        .iter()
        .map(|value| format_complex_scalar_parts(value.re, value.im))
        .collect::<Vec<_>>();
    let p = polynomials
        .p
        .coefficients
        .iter()
        .map(|value| format_complex_scalar_parts(value.re, value.im))
        .collect::<Vec<_>>();
    let row_count = e.len().max(f.len()).max(p.len());

    Table {
        headers: vec![
            "i".to_string(),
            "E(s)".to_string(),
            "F(s)".to_string(),
            "P(s)".to_string(),
        ],
        rows: (0..row_count)
            .map(|index| {
                vec![
                    index.to_string(),
                    e.get(index).cloned().unwrap_or_default(),
                    f.get(index).cloned().unwrap_or_default(),
                    p.get(index).cloned().unwrap_or_default(),
                ]
            })
            .collect(),
    }
}

pub fn format_singularity_table_data(
    spec: &FilterSpec,
    generalized: Option<&GeneralizedChebyshevData>,
) -> crate::error::Result<Table> {
    let reflection_zeros = generalized
        .map(|data| data.f_s.roots())
        .transpose()?
        .unwrap_or_default()
        .into_iter()
        .map(|root| (root.re, root.im))
        .collect::<Vec<_>>();
    let reflection_zeros = sort_complex_pairs(&reflection_zeros)
        .into_iter()
        .map(|(re, im)| format_complex_scalar_parts(re, im))
        .collect::<Vec<_>>();

    let prescribed_zeros = generalized
        .map(|data| data.padded_zeros.clone())
        .unwrap_or_else(|| {
            spec.transmission_zeros
                .iter()
                .filter_map(|zero| validated_transmission_zero(*zero).ok())
                .collect::<Vec<_>>()
        })
        .into_iter()
        .map(|value: f64| {
            if value.is_infinite() {
                format!("j{}", '\u{221E}')
            } else {
                format_complex_scalar_parts(0.0, value)
            }
        })
        .collect::<Vec<_>>();

    let poles = generalized
        .and_then(|data| data.e_stage.as_ref())
        .map(|stage| {
            sort_complex_pairs(
                &stage
                    .e_s_roots
                    .iter()
                    .map(|root| (root.re, root.im))
                    .collect::<Vec<_>>(),
            )
            .into_iter()
            .map(|(re, im)| format_complex_scalar_parts(re, im))
            .collect::<Vec<_>>()
        })
        .unwrap_or_default();

    let row_count = reflection_zeros
        .len()
        .max(prescribed_zeros.len())
        .max(poles.len());

    Ok(Table {
        headers: vec![
            "i".to_string(),
            "Reflection Zeros (Roots of F(s))".to_string(),
            "Transmission Zeros (Prescribed)".to_string(),
            "Transmission/Reflection Poles (Roots of E(s))".to_string(),
        ],
        rows: (0..row_count)
            .map(|index| {
                vec![
                    (index + 1).to_string(),
                    reflection_zeros.get(index).cloned().unwrap_or_default(),
                    prescribed_zeros.get(index).cloned().unwrap_or_default(),
                    poles.get(index).cloned().unwrap_or_default(),
                ]
            })
            .collect(),
    })
}

pub fn format_out_of_band_attenuation_table_data(
    out_of_band_attenuation: &OutOfBandAttenuationSpec,
) -> Table {
    Table {
        headers: vec![
            "window".to_string(),
            "start_freq_hz".to_string(),
            "stop_freq_hz".to_string(),
            "attenuation_db".to_string(),
        ],
        rows: out_of_band_attenuation
            .windows
            .iter()
            .enumerate()
            .map(|(index, window)| {
                vec![
                    (index + 1).to_string(),
                    format_unsigned_number(window.start_freq_hz),
                    format_unsigned_number(window.stop_freq_hz),
                    format_unsigned_number(window.attenuation_db),
                ]
            })
            .collect(),
    }
}

pub fn format_reference_actual_polynomial_table_data(
    reference_e: &[ComplexValue],
    actual_e: &[ComplexValue],
    reference_f: &[ComplexValue],
    actual_f: &[ComplexValue],
    reference_p: &[ComplexValue],
    actual_p: &[ComplexValue],
) -> Table {
    let row_count = reference_e
        .len()
        .max(actual_e.len())
        .max(reference_f.len())
        .max(actual_f.len())
        .max(reference_p.len())
        .max(actual_p.len());

    Table {
        headers: vec![
            "i".to_string(),
            "reference E(s)".to_string(),
            "actual E(s)".to_string(),
            "reference F(s)".to_string(),
            "actual F(s)".to_string(),
            "reference P(s)".to_string(),
            "actual P(s)".to_string(),
        ],
        rows: (0..row_count)
            .map(|index| {
                vec![
                    index.to_string(),
                    reference_e
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    actual_e
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    reference_f
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    actual_f
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    reference_p
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    actual_p
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                ]
            })
            .collect(),
    }
}

pub fn format_root_comparison_table_data(
    reference_reflection_zeros: &[ComplexValue],
    actual_f_roots: &[ComplexValue],
    reference_reflection_poles: &[ComplexValue],
    actual_e_roots: &[ComplexValue],
) -> Table {
    let row_count = reference_reflection_zeros
        .len()
        .max(actual_f_roots.len())
        .max(reference_reflection_poles.len())
        .max(actual_e_roots.len());

    Table {
        headers: vec![
            "i".to_string(),
            "reference reflection zeros".to_string(),
            "actual F(s) roots".to_string(),
            "reference reflection poles".to_string(),
            "actual E(s) roots".to_string(),
        ],
        rows: (0..row_count)
            .map(|index| {
                vec![
                    (index + 1).to_string(),
                    reference_reflection_zeros
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    actual_f_roots
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    reference_reflection_poles
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                    actual_e_roots
                        .get(index)
                        .map(|value| format_complex_scalar_parts(value.re, value.im))
                        .unwrap_or_default(),
                ]
            })
            .collect(),
    }
}

fn header_label(index: usize, side: usize) -> String {
    if index == 0 {
        "S".to_string()
    } else if index + 1 == side {
        "L".to_string()
    } else {
        index.to_string()
    }
}

fn matrix_headers(side: usize) -> impl Iterator<Item = String> {
    (0..side).map(move |index| header_label(index, side))
}

fn sample_indices(len: usize) -> Vec<usize> {
    if len == 0 {
        return Vec::new();
    }

    let mut indices = vec![0, len / 2, len - 1];
    indices.sort_unstable();
    indices.dedup();
    indices
}

fn sort_complex_pairs(values: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut sorted = values.to_vec();
    sorted.sort_by(|lhs, rhs| {
        lhs.0
            .partial_cmp(&rhs.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                lhs.1
                    .partial_cmp(&rhs.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    sorted
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
