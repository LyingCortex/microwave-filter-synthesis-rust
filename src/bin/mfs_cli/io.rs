//! §D.5: CLI I/O — the output-format enum, the CLI error type and the
//! JSON error-fallback, input reading, and all result output (full-context
//! table, per-stage artifact JSON/table, JSON pretty printing).
//!
//! Pure move out of the monolithic `src/bin/mfs_cli.rs` — no behavior
//! change. §D.7: the hand-rolled `10*log10(re^2+im^2)` dB sites now route
//! through `ResponseSample::s11_db()` / `s21_db()` (−300 dB floor).
//!
use std::io::{self, Read};

use serde::Serialize;

use mfs::pipeline::SynthesisContext;

use crate::{Cli, StageName};
/// Supported output formats.
#[derive(Debug, Clone, clap::ValueEnum)]
pub(crate) enum OutputFormat {
    /// JSON output (default).
    Json,
    /// Human-readable table output.
    Table,
}

/// Structured error object written to stderr on failure.
#[derive(Serialize)]
struct CliError {
    error_type: String,
    message: String,
    context: Option<String>,
}

/// Outputs only the artifact produced by the specified stage.
pub(crate) fn output_stage_artifact(
    ctx: &SynthesisContext,
    stage: &StageName,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    match format {
        OutputFormat::Json => {
            let json = serialize_stage_artifact(ctx, stage)?;
            println!("{json}");
        }
        OutputFormat::Table => {
            print_stage_table(ctx, stage);
        }
    }
    Ok(())
}

/// Serializes only the artifact for a specific stage to pretty JSON.
fn serialize_stage_artifact(
    ctx: &SynthesisContext,
    stage: &StageName,
) -> Result<String, CliRunError> {
    let value: serde_json::Value = match stage {
        StageName::Approximation => {
            let polys = ctx.polynomials().ok_or_else(|| {
                CliRunError::Serialization(
                    "approximation stage has not produced output".to_string(),
                )
            })?;
            serde_json::json!({
                "stage": "approximation",
                "order": polys.order,
                "epsilon": polys.eps,
                "epsilon_r": polys.eps_r,
                "transmission_zeros_normalized": polys.transmission_zeros_normalized,
            })
        }
        StageName::MatrixSynthesis => {
            let matrix = ctx.matrix().ok_or_else(|| {
                CliRunError::Serialization(
                    "matrix_synthesis stage has not produced output".to_string(),
                )
            })?;
            serde_json::json!({
                "stage": "matrix_synthesis",
                "order": matrix.order(),
                "side": matrix.side(),
                "topology": format!("{:?}", matrix.topology()),
                "data": matrix.as_slice(),
            })
        }
        StageName::TopologyTransform => {
            let transform = ctx.transform().ok_or_else(|| {
                CliRunError::Serialization(
                    "topology_transform stage has not produced output".to_string(),
                )
            })?;
            serde_json::json!({
                "stage": "topology_transform",
                "topology": format!("{:?}", transform.topology),
                "pattern_verified": transform.report.pattern_verified,
                "notes": transform.report.notes,
            })
        }
        StageName::ResponseEvaluation => {
            let response = ctx.response().ok_or_else(|| {
                CliRunError::Serialization(
                    "response_evaluation stage has not produced output".to_string(),
                )
            })?;
            let samples: Vec<serde_json::Value> = response
                .samples
                .iter()
                .map(|s| {
                    // §D.7: route through ResponseSample::s11_db()/s21_db()
                    // (−300 dB floor) instead of the hand-rolled 10*log10(re^2+im^2).
                    let s11_db = s.s11_db();
                    let s21_db = s.s21_db();
                    serde_json::json!({
                        "frequency_hz": s.frequency_hz,
                        "s11_db": s11_db,
                        "s21_db": s21_db,
                    })
                })
                .collect();
            serde_json::json!({
                "stage": "response_evaluation",
                "sample_count": response.samples.len(),
                "samples": samples,
            })
        }
    };

    serde_json::to_string_pretty(&value).map_err(|e| CliRunError::Serialization(e.to_string()))
}

/// Prints a table-formatted view of a single stage's artifact.
fn print_stage_table(ctx: &SynthesisContext, stage: &StageName) {
    match stage {
        StageName::Approximation => {
            if let Some(polys) = ctx.polynomials() {
                println!("── Approximation Stage Output ──");
                println!("  Order:    {}", polys.order);
                println!("  Epsilon:  {:.6}", polys.eps);
                println!("  Epsilon_R: {:.6}", polys.eps_r);
                if !polys.transmission_zeros_normalized.is_empty() {
                    let tz: Vec<String> = polys
                        .transmission_zeros_normalized
                        .iter()
                        .map(|z| {
                            if z.is_infinite() {
                                "∞".to_string()
                            } else {
                                format!("{:.4}", z)
                            }
                        })
                        .collect();
                    println!("  Norm. Zeros: [{}]", tz.join(", "));
                }
            } else {
                println!("(approximation stage has not produced output)");
            }
        }
        StageName::MatrixSynthesis => {
            if let Some(matrix) = ctx.matrix() {
                println!("── Matrix Synthesis Stage Output ──");
                println!("  Order:    {}", matrix.order());
                println!("  Size:     {}×{}", matrix.side(), matrix.side());
                println!("  Topology: {:?}", matrix.topology());
                let side = matrix.side();
                if side <= 8 {
                    println!("  Data:");
                    for i in 0..side {
                        let row: Vec<String> = (0..side)
                            .map(|j| format!("{:>8.4}", matrix.at(i, j).unwrap_or(0.0)))
                            .collect();
                        println!("    [{}]", row.join(" "));
                    }
                }
            } else {
                println!("(matrix_synthesis stage has not produced output)");
            }
        }
        StageName::TopologyTransform => {
            if let Some(transform) = ctx.transform() {
                println!("── Topology Transform Stage Output ──");
                println!("  Topology:         {:?}", transform.topology);
                println!("  Pattern Verified: {}", transform.report.pattern_verified);
                if !transform.report.notes.is_empty() {
                    for note in &transform.report.notes {
                        println!("  Note: {note}");
                    }
                }
            } else {
                println!("(topology_transform stage has not produced output)");
            }
        }
        StageName::ResponseEvaluation => {
            if let Some(response) = ctx.response() {
                println!("── Response Evaluation Stage Output ──");
                println!("  Samples: {}", response.samples.len());
                let n = response.samples.len();
                let show = 5.min(n);
                println!(
                    "  {:>14} {:>10} {:>10}",
                    "Freq (Hz)", "S11 (dB)", "S21 (dB)"
                );
                println!(
                    "  {:>14} {:>10} {:>10}",
                    "──────────", "────────", "────────"
                );
                for sample in response.samples.iter().take(show) {
                    // §D.7: route through ResponseSample::s11_db()/s21_db()
                    // (−300 dB floor) instead of the hand-rolled 10*log10(re^2+im^2).
                    let s11_db = sample.s11_db();
                    let s21_db = sample.s21_db();
                    println!(
                        "  {:>14.4} {:>10.3} {:>10.3}",
                        sample.frequency_hz, s11_db, s21_db
                    );
                }
                if n > show * 2 {
                    println!("  {:>14}", "...");
                }
                if n > show {
                    for sample in response.samples.iter().skip(n - show) {
                        // §D.7: route through ResponseSample::s11_db()/s21_db()
                        // (−300 dB floor) instead of the hand-rolled 10*log10(re^2+im^2).
                        let s11_db = sample.s11_db();
                        let s21_db = sample.s21_db();
                        println!(
                            "  {:>14.4} {:>10.3} {:>10.3}",
                            sample.frequency_hz, s11_db, s21_db
                        );
                    }
                }
            } else {
                println!("(response_evaluation stage has not produced output)");
            }
        }
    }
}

/// Outputs the full context (used after --resume without --stage).
pub(crate) fn output_context(ctx: &SynthesisContext, format: &OutputFormat) -> Result<(), CliRunError> {
    match format {
        OutputFormat::Json => {
            let json = serde_json::to_string_pretty(ctx)
                .map_err(|e| CliRunError::Serialization(e.to_string()))?;
            println!("{json}");
        }
        OutputFormat::Table => {
            print_table_output(ctx);
        }
    }
    Ok(())
}

pub(crate) fn read_input(cli: &Cli) -> Result<String, CliRunError> {
    // --input flag takes priority over positional argument
    let path = cli.input.as_ref().or(cli.input_file.as_ref());

    match path {
        Some(file_path) => std::fs::read_to_string(file_path).map_err(|e| {
            CliRunError::InputRead(format!("failed to read '{}': {e}", file_path.display()))
        }),
        None => {
            // Read from stdin
            let mut buffer = String::new();
            io::stdin()
                .read_to_string(&mut buffer)
                .map_err(|e| CliRunError::InputRead(format!("failed to read from stdin: {e}")))?;
            Ok(buffer)
        }
    }
}

pub(crate) fn print_table_output(ctx: &mfs::pipeline::SynthesisContext) {
    println!("╔══════════════════════════════════════╗");
    println!("║   MFS Synthesis Pipeline Results     ║");
    println!("╚══════════════════════════════════════╝");
    println!();

    // Spec summary
    println!("── Filter Specification ──");
    println!("  Order:          {}", ctx.spec.order);
    println!("  Return Loss:    {:.2} dB", ctx.spec.return_loss_db);
    if !ctx.spec.transmission_zeros.is_empty() {
        let zeros: Vec<String> = ctx
            .spec
            .transmission_zeros
            .iter()
            .map(|tz| {
                if tz.value.is_infinite() {
                    "∞".to_string()
                } else {
                    format!("{:.4}", tz.value)
                }
            })
            .collect();
        println!("  Zeros:          [{}]", zeros.join(", "));
    }
    if let Some(q) = ctx.spec.unloaded_q {
        println!("  Unloaded Q:     {:.0}", q);
    }
    println!();

    // Metadata
    println!("── Execution Metadata ──");
    println!("  Version:        {}", ctx.metadata.version);
    println!(
        "  Stages:         {}",
        ctx.metadata.stages_executed.join(" → ")
    );
    if !ctx.metadata.stage_timings_ms.is_empty() {
        let timings: Vec<String> = ctx
            .metadata
            .stages_executed
            .iter()
            .zip(ctx.metadata.stage_timings_ms.iter())
            .map(|(name, ms)| format!("{name}: {ms:.2}ms"))
            .collect();
        println!("  Timings:        {}", timings.join(", "));
    }
    if !ctx.metadata.warnings.is_empty() {
        println!("  Warnings:       {}", ctx.metadata.warnings.join("; "));
    }
    println!();

    // Polynomials
    if let Some(polys) = ctx.polynomials() {
        println!("── Approximation ──");
        println!("  Polynomial Order: {}", polys.order);
        println!("  Epsilon:          {:.6}", polys.eps);
        println!("  Epsilon_R:        {:.6}", polys.eps_r);
        if !polys.transmission_zeros_normalized.is_empty() {
            let tz: Vec<String> = polys
                .transmission_zeros_normalized
                .iter()
                .map(|z| {
                    if z.is_infinite() {
                        "∞".to_string()
                    } else {
                        format!("{:.4}", z)
                    }
                })
                .collect();
            println!("  Norm. Zeros:      [{}]", tz.join(", "));
        }
        println!();
    }

    // Matrix
    if let Some(matrix) = ctx.matrix() {
        println!("── Coupling Matrix ──");
        println!("  Order:    {}", matrix.order());
        println!("  Size:     {}×{}", matrix.side(), matrix.side());
        println!("  Topology: {:?}", matrix.topology());
        // Print matrix entries in a compact format
        let side = matrix.side();
        if side <= 8 {
            println!("  Data:");
            for i in 0..side {
                let row: Vec<String> = (0..side)
                    .map(|j| format!("{:>8.4}", matrix.at(i, j).unwrap_or(0.0)))
                    .collect();
                println!("    [{}]", row.join(" "));
            }
        } else {
            println!("  (matrix too large to display, {} entries)", side * side);
        }
        println!();
    }

    // Transform
    if let Some(transform) = ctx.transform() {
        println!("── Topology Transform ──");
        println!("  Topology:         {:?}", transform.topology);
        println!("  Pattern Verified: {}", transform.report.pattern_verified);
        if !transform.report.notes.is_empty() {
            for note in &transform.report.notes {
                println!("  Note: {note}");
            }
        }
        println!();
    }

    // Response
    if let Some(response) = ctx.response() {
        println!("── S-Parameter Response ──");
        println!("  Samples: {}", response.samples.len());
        // Show first few and last few samples
        let n = response.samples.len();
        let show = 5.min(n);
        println!(
            "  {:>14} {:>10} {:>10}",
            "Freq (Hz)", "S11 (dB)", "S21 (dB)"
        );
        println!(
            "  {:>14} {:>10} {:>10}",
            "──────────", "────────", "────────"
        );
        for sample in response.samples.iter().take(show) {
            // §D.7: route through ResponseSample::s11_db()/s21_db()
            // (−300 dB floor) instead of the hand-rolled 10*log10(re^2+im^2).
            let s11_db = sample.s11_db();
            let s21_db = sample.s21_db();
            println!(
                "  {:>14.4} {:>10.3} {:>10.3}",
                sample.frequency_hz, s11_db, s21_db
            );
        }
        if n > show * 2 {
            println!("  {:>14}", "...");
        }
        if n > show {
            for sample in response.samples.iter().skip(n - show) {
                // §D.7: route through ResponseSample::s11_db()/s21_db()
                // (−300 dB floor) instead of the hand-rolled 10*log10(re^2+im^2).
                let s11_db = sample.s11_db();
                let s21_db = sample.s21_db();
                println!(
                    "  {:>14.4} {:>10.3} {:>10.3}",
                    sample.frequency_hz, s11_db, s21_db
                );
            }
        }
        println!();
    }
}

// ---------------------------------------------------------------------------
// Error handling
// ---------------------------------------------------------------------------

/// Internal error type for CLI operations.
#[derive(Debug)]
pub(crate) enum CliRunError {
    /// Failed to read input file or stdin.
    InputRead(String),
    /// Failed to parse input JSON.
    InputParse(String),
    /// A command-line argument is invalid (e.g. an unknown enum value).
    ArgValidation(String),
    /// Pipeline execution error.
    Pipeline(mfs::MfsError),
    /// Serialization error.
    Serialization(String),
    /// Touchstone file read error.
    Touchstone(String),
}

impl From<serde_json::Error> for CliRunError {
    fn from(e: serde_json::Error) -> Self {
        Self::Serialization(e.to_string())
    }
}

impl From<mfs::MfsError> for CliRunError {
    fn from(e: mfs::MfsError) -> Self {
        Self::Pipeline(e)
    }
}

impl std::fmt::Display for CliRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputRead(msg) => write!(f, "{msg}"),
            Self::InputParse(msg) => write!(f, "invalid input JSON: {msg}"),
            Self::ArgValidation(msg) => write!(f, "invalid argument: {msg}"),
            Self::Pipeline(err) => write!(f, "{err}"),
            Self::Serialization(msg) => write!(f, "serialization error: {msg}"),
            Self::Touchstone(msg) => write!(f, "touchstone error: {msg}"),
        }
    }
}

pub(crate) fn error_type_name(err: &CliRunError) -> String {
    match err {
        CliRunError::InputRead(_) => "io_error".to_string(),
        CliRunError::InputParse(_) => "parse_error".to_string(),
        CliRunError::ArgValidation(_) => "invalid_argument".to_string(),
        CliRunError::Pipeline(mfs_err) => match mfs_err {
            mfs::MfsError::InvalidOrder { .. } => "invalid_order".to_string(),
            mfs::MfsError::InvalidReturnLoss { .. } => "invalid_return_loss".to_string(),
            mfs::MfsError::InvalidFrequency(_) => "invalid_frequency".to_string(),
            mfs::MfsError::InvalidGridSize { .. } => "invalid_grid_size".to_string(),
            mfs::MfsError::InvalidTransmissionZero(_) => "invalid_transmission_zero".to_string(),
            mfs::MfsError::DimensionMismatch { .. } => "dimension_mismatch".to_string(),
            mfs::MfsError::NumericalFailure(_) => "numerical_failure".to_string(),
            mfs::MfsError::NotImplemented(_) => "not_implemented".to_string(),
            mfs::MfsError::PreconditionViolation(_) => "precondition_violation".to_string(),
        },
        CliRunError::Serialization(_) => "serialization_error".to_string(),
        CliRunError::Touchstone(_) => "touchstone_error".to_string(),
    }
}
/// §D.5: the structured error JSON written to stderr on failure (moved here
/// from `main`). Falls back to a `serde_json::json!` object (so the message
/// is escaped correctly) and finally to a hardcoded literal, if struct
/// serialization itself fails.
pub(crate) fn format_cli_error(error_type: String, message: String) -> String {
    let cli_error = CliError {
        error_type,
        message: message.clone(),
        context: None,
    };
    serde_json::to_string(&cli_error).unwrap_or_else(|_| {
        // Last-resort fallback: build the JSON object via the
        // serializer so the message is escaped correctly (the
        // previous hand-formatted string could emit invalid JSON
        // when the message contained quotes or backslashes).
        let value = serde_json::json!({
            "error_type": "internal",
            "message": message,
            "context": null,
        });
        serde_json::to_string(&value)
            .unwrap_or_else(|_| r#"{"error_type":"internal","message":"serialization failed"}"#.into())
    })
}

#[cfg(test)]
mod io_tests {
    //! §D.7 regression: stage-artifact JSON must use the −300 dB floor
    //! (`ResponseSample::s11_db()`/`s21_db()`) at zero magnitude — the old
    //! hand-rolled `10.0 * (re.powi(2) + im.powi(2)).log10()` produced
    //! `log10(0) = −inf`, which `serde_json` silently degrades to `null`.
    use super::*;
    use crate::StageName;
    use mfs::pipeline::SynthesisContext;
    use mfs::response::{ResponseSample, SParameterResponse};
    use mfs::FilterSpec;

    #[test]
    fn stage_artifact_db_at_transmission_zero_uses_300db_floor() {
        let mut ctx = SynthesisContext::new(FilterSpec::new(2, 20.0).expect("spec"));
        // A sample exactly at a transmission zero: |S11| = |S21| = 0.
        ctx.response = Some(SParameterResponse {
            samples: vec![ResponseSample {
                frequency_hz: 1.0,
                normalized_omega: 1.0,
                group_delay: None,
                s11_re: 0.0,
                s11_im: 0.0,
                s21_re: 0.0,
                s21_im: 0.0,
            }],
        });
        let json =
            serialize_stage_artifact(&ctx, &StageName::ResponseEvaluation).expect("serialize");
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        let s11_db = value["samples"][0]["s11_db"].as_f64();
        let s21_db = value["samples"][0]["s21_db"].as_f64();
        assert_eq!(s11_db, Some(-300.0), "S11 at a null must be the -300 dB floor, got {s11_db:?}");
        assert_eq!(s21_db, Some(-300.0), "S21 at a null must be the -300 dB floor, got {s21_db:?}");
    }
}
