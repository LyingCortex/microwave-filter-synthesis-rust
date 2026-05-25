//! MFS CLI — command-line interface for the microwave filter synthesis pipeline.
//!
//! Accepts a JSON synthesis request from a file or stdin, runs the full pipeline,
//! and outputs the result as JSON (default) or a human-readable table.
//!
//! # Usage
//!
//! ```text
//! mfs_cli [OPTIONS] [INPUT_FILE]
//!
//! Options:
//!   --input <FILE>       JSON input file (reads stdin if omitted)
//!   --format <FORMAT>    Output format: json (default) | table
//!   --stage <STAGE>      Execute only this stage: approximation | matrix_synthesis |
//!                        topology_transform | response_evaluation
//!   --resume <FILE>      Resume from a previously saved SynthesisContext JSON file
//! ```

use std::io::{self, Read};
use std::path::PathBuf;
use std::process;

use clap::Parser;
use serde::Serialize;

use mfs::pipeline::{
    load_context, run_from_json, run_full_pipeline, run_stage, PipelineOptions, SynthesisContext,
    SynthesisRequest,
};

/// MFS — Microwave Filter Synthesis CLI
///
/// Runs the full synthesis pipeline from a JSON request and outputs the result.
#[derive(Parser, Debug)]
#[command(name = "mfs_cli", version, about)]
struct Cli {
    /// JSON input file path. Reads from stdin if omitted.
    #[arg(short, long, value_name = "FILE")]
    input: Option<PathBuf>,

    /// Output format: json (default) or table.
    #[arg(short, long, value_name = "FORMAT", default_value = "json")]
    format: OutputFormat,

    /// Execute only this pipeline stage. Outputs only that stage's artifact.
    /// Valid values: approximation, matrix_synthesis, topology_transform, response_evaluation
    #[arg(short, long, value_name = "STAGE")]
    stage: Option<StageName>,

    /// Resume from a previously saved SynthesisContext JSON file.
    /// Loads the context and continues from the last completed stage.
    /// Can be combined with --stage to run just one more stage from the saved state.
    #[arg(short, long, value_name = "FILE")]
    resume: Option<PathBuf>,

    /// Positional input file (alternative to --input).
    #[arg(value_name = "INPUT_FILE")]
    input_file: Option<PathBuf>,
}

/// Valid pipeline stage names for the --stage flag.
#[derive(Debug, Clone, clap::ValueEnum)]
enum StageName {
    /// Polynomial approximation stage.
    Approximation,
    /// Coupling matrix synthesis stage.
    MatrixSynthesis,
    /// Topology transformation stage.
    TopologyTransform,
    /// S-parameter response evaluation stage.
    ResponseEvaluation,
}

impl StageName {
    /// Returns the stage name string used by the pipeline execution module.
    fn as_str(&self) -> &'static str {
        match self {
            Self::Approximation => "approximation",
            Self::MatrixSynthesis => "matrix_synthesis",
            Self::TopologyTransform => "topology_transform",
            Self::ResponseEvaluation => "response_evaluation",
        }
    }
}

/// Supported output formats.
#[derive(Debug, Clone, clap::ValueEnum)]
enum OutputFormat {
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

fn main() {
    let cli = Cli::parse();

    match run(&cli) {
        Ok(()) => {}
        Err(e) => {
            let cli_error = CliError {
                error_type: error_type_name(&e),
                message: e.to_string(),
                context: None,
            };
            let json = serde_json::to_string(&cli_error).unwrap_or_else(|_| {
                format!(r#"{{"error_type":"internal","message":"{}"}}"#, e)
            });
            eprintln!("{json}");
            process::exit(1);
        }
    }
}

fn run(cli: &Cli) -> Result<(), CliRunError> {
    // Determine execution mode based on flags
    match (&cli.resume, &cli.stage) {
        // --resume with optional --stage: load context and continue
        (Some(resume_path), stage_opt) => {
            let mut ctx = load_context(resume_path).map_err(CliRunError::Pipeline)?;

            match stage_opt {
                // --resume + --stage: run just the specified stage from saved state
                Some(stage_name) => {
                    run_stage(&mut ctx, stage_name.as_str()).map_err(CliRunError::Pipeline)?;
                    output_stage_artifact(&ctx, stage_name, &cli.format)?;
                }
                // --resume only: run all remaining stages
                None => {
                    run_remaining_stages(&mut ctx).map_err(CliRunError::Pipeline)?;
                    output_context(&ctx, &cli.format)?;
                }
            }
        }
        // --stage without --resume: create fresh context from input and run single stage
        (None, Some(stage_name)) => {
            let json_input = read_input(cli)?;
            let request: SynthesisRequest = serde_json::from_str(&json_input)
                .map_err(|e| CliRunError::InputParse(e.to_string()))?;
            let mut ctx = build_context_from_request(request).map_err(CliRunError::Pipeline)?;

            run_stage(&mut ctx, stage_name.as_str()).map_err(CliRunError::Pipeline)?;
            output_stage_artifact(&ctx, stage_name, &cli.format)?;
        }
        // No --stage, no --resume: full pipeline (original behavior)
        (None, None) => {
            let json_input = read_input(cli)?;

            match cli.format {
                OutputFormat::Json => {
                    let output = run_from_json(&json_input).map_err(CliRunError::Pipeline)?;
                    let value: serde_json::Value = serde_json::from_str(&output)
                        .map_err(|e| CliRunError::Serialization(e.to_string()))?;
                    let pretty = serde_json::to_string_pretty(&value)
                        .map_err(|e| CliRunError::Serialization(e.to_string()))?;
                    println!("{pretty}");
                }
                OutputFormat::Table => {
                    let request: SynthesisRequest = serde_json::from_str(&json_input)
                        .map_err(|e| CliRunError::InputParse(e.to_string()))?;
                    let ctx = run_full_pipeline(request).map_err(CliRunError::Pipeline)?;
                    print_table_output(&ctx);
                }
            }
        }
    }

    Ok(())
}

/// Builds a `SynthesisContext` from a `SynthesisRequest` without running any stages.
/// This prepares the context for incremental stage execution via `run_stage`.
fn build_context_from_request(request: SynthesisRequest) -> mfs::Result<SynthesisContext> {
    // Validate the request
    request.validate().map_err(|ve| {
        mfs::MfsError::PreconditionViolation(ve.to_string())
    })?;

    // Build FilterSpec from request
    let mut spec = mfs::FilterSpec::new(request.order, request.return_loss_db)?;
    if !request.transmission_zeros.is_empty() {
        spec = spec.with_normalized_transmission_zeros(request.transmission_zeros.iter().copied());
    }
    if let Some(q) = request.unloaded_q {
        spec = spec.with_unloaded_q(q);
    }

    // Build PipelineOptions from request
    let options = PipelineOptions {
        topology: request.topology,
        grid: request.grid,
        mapping: request.mapping,
        response_tolerance: None,
    };

    Ok(SynthesisContext::with_options(spec, options))
}

/// Determines which stages haven't been executed yet and runs them in order.
fn run_remaining_stages(ctx: &mut SynthesisContext) -> mfs::Result<()> {
    let all_stages = ["approximation", "matrix_synthesis", "topology_transform", "response_evaluation"];

    for stage_name in &all_stages {
        if ctx.metadata.stages_executed.contains(&stage_name.to_string()) {
            continue;
        }
        // Skip optional stages that lack required configuration
        match *stage_name {
            "topology_transform" => {
                if ctx.options.topology.is_none() {
                    continue;
                }
            }
            "response_evaluation" => {
                // Always attempt response evaluation if matrix is available
            }
            _ => {}
        }
        run_stage(ctx, stage_name)?;
    }

    Ok(())
}

/// Outputs only the artifact produced by the specified stage.
fn output_stage_artifact(
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
                CliRunError::Serialization("approximation stage has not produced output".to_string())
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
                    let s11_db = 10.0 * (s.s11_re.powi(2) + s.s11_im.powi(2)).log10();
                    let s21_db = 10.0 * (s.s21_re.powi(2) + s.s21_im.powi(2)).log10();
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
                println!("  {:>14} {:>10} {:>10}", "Freq (Hz)", "S11 (dB)", "S21 (dB)");
                println!("  {:>14} {:>10} {:>10}", "──────────", "────────", "────────");
                for sample in response.samples.iter().take(show) {
                    let s11_db = 10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
                    let s21_db = 10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
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
                        let s11_db =
                            10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
                        let s21_db =
                            10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
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
fn output_context(ctx: &SynthesisContext, format: &OutputFormat) -> Result<(), CliRunError> {
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

fn read_input(cli: &Cli) -> Result<String, CliRunError> {
    // --input flag takes priority over positional argument
    let path = cli.input.as_ref().or(cli.input_file.as_ref());

    match path {
        Some(file_path) => {
            std::fs::read_to_string(file_path).map_err(|e| {
                CliRunError::InputRead(format!(
                    "failed to read '{}': {e}",
                    file_path.display()
                ))
            })
        }
        None => {
            // Read from stdin
            let mut buffer = String::new();
            io::stdin().read_to_string(&mut buffer).map_err(|e| {
                CliRunError::InputRead(format!("failed to read from stdin: {e}"))
            })?;
            Ok(buffer)
        }
    }
}

fn print_table_output(ctx: &mfs::pipeline::SynthesisContext) {
    println!("╔══════════════════════════════════════╗");
    println!("║   MFS Synthesis Pipeline Results     ║");
    println!("╚══════════════════════════════════════╝");
    println!();

    // Spec summary
    println!("── Filter Specification ──");
    println!("  Order:          {}", ctx.spec.order);
    println!("  Return Loss:    {:.2} dB", ctx.spec.return_loss_db);
    if !ctx.spec.transmission_zeros.is_empty() {
        let zeros: Vec<String> = ctx.spec.transmission_zeros.iter().map(|tz| {
            if tz.value.is_infinite() {
                "∞".to_string()
            } else {
                format!("{:.4}", tz.value)
            }
        }).collect();
        println!("  Zeros:          [{}]", zeros.join(", "));
    }
    if let Some(q) = ctx.spec.unloaded_q {
        println!("  Unloaded Q:     {:.0}", q);
    }
    println!();

    // Metadata
    println!("── Execution Metadata ──");
    println!("  Version:        {}", ctx.metadata.version);
    println!("  Stages:         {}", ctx.metadata.stages_executed.join(" → "));
    if !ctx.metadata.stage_timings_ms.is_empty() {
        let timings: Vec<String> = ctx.metadata.stages_executed.iter()
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
            let tz: Vec<String> = polys.transmission_zeros_normalized.iter()
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
        println!("  {:>14} {:>10} {:>10}", "Freq (Hz)", "S11 (dB)", "S21 (dB)");
        println!("  {:>14} {:>10} {:>10}", "──────────", "────────", "────────");
        for sample in response.samples.iter().take(show) {
            let s11_db = 10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
            let s21_db = 10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
            println!("  {:>14.4} {:>10.3} {:>10.3}", sample.frequency_hz, s11_db, s21_db);
        }
        if n > show * 2 {
            println!("  {:>14}", "...");
        }
        if n > show {
            for sample in response.samples.iter().skip(n - show) {
                let s11_db = 10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
                let s21_db = 10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
                println!("  {:>14.4} {:>10.3} {:>10.3}", sample.frequency_hz, s11_db, s21_db);
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
enum CliRunError {
    /// Failed to read input file or stdin.
    InputRead(String),
    /// Failed to parse input JSON.
    InputParse(String),
    /// Pipeline execution error.
    Pipeline(mfs::MfsError),
    /// Serialization error.
    Serialization(String),
}

impl std::fmt::Display for CliRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputRead(msg) => write!(f, "{msg}"),
            Self::InputParse(msg) => write!(f, "invalid input JSON: {msg}"),
            Self::Pipeline(err) => write!(f, "{err}"),
            Self::Serialization(msg) => write!(f, "serialization error: {msg}"),
        }
    }
}

fn error_type_name(err: &CliRunError) -> String {
    match err {
        CliRunError::InputRead(_) => "io_error".to_string(),
        CliRunError::InputParse(_) => "parse_error".to_string(),
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
    }
}
