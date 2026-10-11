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

use std::path::PathBuf;
use std::process;

use clap::Parser;

use mfs::pipeline::{
    context_from_request, load_context, run_from_json, run_full_pipeline, run_remaining_stages,
    run_stage, SynthesisRequest,
};

// §D.5: one module per subcommand group — a pure move out of the monolithic
// `src/bin/mfs_cli.rs` (Cargo.toml `[[bin]] path` now points at this file).
mod design;
mod io;
mod multiplexer;
mod optimize;
mod target;

use crate::design::{run_design, DesignArgs};
use crate::multiplexer::run_multiplexer;
use crate::io::{
    CliRunError, OutputFormat, error_type_name, format_cli_error, output_context,
    output_stage_artifact, print_table_output, read_input,
};
use crate::optimize::{run_optimize, run_extract};

/// MFS — Microwave Filter Synthesis CLI
///
/// Runs the full synthesis pipeline from a JSON request and outputs the result.
/// Also supports a quick `design` subcommand for direct command-line usage.
#[derive(Parser, Debug)]
#[command(name = "mfs_cli", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// JSON input file path. Reads from stdin if omitted.
    #[arg(short, long, value_name = "FILE", global = true)]
    input: Option<PathBuf>,

    /// Output format: json (default) or table.
    #[arg(
        short,
        long,
        value_name = "FORMAT",
        default_value = "json",
        global = true
    )]
    format: OutputFormat,

    /// Execute only this pipeline stage.
    #[arg(short, long, value_name = "STAGE", global = true)]
    stage: Option<StageName>,

    /// Resume from a previously saved SynthesisContext JSON file.
    #[arg(long, value_name = "FILE", global = true)]
    resume: Option<PathBuf>,

    /// Positional input file (alternative to --input).
    #[arg(value_name = "INPUT_FILE")]
    input_file: Option<PathBuf>,
}

/// Subcommands for quick usage without JSON.
#[derive(Debug, clap::Subcommand)]
enum Command {
    /// Quick filter design from command-line parameters.
    Design {
        /// Filter order (number of resonators).
        #[arg(short = 'n', long)]
        order: usize,

        /// Return loss in dB.
        #[arg(short = 'r', long, alias = "rl")]
        return_loss: f64,

        /// Center frequency in Hz (enables bandpass mode).
        #[arg(short, long)]
        center: Option<f64>,

        /// Bandwidth in Hz.
        #[arg(short, long)]
        bandwidth: Option<f64>,

        /// Lower passband edge in Hz (alternative to center+bandwidth).
        #[arg(long, alias = "f1")]
        start: Option<f64>,

        /// Upper passband edge in Hz (alternative to center+bandwidth).
        #[arg(long, alias = "f2")]
        stop: Option<f64>,

        /// Transmission zeros (comma-separated, in Hz if bandpass, normalized otherwise).
        #[arg(short, long, value_delimiter = ',', allow_hyphen_values = true)]
        zeros: Option<Vec<f64>>,

        /// Unloaded Q factor.
        #[arg(short, long)]
        q: Option<f64>,

        /// Output topology: transversal, folded, arrow.
        #[arg(short, long, default_value = "folded")]
        topology: String,

        /// Number of frequency points for response.
        #[arg(long, default_value = "201")]
        points: usize,

        /// Response frequency range start (Hz). Defaults to 0.5·bw below the outermost passband lower edge.
        #[arg(long)]
        resp_start: Option<f64>,

        /// Response frequency range stop (Hz). Defaults to 0.5·bw above the outermost passband upper edge.
        #[arg(long)]
        resp_stop: Option<f64>,

        /// Export Touchstone file to this path.
        #[arg(long, value_name = "FILE")]
        s2p: Option<PathBuf>,
    },

    /// Synthesize an N-channel multiplexer from repeated channel specs.
    Multiplexer {
        /// Channel spec as `order:return_loss_db:center_hz:bandwidth_hz`.
        /// Repeat for each channel (at least two required).
        #[arg(long, value_name = "CHANNEL", required = true)]
        channel: Vec<String>,

        /// Maximum number of iterative convergence steps.
        #[arg(long, default_value = "40")]
        max_iterations: usize,

        /// Convergence tolerance on root displacement.
        #[arg(long, default_value = "1e-8")]
        tol: f64,

        /// Response frequency range start (Hz). Defaults to the lowest channel
        /// center minus 1.5× its bandwidth.
        #[arg(long)]
        resp_start: Option<f64>,

        /// Response frequency range stop (Hz). Defaults to the highest channel
        /// center plus 1.5× its bandwidth.
        #[arg(long)]
        resp_stop: Option<f64>,

        /// Number of frequency points for the response.
        #[arg(long, default_value = "401")]
        points: usize,
    },

    /// Optimize a coupling matrix against a target response.
    ///
    /// The input matrix and the target response are both read from the JSON
    /// file named by --input / stdin (matrix may be flat `data` or nested
    /// rows; target may be a Touchstone .s2p file named by --target).
    ///
    /// §C.7: the target's `frequency_hz` samples are used **directly as the
    /// normalized grid ω** (the optimizer evaluates with `evaluate_normalized`);
    /// there is no `--normalized` flag — the grid is always normalized. A
    /// target whose frequency range spans > 1e6 (physical Hz) is rejected.
    Optimize {
        /// Filter order (number of resonators).
        #[arg(short = 'n', long)]
        order: usize,

        /// Target response file (Touchstone .s2p). If omitted, the target
        /// is read from the JSON input.
        #[arg(long, value_name = "FILE")]
        target: Option<PathBuf>,

        /// Maximum number of Levenberg-Marquardt iterations.
        #[arg(long, default_value = "200")]
        max_iterations: usize,

        /// Convergence tolerance on the cost change.
        #[arg(long, default_value = "1e-8")]
        tolerance: f64,

        /// Initial Levenberg-Marquardt damping factor.
        #[arg(long, default_value = "1e-3")]
        lambda: f64,

        /// Unloaded Q (finite = lossy, omit for lossless).
        #[arg(long)]
        q: Option<f64>,

        /// Export the optimized matrix as Touchstone-free JSON to this path.
        #[arg(long, value_name = "FILE")]
        output: Option<PathBuf>,
    },

    /// Extract a coupling matrix from a measured response.
    ///
    /// Builds a crude inline-chain guess for the given order and optimizes
    /// it against the measured response (Touchstone .s2p file or JSON input).
    ///
    /// §C.7: as with `optimize`, the measured response's `frequency_hz`
    /// samples are used **directly as the normalized grid ω** (no
    /// `--normalized` flag); a frequency range spanning > 1e6 (physical Hz)
    /// is rejected.
    Extract {
        /// Filter order (number of resonators).
        #[arg(short = 'n', long)]
        order: usize,

        /// Measured response file (Touchstone .s2p or JSON). If omitted,
        /// the response is read from the JSON input / stdin.
        #[arg(long, value_name = "FILE")]
        target: Option<PathBuf>,

        /// Export the extracted matrix JSON to this path.
        #[arg(long, value_name = "FILE")]
        output: Option<PathBuf>,
    },
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

fn main() {
    let cli = Cli::parse();

    match run(&cli) {
        Ok(()) => {}
        Err(e) => {
            // §D.5: structured error JSON + fallback rendered in `io`.
            let json = format_cli_error(error_type_name(&e), e.to_string());
            eprintln!("{json}");
            process::exit(1);
        }
    }
}

fn run(cli: &Cli) -> Result<(), CliRunError> {
    // Handle subcommands first
    if let Some(Command::Design {
        order,
        return_loss,
        center,
        bandwidth,
        start,
        stop,
        zeros,
        q,
        topology,
        points,
        resp_start,
        resp_stop,
        s2p,
    }) = &cli.command
    {
        let args = DesignArgs {
            order: *order,
            return_loss: *return_loss,
            center: *center,
            bandwidth: *bandwidth,
            start: *start,
            stop: *stop,
            zeros: zeros.as_deref(),
            q: *q,
            topology,
            points: *points,
            resp_start: *resp_start,
            resp_stop: *resp_stop,
            s2p: s2p.as_deref(),
        };
        return run_design(&args, &cli.format);
    }

    if let Some(Command::Multiplexer {
        channel,
        max_iterations,
        tol,
        resp_start,
        resp_stop,
        points,
    }) = &cli.command
    {
        return run_multiplexer(
            channel,
            *max_iterations,
            *tol,
            *resp_start,
            *resp_stop,
            *points,
            &cli.format,
        );
    }

    if let Some(Command::Optimize {
        order,
        target,
        max_iterations,
        tolerance,
        lambda,
        q,
        output,
    }) = &cli.command
    {
        use mfs::optimize::OptimizeConfig;
        let config = OptimizeConfig {
            max_iterations: *max_iterations,
            tolerance: *tolerance,
            lambda: *lambda,
            unloaded_q: q.filter(|q| q.is_finite() && *q > 0.0).unwrap_or(f64::INFINITY),
            ..Default::default()
        };
        return run_optimize(cli, *order, target.as_deref(), &config, output.as_deref(), &cli.format);
    }

    if let Some(Command::Extract { order, target, output }) = &cli.command {
        return run_extract(cli, *order, target.as_deref(), output.as_deref(), &cli.format);
    }

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
            let mut ctx = context_from_request(request).map_err(CliRunError::Pipeline)?;

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
