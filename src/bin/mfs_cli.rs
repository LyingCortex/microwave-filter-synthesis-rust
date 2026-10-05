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
    context_from_request, load_context, run_from_json, run_full_pipeline, run_stage,
    SynthesisContext, SynthesisRequest,
};

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

        /// Response frequency range start (Hz). Defaults to center - 1.5*bw.
        #[arg(long)]
        resp_start: Option<f64>,

        /// Response frequency range stop (Hz). Defaults to center + 1.5*bw.
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
            let json = serde_json::to_string(&cli_error)
                .unwrap_or_else(|_| format!(r#"{{"error_type":"internal","message":"{}"}}"#, e));
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
        return run_design(
            *order,
            *return_loss,
            *center,
            *bandwidth,
            *start,
            *stop,
            zeros.as_deref(),
            *q,
            topology,
            *points,
            *resp_start,
            *resp_stop,
            s2p.as_deref(),
            &cli.format,
        );
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

/// Runs the quick `design` subcommand.
fn run_design(
    order: usize,
    return_loss: f64,
    center: Option<f64>,
    bandwidth: Option<f64>,
    start: Option<f64>,
    stop: Option<f64>,
    zeros: Option<&[f64]>,
    q: Option<f64>,
    topology: &str,
    points: usize,
    resp_start: Option<f64>,
    resp_stop: Option<f64>,
    s2p: Option<&std::path::Path>,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::design::FilterDesign;
    use mfs::matrix::MatrixTopology;
    use mfs::touchstone::{self, DataFormat, FreqUnit, TouchstoneConfig, TouchstoneVersion};

    // Resolve center/bandwidth from either (center, bw) or (start, stop) notation
    let (resolved_center, resolved_bw) = match (center, bandwidth, start, stop) {
        (Some(c), Some(bw), _, _) => (Some(c), Some(bw)),
        (None, None, Some(f1), Some(f2)) => {
            if f2 <= f1 {
                return Err(CliRunError::InputParse(
                    "--stop must be greater than --start".into(),
                ));
            }
            // Geometric-mean center, matching FilterDesign::bandpass_edges so the
            // Ω(f1)=-1/Ω(f2)=+1 convention holds (single source in freq.rs).
            let (c, bw) = mfs::freq::edges_to_center_bandwidth(f1, f2);
            (Some(c), Some(bw))
        }
        (Some(c), None, Some(f1), _) => {
            // center + start → derive bandwidth = 2*(center - start)
            let bw = 2.0 * (c - f1);
            if bw <= 0.0 {
                return Err(CliRunError::InputParse(
                    "--start must be less than --center".into(),
                ));
            }
            (Some(c), Some(bw))
        }
        _ => {
            // No resolvable (center, bandwidth) pair. A pure prototype is only
            // legitimate when NO band-pass parameter was supplied; otherwise a
            // parameter would be silently ignored (review P2-2). Arm 1 above
            // already handles (center, bandwidth), so reaching here with any
            // band-pass parameter set means the combination is incomplete.
            if center.is_none() && bandwidth.is_none() && start.is_none() && stop.is_none() {
                (None, None) // pure normalized prototype: no band-pass params
            } else {
                let parts: Vec<String> = [
                    center.map(|c| format!("--center {c}")),
                    bandwidth.map(|b| format!("--bandwidth {b}")),
                    start.map(|f1| format!("--start {f1}")),
                    stop.map(|f2| format!("--stop {f2}")),
                ]
                .into_iter()
                .flatten()
                .collect();
                return Err(CliRunError::InputParse(format!(
                    "incomplete band-pass parameters [{}]: pair --center with \
                     --bandwidth or --start, or use --start + --stop",
                    parts.join(", ")
                )));
            }
        }
    };

    // Build the design
    let design = if let (Some(c), Some(bw)) = (resolved_center, resolved_bw) {
        let mut builder = FilterDesign::bandpass(order, return_loss, c, bw);
        if let Some(z) = zeros {
            builder = builder.zeros_hz(z.iter().copied());
        }
        if let Some(q_val) = q {
            builder = builder.unloaded_q(q_val);
        }
        builder.synthesize().map_err(CliRunError::Pipeline)?
    } else {
        let mut builder = FilterDesign::prototype(order, return_loss);
        if let Some(z) = zeros {
            builder = builder.zeros(z.iter().copied());
        }
        if let Some(q_val) = q {
            builder = builder.unloaded_q(q_val);
        }
        builder.synthesize().map_err(CliRunError::Pipeline)?
    };

    // Get topology matrix
    let topo = match topology.to_lowercase().as_str() {
        "arrow" => MatrixTopology::Arrow,
        "transversal" => MatrixTopology::Transversal,
        _ => MatrixTopology::Folded,
    };
    let matrix = design.to_topology(topo).map_err(CliRunError::Pipeline)?;

    // Get response
    let response = if let (Some(c), Some(bw)) = (resolved_center, resolved_bw) {
        // Default response range: center ± 1.5*bw (covers passband + nearby stopband)
        let r_start = resp_start.unwrap_or(c - 1.5 * bw);
        let r_stop = resp_stop.unwrap_or(c + 1.5 * bw);
        if let Some(q_val) = q {
            design
                .response_lossy(r_start, r_stop, points, q_val)
                .map_err(CliRunError::Pipeline)?
        } else {
            design
                .response(r_start, r_stop, points)
                .map_err(CliRunError::Pipeline)?
        }
    } else {
        let r_start = resp_start.unwrap_or(-3.0);
        let r_stop = resp_stop.unwrap_or(3.0);
        if let Some(q_val) = q {
            design
                .response_lossy_normalized(r_start, r_stop, points, q_val)
                .map_err(CliRunError::Pipeline)?
        } else {
            design
                .response_normalized(r_start, r_stop, points)
                .map_err(CliRunError::Pipeline)?
        }
    };

    // Export Touchstone if requested
    if let Some(s2p_path) = s2p {
        let config = TouchstoneConfig {
            freq_unit: FreqUnit::GHz,
            format: DataFormat::RI,
            impedance: 50.0,
            version: TouchstoneVersion::V1,
            comments: design.auto_comments(),
        };
        touchstone::write_touchstone(&response, &config, s2p_path)
            .map_err(CliRunError::Pipeline)?;
        eprintln!("Touchstone file written to: {}", s2p_path.display());
    }

    // Output
    match format {
        OutputFormat::Json => {
            let side = matrix.side();
            let matrix_data: Vec<Vec<f64>> = (0..side)
                .map(|r| (0..side).map(|c| matrix.at(r, c).unwrap_or(0.0)).collect())
                .collect();
            let samples: Vec<serde_json::Value> = response
                .samples
                .iter()
                .map(|s| {
                    serde_json::json!({
                        "frequency_hz": s.frequency_hz,
                        "s11_db": s.s11_db(),
                        "s21_db": s.s21_db(),
                        "s11_phase_deg": s.s11_phase_deg(),
                        "s21_phase_deg": s.s21_phase_deg(),
                    })
                })
                .collect();

            let output = serde_json::json!({
                "order": order,
                "return_loss_db": return_loss,
                "topology": topology,
                "matrix": matrix_data,
                "response": {
                    "points": response.samples.len(),
                    "samples": samples,
                },
            });
            println!(
                "{}",
                serde_json::to_string_pretty(&output)
                    .map_err(|e| CliRunError::Serialization(e.to_string()))?
            );
        }
        OutputFormat::Table => {
            println!("── Filter Design ──");
            println!("  Order:       {order}");
            println!("  Return Loss: {return_loss:.1} dB");
            if let (Some(c), Some(bw)) = (resolved_center, resolved_bw) {
                println!("  Center:      {:.6} GHz", c / 1e9);
                println!("  Bandwidth:   {:.3} MHz", bw / 1e6);
                println!(
                    "  Passband:    {:.6} – {:.6} GHz",
                    (c - bw / 2.0) / 1e9,
                    (c + bw / 2.0) / 1e9
                );
            }
            if let Some(z) = zeros {
                let zs: Vec<String> = z.iter().map(|v| format!("{v:.4}")).collect();
                println!("  Zeros:       [{}]", zs.join(", "));
            }
            if let Some(q_val) = q {
                println!("  Unloaded Q:  {q_val:.0}");
            }
            println!("  Topology:    {topology}");
            println!();

            // Matrix
            let side = matrix.side();
            if side <= 10 {
                println!("── Coupling Matrix ({topology}) ──");
                for i in 0..side {
                    let row: Vec<String> = (0..side)
                        .map(|j| format!("{:>8.4}", matrix.at(i, j).unwrap_or(0.0)))
                        .collect();
                    println!("  [{}]", row.join(" "));
                }
                println!();
            }

            // Response summary
            println!(
                "── S-Parameter Response ({} points) ──",
                response.samples.len()
            );
            println!("  {:>12} {:>10} {:>10}", "Freq", "S21 (dB)", "S11 (dB)");
            println!("  {:>12} {:>10} {:>10}", "────────", "────────", "────────");
            let n = response.samples.len();
            let show = 5.min(n);
            for s in response.samples.iter().take(show) {
                let freq_label = if resolved_center.is_some() {
                    format!("{:.4} GHz", s.frequency_hz / 1e9)
                } else {
                    format!("{:.4}", s.frequency_hz)
                };
                println!(
                    "  {:>12} {:>10.3} {:>10.3}",
                    freq_label,
                    s.s21_db(),
                    s.s11_db()
                );
            }
            if n > show * 2 {
                println!("  {:>12}", "...");
            }
            if n > show {
                for s in response.samples.iter().skip(n - show) {
                    let freq_label = if resolved_center.is_some() {
                        format!("{:.4} GHz", s.frequency_hz / 1e9)
                    } else {
                        format!("{:.4}", s.frequency_hz)
                    };
                    println!(
                        "  {:>12} {:>10.3} {:>10.3}",
                        freq_label,
                        s.s21_db(),
                        s.s11_db()
                    );
                }
            }
        }
    }

    Ok(())
}

/// Runs the quick `multiplexer` subcommand: N-channel synthesis + analysis.
fn run_multiplexer(
    channel: &[String],
    max_iterations: usize,
    tol: f64,
    resp_start: Option<f64>,
    resp_stop: Option<f64>,
    points: usize,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::multiplexer::{Multiplexer, MultiplexerSynthesis};

    // Parse each `order:rl:center:bw` string into a ChannelSpec.
    let mut specs: Vec<mfs::ChannelSpec> = Vec::with_capacity(channel.len());
    for raw in channel {
        let parts: Vec<&str> = raw.split(':').collect();
        if parts.len() != 4 {
            return Err(CliRunError::InputParse(format!(
                "channel spec '{raw}' must be formatted as order:return_loss_db:center_hz:bandwidth_hz"
            )));
        }
        let order = parts[0].trim().parse::<usize>().map_err(|_| {
            CliRunError::InputParse(format!(
                "invalid order '{}' in channel spec '{raw}'",
                parts[0]
            ))
        })?;
        let return_loss_db = parts[1].trim().parse::<f64>().map_err(|_| {
            CliRunError::InputParse(format!(
                "invalid return loss '{}' in channel spec '{raw}'",
                parts[1]
            ))
        })?;
        let center_hz = parts[2].trim().parse::<f64>().map_err(|_| {
            CliRunError::InputParse(format!(
                "invalid center frequency '{}' in channel spec '{raw}'",
                parts[2]
            ))
        })?;
        let bandwidth_hz = parts[3].trim().parse::<f64>().map_err(|_| {
            CliRunError::InputParse(format!(
                "invalid bandwidth '{}' in channel spec '{raw}'",
                parts[3]
            ))
        })?;
        specs.push(mfs::ChannelSpec {
            order,
            return_loss_db,
            zeros: Vec::new(),
            center_hz,
            bandwidth_hz,
        });
    }

    // Run iterative N-channel synthesis.
    let mut synthesis = MultiplexerSynthesis::new()
        .channels(specs.clone())
        .max_iterations(max_iterations)
        .convergence_tol(tol);
    if specs.len() == 2 {
        // Mirror the diplexer's default ideal T-junction.
        synthesis = synthesis.junction(mfs::JunctionType::TypeI { n: 1.0, b0: 0.0 });
    }
    let result = synthesis.synthesize().map_err(CliRunError::Pipeline)?;

    // Build the multiplexer and compute a response.
    let mux = Multiplexer::new(result.channels.clone()).map_err(CliRunError::Pipeline)?;

    let (r_start, r_stop) = resolve_multiplexer_range(&result.channels, resp_start, resp_stop);
    let response = mux
        .response(r_start, r_stop, points)
        .map_err(CliRunError::Pipeline)?;

    // Output.
    match format {
        OutputFormat::Json => {
            let channels: Vec<serde_json::Value> = result
                .channels
                .iter()
                .map(|ch| {
                    serde_json::json!({
                        "order": ch.order(),
                        "center_hz": ch.center_hz(),
                        "bandwidth_hz": ch.bandwidth_hz(),
                    })
                })
                .collect();
            let output = serde_json::json!({
                "channels": channels,
                "iterations": result.iterations,
                "convergence_error": result.convergence_error,
                "junction": format!("{:?}", result.junction),
                "response": {
                    "frequencies": response.frequencies,
                    "antenna_s11_db": response.antenna_s11_db,
                    "channel_s21_db": response.channel_s21_db,
                    "channel_s11_db": response.channel_s11_db,
                    "isolation_db": response.isolation_db,
                },
            });
            println!(
                "{}",
                serde_json::to_string_pretty(&output)
                    .map_err(|e| CliRunError::Serialization(e.to_string()))?
            );
        }
        OutputFormat::Table => {
            println!("── Multiplexer Synthesis ──");
            println!("  Channels:      {}", result.channels.len());
            println!("  Iterations:    {}", result.iterations);
            println!("  Conv. Error:   {:.3e}", result.convergence_error);
            println!("  Junction:      {:?}", result.junction);
            println!();
            for (i, ch) in result.channels.iter().enumerate() {
                println!(
                    "  Channel {}: order={}, center={:.6} GHz, bw={:.3} MHz",
                    i,
                    ch.order(),
                    ch.center_hz().unwrap_or(0.0) / 1e9,
                    ch.bandwidth_hz().unwrap_or(0.0) / 1e6,
                );
                match mux.channel_folded(i) {
                    Ok(matrix) => {
                        let side = matrix.side();
                        if side <= 10 {
                            println!("    Folded coupling matrix:");
                            for r in 0..side {
                                let row: Vec<String> = (0..side)
                                    .map(|c| format!("{:>8.4}", matrix.at(r, c).unwrap_or(0.0)))
                                    .collect();
                                println!("    [{}]", row.join(" "));
                            }
                        }
                    }
                    Err(e) => println!("    (folded matrix unavailable: {e})"),
                }
                println!();
            }
            println!("── Response ({} points) ──", response.frequencies.len());
            println!("  Freq: {:.4} GHz – {:.4} GHz", r_start / 1e9, r_stop / 1e9);
            println!(
                "  Antenna S11:   min {:.1} dB",
                min_db(&response.antenna_s11_db)
            );
            for (i, s21) in response.channel_s21_db.iter().enumerate() {
                println!("  Channel {i} S21: peak {:.1} dB", max_db(s21));
            }
        }
    }

    Ok(())
}
// ---------------------------------------------------------------------------
// Optimizer subcommands (D3-opt): `optimize` and `extract`
// ---------------------------------------------------------------------------

/// Parses a measured / target S-parameter response from JSON.
///
/// Accepts `{"samples":[...]}` or a bare array, each sample carrying
/// `frequency_hz` plus any of `s11_re/s11_im` and `s21_re/s21_im` (missing
/// channels default to 0). Returns a `PreconditionViolation` on malformed
/// input so the CLI reports it as a precondition error, not a crash.
fn parse_response_json(json: &str) -> mfs::Result<mfs::response::SParameterResponse> {
    use mfs::response::{ResponseSample, SParameterResponse};
    let v: serde_json::Value = serde_json::from_str(json)
        .map_err(|e| mfs::error::MfsError::PreconditionViolation(format!("invalid target JSON: {e}")))?;
    let arr = v
        .get("samples")
        .and_then(|s| s.as_array())
        .or_else(|| v.as_array())
        .ok_or_else(|| {
            mfs::error::MfsError::PreconditionViolation(
                "target JSON must contain a `samples` array".into(),
            )
        })?;
    let samples = arr
        .iter()
        .map(|s| {
            let get = |k: &str| s.get(k).and_then(|x| x.as_f64()).unwrap_or(0.0);
            Ok(ResponseSample {
                frequency_hz: get("frequency_hz"),
                normalized_omega: 0.0,
                group_delay: None,
                s11_re: get("s11_re"),
                s11_im: get("s11_im"),
                s21_re: get("s21_re"),
                s21_im: get("s21_im"),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(SParameterResponse { samples })
}

/// Parses a coupling matrix from JSON.
///
/// Accepts `{"order": N, "data": [row-major]}` / `{"data": [[rows]]}` /
/// `{"order": N, "matrix": [flat]}` / `{"matrix": [[rows]]}`. The topology is
/// taken from the input if present (one of `transversal`/`folded`/`arrow`),
/// else defaulted to transversal.
fn parse_matrix_json(json: &str) -> mfs::Result<mfs::matrix::CouplingMatrix> {
    use mfs::error::MfsError;
    use mfs::matrix::{CouplingMatrix, MatrixTopology};
    let v: serde_json::Value = serde_json::from_str(json)
        .map_err(|e| MfsError::PreconditionViolation(format!("invalid matrix JSON: {e}")))?;

    // Locate the raw matrix value: an object with "data"/"matrix", a bare
    // object, or a bare array.
    let (raw, order_hint, topo_hint) = if v.is_object() {
        (
            v.get("data")
                .cloned()
                .or_else(|| v.get("matrix").cloned()),
            v.get("order").and_then(|x| x.as_u64()).map(|n| n as usize),
            v.get("topology").and_then(|t| t.as_str()),
        )
    } else {
        (Some(v.clone()), None, None)
    };
    let raw = raw.ok_or_else(|| {
        MfsError::PreconditionViolation(
            "matrix JSON must be {{order, data}} or an array (flat row-major or nested rows)"
                .into(),
        )
    })?;

    // Accept either a flat row-major vector or nested rows.
    let rows: Vec<Vec<f64>> = if let Some(a) = raw.as_array() {
        match a.first() {
            Some(inner) if inner.is_array() => {
                // Nested rows: array of arrays.
                a.iter()
                    .map(|r| {
                        r.as_array().ok_or_else(|| {
                            MfsError::PreconditionViolation("matrix rows must be arrays".into())
                        })?
                        .iter()
                        .map(|x| x.as_f64())
                        .collect::<Option<Vec<f64>>>()
                        .ok_or_else(|| {
                            MfsError::PreconditionViolation("matrix entries must be numbers".into())
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            _ => {
                // Flat row-major vector: reshape into a square matrix.
                let flat: Vec<f64> = a
                    .iter()
                    .map(|x| x.as_f64())
                    .collect::<Option<Vec<f64>>>()
                    .ok_or_else(|| {
                        MfsError::PreconditionViolation(
                            "matrix `data` must be a flat array of numbers".into(),
                        )
                    })?;
                let side = (flat.len() as f64).sqrt().round() as usize;
                if side * side != flat.len() || side < 3 {
                    return Err(MfsError::DimensionMismatch {
                        expected: side * side,
                        actual: flat.len(),
                    });
                }
                flat.chunks_exact(side).map(|c| c.to_vec()).collect()
            }
        }
    } else {
        return Err(MfsError::PreconditionViolation(
            "matrix JSON must be an array or an object with `data`/`matrix`".into(),
        ));
    };

    let side = rows.len();
    for r in &rows {
        if r.len() != side {
            return Err(MfsError::DimensionMismatch {
                expected: side,
                actual: r.len(),
            });
        }
    }
    if side < 3 {
        return Err(MfsError::PreconditionViolation(
            "matrix side must be >= 3 (order >= 1)".into(),
        ));
    }

    let order = match order_hint {
        Some(o) => {
            if o + 2 != side {
                return Err(MfsError::DimensionMismatch {
                    expected: (o + 2) * (o + 2),
                    actual: side * side,
                });
            }
            o
        }
        None => side - 2,
    };
    let topology = match topo_hint {
        Some("folded") => MatrixTopology::Folded,
        Some("arrow") => MatrixTopology::Arrow,
        Some("transversal") | None => MatrixTopology::Transversal,
        Some(other) => {
            return Err(MfsError::PreconditionViolation(format!(
                "unknown matrix topology '{other}' (expected transversal|folded|arrow)"
            )));
        }
    };
    CouplingMatrix::new_with_topology(order, topology, rows.into_iter().flatten().collect())
}

/// Loads a target response from a Touchstone file (`--target` pointing at a
/// `.s2p`) or from the JSON string otherwise.
fn load_target(
    target: Option<&std::path::Path>,
    json: &str,
) -> Result<mfs::response::SParameterResponse, CliRunError> {
    match target {
        Some(path) if path.extension().and_then(|e| e.to_str()) == Some("s2p") => {
            mfs::touchstone::read_touchstone(path).map_err(|e| CliRunError::Touchstone(e.to_string()))
        }
        Some(path) => Err(CliRunError::InputParse(format!(
            "--target must be a Touchstone (.s2p) file, got '{}'",
            path.display()
        ))),
        None => parse_response_json(json).map_err(CliRunError::Pipeline),
    }
}

/// Re-evaluates a matrix on the grid defined by the target's frequency samples,
/// for display next to the target.
fn recovered_response(
    matrix: &mfs::matrix::CouplingMatrix,
    target: &mfs::response::SParameterResponse,
) -> mfs::Result<mfs::response::SParameterResponse> {
    use mfs::freq::FrequencyGrid;
    use mfs::response::ResponseSolver;
    let n = target.samples.len();
    let first = target.samples.first().map(|s| s.frequency_hz).unwrap_or(-2.0);
    let last = target.samples.last().map(|s| s.frequency_hz).unwrap_or(2.0);
    let grid = FrequencyGrid::linspace(first, last, n.max(2))?;
    ResponseSolver.evaluate_normalized(matrix, &grid)
}

/// Renders the optimized / extracted matrix as JSON.
fn optimize_result_json(
    order: usize,
    result: &mfs::optimize::OptimizeResult,
    target: &mfs::response::SParameterResponse,
) -> Result<serde_json::Value, CliRunError> {
    let side = result.matrix.side();
    let rows: Vec<Vec<f64>> = (0..side)
        .map(|r| (0..side).map(|c| result.matrix.at(r, c).unwrap_or(0.0)).collect())
        .collect();
    let mut value = serde_json::json!({
        "order": order,
        "topology": format!("{:?}", result.matrix.topology()),
        "matrix": rows,
        "data": rows.iter().flatten().copied().collect::<Vec<f64>>(),
        "cost": result.cost,
        "iterations": result.iterations,
        "converged": result.converged,
    });
    if let Ok(rec) = recovered_response(&result.matrix, target) {
        let samples: Vec<serde_json::Value> = rec
            .samples
            .iter()
            .map(|s| {
                serde_json::json!({
                    "frequency_hz": s.frequency_hz,
                    "s11_re": s.s11_re,
                    "s11_im": s.s11_im,
                    "s21_re": s.s21_re,
                    "s21_im": s.s21_im,
                })
            })
            .collect();
        value["recovered"] = serde_json::json!({
            "samples": samples,
        });
    }
    Ok(value)
}

/// Runs the `optimize` subcommand.
fn run_optimize(
    cli: &Cli,
    order: usize,
    target: Option<&std::path::Path>,
    config: &mfs::optimize::OptimizeConfig,
    output: Option<&std::path::Path>,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::optimize::optimize_matrix;
    // The matrix always comes from --input/stdin (JSON); only the *target*
    // may be a .s2p. The documented "matrix from input + Touchstone target"
    // usage requires reading the matrix even when the target is a .s2p
    // (mirrors `run_extract`, review P2-1).
    let json_input = read_input(cli)?;
    let target_resp = load_target(target, &json_input)?;
    let matrix = parse_matrix_json(&json_input)?;
    if matrix.order() != order {
        return Err(CliRunError::InputParse(format!(
            "--order {order} does not match the input matrix order {}",
            matrix.order()
        )));
    }
    let result = optimize_matrix(&matrix, &target_resp, &response_grid(&target_resp)?, config)
        .map_err(CliRunError::Pipeline)?;
    let value = optimize_result_json(order, &result, &target_resp)?;
    if let Some(out) = output {
        std::fs::write(out, serde_json::to_string_pretty(&value)?)
            .map_err(|e| CliRunError::InputRead(format!("failed to write '{}': {e}", out.display())))?;
        eprintln!("Optimized matrix written to: {}", out.display());
    }
    match format {
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Table => {
            println!("── Optimized Matrix ──");
            println!("  Order:     {order}");
            println!("  Topology:  {:?}", result.matrix.topology());
            println!("  Cost:      {:.3e}", result.cost);
            println!("  Iters:     {}", result.iterations);
            println!("  Converged: {}", result.converged);
            let side = result.matrix.side();
            if side <= 10 {
                println!("  Matrix:");
                for i in 0..side {
                    let row: Vec<String> = (0..side)
                        .map(|j| format!("{:>8.4}", result.matrix.at(i, j).unwrap_or(0.0)))
                        .collect();
                    println!("    [{}]", row.join(" "));
                }
            }
        }
    }
    Ok(())
}

/// Builds the grid the optimizer evaluates on, from the target response.
fn response_grid(target: &mfs::response::SParameterResponse) -> Result<mfs::freq::FrequencyGrid, CliRunError> {
    use mfs::freq::FrequencyGrid;
    let n = target.samples.len().max(2);
    let first = target.samples.first().map(|s| s.frequency_hz).unwrap_or(-2.0);
    let last = target.samples.last().map(|s| s.frequency_hz).unwrap_or(2.0);
    FrequencyGrid::linspace(first, last, n).map_err(CliRunError::Pipeline)
}

/// Runs the `extract` subcommand.
fn run_extract(
    cli: &Cli,
    order: usize,
    target: Option<&std::path::Path>,
    output: Option<&std::path::Path>,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::optimize::extract_matrix;
    let json_input = read_input(cli)?;
    let target_resp = load_target(target, &json_input)?;
    let result = extract_matrix(&target_resp, &response_grid(&target_resp)?, order)
        .map_err(CliRunError::Pipeline)?;
    let value = optimize_result_json(order, &result, &target_resp)?;
    if let Some(out) = output {
        std::fs::write(out, serde_json::to_string_pretty(&value)?)
            .map_err(|e| CliRunError::InputRead(format!("failed to write '{}': {e}", out.display())))?;
        eprintln!("Extracted matrix written to: {}", out.display());
    }
    match format {
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Table => {
            println!("── Extracted Matrix ──");
            println!("  Order:     {order}");
            println!("  Topology:  {:?}", result.matrix.topology());
            println!("  Cost:      {:.3e}", result.cost);
            println!("  Iters:     {}", result.iterations);
            println!("  Converged: {}", result.converged);
            let side = result.matrix.side();
            if side <= 10 {
                println!("  Matrix:");
                for i in 0..side {
                    let row: Vec<String> = (0..side)
                        .map(|j| format!("{:>8.4}", result.matrix.at(i, j).unwrap_or(0.0)))
                        .collect();
                    println!("    [{}]", row.join(" "));
                }
            }
        }
    }
    Ok(())
}
/// Chooses a response frequency range, defaulting to the outer channel bands.
fn resolve_multiplexer_range(
    channels: &[mfs::design::FilterDesign],
    resp_start: Option<f64>,
    resp_stop: Option<f64>,
) -> (f64, f64) {
    let lo = channels
        .iter()
        .map(|ch| ch.center_hz().unwrap_or(0.0) - 1.5 * ch.bandwidth_hz().unwrap_or(0.0))
        .fold(f64::INFINITY, f64::min);
    let hi = channels
        .iter()
        .map(|ch| ch.center_hz().unwrap_or(0.0) + 1.5 * ch.bandwidth_hz().unwrap_or(0.0))
        .fold(f64::NEG_INFINITY, f64::max);
    (resp_start.unwrap_or(lo), resp_stop.unwrap_or(hi))
}

fn min_db(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

fn max_db(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Determines which stages haven't been executed yet and runs them in order.
fn run_remaining_stages(ctx: &mut SynthesisContext) -> mfs::Result<()> {
    let all_stages = [
        "approximation",
        "matrix_synthesis",
        "topology_transform",
        "response_evaluation",
    ];

    for stage_name in &all_stages {
        if ctx
            .metadata
            .stages_executed
            .contains(&stage_name.to_string())
        {
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
                println!(
                    "  {:>14} {:>10} {:>10}",
                    "Freq (Hz)", "S11 (dB)", "S21 (dB)"
                );
                println!(
                    "  {:>14} {:>10} {:>10}",
                    "──────────", "────────", "────────"
                );
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
                        let s11_db = 10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
                        let s21_db = 10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
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
                let s11_db = 10.0 * (sample.s11_re.powi(2) + sample.s11_im.powi(2)).log10();
                let s21_db = 10.0 * (sample.s21_re.powi(2) + sample.s21_im.powi(2)).log10();
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
enum CliRunError {
    /// Failed to read input file or stdin.
    InputRead(String),
    /// Failed to parse input JSON.
    InputParse(String),
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
            Self::Pipeline(err) => write!(f, "{err}"),
            Self::Serialization(msg) => write!(f, "serialization error: {msg}"),
            Self::Touchstone(msg) => write!(f, "touchstone error: {msg}"),
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
        CliRunError::Touchstone(_) => "touchstone_error".to_string(),
    }
}

#[cfg(test)]
mod optimizer_tests {
    use super::*;
    use mfs::design::FilterDesign;
    use mfs::freq::FrequencyGrid;
    use mfs::response::ResponseSolver;

    /// Scales the first `count` non-zero inter-resonator couplings (upper
    /// triangle, off-diagonal) of `matrix` by `factor`, mirroring the
    /// integration-test perturbation helper.
    fn perturb(matrix: &mfs::matrix::CouplingMatrix, count: usize, factor: f64) -> mfs::matrix::CouplingMatrix {
        let side = matrix.side();
        let order = matrix.order();
        let mut data = matrix.as_slice().to_vec();
        let mut done = 0;
        for row in 1..=order {
            for col in (row + 1)..=order {
                let idx = row * side + col;
                if data[idx].abs() > 1e-9 {
                    data[idx] *= factor;
                    data[col * side + row] *= factor;
                    done += 1;
                    if done == count {
                        return mfs::matrix::CouplingMatrix::new_with_topology(
                            order,
                            matrix.topology(),
                            data,
                        )
                        .expect("rebuild");
                    }
                }
            }
        }
        mfs::matrix::CouplingMatrix::new_with_topology(order, matrix.topology(), data).expect("rebuild")
    }

    /// Worst relative deviation between two matrices over entries with |orig| > 1e-9.
    fn worst_rel(a: &mfs::matrix::CouplingMatrix, b: &mfs::matrix::CouplingMatrix) -> f64 {
        let side = a.side();
        let mut w = 0.0_f64;
        for i in 0..side {
            for j in i..side {
                let o = a.at(i, j).unwrap_or(0.0);
                let r = b.at(i, j).unwrap_or(0.0);
                if o.abs() > 1e-9 {
                    w = w.max((o - r).abs() / o.abs());
                }
            }
        }
        w
    }

    #[test]
    fn parse_matrix_json_flat_and_nested() {
        // Flat row-major with order + topology.
        let m = parse_matrix_json(r#"{"order":1,"data":[0,0.7,0.7,0,0.7,0,0,0,0]}"#).unwrap();
        assert_eq!(m.order(), 1);
        assert!((m.at(0, 1).unwrap() - 0.7).abs() < 1e-12);
        assert_eq!(m.topology(), mfs::matrix::MatrixTopology::Transversal);

        // Folded topology honored.
        let m = parse_matrix_json(r#"{"order":1,"topology":"folded","data":[0,0.7,0.7,0,0.7,0,0,0,0]}"#)
            .unwrap();
        assert_eq!(m.topology(), mfs::matrix::MatrixTopology::Folded);

        // Nested rows, order inferred (side 4 => order 2).
        let m = parse_matrix_json(r#"[ [0,0,0,0],[0,1,1,0],[0,1,2,0],[0,0,0,0] ]"#).unwrap();
        assert_eq!(m.order(), 2);
        assert!((m.at(1, 2).unwrap() - 1.0).abs() < 1e-12);

        // Order / side mismatch is a typed dimension error.
        let err = parse_matrix_json(r#"{"order":2,"data":[0,0.7,0.7,0,0.7,0,0,0,0]}"#).unwrap_err();
        assert!(matches!(err, mfs::error::MfsError::DimensionMismatch { .. }));
    }

    #[test]
    fn parse_response_json_samples() {
        let r = parse_response_json(
            r#"{"samples":[
               {"frequency_hz":-1.0,"s11_re":0.1,"s11_im":-0.2,"s21_re":0.9,"s21_im":0.0},
               {"frequency_hz":1.0,"s11_re":0.1,"s11_im":0.2,"s21_re":0.9,"s21_im":0.0}
            ]}"#,
        )
        .unwrap();
        assert_eq!(r.samples.len(), 2);
        assert!((r.samples[0].frequency_hz - (-1.0)).abs() < 1e-12);
        assert!((r.samples[1].s21_im - 0.0).abs() < 1e-12);
    }

    #[test]
    fn optimize_subcommand_recovers_perturbed_matrix() {
        // Order-6 folded, perturb 4 couplings x1.07, recover via the CLI core.
        let original = FilterDesign::prototype(6, 20.0)
            .synthesize()
            .unwrap()
            .to_folded()
            .unwrap();
        let g = FrequencyGrid::linspace(-2.0, 2.0, 41).unwrap();
        let target = ResponseSolver.evaluate_normalized(&original, &g).unwrap();
        let perturbed = perturb(&original, 4, 1.07);

        let config = mfs::optimize::OptimizeConfig {
            max_iterations: 500,
            tolerance: 1e-10,
            lambda: 1e-2,
            ..Default::default()
        };
        let result = mfs::optimize::optimize_matrix(&perturbed, &target, &g, &config).unwrap();
        assert!(result.converged);
        assert!(worst_rel(&original, &result.matrix) < 1e-4);

        // The JSON round-trip through parse_matrix_json must preserve the
        // perturbed matrix so a real CLI invocation sees the same start.
        let side = perturbed.side();
        let nested: Vec<Vec<f64>> = (0..side)
            .map(|r| (0..side).map(|c| perturbed.at(r, c).unwrap_or(0.0)).collect())
            .collect();
        let json = format!(
            r#"{{"order": 6, "topology": "folded", "data": {}}}"#,
            serde_json::to_string(&nested).unwrap()
        );
        let parsed = parse_matrix_json(&json).unwrap();
        assert_eq!(parsed.topology(), mfs::matrix::MatrixTopology::Folded);
        assert!(worst_rel(&perturbed, &parsed) < 1e-12);

    }

    #[test]
    fn extract_subcommand_recovers_transversal() {
        // Order-6 transversal target; extract from a crude guess via the CLI core.
        let original = FilterDesign::prototype(6, 20.0).synthesize().unwrap().matrix().clone();
        let g = FrequencyGrid::linspace(-2.0, 2.0, 41).unwrap();
        let target = ResponseSolver.evaluate_normalized(&original, &g).unwrap();
        let result = mfs::optimize::extract_matrix(&target, &g, 6).unwrap();
        assert!(result.converged);

        let rec = ResponseSolver.evaluate_normalized(&result.matrix, &g).unwrap();
        let mut max_dev = 0.0_f64;
        for (a, b) in rec.samples.iter().zip(target.samples.iter()) {
            max_dev = max_dev
                .max((a.s11_re - b.s11_re).abs())
                .max((a.s11_im - b.s11_im).abs())
                .max((a.s21_re - b.s21_re).abs())
                .max((a.s21_im - b.s21_im).abs());
        }
        assert!(max_dev < 1e-3);
    }
}
