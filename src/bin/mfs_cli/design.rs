//! §D.5: the `design` subcommand (incl. response + Touchstone export).
//!
//! Pure move out of the monolithic `src/bin/mfs_cli.rs` — no behavior
//! change; `main::run` dispatches to `run_design`.
//!
use crate::io::{CliRunError, OutputFormat};
/// Parameters for the `design` subcommand, grouped into a single value so
/// [`run_design`] stays under clippy's argument-count limit.
pub(crate) struct DesignArgs<'a> {
    pub(crate) order: usize,
    pub(crate) return_loss: f64,
    pub(crate) center: Option<f64>,
    pub(crate) bandwidth: Option<f64>,
    pub(crate) start: Option<f64>,
    pub(crate) stop: Option<f64>,
    pub(crate) zeros: Option<&'a [f64]>,
    pub(crate) q: Option<f64>,
    pub(crate) topology: &'a str,
    pub(crate) points: usize,
    pub(crate) resp_start: Option<f64>,
    pub(crate) resp_stop: Option<f64>,
    pub(crate) s2p: Option<&'a std::path::Path>,
}

/// Runs the quick `design` subcommand.
pub(crate) fn run_design(args: &DesignArgs<'_>, format: &OutputFormat) -> Result<(), CliRunError> {
    use mfs::design::FilterDesign;
    use mfs::matrix::MatrixTopology;
    use mfs::touchstone::{self, DataFormat, FreqUnit, TouchstoneConfig, TouchstoneVersion};

    let order = args.order;
    let return_loss = args.return_loss;
    let center = args.center;
    let bandwidth = args.bandwidth;
    let start = args.start;
    let stop = args.stop;
    let zeros = args.zeros;
    let q = args.q;
    let topology = args.topology;
    let points = args.points;
    let resp_start = args.resp_start;
    let resp_stop = args.resp_stop;
    let s2p = args.s2p;

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
        "folded" => MatrixTopology::Folded,
        // P2-13: reject unknown values instead of silently degrading to Folded
        // (the no-silent-degradation policy; `parse_matrix_json` already does
        // this for topology in the JSON path).
        other => {
            return Err(CliRunError::ArgValidation(format!(
                "unknown topology '{other}'; expected 'transversal', 'folded', or 'arrow'"
            )))
        }
    };
    let matrix = design.to_topology(topo).map_err(CliRunError::Pipeline)?;

    // Get response
    let response = if let (Some(c), Some(bw)) = (resolved_center, resolved_bw) {
        // Default response range: the geometric passband edges padded by
        // 0.5·bw — stays positive for wideband designs where the old
        // `c ± 1.5·bw` window went below 0 Hz (P2-3).
        let (lower_edge, upper_edge) = mfs::freq::band_edge_pair(c, bw);
        let r_start = resp_start.unwrap_or((lower_edge - 0.5 * bw).max(1e-9 * c));
        let r_stop = resp_stop.unwrap_or(upper_edge + 0.5 * bw);
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
