//! §D.5: the `multiplexer` subcommand: N-channel synthesis + analysis.
//!
//! Pure move out of the monolithic `src/bin/mfs_cli.rs` — no behavior
//! change; `main::run` dispatches to `run_multiplexer`.
//!
use crate::io::{CliRunError, OutputFormat};
/// Runs the quick `multiplexer` subcommand: N-channel synthesis + analysis.
pub(crate) fn run_multiplexer(
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

/// Chooses a response frequency range, defaulting to the outer channel bands.
///
/// The default window spans from `0.5·bw` below the outermost passband lower
/// edge to `0.5·bw` above the outermost passband upper edge, where the edges
/// are the geometric-mean `band_edge_pair` values (review P2-3, matching
/// `FilterDesign::default_response`). This mirrors the former `center ∓
/// 1.5·bw` window for narrowband designs but is derived from the *edges*, so a
/// wideband channel (fbw > 2/3, where `center − 1.5·bw < 0`) still yields a
/// positive window instead of one that fails in `map_hz_to_normalized`.
pub(crate) fn resolve_multiplexer_range(
    channels: &[mfs::design::FilterDesign],
    resp_start: Option<f64>,
    resp_stop: Option<f64>,
) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for ch in channels.iter() {
        if let (Some(center), Some(bw)) = (ch.center_hz(), ch.bandwidth_hz()) {
            let (lower, upper) = mfs::freq::band_edge_pair(center, bw);
            let pad = 0.5 * bw;
            // Each channel's default window (per `default_response`): passband
            // edges padded by half a bandwidth, floored so wideband channels
            // (fbw > 2/3) stay positive.
            lo = lo.min((lower - pad).max(1e-9 * center));
            hi = hi.max(upper + pad);
        }
    }
    (resp_start.unwrap_or(lo), resp_stop.unwrap_or(hi))
}

/// §D.7: min over a *raw* precomputed dB array (`Multiplexer::response`
/// `antenna_s11_db` / `channel_s21_db` are magnitudes already in dB, not
/// `ResponseSample`s) — no `ResponseSample::s*_db()` to route through.
pub(crate) fn min_db(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

/// §D.7: max over a raw precomputed dB array (see `min_db`).
pub(crate) fn max_db(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}


#[cfg(test)]
mod multiplexer_tests {
    //! §D.5: moved out of `optimizer_tests` with its code.
    use super::*;
    use mfs::design::FilterDesign;

    /// The default multiplexer response window is derived from the passband
    /// *edges* (geometric-mean `band_edge_pair`) padded by 0.5·bw, floored at
    /// `1e-9·center` — not the old arithmetic `center ∓ 1.5·bw`, which goes
    /// negative for wideband channels (fbw > 2/3) and fails in
    /// `map_hz_to_normalized` (review round-2, matching P2-3).
    #[test]
    fn multiplexer_default_range_is_edge_based_and_positive() {
        // Two bandpass channels: the default window is the geometric-mean
        // `band_edge_pair` edges padded by 0.5·bw on each side (the same
        // formula as `FilterDesign::default_response`).
        let ch0 = FilterDesign::bandpass(2, 20.0, 6.75e9, 0.3e9).synthesize().unwrap();
        let ch1 = FilterDesign::bandpass(2, 20.0, 6.95e9, 0.3e9).synthesize().unwrap();
        let (lo, hi) = resolve_multiplexer_range(&[ch0.clone(), ch1.clone()], None, None);

        // The outer edges come from the outermost channels (ch0's lower edge,
        // ch1's upper edge); pad is each channel's own 0.5·bw.
        let (l0, _u0) = mfs::freq::band_edge_pair(
            ch0.center_hz().unwrap(),
            ch0.bandwidth_hz().unwrap(),
        );
        let (_l1, u1) = mfs::freq::band_edge_pair(
            ch1.center_hz().unwrap(),
            ch1.bandwidth_hz().unwrap(),
        );
        let expect_lo = l0 - 0.5 * ch0.bandwidth_hz().unwrap();
        let expect_hi = u1 + 0.5 * ch1.bandwidth_hz().unwrap();
        assert!((lo - expect_lo).abs() < 1e-6 * lo, "lo {lo} vs {expect_lo}");
        assert!((hi - expect_hi).abs() < 1e-6 * hi, "hi {hi} vs {expect_hi}");
        assert!(lo > 0.0);

        // Wideband channel (fbw = 1.2 > 2/3): the old arithmetic window
        // `center - 1.5·bw` is negative and would fail in
        // `map_hz_to_normalized`; the edge-based floor keeps the start positive.
        let wide = FilterDesign::bandpass(2, 20.0, 1.0e9, 1.2e9).synthesize().unwrap();
        let (wlo, _whi) = resolve_multiplexer_range(std::slice::from_ref(&wide), None, None);
        let wc = wide.center_hz().unwrap();
        let wbw = wide.bandwidth_hz().unwrap();
        let (wlower, _wupper) = mfs::freq::band_edge_pair(wc, wbw);
        assert!(
            wlower - 0.5 * wbw < 1e-9 * wc,
            "wideband lower edge minus pad must need the positive floor"
        );
        assert!(wlo > 0.0, "wideband window start must stay positive, got {wlo}");
        assert!((wlo - 1e-9 * wc).abs() < 1e-3 * wc, "lo {wlo} vs floor {}", 1e-9 * wc);
    }
}
