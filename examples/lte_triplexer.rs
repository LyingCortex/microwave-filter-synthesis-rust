//! LTE triplexer: B1 + B3 + B7.
//!
//! Synthesizes an N = 3 multiplexer over three real 3GPP bands:
//!
//! | band | downlink      | width |
//! |------|---------------|-------|
//! | B3   | 1805–1890 MHz | 85 MHz |
//! | B1   | 1920–2170 MHz | 250 MHz |
//! | B7   | 2500–2690 MHz | 190 MHz |
//!
//! B3/B1 have a 30 MHz guard band; B1/B7 are well separated.
//! Per-band order: B3 → 4, B1 → 7, B7 → 5 (all 20 dB RL).
 //! One finite transmission zero each on B3 (at 1925 MHz, the B1 band's
 //! lower edge) and on B1 (at 1885 MHz, the B3 band's upper edge)
 //! suppresses the 30-MHz-guard-band cross leakage. Each is a |omega| > 1
 //! stopband zero in the channel's own normalized prototype domain,
 //! entering the common-domain numerator. Zeros at the *neighboring*
 //! band edge are the worst case for the full-path spectral
 //! factorization (a root cluster of the junction R polynomial sits
 //! within ΔΩ ≈ 0.01 of the jω axis); they are handled by the Newton
 //! polish in `companion_spectral_factor`.
//! Runs the full Macchiarella–Tamiazzo spectral path
//! (`full_iteration(true)`, ideal Type-I junction).
//!
//! Observables (per the design docs):
//! - the *junction-level* rational response (`junction_response` on
//!   `junction_common`) is the physical (N+1)-port of the full path:
//!   0 dB passband at each band center, antenna RL pinned near the
//!   outer-edge spec, lossless on the passband axis;
//! - the *recovered channel matrices* are the per-channel coupling
//!   filters (each keeps a real standalone passband); composing them at
//!   the common node (`Multiplexer::response`) is an admittance model
//!   that is exact for the simplified path, so it is reported for
//!   isolation/unitarity but not asserted to match the junction 0 dB.
//!
//! Also exports, next to the binary (cwd-relative):
//! - `lte_triplexer.s4p` — the 4-port Touchstone (port 0 = antenna,
//!   ports 1..3 = B3/B1/B7), built from the recovered-matrix S-matrix;
//! - `lte_triplexer_curves.svg` — a 3-panel plot (antenna S11, per-band
//!   S21, cross-band isolation) drawn from the junction-level response,
//!   with the three passbands shaded.

use mfs::prelude::*;
use std::fs::File;
use std::io::Write;

/// Builds a 3-panel SVG from the junction-level response: (1) antenna S11,
/// (2) per-band S21, (3) cross-band isolation, with passbands shaded.
fn build_svg(
    resp: &mfs::synthesis::JunctionResponse,
    bands: &[(&str, usize, f64, f64, &[f64])],
    f_lo: f64,
    f_hi: f64,
) -> String {
    const W: f64 = 1000.0;
    const H: f64 = 760.0;
    const ML: f64 = 70.0; // left margin
    const MR: f64 = 130.0; // right margin (legend)
    const MT: f64 = 40.0;
    const PG: f64 = 45.0; // gap between panels
    let plot_w = W - ML - MR;
    let panel_h = (H - MT - 50.0 - 2.0 * PG) / 3.0;

    let y_db = |db: f64, top: f64| -> f64 {
        // y range: 0 dB at top, -60 dB at bottom, clamped
        let c = db.clamp(-60.0, 0.0);
        top + (-c / 60.0) * panel_h
    };
    let x_hz = |f: f64| ML + (f - f_lo) / (f_hi - f_lo) * plot_w;
    // downsample for the SVG: ~2000 points is plenty for a 1000 px plot
    let stride = (resp.frequencies.len() / 2000).max(1);
    let mut s = String::with_capacity(1 << 20);
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {W:.0} {H:.0}\" font-family=\"sans-serif\">\n"
    ));
    s.push_str(&format!(
        "<text x=\"{:.0}\" y=\"24\" font-size=\"17\" font-weight=\"bold\">LTE Triplexer B1/B3/B7 — junction (4-port) response</text>\n",
        ML
    ));

    let colors = ["#e6194b", "#3cb44b", "#4363d8"];

    // Panel titles
    let titles = ["Antenna reflection S11 (dB)", "Channel transmission S21 (dB)", "Cross-band leakage: antenna→port i at other bands' passband (dB)"];

    for panel in 0..3 {
        let top = MT + panel as f64 * (panel_h + PG);
        // frame
        s.push_str(&format!(
            "<rect x=\"{ML:.0}\" y=\"{top:.0}\" width=\"{plot_w:.0}\" height=\"{panel_h:.0}\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1\" opacity=\"0.5\"/>\n"
        ));
        s.push_str(&format!(
            "<text x=\"{ML:.0}\" y=\"{:.0}\" font-size=\"13\" font-weight=\"bold\">{}</text>\n",
            top + 16.0,
            titles[panel]
        ));
        // y grid lines every 10 dB
        for db in 0..=6 {
            let val = -(db as f64 * 10.0);
            let y = y_db(val, top);
            s.push_str(&format!(
                "<line x1=\"{ML:.0}\" y1=\"{y:.1}\" x2=\"{:.0}\" y2=\"{y:.1}\" stroke=\"currentColor\" stroke-width=\"0.5\" opacity=\"0.15\"/>\n",
                ML + plot_w
            ));
            s.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"end\">{:.0}</text>\n",
                ML - 6.0,
                y + 3.5,
                val
            ));
        }
        // x tick labels (GHz)
        for ghz in (1800u32..=2800).step_by(50) {
            let f = ghz as f64 * 1e9;
            if f < f_lo || f > f_hi {
                continue;
            }
            let x = x_hz(f);
            s.push_str(&format!(
                "<line x1=\"{x:.1}\" y1=\"{:.0}\" x2=\"{x:.1}\" y2=\"{:.0}\" stroke=\"currentColor\" stroke-width=\"0.5\" opacity=\"0.15\"/>\n",
                top,
                top + panel_h
            ));
            if panel == 2 || ghz % 100 == 0 {
                s.push_str(&format!(
                    "<text x=\"{x:.1}\" y=\"{:.0}\" font-size=\"10\" text-anchor=\"middle\">{:.2}</text>\n",
                    top + panel_h + 14.0,
                    ghz as f64 / 1000.0
                ));
            }
        }

        // shaded passbands (panels 0 & 1)
        if panel != 2 {
            for (i, (_, _, lo, hi, _)) in bands.iter().enumerate() {
                let x1 = x_hz(*lo * 1e6);
                let x2 = x_hz(*hi * 1e6);
                s.push_str(&format!(
                    "<rect x=\"{x1:.1}\" y=\"{top:.0}\" width=\"{:.1}\" height=\"{panel_h:.0}\" fill=\"{}\" opacity=\"0.07\"/>\n",
                    x2 - x1,
                    colors[i]
                ));
            }
        }

        // data paths
        let mut paths: Vec<String> = Vec::new();
        match panel {
            0 => {
                let last = resp.frequencies.len() - 1;
                let d = resp
                    .antenna_s11_db
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| *k % stride == 0 || *k == last)
                    .map(|(k, &db)| {
                        format!("{:.1},{:.1}", x_hz(resp.frequencies[k]), y_db(db, top))
                    })
                    .collect::<Vec<_>>()
                    .join(" ");
                paths.push(format!(
                    "<polyline points=\"{d}\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.8\" opacity=\"0.9\"/>"
                ));
            }
            1 | 2 => {
                if panel == 1 {
                    // per-band transmission S21 over the full sweep
                    let last = resp.frequencies.len() - 1;
                    for i in 0..bands.len() {
                        let d = resp
                            .channel_s21_db[i]
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| *k % stride == 0 || *k == last)
                            .map(|(k, &db)| {
                                format!("{:.1},{:.1}", x_hz(resp.frequencies[k]), y_db(db, top))
                            })
                            .collect::<Vec<_>>()
                            .join(" ");
                        paths.push(format!(
                            "<polyline points=\"{d}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.6\" opacity=\"0.95\"/>",
                            colors[i]
                        ));
                    }
                } else {
                    // cross-band isolation: leakage into port i at frequencies
                    // where another band j is in its passband (segments of
                    // channel_s21_db[i] clipped to each foreign band)
                    for i in 0..bands.len() {
                    for (j, (_, _, lo, hi, _)) in bands.iter().enumerate() {
                            if i == j {
                                continue;
                            }
                            let (lo, hi) = (*lo * 1e6, *hi * 1e6);
                            let d: Vec<String> = (0..resp.frequencies.len())
                                .step_by(stride)
                                .filter(|&k| {
                                    let f = resp.frequencies[k];
                                    f >= lo && f <= hi
                                })
                                .map(|k| {
                                    format!("{:.1},{:.1}", x_hz(resp.frequencies[k]), y_db(resp.channel_s21_db[i][k], top))
                                })
                                .collect();
                            if d.len() > 1 {
                                paths.push(format!(
                                    "<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.6\" opacity=\"0.95\" stroke-dasharray=\"4 2\"/>",
                                    d.join(" "),
                                    colors[i]
                                ));
                            }
                        }
                    }
                }
            }
            _ => {}
        }
        s.push_str(&paths.join("\n"));
        s.push('\n');
    }

    // legend (right side): one row per band, then a gray row for S11
    let legend_x = ML + plot_w + 20.0;
    let mut ly = MT + 2.0 * (panel_h + PG) + 14.0;
    let entries: Vec<(&str, &str)> = bands
        .iter()
        .enumerate()
        .map(|(i, (n, _, _, _, _))| (*n, colors[i]))
        .chain(std::iter::once(("Antenna S11", "currentColor")))
        .collect();
    for (label, color) in entries {
        s.push_str(&format!(
            "<line x1=\"{legend_x:.0}\" y1=\"{ly:.1}\" x2=\"{:.0}\" y2=\"{ly:.1}\" stroke=\"{color}\" stroke-width=\"2.5\"/>\n",
            legend_x + 22.0
        ));
        s.push_str(&format!(
            "<text x=\"{:.0}\" y=\"{:.1}\" font-size=\"11\">{label}</text>\n",
            legend_x + 28.0,
            ly + 3.5
        ));
        ly += 20.0;
    }
    s.push_str("</svg>\n");
    s
}
fn band(order: usize, rl: f64, name: &str, lo_mhz: f64, hi_mhz: f64, tz_mhz: &[f64]) -> ChannelSpec {
    let (center, bw) = mfs::freq::edges_to_center_bandwidth(lo_mhz * 1e6, hi_mhz * 1e6);
    // Finite transmission zeros in the channel's own normalized prototype domain
    // (stopband => |omega| > 1). In the full path these enter the common-domain
    // numerator P_i and the generalized reflection F polynomial.
    let mapping = mfs::freq::BandPassMapping::new(center, bw).expect("channel mapping");
    let zeros = tz_mhz
        .iter()
        .map(|f| mapping.map_hz_to_normalized(*f * 1e6).expect("tz mapping").omega)
        .collect::<Vec<f64>>();
    println!(
        "  {name}: {lo_mhz:.0}-{hi_mhz:.0} MHz (center {:.4} GHz, bw {:.1} MHz, order {order}, RL {rl:.0} dB, TZ {tz_mhz:?} MHz)",
        center / 1e9,
        bw / 1e6
    );
    ChannelSpec { order, return_loss_db: rl, zeros, center_hz: center, bandwidth_hz: bw }
}

fn closest(freqs: &[f64], f: f64) -> usize {
    let mut best = 0;
    for (i, &x) in freqs.iter().enumerate().skip(1) {
        if (x - f).abs() < (freqs[best] - f).abs() {
            best = i;
        }
    }
    best
}

fn db(mag: f64) -> f64 {
    if mag > 1e-15 { 20.0 * mag.log10() } else { f64::NEG_INFINITY }
}

fn main() -> mfs::error::Result<()> {
    let bands: [(&str, usize, f64, f64, &[f64]); 3] =
         [("B3", 4, 1805.0, 1890.0, &[1925.0]),
          ("B1", 7, 1920.0, 2170.0, &[1885.0]),
          ("B7", 5, 2500.0, 2690.0, &[])];
    let specs: Vec<ChannelSpec> = bands
        .iter()
        .map(|(n, ord, lo, hi, tz)| band(*ord, 20.0, n, *lo, *hi, *tz))
        .collect();

    // ── Synthesis: full spectral path ────────────────────────────────────
    println!("\nSynthesizing (full_iteration = true, Type-I junction) ...");
    let t = std::time::Instant::now();
    let result = MultiplexerSynthesis::new()
        .channels(specs.clone())
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(80)
        .convergence_tol(1e-8)
        .full_iteration(true)
        .synthesize()?;
    let dt = t.elapsed();

    println!(
        "converged={}  iterations={}  conv_err={:.3e}  ({:.2?})",
        result.converged,
        result.iterations,
        result.convergence_error,
        dt
    );
    assert!(
        result.converged,
        "full path did not stop on its own criterion (conv_err {:.3e})",
        result.convergence_error
    );

    // ── Junction-level response (the physical observable, full path) ────
    let (f_lo, f_hi) = (1750e6, 2750e6);
    let common = result.junction_common.as_ref().expect("full path must capture junction state");
    let jresp = mfs::synthesis::junction_response(common, 3, f_lo, f_hi, 20001)?;

    println!("\nJunction-level (4-port) response, dB:");
    println!("  {:<4} {:>12} {:>12} {:>12} {:>8}", "band", "S21@center", "S11@lo-edge", "S11@hi-edge", "spec RL");
    for (i, (name, _, lo, hi, _)) in bands.iter().enumerate() {
        let ci = closest(&jresp.frequencies, specs[i].center_hz);
        let li = closest(&jresp.frequencies, *lo * 1e6);
        let hi_ = closest(&jresp.frequencies, *hi * 1e6);
        let s21 = jresp.channel_s21_db[i][ci];
        println!(
            "  {name:<4} {:>12.3} {:>12.3} {:>12.3} {:>8.1}",
            s21,
            jresp.antenna_s11_db[li],
            jresp.antenna_s11_db[hi_],
            specs[i].return_loss_db
        );
        assert!(
            s21 > -0.5,
            "{name} junction passband S21 {s21:.3} dB (want > -0.5)"
        );
    }
    assert!(
        jresp.antenna_s11_db.iter().all(|&x| x.is_finite()),
        "antenna S11 non-finite"
    );

    // Junction-level cross-band isolation: worst (max) channel-i S21 while
    // the antenna is driven inside channel j's passband — the physical
    // inter-channel isolation of the (N+1)-port.
    println!("\nJunction-level cross-band isolation (worst within the other band), dB:");
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                continue;
            }
            let (lo, hi) = (bands[j].2 * 1e6, bands[j].3 * 1e6);
            let worst = jresp
                .frequencies
                .iter()
                .enumerate()
                .filter(|&(_, &f)| f >= lo && f <= hi)
                .map(|(k, _)| jresp.channel_s21_db[i][k])
                .fold(f64::NEG_INFINITY, |a, b| a.max(b));
            println!(
                "  {} driven in {} band: {worst:.3} {}",
                bands[i].0,
                bands[j].0,
                if worst < -20.0 { "OK" } else { "⚠" }
            );
        }
    }

    // ── Recovered channel matrices (per-channel filters) ────────────────
    let mux = Multiplexer::new(result.channels)?;

    println!("\nRecovered channel filters (standalone peak S21):");
    for (i, (name, _, _, _, _)) in bands.iter().enumerate() {
        let ch = mux.channel(i);
        let center = ch.center_hz().expect("center present");
        let bw = ch.bandwidth_hz().expect("bandwidth present");
        let resp = ch.response_bandpass(center, bw, center - bw, center + bw, 101)?;
        let peak = resp.samples.iter().map(|s| s.s21_mag()).fold(0.0_f64, f64::max);
        let peak_db = db(peak);
        println!("  {name}: peak S21 {peak_db:.3} dB, order {}", ch.order());
        assert!(
            peak_db > -10.0,
            "{name} recovered channel degenerate (standalone peak {peak_db} dB)"
        );
    }

    // Inter-channel isolation from the recovered matrices (admittance
    // model): worst (max) |S| over the full sweep.
    let resp = mux.response(f_lo, f_hi, 2001)?;
    println!("\nInter-channel isolation (recovered-matrix admittance model, worst over 1.75–2.75 GHz), dB:");
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                continue;
            }
            let worst = resp.isolation_db[i][j].iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            println!(
                "  S{}{} (drive {}, measure {}): {worst:.3} {}",
                i + 2,
                j + 2,
                bands[j].0,
                bands[i].0,
                if worst < -40.0 { "OK" } else { "⚠" }
            );
        }
    }

    // Antenna-row power conservation of the composed (N+1)-port.
    let smat = mux.response_matrix(f_lo, f_hi, 401)?;
    let p = smat.num_ports;
    let mut max_dev = 0.0f64;
    for row in &smat.s {
        let sum = (0..p).map(|j| row[j * p].norm_sqr()).sum::<f64>();
        max_dev = max_dev.max((sum - 1.0).abs());
    }
    println!("\nAntenna-row power conservation: max|1 − (|S00|²+Σ|S0i|²)| = {max_dev:.3e}");
    assert!(max_dev < 1e-8, "antenna row not power-conserving");

    // ── Export: 4-port Touchstone (.s4p) from the recovered matrices ────
    let smat_full = mux.response_matrix(f_lo, f_hi, 1601)?;
    let p = smat_full.num_ports;
    let mut builder = touchstone::NetworkBuilder::new("lte_triplexer.s4p", p)
        .frequency_unit("MHz")
         .comment("ports: 0=antenna 1=B3 2=B1 3=B7 (recovered-matrix node model); TZ: B3@1925MHz B1@1885MHz");
    for (q, &f) in smat_full.frequencies.iter().enumerate() {
        let data = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        let c = smat_full.s[q][i * p + j];
                        touchstone::Complex { re: c.re, im: c.im }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        builder = builder.point(f, touchstone::SMatrix { rank: p, data });
    }
    let net = builder.build().expect("build 4-port network");
    let s4p_path = "lte_triplexer.s4p";
    net.save(s4p_path).expect("write .s4p");
    println!("\nWrote {s4p_path} ({p}-port, {} freq points).", smat_full.frequencies.len());

    // ── Export: SVG curves from the junction-level response ─────────────
    let svg = build_svg(&jresp, &bands, f_lo, f_hi);
    let svg_path = "lte_triplexer_curves.svg";
    File::create(svg_path)
        .expect("open svg")
        .write_all(svg.as_bytes())
        .expect("write svg");
    println!("Wrote {svg_path}.");

    println!("\nB1/B3/B7 triplexer synthesized successfully.");
    Ok(())
}
