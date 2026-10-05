//! C3 example: render the synthesized multi-band (dual-band) response — the
//! companion to `tests/multi_band.rs`. Synthesizes the Cameron 2011 Fig. 25
//! design (route-α, `docs/internal/multi-band-synthesis-design.md`) via
//! `multi_band_chebyshev_polynomials` + `MatrixSynthesisEngine`, evaluates it
//! on a normalized sweep, prints the per-band / gap / stopband structure, and
//! writes a Touchstone .s2p of the full response.
use mfs::approx::multi_band_chebyshev_polynomials;
use mfs::freq::FrequencyGrid;
use mfs::response::ResponseSolver;
use mfs::synthesis::MatrixSynthesisEngine;
use mfs::touchstone::{write_touchstone, DataFormat, FreqUnit, TouchstoneConfig, TouchstoneVersion};
use mfs::MfsError;
use mfs::Result;

/// Cameron 2011 Fig. 25 dual-band symmetric spec (normalized prototype domain).
const BANDS: &[(f64, f64)] = &[(-1.0, -0.35), (0.35, 1.0)];
const ORDER: usize = 10;
const RL: f64 = 22.0;
const TZ: &[f64] = &[0.0, 0.25420871, -0.25420871, 1.15728568, -1.15728568];

fn main() -> Result<()> {
    // Synthesize: approximation (C2b) → engine (C2c) → transversal N+2 matrix.
    let set = multi_band_chebyshev_polynomials(BANDS, ORDER, RL, TZ)?;
    let matrix = MatrixSynthesisEngine {}.synthesize(&set)?;
    println!("matrix shape: {:?}", matrix.shape());
    println!(
        "polynomial degrees: E={} F={} P={}",
        set.e.degree(),
        set.f.degree(),
        set.p.degree()
    );

    // Evaluate on a normalized sweep covering both bands, the gap, and the
    // outer stopbands.
    let grid = FrequencyGrid::linspace(-2.0, 2.0, 401)?;
    let omegas = grid.as_slice().to_vec();
    let response = ResponseSolver.evaluate_normalized(&matrix, &grid)?;

    println!("\n|w|    |S21| dB   |S11| dB   region");
    println!("----  ---------  ----------  -------------------------");
    for w in [-1.0, -0.675, -0.35, -0.25, 0.0, 0.25, 0.35, 0.675, 1.0, 1.2, 1.5, 2.0] {
        let i = omegas
            .iter()
            .enumerate()
            .min_by(|x, y| (x.1 - w).abs().total_cmp(&((y.1 - w).abs())))
            .unwrap()
            .0;
        let region = if w.abs() >= 0.35 && w.abs() <= 1.0 {
            "passband"
        } else if w.abs() < 0.35 {
            "inter-band gap"
        } else {
            "outer stopband"
        };
        println!(
            "{:+4.2}  {:>9.2}  {:>10.2}  {region}",
            omegas[i],
            response.samples[i].s21_db(),
            response.samples[i].s11_db()
        );
    }

    // Summary of the C3 gates (the integration test asserts these).
    let mut dunit = 0.0_f64;
    for s in &response.samples {
        let a = s.s11_mag();
        let b = s.s21_mag();
        dunit = dunit.max((a * a + b * b - 1.0).abs());
    }
    let worst_in_band = [0.35_f64, -1.0_f64]
        .iter()
        .map(|&lo| {
            omegas
                .iter()
                .enumerate()
                .filter(|(_, w)| **w >= lo && **w <= lo + 0.65)
                .map(|(i, _)| response.samples[i].s11_db())
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .fold(f64::NEG_INFINITY, f64::max);
    println!(
        "\nworst in-band |S11| = {worst_in_band:.2} dB (spec -{RL} dB RL)"
    );
    println!("max unitarity residual |S11|²+|S21|²-1| = {dunit:.2e}");
    println!(
        "TZ nulls at w=0, ±0.2542, ±1.1573 (|S21| ≈ 0 at each)"
    );

    // Render the full response as a Touchstone file (normalized axis, 50 Ω).
    let out = std::env::current_dir()
        .map_err(|e| MfsError::PreconditionViolation(format!("cannot read current dir: {e}")))?
        .join("target/dual_band_fig25.s2p");
    let config = TouchstoneConfig {
        freq_unit: FreqUnit::GHz,
        format: DataFormat::DB,
        impedance: 50.0,
        version: TouchstoneVersion::V1,
        comments: vec![
            "Cameron 2011 Fig.25 dual-band symmetric".into(),
            "route-α multi-band Chebyshev (mfs)".into(),
        ],
    };
    write_touchstone(&response, &config, &out).map_err(|e| {
        MfsError::PreconditionViolation(format!("failed to write {}: {e}", out.display()))
    })?;
    println!("touchstone: {}", out.display());
    Ok(())
}
