//! Verifies MFS synthesis results against Cameron literature reference data.

use mfs::prelude::*;
use mfs::approx::generalized_chebyshev_polynomials;
use mfs::spec::FilterSpec;

fn approx_eq(a: f64, b: f64, tol: f64, label: &str) {
    let diff = (a - b).abs();
    if diff > tol {
        eprintln!("  MISMATCH {label}: got {a:.6}, expected {b:.6}, diff={diff:.2e}");
    }
}

fn main() -> mfs::Result<()> {
    println!("=== Literature Verification ===\n");

    // ─── Cameron 1999/2011: 4-2 asymmetric, 22 dB RL, TZs at +j1.3217, +j1.8082 ───
    println!("--- Case: 4-2 Asymmetric (Cameron 1999 Table I / 2011 Table 1) ---");
    let spec = FilterSpec::new(4, 22.0)?
        .with_normalized_transmission_zeros(vec![1.3217, 1.8082]);
    let polys = generalized_chebyshev_polynomials(&spec)?;

    println!("  Order: {}", polys.order);
    println!("  Epsilon: {:.4} (expected: 1.1548)", polys.eps);
    println!("  Epsilon_R: {:.4} (expected: 1.0)", polys.eps_r);
    approx_eq(polys.eps, 1.1548, 0.01, "epsilon");
    approx_eq(polys.eps_r, 1.0, 0.01, "epsilon_r");

    // Check E(s) polynomial coefficients
    println!("\n  E(s) coefficients (ascending power):");
    let e_ref = [
        (-0.1268, -2.0658),
        (2.4874, -3.6255),
        (3.6706, -2.1950),
        (2.4015, -0.7591),
        (1.0, 0.0),
    ];
    for (i, coeff) in polys.e.coefficients.iter().enumerate() {
        if i < e_ref.len() {
            let (ref_re, ref_im) = e_ref[i];
            let re_err = (coeff.re - ref_re).abs();
            let im_err = (coeff.im - ref_im).abs();
            let status = if re_err < 0.01 && im_err < 0.01 { "OK" } else { "DIFF" };
            println!("    s^{i}: ({:.4}, {:.4}) ref=({ref_re:.4}, {ref_im:.4}) [{status}]",
                coeff.re, coeff.im);
        }
    }

    // Check F(s) polynomial coefficients
    println!("\n  F(s) coefficients:");
    let f_ref = [
        (0.0208, 0.0),
        (0.0, -0.5432),
        (0.7869, 0.0),
        (0.0, -0.7591),
        (1.0, 0.0),
    ];
    for (i, coeff) in polys.f.coefficients.iter().enumerate() {
        if i < f_ref.len() {
            let (ref_re, ref_im) = f_ref[i];
            let re_err = (coeff.re - ref_re).abs();
            let im_err = (coeff.im - ref_im).abs();
            let status = if re_err < 0.01 && im_err < 0.01 { "OK" } else { "DIFF" };
            println!("    s^{i}: ({:.4}, {:.4}) ref=({ref_re:.4}, {ref_im:.4}) [{status}]",
                coeff.re, coeff.im);
        }
    }

    // Synthesize matrix and check folded topology
    let design = FilterDesign::from_spec(&spec)?;
    let folded = design.to_folded()?;
    println!("\n  Folded matrix source coupling: {:.4}", folded.source_coupling());
    println!("  Folded matrix load coupling: {:.4}", folded.load_coupling());
    println!("  Chain couplings: {:?}", folded.chain_couplings());

    // Check S-parameter response at key points
    let resp = design.response_normalized(-3.0, 3.0, 61)?;
    // At omega=0 (center), S21 should be near 0 dB (passband)
    let center = &resp.samples[30];
    println!("\n  Response at center (omega=0): S21={:.2} dB, S11={:.2} dB",
        center.s21_db(), center.s11_db());
    // S11 at center should be < -22 dB (return loss)
    if center.s11_db() > -20.0 {
        eprintln!("  WARNING: S11 at center = {:.2} dB, expected < -22 dB", center.s11_db());
    }

    // At omega=1.3217 (first TZ), S21 should be very low
    // Find closest point to omega=1.3217
    let tz1_idx = resp.samples.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (a.normalized_omega - 1.3217).abs()
                .partial_cmp(&(b.normalized_omega - 1.3217).abs())
                .unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    let tz1 = &resp.samples[tz1_idx];
    println!("  Response near TZ1 (omega={:.3}): S21={:.1} dB",
        tz1.normalized_omega, tz1.s21_db());

    // ─── All-pole 4th order, 20 dB (basic sanity check) ───
    println!("\n--- Case: 4th order all-pole, 20 dB RL ---");
    let design2 = FilterDesign::prototype(4, 20.0).synthesize()?;
    let resp2 = design2.response_normalized(-2.0, 2.0, 41)?;
    let c2 = &resp2.samples[20];
    println!("  S21 at center: {:.4} dB (should be ~0)", c2.s21_db());
    println!("  S11 at center: {:.2} dB (should be < -20)", c2.s11_db());
    let power = c2.s11_mag().powi(2) + c2.s21_mag().powi(2);
    println!("  Power conservation: |S11|²+|S21|² = {:.10} (should be 1.0)", power);

    // ─── 6th order with 3 TZs (higher order test) ───
    println!("\n--- Case: 6th order, 23 dB RL, 3 TZs at -2.0, -1.2, 1.5 ---");
    let design3 = FilterDesign::prototype(6, 23.0)
        .zeros([-2.0, -1.2, 1.5])
        .synthesize()?;
    println!("  Order: {}", design3.order());
    println!("  Epsilon: {:.6}", design3.polynomials().eps);
    println!("  Epsilon_R: {:.6}", design3.polynomials().eps_r);

    let folded3 = design3.to_folded()?;
    println!("  Folded source coupling: {:.4}", folded3.source_coupling());
    println!("  Folded load coupling: {:.4}", folded3.load_coupling());
    println!("  Chain couplings: {:?}", folded3.chain_couplings());

    // Verify response
    let resp3 = design3.response_normalized(-3.0, 3.0, 61)?;
    let c3 = &resp3.samples[30];
    println!("  S21 at center: {:.4} dB", c3.s21_db());
    println!("  S11 at center: {:.2} dB (should be < -23)", c3.s11_db());
    let power3 = c3.s11_mag().powi(2) + c3.s21_mag().powi(2);
    println!("  Power conservation: {:.10}", power3);

    // Check TZ locations - S21 should be very low near zeros
    for &tz in &[-2.0, -1.2, 1.5] {
        let idx = resp3.samples.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (a.normalized_omega - tz).abs()
                    .partial_cmp(&(b.normalized_omega - tz).abs())
                    .unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        let s = &resp3.samples[idx];
        println!("  S21 near TZ={tz:.1}: {:.1} dB (omega={:.3})",
            s.s21_db(), s.normalized_omega);
    }

    // ─── 8th order all-pole, 20 dB (high order stability) ───
    println!("\n--- Case: 8th order all-pole, 20 dB RL ---");
    let design4 = FilterDesign::prototype(8, 20.0).synthesize()?;
    let folded4 = design4.to_folded()?;
    let arrow4 = design4.to_arrow()?;
    println!("  Folded matrix order: {}", folded4.order());
    println!("  Arrow matrix order: {}", arrow4.order());

    // Verify folded and arrow give same response
    let grid4 = mfs::FrequencyGrid::linspace(-2.0, 2.0, 21)?;
    let resp_f = mfs::ResponseSolver.evaluate_normalized(&folded4, &grid4)?;
    let resp_a = mfs::ResponseSolver.evaluate_normalized(&arrow4, &grid4)?;
    let mut max_diff = 0.0_f64;
    for (f, a) in resp_f.samples.iter().zip(resp_a.samples.iter()) {
        let diff = (f.s21_mag() - a.s21_mag()).abs();
        max_diff = max_diff.max(diff);
    }
    println!("  Folded vs Arrow max |S21| diff: {:.2e} (should be < 1e-10)", max_diff);

    // ─── 12th order with 6 TZs (stress test) ───
    println!("\n--- Case: 12th order, 20 dB RL, 6 TZs ---");
    let design5 = FilterDesign::prototype(12, 20.0)
        .zeros([-1.3, 1.3, -1.8, 1.8, -2.5, 2.5])
        .synthesize()?;
    println!("  Order: {}", design5.order());
    let resp5 = design5.response_normalized(-3.0, 3.0, 101)?;
    let c5 = &resp5.samples[50];
    let power5 = c5.s11_mag().powi(2) + c5.s21_mag().powi(2);
    println!("  S21 at center: {:.4} dB", c5.s21_db());
    println!("  Power conservation: {:.10}", power5);
    println!("  Folded matrix symmetric: {}", {
        let f5 = design5.to_folded()?;
        let side = f5.side();
        let mut max_asym = 0.0_f64;
        for r in 0..side {
            for c in (r+1)..side {
                let diff = (f5.at(r,c).unwrap() - f5.at(c,r).unwrap()).abs();
                max_asym = max_asym.max(diff);
            }
        }
        format!("{:.2e} (should be < 1e-14)", max_asym)
    });

    println!("\n=== Done ===");
    Ok(())
}
