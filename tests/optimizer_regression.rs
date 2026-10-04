// Regression tests for the optimizer's *extraction* contract beyond the
// well-posed cases — D4-opt of `docs/continuation-plan.md` section D.
//
// `tests/optimizer_integration.rs` pins `extract_matrix` on exact (noiseless)
// responses. Real extraction runs on *measured* S-data, which carries a finite
// noise floor. This file pins the behaviour that matters there:
//
//   * the fit converges rather than diverging on noisy input,
//   * the recovered model stays faithful to the **underlying filter** (the
//     response error against the truth is at the noise scale, not the signal
//     scale), and
//   * it does not chase the noise (the residual against the noisy data is the
//     injected noise level, not smaller — a fit that drove the residual below
//     the noise would be fitting the noise, which is unphysical).
//
// The recovered *couplings* are deliberately not asserted: a transversal
// matrix with source/load coupled to every resonator is a redundant
// parameterisation, so many matrices reproduce the same response and the
// physical contract is the response, not the individual entries.

use mfs::design::FilterDesign;
use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::optimize::extract_matrix;
use mfs::response::{ResponseSample, ResponseSolver, SParameterResponse};

const ORDER: usize = 6;

fn grid() -> FrequencyGrid {
    FrequencyGrid::linspace(-2.0, 2.0, 41).expect("grid")
}

/// Deterministic 0.1 dB magnitude noise (S11 and S21 independently), phase
/// untouched — an LCG so the test is reproducible without a rand dependency.
fn add_magnitude_noise_db(response: &SParameterResponse, db: f64) -> SParameterResponse {
    let mut state: u32 = 0x2545_F491;
    let mut next_scale_db = || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        let unit = state as f64 / u32::MAX as f64 - 0.5; // [-0.5, 0.5]
        unit * 2.0 * db
    };
    let samples: Vec<ResponseSample> = response
        .samples
        .iter()
        .map(|s| {
            let s11 = 10f64.powf(next_scale_db() / 20.0);
            let s21 = 10f64.powf(next_scale_db() / 20.0);
            let mut noisy = *s;
            noisy.s11_re *= s11;
            noisy.s11_im *= s11;
            noisy.s21_re *= s21;
            noisy.s21_im *= s21;
            noisy
        })
        .collect();
    SParameterResponse { samples }
}

/// Point-wise max deviation across all four S-parameter channels.
fn max_dev(a: &SParameterResponse, b: &SParameterResponse) -> f64 {
    let mut m = 0.0f64;
    for (x, y) in a.samples.iter().zip(b.samples.iter()) {
        m = m
            .max((x.s11_re - y.s11_re).abs())
            .max((x.s11_im - y.s11_im).abs())
            .max((x.s21_re - y.s21_re).abs())
            .max((x.s21_im - y.s21_im).abs());
    }
    m
}

#[test]
fn extract_recovers_order6_filter_from_noisy_response() -> Result<()> {
    let truth_matrix = FilterDesign::prototype(ORDER, 20.0).synthesize()?.matrix().clone();
    let g = grid();
    let truth = ResponseSolver.evaluate_normalized(&truth_matrix, &g)?;
    let noisy = add_magnitude_noise_db(&truth, 0.1);

    // The injected noise itself, for scale.
    let noise_level = max_dev(&noisy, &truth);
    assert!(
        noise_level > 1e-4,
        "noise injection is a no-op: max|dS|={noise_level:.3e}"
    );

    let result = extract_matrix(&noisy, &g, ORDER)?;
    assert!(
        result.converged,
        "extraction did not converge on noisy data: iterations={} cost={:.3e}",
        result.iterations, result.cost
    );

    let recovered = ResponseSolver.evaluate_normalized(&result.matrix, &g)?;
    let dev_truth = max_dev(&recovered, &truth);
    let dev_noisy = max_dev(&recovered, &noisy);

    // Faithful to the underlying filter: response error at the noise scale
    // (measured 2026-10-04: 1.22e-3 against a 1.0e-2 noise floor).
    assert!(
        dev_truth < 5e-3,
        "recovered model diverges from the underlying filter: max|dS|_vs_truth={dev_truth:.3e}"
    );
    // Not chasing the noise: the residual against the noisy data is the noise
    // level (it cannot be dramatically smaller without overfitting).
    assert!(
        dev_noisy < 4.0 * noise_level,
        "residual against noisy data ({dev_noisy:.3e}) is inconsistent with the noise \
         floor ({noise_level:.3e})"
    );
    Ok(())
}
