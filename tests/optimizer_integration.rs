//! Integration coverage for the LM optimizer (`optimize::tune_matrix` /
//! `extract_matrix`) — the RED→GREEN gate for section D of
//! `docs/continuation-plan.md` ("Optimization — from 'works at order 3' to a
//! supported feature").
//!
//! The plan's D1-opt contract: **synthesize → perturb → extract**, then assert
//! the recovered couplings are within 1e-4 *relative* of the originals. The
//! existing lib unit tests are all order 3–4 and well-posed; this file adds the
//! order-6 integration evidence the plan calls for, and it is the gate that
//! D2-opt (the scale-aware Jacobian step) is meant to turn green.
//!
//! Two scenarios:
//!
//! 1. **Good-guess recovery (the core gate).** Synthesize an order-6 folded
//!    matrix, record its response, perturb a few inter-resonator couplings by
//!    5–10 %, and re-optimize from the perturbed matrix against the recorded
//!    response. The perturbed matrix is a *good* initial guess (the real filter
//!    plus small noise — the realistic extraction scenario), and the solution
//!    is locally unique, so each recovered coupling must land within 1e-4
//!    relative of the original. This is what the fixed-δ Jacobian currently
//!    struggles with: a 1e-7 FD step is only scale-safe when couplings are O(1).
//!
//! 2. **Crude-guess extraction.** Feed `extract_matrix` the order-6 response
//!    with no help (it builds its own crude uniform transversal guess) and
//!    require the recovered matrix to reproduce the target response. This pins
//!    the `extract_matrix` entry point end-to-end, which has no integration
//!    coverage today.

use mfs::design::FilterDesign;
use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::optimize::{extract_matrix, optimize_matrix, OptimizeConfig};
use mfs::response::ResponseSolver;

const ORDER: usize = 6;
const REL_TOL: f64 = 1e-4;

/// A 41-point normalized sweep, the same band the optimizer's own lib tests use
/// (the `-2..2` prototype domain).
fn grid() -> FrequencyGrid {
    FrequencyGrid::linspace(-2.0, 2.0, 41).expect("grid")
}
fn synthesize_order6(topology: mfs::matrix::MatrixTopology) -> Result<CouplingMatrix> {
    let design = FilterDesign::prototype(ORDER, 20.0).synthesize()?;
    if topology == mfs::matrix::MatrixTopology::Folded {
        design.to_folded()
    } else {
        Ok(design.matrix().clone())
    }
}

/// Build a perturbed copy of `matrix`: the first `count` non-zero
/// inter-resonator couplings (upper triangle, off-diagonal) are scaled by
/// `factor` (5–10 %), everything else is left untouched.
fn perturb_inter_resonator(matrix: &CouplingMatrix, count: usize, factor: f64) -> CouplingMatrix {
    let side = matrix.side();
    let order = matrix.order();
    let mut data = matrix.as_slice().to_vec();
    let mut done = 0usize;
    for row in 1..=order {
        for col in (row + 1)..=order {
            let idx = row * side + col;
            if data[idx].abs() > 1e-9 {
                data[idx] *= factor;
                data[col * side + row] *= factor; // keep symmetry
                done += 1;
                if done == count {
                    return CouplingMatrix::new_with_topology(order, matrix.topology(), data)
                        .expect("rebuild perturbed matrix");
                }
            }
        }
    }
    CouplingMatrix::new_with_topology(order, matrix.topology(), data).expect("rebuild")
}

/// Assert every significant entry of `original` is reproduced by `recovered`
/// within `REL_TOL` relative. Couplings (off-diagonal) are the meaningful
/// contract; zero entries are skipped on both sides.
fn assert_relative_recovery(original: &CouplingMatrix, recovered: &CouplingMatrix) {
    let side = original.side();
    let mut checked = 0usize;
    for r in 0..side {
        for c in 0..side {
            let o = original.at(r, c).unwrap_or(0.0);
            let rv = recovered.at(r, c).unwrap_or(0.0);
            if o.abs() > 1e-9 {
                let rel = (rv - o).abs() / o.abs();
                assert!(
                    rel < REL_TOL,
                    "entry ({r},{c}): original={o:.6e} recovered={rv:.6e} rel={rel:.3e} \
                     (>= {REL_TOL:.1e})"
                );
                checked += 1;
            } else {
                assert!(
                    rv.abs() < 1e-6,
                    "entry ({r},{c}): original≈0 but recovered={rv:.3e}"
                );
            }
        }
    }
    assert!(checked >= 8, "expected several couplings to check, got {checked}");
}

/// Scenario 1: synthesize order-6 folded → perturb a few couplings 7 % →
/// re-optimize from the perturbed (good-guess) matrix against the recorded
/// response. The recovered couplings must be within 1e-4 relative of the
/// originals. This is the core D1-opt RED→GREEN gate.
///
/// **Status (2026-09-29): RED.** The fixed `δ = 1e-7` finite-difference
/// Jacobian is not scale-safe for order-6 couplings (some entries are O(0.1)),
/// so the optimizer converges to ~1e-2 relative instead of 1e-4. The
/// scale-aware step (section D, D2-opt) is the fix; un-ignore this test when
/// D2-opt lands.
#[test]
#[ignore = "D1-opt RED gate: fixed-δ Jacobian not scale-safe at order 6 (D2-opt is the fix)"]
fn optimizer_recovers_perturbed_order6_folded_from_good_guess() -> Result<()> {
    let original = synthesize_order6(mfs::matrix::MatrixTopology::Folded)?;
    let g = grid();
    let target = ResponseSolver.evaluate_normalized(&original, &g)?;

    let perturbed = perturb_inter_resonator(&original, 4, 1.07);
    // Sanity: the perturbation actually moved the matrix.
    assert_ne!(perturbed.as_slice(), original.as_slice(), "perturbation is a no-op");

    let config = OptimizeConfig {
        max_iterations: 500,
        tolerance: 1e-10,
        lambda: 1e-2,
        ..Default::default()
    };
    let result = optimize_matrix(&perturbed, &target, &g, &config)?;

    assert!(
        result.converged,
        "optimizer did not report convergence: iterations={} cost={:.3e}",
        result.iterations,
        result.cost
    );
    assert_relative_recovery(&original, &result.matrix);
    Ok(())
}

/// Scenario 2: feed `extract_matrix` the order-6 response with no help (crude
/// uniform transversal guess) and require the recovered matrix to reproduce the
/// target response to a tight point-wise bound. Pins the `extract_matrix`
/// entry point end-to-end.
///
/// **Status (2026-09-29): RED.** From the crude uniform guess the LM path
/// wanders into singular response matrices ("response matrix became singular
/// during solve") before recovering the response. Same root-cause family as
/// scenario 1 (fixed-δ Jacobian + crude starting point); D2-opt is the fix.
/// Un-ignore when D2-opt lands.
#[test]
#[ignore = "D1-opt RED gate: extract_matrix from crude guess hits a singular solve (D2-opt is the fix)"]
fn optimizer_extracts_order6_from_crude_guess() -> Result<()> {
    let original = synthesize_order6(mfs::matrix::MatrixTopology::Transversal)?;
    let g = grid();
    let target = ResponseSolver.evaluate_normalized(&original, &g)?;

    let result = extract_matrix(&target, &g, ORDER)?;

    // The recovered matrix must reproduce the target response, not merely
    // return some matrix. Compare all four S-parameter channels point-wise.
    let recovered = ResponseSolver.evaluate_normalized(&result.matrix, &g)?;
    let mut max_dev = 0.0f64;
    for (a, b) in recovered.samples.iter().zip(target.samples.iter()) {
        max_dev = max_dev
            .max((a.s11_re - b.s11_re).abs())
            .max((a.s11_im - b.s11_im).abs())
            .max((a.s21_re - b.s21_re).abs())
            .max((a.s21_im - b.s21_im).abs());
    }
    assert!(
        max_dev < 1e-3,
        "extracted matrix does not reproduce the target response: max|ΔS|={max_dev:.3e}, \
         converged={} iterations={} cost={:.3e}",
        result.converged,
        result.iterations,
        result.cost
    );
    Ok(())
}
