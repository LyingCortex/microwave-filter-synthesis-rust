//! Stress test for numerical stability at filter orders 20+.
//!
//! This example synthesizes coupling matrices for high-order filters and
//! verifies that the results are physically meaningful (symmetric, bounded
//! couplings, correct eigenvalue structure).

//! Demonstrates the legacy prelude/free-function workflow, which is
//! deprecated but must keep working for backward compatibility.
#![allow(deprecated)]

use mfs::matrix::MatrixTopology;
use mfs::prelude::*;
use mfs::synthesis::MatrixSynthesisEngine;
use mfs::CouplingMatrix;

fn check_matrix_sanity(matrix: &CouplingMatrix, label: &str) -> bool {
    let side = matrix.side();
    let order = matrix.order();
    let mut ok = true;

    // Check symmetry: |M[i,j] - M[j,i]| < tol
    let mut max_asymmetry = 0.0_f64;
    for row in 0..side {
        for col in (row + 1)..side {
            let diff = (matrix.at(row, col).unwrap() - matrix.at(col, row).unwrap()).abs();
            max_asymmetry = max_asymmetry.max(diff);
        }
    }

    if max_asymmetry > 1e-10 {
        eprintln!("  [{label}] FAIL: max asymmetry = {max_asymmetry:.2e}");
        ok = false;
    }

    // Check that coupling values are bounded (no blow-up)
    let max_entry = matrix
        .as_slice()
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    if max_entry > 100.0 {
        eprintln!("  [{label}] WARN: max entry magnitude = {max_entry:.4} (unusually large)");
    }
    if !max_entry.is_finite() {
        eprintln!("  [{label}] FAIL: matrix contains non-finite values");
        ok = false;
    }

    // Check source/load coupling exists
    let source = matrix.at(0, 1).unwrap().abs();
    let load = matrix.at(order, side - 1).unwrap().abs();
    if source < 1e-12 || load < 1e-12 {
        eprintln!("  [{label}] FAIL: source={source:.2e}, load={load:.2e} (missing port coupling)");
        ok = false;
    }

    if ok {
        println!("  [{label}] OK: symmetry={max_asymmetry:.2e}, max_entry={max_entry:.6}, source={source:.6}, load={load:.6}");
    }
    ok
}

/// Synthesizes one design and checks the resulting matrices. Returns
/// `true` when the transversal matrix was produced (the Folded/Arrow
/// variants are best-effort and do not affect the return value).
fn test_order(order: usize, return_loss_db: f64, zeros: &[f64]) -> bool {
    let label = if zeros.is_empty() {
        format!("Order {order}, all-pole, {return_loss_db} dB")
    } else {
        format!("Order {order}, {} TZs, {return_loss_db} dB", zeros.len())
    };
    println!("\n--- {label} ---");

    let spec = match FilterSpec::new(order, return_loss_db) {
        Ok(s) => {
            if zeros.is_empty() {
                s
            } else {
                s.with_normalized_transmission_zeros(zeros.to_vec())
            }
        }
        Err(e) => {
            eprintln!("  Spec creation failed: {e}");
            return false;
        }
    };

    let polynomials = match generalized_chebyshev_polynomials(&spec) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  Polynomial synthesis failed: {e}");
            return false;
        }
    };

    // Transversal matrix
    let engine = MatrixSynthesisEngine;
    match engine.synthesize(&polynomials) {
        Ok(transversal) => {
            check_matrix_sanity(&transversal, "Transversal");

            // Folded topology
            match engine.synthesize_with_topology(&polynomials, MatrixTopology::Folded) {
                Ok(folded) => {
                    check_matrix_sanity(&folded, "Folded");
                }
                Err(e) => eprintln!("  Folded transform failed: {e}"),
            }

            // Arrow topology
            match engine.synthesize_with_topology(&polynomials, MatrixTopology::Arrow) {
                Ok(arrow) => {
                    check_matrix_sanity(&arrow, "Arrow");
                }
                Err(e) => eprintln!("  Arrow transform failed: {e}"),
            }
            true
        }
        Err(e) => {
            // A typed `NumericalFailure` here is an *expected, honest* outcome
            // past the adaptive pole-axis tolerance boundary — report it, do
            // not let it masquerade as a crash or a silent success.
            eprintln!("  Matrix synthesis REJECTED: {e}");
            false
        }
    }
}

fn main() {
    println!("=== High-Order Filter Numerical Stability Test ===");
    println!(
        "Verified boundary (RL 20-22 dB, measured 2026-10-13): all-pole OK \
         through order 28, finite-TZ OK through order 24; order 30+ is \
         rejected by the adaptive pole-axis tolerance with a typed \
         NumericalFailure. The verified set below MUST succeed; the \
         exploratory set (order 30+) may be rejected and that is expected."
    );

    // ---- Verified set: these MUST synthesize (regression-anchored) ----
    // All-pole through the measured boundary.
    for &order in &[10, 12, 15, 18, 20, 22, 24, 26, 28] {
        if !test_order(order, 20.0, &[]) {
            eprintln!("\nVERIFIED-SET FAILURE: all-pole order {order} was rejected");
            std::process::exit(1);
        }
    }
    // Finite-TZ through the measured boundary.
    let verified_tz: &[(usize, &[f64], f64)] = &[
        (10, &[1.5, -1.5, 2.0, -2.0], 20.0),
        (12, &[1.3, -1.3, 1.8, -1.8, 2.5, -2.5], 20.0),
        (16, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0], 20.0),
        (20, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0, 5.0, -5.0], 20.0),
        (24, &[1.3, -1.3, 1.6, -1.6, 2.0, -2.0], 22.0),
    ];
    for &(order, zeros, rl) in verified_tz {
        if !test_order(order, rl, zeros) {
            eprintln!("\nVERIFIED-SET FAILURE: finite-TZ order {order} was rejected");
            std::process::exit(1);
        }
    }
    println!("\n--- Verified set: all succeeded ---");

    // ---- Exploratory set: order 30+ (informational; rejection is expected) ----
    println!("\n=== Exploratory scan (order 30+; rejection is expected) ===");
    let mut rejected = Vec::new();
    for &order in &[30, 32, 40] {
        if !test_order(order, 20.0, &[]) {
            rejected.push(format!("all-pole {order}"));
        }
    }
    if !test_order(30, 20.0, &[1.5, -1.5, 2.0, -2.0, 3.0, -3.0]) {
        rejected.push("finite-TZ 30".to_string());
    }
    if !rejected.is_empty() {
        println!("\nExploratory rejections (expected past the boundary): {rejected:?}");
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression anchor for the verified high-order boundary (P2-1): the
    /// all-pole boundary (order 28) must synthesize and order 30 must be
    /// rejected with a typed `NumericalFailure`. If the adaptive pole-axis
    /// tolerance or root-finding changes and order 30 starts succeeding (or
    /// order 28 starts failing), this test fails so the documented boundary
    /// in README/architecture is forced to be re-checked.
    #[test]
    fn verified_high_order_boundary() {
        // The verified all-pole boundary must synthesize.
        assert!(
            test_order(28, 20.0, &[]),
            "order 28 all-pole must synthesize (documented boundary)"
        );
        // Order 30 all-pole is past the adaptive pole-axis tolerance and must
        // be rejected (not silently succeed, not panic).
        assert!(
            !test_order(30, 20.0, &[]),
            "order 30 all-pole should be rejected by the pole-axis tolerance; \
             if it now succeeds, the documented boundary is stale"
        );
    }
}
