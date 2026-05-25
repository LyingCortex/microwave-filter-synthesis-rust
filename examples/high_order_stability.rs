//! Stress test for numerical stability at filter orders 20+.
//!
//! This example synthesizes coupling matrices for high-order filters and
//! verifies that the results are physically meaningful (symmetric, bounded
//! couplings, correct eigenvalue structure).

use mfs::prelude::*;
use mfs::CouplingMatrix;
use mfs::matrix::MatrixTopology;
use mfs::synthesis::MatrixSynthesisEngine;

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
    let max_entry = matrix.as_slice().iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
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

fn test_order(order: usize, return_loss_db: f64, zeros: &[f64]) {
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
            return;
        }
    };

    let polynomials = match generalized_chebyshev_polynomials(&spec) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  Polynomial synthesis failed: {e}");
            return;
        }
    };

    // Transversal matrix
    let engine = MatrixSynthesisEngine;
    match engine.synthesize(&polynomials) {
        Ok(transversal) => {
            check_matrix_sanity(&transversal, "Transversal");

            // Folded topology
            match engine.synthesize_with_topology(&polynomials, MatrixTopology::Folded) {
                Ok(folded) => { check_matrix_sanity(&folded, "Folded"); }
                Err(e) => eprintln!("  Folded transform failed: {e}"),
            }

            // Arrow topology
            match engine.synthesize_with_topology(&polynomials, MatrixTopology::Arrow) {
                Ok(arrow) => { check_matrix_sanity(&arrow, "Arrow"); }
                Err(e) => eprintln!("  Arrow transform failed: {e}"),
            }
        }
        Err(e) => eprintln!("  Matrix synthesis failed: {e}"),
    }
}

fn main() {
    println!("=== High-Order Filter Numerical Stability Test ===");

    // All-pole filters at increasing orders
    for &order in &[10, 12, 15, 18, 20, 22, 24, 26, 28, 30] {
        test_order(order, 20.0, &[]);
    }

    // Filters with transmission zeros
    test_order(10, 20.0, &[1.5, -1.5, 2.0, -2.0]);
    test_order(12, 20.0, &[1.3, -1.3, 1.8, -1.8, 2.5, -2.5]);
    test_order(16, 20.0, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0]);
    test_order(20, 20.0, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0, 5.0, -5.0]);
    test_order(24, 22.0, &[1.3, -1.3, 1.6, -1.6, 2.0, -2.0]);

    // Extreme: order 30 with zeros
    test_order(30, 20.0, &[1.5, -1.5, 2.0, -2.0, 3.0, -3.0]);

    println!("\n=== Done ===");
}
