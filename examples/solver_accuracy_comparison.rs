//! Compares root solver accuracy across different orders.
//!
//! Tests each solver on the actual polynomials produced by the filter synthesis
//! pipeline, measuring residual |p(z_k)| at each computed root.

use mfs::prelude::*;
use mfs::approx::{
    ComplexRootSolver, DurandKernerRootSolver, AberthRootSolver,
    CompanionMatrixRootSolver, AdaptiveRootSolver,
    generalized_chebyshev_polynomials,
};

fn test_solver_on_filter(
    order: usize,
    zeros: &[f64],
    rl: f64,
    label: &str,
) {
    let spec = match FilterSpec::new(order, rl) {
        Ok(s) => {
            if zeros.is_empty() { s }
            else { s.with_normalized_transmission_zeros(zeros.to_vec()) }
        }
        Err(e) => { eprintln!("  [{label}] spec error: {e}"); return; }
    };

    let polys = match generalized_chebyshev_polynomials(&spec) {
        Ok(p) => p,
        Err(e) => { eprintln!("  [{label}] poly error: {e}"); return; }
    };

    // Get the E(s) polynomial — this is the one we need to root-find
    let e_poly = &polys.e;
    let degree = e_poly.degree();

    println!("\n  [{label}] degree={degree}, coeff range: {:.2e} to {:.2e}",
        e_poly.coefficients.iter().map(|c| c.norm()).filter(|&n| n > 1e-20).fold(f64::MAX, f64::min),
        e_poly.coefficients.iter().map(|c| c.norm()).fold(0.0_f64, f64::max),
    );

    // Test each solver
    for (name, result) in [
        ("DurandKerner", DurandKernerRootSolver.roots_of(e_poly)),
        ("Aberth", AberthRootSolver.roots_of(e_poly)),
        ("Companion", CompanionMatrixRootSolver.roots_of(e_poly)),
        ("Adaptive", AdaptiveRootSolver.roots_of(e_poly)),
    ] {
        match result {
            Ok(roots) => {
                let max_residual = roots.iter()
                    .map(|&r| e_poly.evaluate(r).norm())
                    .fold(0.0_f64, f64::max);
                let avg_residual = roots.iter()
                    .map(|&r| e_poly.evaluate(r).norm())
                    .sum::<f64>() / roots.len() as f64;
                println!("    {name:15} max_residual={max_residual:.2e}  avg={avg_residual:.2e}");
            }
            Err(e) => {
                println!("    {name:15} FAILED: {e}");
            }
        }
    }
}

fn main() {
    println!("=== Root Solver Accuracy Comparison ===");

    // All-pole filters at increasing orders
    for &order in &[4, 8, 12, 16, 20, 24, 28, 30, 35, 40] {
        test_solver_on_filter(order, &[], 20.0,
            &format!("Order {order} all-pole"));
    }

    // With transmission zeros
    test_solver_on_filter(10, &[1.5, -1.5, 2.0, -2.0], 20.0,
        "Order 10, 4 TZ");
    test_solver_on_filter(16, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0], 20.0,
        "Order 16, 6 TZ");
    test_solver_on_filter(20, &[1.2, -1.2, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0], 20.0,
        "Order 20, 8 TZ");
    test_solver_on_filter(30, &[1.5, -1.5, 2.0, -2.0, 3.0, -3.0], 20.0,
        "Order 30, 6 TZ");

    println!("\n=== Done ===");
}
