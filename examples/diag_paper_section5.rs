//! Throwaway: run the WIP common-domain `full_iteration(true)` on the
//! Macchiarella–Tamiazzo Section-V spec (7/7, 20 dB, n=1.47, b0=-0.171)
//! and the D3 band (order 4, 6). Reports convergence, iteration count,
//! and the recovered channel fidelity against the independent bound.
use mfs::duplexer::{DuplexerSynthesis, JunctionType};
fn run(order: usize) {
    let s = DuplexerSynthesis::new()
        .tx(order, 20.0, 14.9e9, 15.1e9)
        .rx(order, 20.0, 15.15e9, 15.35e9)
        .junction(JunctionType::TypeI { n: 1.47, b0: -0.171 })
        .max_iterations(60)
        .convergence_tol(1e-8)
        .full_iteration(true);
    match s.synthesize() {
        Ok(r) => {
            let ok = r.iterations > 0 && r.convergence_error < 1e-6;
            println!(
                "[order {}] converged={} err={:.3e} iters={}  (satisfied tol? {})",
                order,
                ok,
                r.convergence_error,
                r.iterations,
                r.convergence_error < 1e-6,
            );
        }
        Err(e) => println!("[order {}] ERROR: {e}", order),
    }
}

fn main() {
    run(2);   // D1 RED-gate setup (currently #[ignore]d)
    run(4);   // D3 band
    run(6);   // D3 band
    run(7);   // Section-V spec
}
