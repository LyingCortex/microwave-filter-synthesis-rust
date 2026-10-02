//! Throwaway: test whether the full Macchiarella path converges when the
//! bands follow the paper's convention (RX = lower 14.9-15.1, TX = upper
//! 15.15-15.35): `assign_roots` assumes RX is the lower band.
use mfs::duplexer::{DuplexerSynthesis, JunctionType};

fn run(label: &str, order: usize, tx_lo: f64, tx_hi: f64, rx_lo: f64, rx_hi: f64) {
    let s = DuplexerSynthesis::new()
        .tx(order, 20.0, tx_lo, tx_hi)
        .rx(order, 20.0, rx_lo, rx_hi)
        .junction(JunctionType::TypeI { n: 1.47, b0: -0.171 })
        .max_iterations(60)
        .convergence_tol(1e-8)
        .full_iteration(true);
    match s.synthesize() {
        Ok(r) => println!(
            "[{label} order {order}] converged={} err={:.3e} iters={}",
            r.convergence_error < 1e-6,
            r.convergence_error,
            r.iterations
        ),
        Err(e) => println!("[{label} order {order}] ERROR: {e}"),
    }
}

fn main() {
    // Paper Section V convention: RX lower, TX upper.
    run("paper RX-low", 2, 15.15e9, 15.35e9, 14.9e9, 15.1e9);
    run("paper RX-low", 4, 15.15e9, 15.35e9, 14.9e9, 15.1e9);
    run("paper RX-low", 7, 15.15e9, 15.35e9, 14.9e9, 15.1e9);
    // diag1 convention: TX lower, RX upper.
    run("code TX-low", 2, 14.9e9, 15.1e9, 15.15e9, 15.35e9);
}
