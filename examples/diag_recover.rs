//! Throwaway: inspect the recovered TX/RX channel polynomials and passband
//! response from the full Macchiarella path vs the simplified path.
use mfs::duplexer::{Duplexer, DuplexerSynthesis, JunctionType};

fn run(label: &str, full: bool) {
    let r = DuplexerSynthesis::new()
        .tx(2, 22.0, 14.9e9, 15.1e9)
        .rx(2, 22.0, 15.15e9, 15.35e9)
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(40)
        .convergence_tol(1e-9)
        .full_iteration(full)
        .synthesize()
        .unwrap();
    let resp = Duplexer::new(r.tx.clone(), r.rx.clone())
        .unwrap()
        .response(14.8e9, 15.45e9, 651)
        .unwrap();
    let best_tx = resp
        .frequencies
        .iter()
        .enumerate()
        .filter(|(_, f)| **f > 14.9e9 && **f < 15.1e9)
        .map(|(i, _)| resp.tx_s21_db[i])
        .fold(f64::NEG_INFINITY, f64::max);
    let best_rx = resp
        .frequencies
        .iter()
        .enumerate()
        .filter(|(_, f)| **f > 15.15e9 && **f < 15.35e9)
        .map(|(i, _)| resp.rx_s21_db[i])
        .fold(f64::NEG_INFINITY, f64::max);
    println!("[{label}] err={:.2e} iters={} best_tx_s21={best_tx:.3} dB best_rx_s21={best_rx:.3} dB",
        r.convergence_error, r.iterations);
    let pol = r.tx.polynomials();
    println!("  tx E={:?}", pol.e.coefficients.iter().map(|c| (c.re, c.im)).collect::<Vec<_>>());
    println!("  tx F={:?}", pol.f.coefficients.iter().map(|c| (c.re, c.im)).collect::<Vec<_>>());
    println!("  tx P={:?}", pol.p.coefficients.iter().map(|c| (c.re, c.im)).collect::<Vec<_>>());
}

fn main() {
    run("simplified", false);
    run("full", true);
}
