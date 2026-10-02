use mfs::duplexer::{DuplexerSynthesis, JunctionType};
fn main() {
    let s = DuplexerSynthesis::new()
        .tx(2, 22.0, 14.9e9, 15.1e9)
        .rx(2, 22.0, 15.15e9, 15.35e9)
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(60).convergence_tol(1e-8).full_iteration(true);
    match s.synthesize() {
        Ok(r) => println!("OK err={:.3e} iters={}", r.convergence_error, r.iterations),
        Err(e) => println!("ERR {e}"),
    }
}
