//! Throwaway: compare independent vs recovered channel standalone response.
use mfs::duplexer::{DuplexerSynthesis, JunctionType};
use mfs::design::FilterDesign;

fn peak(tag: &str, d: &FilterDesign, lo: f64, hi: f64) {
    let resp = d.response(lo, hi, 201).unwrap();
    let (mut best, mut bf, mut bs11) = (f64::NEG_INFINITY, 0.0, 0.0);
    for s in resp.samples.iter() {
        if s.s21_db() > best { best = s.s21_db(); bf = s.frequency_hz; bs11 = s.s11_db(); }
    }
    println!("{tag}: peak S21 = {best:.4} dB at {:.4} GHz (S11 there {bs11:.2} dB)", bf / 1e9);
}

fn main() {
    let ind = FilterDesign::bandpass_edges(2, 22.0, 14.9e9, 15.1e9).synthesize().unwrap();
    println!("ind E={:?}", ind.polynomials().e.coefficients);
    println!("ind F={:?}", ind.polynomials().f.coefficients);
    println!("ind P={:?}", ind.polynomials().p.coefficients);
    peak("independent", &ind, 14.9e9, 15.1e9);

    let r = DuplexerSynthesis::new()
        .tx(2, 22.0, 14.9e9, 15.1e9)
        .rx(2, 22.0, 15.15e9, 15.35e9)
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(40)
        .convergence_tol(1e-9)
        .full_iteration(true)
        .synthesize()
        .unwrap();
    println!("rec E={:?}", r.tx.polynomials().e.coefficients);
    println!("rec F={:?}", r.tx.polynomials().f.coefficients);
    println!("rec P={:?}", r.tx.polynomials().p.coefficients);
    peak("recovered", &r.tx, 14.9e9, 15.1e9);
}
