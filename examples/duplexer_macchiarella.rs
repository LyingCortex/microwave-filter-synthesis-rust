//! Macchiarella-Tamiazzo 2006 Duplexer Synthesis Example
//!
//! Reproduces Example 1 from:
//! G. Macchiarella and S. Tamiazzo, "Novel Approach to the Synthesis of Microwave
//! Diplexers," IEEE Trans. MTT, vol. 54, no. 12, Dec 2006.
//!
//! Waveguide diplexer:
//!   - RX filter: 14.9-15.1 GHz, 7-pole, 22 dB RL
//!   - TX filter: 15.15-15.35 GHz, 7-pole, 22 dB RL
//!   - Type-I junction: n = 1, b₀ = 0 (ideal T-junction)
//!
//! NOTE: The current simplified perturbation approach demonstrates the
//! duplexer synthesis API but does not yet achieve the full quasi-equiripple
//! antenna return loss of the paper. The spectral factorization infrastructure
//! is in place for future numerical improvements.
//!
//! Run: cargo run --example duplexer_macchiarella

use mfs::duplexer::{Duplexer, DuplexerSynthesis, DuplexerSynthesisResult, JunctionType};
use mfs::prelude::*;

fn main() -> Result<()> {
    println!("═══════════════════════════════════════════════════════════════");
    println!("  Macchiarella-Tamiazzo 2006 - Diplexer Synthesis");
    println!("  Example 1: Waveguide Diplexer (14.9-15.35 GHz)");
    println!("═══════════════════════════════════════════════════════════════\n");

    // ─── Paper Example 1: 7-pole waveguide diplexer ──────────────────────
    println!("─── Example 1: All-pole 7th-order waveguide diplexer ───");
    println!("  RX: 14.9 - 15.1 GHz (BW = 200 MHz), 7-pole, RL = 22 dB");
    println!("  TX: 15.15 - 15.35 GHz (BW = 200 MHz), 7-pole, RL = 22 dB");
    println!("  Junction: Type-I (n=1, b₀=0)\n");

    let result = DuplexerSynthesis::new()
        .rx(7, 22.0, 14.9e9, 15.1e9)
        .tx(7, 22.0, 15.15e9, 15.35e9)
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(50)
        .convergence_tol(1e-10)
        .synthesize()?;

    print_synthesis_result(&result, "Example 1");

    // ─── Individual filter responses ─────────────────────────────────────
    println!("\n─── RX Filter Response (14.9-15.1 GHz passband) ───\n");
    let rx_resp = result.rx.response(14.8e9, 15.2e9, 21)?;
    println!("  {:>10}  {:>10}  {:>10}", "Freq(GHz)", "S21(dB)", "S11(dB)");
    for i in (0..21).step_by(2) {
        println!("  {:>10.4}  {:>10.3}  {:>10.3}",
            rx_resp.samples[i].frequency_hz / 1e9,
            rx_resp.samples[i].s21_db(),
            rx_resp.samples[i].s11_db());
    }

    println!("\n─── TX Filter Response (15.15-15.35 GHz passband) ───\n");
    let tx_resp = result.tx.response(15.05e9, 15.45e9, 21)?;
    println!("  {:>10}  {:>10}  {:>10}", "Freq(GHz)", "S21(dB)", "S11(dB)");
    for i in (0..21).step_by(2) {
        println!("  {:>10.4}  {:>10.3}  {:>10.3}",
            tx_resp.samples[i].frequency_hz / 1e9,
            tx_resp.samples[i].s21_db(),
            tx_resp.samples[i].s11_db());
    }

    // ─── Combined duplexer using Duplexer struct ─────────────────────────
    println!("\n─── Combined Duplexer (using simple T-junction model) ───\n");
    let duplexer = Duplexer::new(result.tx.clone(), result.rx.clone())?;
    let combined = duplexer.response(14.7e9, 15.5e9, 21)?;
    println!("  {:>10}  {:>10}  {:>10}  {:>10}", "Freq(GHz)", "TX(dB)", "RX(dB)", "Iso(dB)");
    for i in 0..21 {
        println!("  {:>10.4}  {:>10.2}  {:>10.2}  {:>10.2}",
            combined.frequencies[i] / 1e9,
            combined.tx_s21_db[i],
            combined.rx_s21_db[i],
            combined.isolation_db[i]);
    }

    // ─── Example with transmission zeros ─────────────────────────────────
    println!("\n\n─── Example 2: 6th-order with transmission zeros ───\n");

    let result2 = DuplexerSynthesis::new()
        .rx(6, 22.0, 14.9e9, 15.1e9)
        .rx_zeros([-2.5, 2.5])
        .tx(6, 22.0, 15.15e9, 15.35e9)
        .tx_zeros([-2.5, 2.5])
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(50)
        .synthesize()?;

    print_synthesis_result(&result2, "Example 2 (with TZ)");

    println!("\n═══════════════════════════════════════════════════════════════");
    println!("  Status: The simplified perturbation approach produces valid");
    println!("  filter designs. Full quasi-equiripple antenna S11 matching");
    println!("  the paper's Fig. 3-4 requires robust spectral factorization");
    println!("  of high-degree polynomials (future enhancement).");
    println!("═══════════════════════════════════════════════════════════════");

    Ok(())
}

fn print_synthesis_result(result: &DuplexerSynthesisResult, name: &str) {
    println!("  {} synthesis result:", name);
    println!("    Iterations: {}", result.iterations);
    println!("    Convergence error: {:.2e}", result.convergence_error);
    println!("    TX order: {}, RX order: {}", result.tx.order(), result.rx.order());
    println!("    Junction: {:?}", result.junction);

    // Print coupling matrix diagonal (resonant frequencies)
    let m = result.tx.matrix();
    println!("\n    TX coupling matrix (transversal):");
    println!("      Input coupling M_S1 = {:.6}", m.at(0, 1).unwrap_or(0.0));
    print!("      Diagonal: [");
    for i in 1..=result.tx.order() {
        print!("{:.4}", m.at(i, i).unwrap_or(0.0));
        if i < result.tx.order() { print!(", "); }
    }
    println!("]");
    println!("      Output coupling M_NL = {:.6}", m.at(result.tx.order(), result.tx.order() + 1).unwrap_or(0.0));

    let m = result.rx.matrix();
    println!("\n    RX coupling matrix (transversal):");
    println!("      Input coupling M_S1 = {:.6}", m.at(0, 1).unwrap_or(0.0));
    print!("      Diagonal: [");
    for i in 1..=result.rx.order() {
        print!("{:.4}", m.at(i, i).unwrap_or(0.0));
        if i < result.rx.order() { print!(", "); }
    }
    println!("]");
    println!("      Output coupling M_NL = {:.6}", m.at(result.rx.order(), result.rx.order() + 1).unwrap_or(0.0));
}
