//! Throwaway: measure the para-identity floor (D·D̃ vs R on-axis) for the
//! N=3 triplexer and N=4 quadruplexer full-path cases, to sanity-check the
//! degree-keyed para-identity test bounds (1e-4 for d_deg ≤ 14, 1e-3 above).
use mfs::prelude::*;

fn para_max(common: &mfs::synthesis::JunctionCommon) -> f64 {
    let d_tilde = common.d.alternating_conjugate().unwrap();
    let n_tilde = common.n.alternating_conjugate().unwrap();
    let mut r_poly = common.n.clone().multiply(&n_tilde).unwrap();
    for (pt, sc) in common.pt.iter().zip(common.p0_sq.iter()) {
        let pt_conj = pt.alternating_conjugate().unwrap();
        let ptpt = pt.multiply(&pt_conj).unwrap();
        r_poly = r_poly.add(&ptpt.scale(mfs::approx::ComplexCoefficient::new(*sc, 0.0)).unwrap()).unwrap();
    }
    let dd = common.d.clone().multiply(&d_tilde).unwrap();
    let mut m: f64 = 0.0;
    for k in 0..257 {
        let w = -1.5 + 3.0 * k as f64 / 256.0;
        let s = mfs::approx::ComplexCoefficient::new(0.0, w);
        let rv = r_poly.evaluate(s);
        let dv = dd.evaluate(s);
        let rel = if rv.norm() > 1e-30 { (dv - rv).norm() / rv.norm() } else { dv.norm() };
        m = m.max(rel);
    }
    m
}

fn run(label: &str, channels: Vec<ChannelSpec>) {
    let result = MultiplexerSynthesis::new()
        .channels(channels)
        .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
        .max_iterations(80)
        .convergence_tol(1e-8)
        .full_iteration(true)
        .synthesize()
        .expect(label);
    let common = result.junction_common.as_ref().unwrap();
    println!(
        "{label}: conv_err={:.3e} d.degree={} para_max={:.3e}",
        result.convergence_error,
        common.d.degree(),
        para_max(common)
    );
}

fn edge_channel(order: usize, rl: f64, lo: f64, hi: f64) -> ChannelSpec {
    let (center, bw) = mfs::freq::edges_to_center_bandwidth(lo, hi);
    ChannelSpec { order, return_loss_db: rl, zeros: Vec::new(), center_hz: center, bandwidth_hz: bw }
}

fn main() {
    run("triplexer N=3 (no TZ)", vec![
        edge_channel(4, 22.0, 14.90e9, 15.00e9),
        edge_channel(4, 22.0, 15.10e9, 15.20e9),
        edge_channel(4, 22.0, 15.30e9, 15.40e9),
    ]);
    run("quadruplexer N=4 (no TZ)", vec![
        edge_channel(4, 22.0, 14.90e9, 15.00e9),
        edge_channel(4, 22.0, 15.10e9, 15.20e9),
        edge_channel(4, 22.0, 15.30e9, 15.40e9),
        edge_channel(4, 22.0, 15.50e9, 15.60e9),
    ]);
}
