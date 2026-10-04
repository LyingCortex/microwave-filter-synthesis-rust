//! Criterion benchmark for the LM optimizer (`optimize::optimize_matrix` /
//! `optimize::extract_matrix`).
//!
//! **Deliberately not part of the cross-language benchmark standard** (the same
//! stance as the program-generated `SWEEP_CASES` in `benches/synthesis.rs`):
//! the standard's `gates.json` defines no `optimizer` family, and inventing
//! case ids or tolerances without provenance is forbidden. This is instead a
//! program-generated order sweep measuring, alongside wall time, the two
//! quantities the optimizer's contract is written in:
//!
//! * **iterations** — how many LM steps the run needed,
//! * **cost reduction** — final cost vs. the initial guess's cost,
//!
//! plus the physical quality of the result (worst *relative* coupling error for
//! tuning; max|ΔS| against the target for extraction).
//!
//! Two scenarios, mirroring `tests/optimizer_integration.rs`:
//!
//! * `tune` — synthesize a order-N folded prototype, perturb 4 inter-resonator
//!   couplings by 7 %, and re-optimize from the perturbed (good-guess) matrix.
//!   The realistic extraction setting: the guess is the real filter plus small
//!   error, and the solution is locally unique.
//! * `extract` — feed `extract_matrix` the response of the synthesized
//!   transversal matrix with no help (it builds its own crude inline-chain
//!   guess).
//!
//! Order-8 `extract` is **not timed**: from the cold chain guess it converges
//! to a local minimum far from the target (measured 2026-10-04: cost 39.2,
//! max|ΔS| ≈ 1.0 after 435 iterations / ~85 s), so a Criterion timing of it
//! would measure noise, not the algorithm. The limitation is documented in
//! `docs/numerical-optimization-report.md` §6.
//!
//! Run: `cargo bench --bench optimizer`

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::design::FilterDesign;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::optimize::{extract_matrix, optimize_matrix, OptimizeConfig};
use mfs::response::{ResponseSolver, SParameterResponse};

/// Order sweep. `tune` covers all three; `extract` covers 4 and 6 (see module
/// docs for why 8 is excluded).
const ORDERS: &[usize] = &[4, 6, 8];
const EXTRACT_ORDERS: &[usize] = &[4, 6];

fn grid() -> FrequencyGrid {
    FrequencyGrid::linspace(-2.0, 2.0, 41).expect("grid")
}

/// A perturbed copy of `matrix`: the first `count` non-zero inter-resonator
/// couplings (upper triangle) scaled by `factor`, everything else untouched.
fn perturb_inter_resonator(matrix: &CouplingMatrix, count: usize, factor: f64) -> CouplingMatrix {
    let side = matrix.side();
    let order = matrix.order();
    let mut data = matrix.as_slice().to_vec();
    let mut done = 0usize;
    'outer: for row in 1..=order {
        for col in (row + 1)..=order {
            let idx = row * side + col;
            if data[idx].abs() > 1e-9 {
                data[idx] *= factor;
                data[col * side + row] *= factor;
                done += 1;
                if done == count {
                    break 'outer;
                }
            }
        }
    }
    CouplingMatrix::new_with_topology(order, matrix.topology(), data).expect("perturbed")
}

/// Worst relative error over the original's non-negligible upper-triangle entries.
fn worst_rel(original: &CouplingMatrix, recovered: &CouplingMatrix) -> f64 {
    let side = original.side();
    let mut worst = 0.0f64;
    for i in 0..side {
        for j in i..side {
            let o = original.at(i, j).unwrap_or(0.0);
            let r = recovered.at(i, j).unwrap_or(0.0);
            if o.abs() > 1e-9 {
                worst = worst.max((o - r).abs() / o.abs());
            }
        }
    }
    worst
}

/// Point-wise max deviation across all four S-parameter channels.
fn max_dev(a: &SParameterResponse, b: &SParameterResponse) -> f64 {
    let mut m = 0.0f64;
    for (x, y) in a.samples.iter().zip(b.samples.iter()) {
        m = m
            .max((x.s11_re - y.s11_re).abs())
            .max((x.s11_im - y.s11_im).abs())
            .max((x.s21_re - y.s21_re).abs())
            .max((x.s21_im - y.s21_im).abs());
    }
    m
}

fn bench_tune(c: &mut Criterion) {
    let g = grid();
    let config = OptimizeConfig::default();

    for &order in ORDERS {
        let original = FilterDesign::prototype(order, 20.0)
            .synthesize()
            .expect("synthesize")
            .to_folded()
            .expect("fold");
        let target = ResponseSolver.evaluate_normalized(&original, &g).expect("target");
        let start = perturb_inter_resonator(&original, (order - 1).min(4), 1.07);
        let start_cost: f64 = {
            let r = ResponseSolver
                .evaluate_normalized(&start, &g)
                .expect("start response");
            max_dev(&r, &target)
        };

        // Untimed run for the contract metrics (iterations + cost reduction).
        let r = optimize_matrix(&start, &target, &g, &config).expect("optimize");
        println!(
            "tune   order={order}: iters={} converged={} cost={:.3e} start_dev={:.3e} \
             worst_rel={:.3e}",
            r.iterations,
            r.converged,
            r.cost,
            start_cost,
            worst_rel(&original, &r.matrix),
        );

        c.bench_function(&format!("optimizer/tune/o{order}"), |b| {
            b.iter(|| {
                black_box(
                    optimize_matrix(
                        black_box(&start),
                        black_box(&target),
                        black_box(&g),
                        black_box(&config),
                    )
                    .expect("optimize"),
                )
            })
        });
    }
}

fn bench_extract(c: &mut Criterion) {
    let g = grid();

    for &order in EXTRACT_ORDERS {
        let target_matrix = FilterDesign::prototype(order, 20.0)
            .synthesize()
            .expect("synthesize")
            .matrix()
            .clone();
        let target =
            ResponseSolver.evaluate_normalized(&target_matrix, &g).expect("target");

        // Untimed run for the contract metrics.
        let r = extract_matrix(&target, &g, order).expect("extract");
        let recovered = ResponseSolver
            .evaluate_normalized(&r.matrix, &g)
            .expect("recovered response");
        println!(
            "extract order={order}: iters={} converged={} cost={:.3e} max|dS|={:.3e}",
            r.iterations,
            r.converged,
            r.cost,
            max_dev(&recovered, &target),
        );

        c.bench_function(&format!("optimizer/extract/o{order}"), |b| {
            b.iter(|| black_box(extract_matrix(black_box(&target), black_box(&g), order).expect("extract")))
        });
    }
}

/// Criterion with a small sample size: one optimizer run is 0.1–0.8 s, so the
/// default 100 samples would turn a sweep into a several-minute run.
fn optimizer_criterion() -> Criterion {
    Criterion::default().sample_size(10)
}

criterion_group! {
    name = benches;
    config = optimizer_criterion();
    targets = bench_tune, bench_extract
}
criterion_main!(benches);
