//! Native Rust (Criterion) synthesis benchmark.
//!
//! Aligned with `fuqtgui/src/scripts/bench_synthesis.py` so the native Rust numbers
//! are directly comparable to the Python (numpy) and Python-via-FFI (`mfs` module)
//! numbers. The same three stages are benchmarked per case:
//!
//!   - synthesis : generalized Chebyshev polynomial + transversal matrix build
//!   - folding   : reduce to folded canonical form (Givens rotations)
//!   - response  : S-parameter sweep over a normalized prototype grid
//!
//! See `fuqtgui/docs/analysis/benchmark-python-vs-rust-2026-09.md` for the methodology.
//!
//! Run: `cargo bench --bench synthesis`

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::design::FilterDesign;
use mfs::freq::FrequencyGrid;
use mfs::response::ResponseSolver;

// Response sweep points, matching the Python side's N_RESPONSE_POINTS.
const N_POINTS: usize = 201;

/// A synthesis case. `zeros` are normalized prototype-domain transmission zeros (Ω).
struct Case {
    name: &'static str,
    order: usize,
    rl: f64,
    zeros: &'static [f64],
}

// Same case matrix as fuqtgui/src/scripts/bench_synthesis.py.
// golden: from Cameron's book golden benchmarks (tests/cameron_golden.py).
// sweep:   order / TZ scaling.
const GOLDEN_CASES: &[Case] = &[
    Case { name: "4-4 fully-canonical asym", order: 4, rl: 22.0, zeros: &[3.7431, 1.8051, -1.5699, -6.1910] },
    Case { name: "6-3 symmetric lattice", order: 6, rl: 22.0, zeros: &[1.0900, 1.1691, 1.5057] },
    Case { name: "4-2 asymmetric", order: 4, rl: 22.0, zeros: &[1.3217, 1.8082] },
    Case { name: "6-2 symmetric", order: 6, rl: 23.0, zeros: &[1.5522, -1.5522] },
    Case { name: "8-3 asymmetric", order: 8, rl: 23.0, zeros: &[-1.3553, 1.1093, 1.2180] },
    Case { name: "4-1 trisection", order: 4, rl: 25.0, zeros: &[2.3940] },
];

const SWEEP_CASES: &[Case] = &[
    Case { name: "order4 tz0", order: 4, rl: 22.0, zeros: &[] },
    Case { name: "order4 tz2", order: 4, rl: 22.0, zeros: &[1.8, -1.8] },
    Case { name: "order6 tz2", order: 6, rl: 22.0, zeros: &[1.8, -1.8] },
    Case { name: "order8 tz4", order: 8, rl: 23.0, zeros: &[1.5, -1.5, 2.2, -2.2] },
    Case { name: "order10 tz4", order: 10, rl: 22.0, zeros: &[1.5, -1.5, 2.2, -2.2] },
    Case { name: "order12 tz6", order: 12, rl: 22.0, zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5] },
    Case { name: "order16 tz8", order: 16, rl: 22.0, zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5, 3.0, -3.0] },
    Case { name: "order20 tz8", order: 20, rl: 22.0, zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5, 3.0, -3.0] },
];

fn build_design(c: &Case) -> mfs::design::FilterDesign {
    FilterDesign::prototype(c.order, c.rl)
        .zeros(c.zeros.iter().copied())
        .synthesize()
        .unwrap()
}

fn bench_synthesis(c: &mut Criterion) {
    for case in GOLDEN_CASES.iter().chain(SWEEP_CASES.iter()) {
        c.bench_function(
            &format!("synthesis/{}", case.name),
            |b| b.iter(|| {
                FilterDesign::prototype(case.order, case.rl)
                    .zeros(case.zeros.iter().copied())
                    .synthesize()
                    .unwrap()
            }),
        );
    }
}

fn bench_folding(c: &mut Criterion) {
    for case in GOLDEN_CASES.iter().chain(SWEEP_CASES.iter()) {
        let design = build_design(case);
        c.bench_function(
            &format!("folding/{}", case.name),
            |b| b.iter(|| black_box(design.to_folded().unwrap())),
        );
    }
}

fn bench_response(c: &mut Criterion) {
    let grid = FrequencyGrid::linspace(-1.3, 1.3, N_POINTS).unwrap();
    for case in GOLDEN_CASES.iter().chain(SWEEP_CASES.iter()) {
        let design = build_design(case);
        c.bench_function(
            &format!("response/{}", case.name),
            |b| b.iter(|| {
                black_box(ResponseSolver.evaluate_normalized(design.matrix(), &grid).unwrap())
            }),
        );
    }
}

criterion_group!(benches, bench_synthesis, bench_folding, bench_response);
criterion_main!(benches);