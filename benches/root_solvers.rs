//! Criterion timing benchmark for the four complex root solvers.
//!
//! Each solver runs on a fixed, deterministic polynomial (known roots at
//! distinct angles/radii) at degrees {10, 20, 30}. The polynomial is built
//! once per degree so only `roots_of` is timed.
//!
//! Run: `cargo bench --bench root_solvers`

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::approx::{
    AberthRootSolver, CompanionMatrixRootSolver, ComplexCoefficient, ComplexPolynomial,
    ComplexRootSolver, DurandKernerRootSolver, FaerRootSolver,
};

const DEGREES: &[usize] = &[10, 20, 30];

fn fixed_complex_roots(degree: usize) -> Vec<ComplexCoefficient> {
    (0..degree)
        .map(|k| {
            let radius = 0.5 + 0.25 * (((k * 7 + 3) % 5) as f64);
            let angle = 2.0 * std::f64::consts::PI * k as f64 / degree as f64 + 0.4;
            ComplexCoefficient::new(radius * angle.cos(), radius * angle.sin())
        })
        .collect()
}

fn bench_root_solvers(c: &mut Criterion) {
    for &degree in DEGREES {
        let polynomial =
            ComplexPolynomial::from_complex_roots(&fixed_complex_roots(degree)).unwrap();

        c.bench_function(&format!("root_solvers/durand_kerner/d{degree}"), |b| {
            b.iter(|| {
                black_box(
                    DurandKernerRootSolver
                        .roots_of(black_box(&polynomial))
                        .unwrap(),
                )
            })
        });
        c.bench_function(&format!("root_solvers/aberth/d{degree}"), |b| {
            b.iter(|| black_box(AberthRootSolver.roots_of(black_box(&polynomial)).unwrap()))
        });
        c.bench_function(&format!("root_solvers/companion_matrix/d{degree}"), |b| {
            b.iter(|| {
                black_box(
                    CompanionMatrixRootSolver
                        .roots_of(black_box(&polynomial))
                        .unwrap(),
                )
            })
        });
        c.bench_function(&format!("root_solvers/faer/d{degree}"), |b| {
            b.iter(|| black_box(FaerRootSolver.roots_of(black_box(&polynomial)).unwrap()))
        });
    }
}

criterion_group!(benches, bench_root_solvers);
criterion_main!(benches);
