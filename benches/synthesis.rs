use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::design::FilterDesign;
use mfs::freq::FrequencyGrid;
use mfs::response::ResponseSolver;

fn bench_synthesis(c: &mut Criterion) {
    c.bench_function("synthesize_order4", |b| {
        b.iter(|| {
            FilterDesign::prototype(4, 20.0)
                .zeros([-1.5, 1.5])
                .synthesize()
                .unwrap()
        })
    });

    c.bench_function("synthesize_order8", |b| {
        b.iter(|| {
            FilterDesign::prototype(8, 20.0)
                .zeros([-1.3, 1.3, -2.0, 2.0])
                .synthesize()
                .unwrap()
        })
    });

    c.bench_function("synthesize_order20", |b| {
        b.iter(|| {
            FilterDesign::prototype(20, 20.0)
                .zeros([-1.2, 1.2, -1.5, 1.5, -2.0, 2.0])
                .synthesize()
                .unwrap()
        })
    });
}

fn bench_topology(c: &mut Criterion) {
    let design = FilterDesign::prototype(8, 20.0)
        .zeros([-1.3, 1.3, -2.0, 2.0])
        .synthesize()
        .unwrap();

    c.bench_function("to_folded_order8", |b| {
        b.iter(|| black_box(design.to_folded().unwrap()))
    });

    c.bench_function("to_arrow_order8", |b| {
        b.iter(|| black_box(design.to_arrow().unwrap()))
    });
}

fn bench_response(c: &mut Criterion) {
    let design = FilterDesign::prototype(8, 20.0)
        .zeros([-1.3, 1.3, -2.0, 2.0])
        .synthesize()
        .unwrap();
    let grid = FrequencyGrid::linspace(-3.0, 3.0, 201).unwrap();

    c.bench_function("response_lu_201pts_order8", |b| {
        b.iter(|| {
            black_box(ResponseSolver.evaluate_normalized(design.matrix(), &grid).unwrap())
        })
    });

    c.bench_function("response_pole_201pts_order8", |b| {
        b.iter(|| {
            black_box(design.response_normalized(-3.0, 3.0, 201).unwrap())
        })
    });

    let grid_2001 = FrequencyGrid::linspace(-3.0, 3.0, 2001).unwrap();
    c.bench_function("response_pole_2001pts_order8", |b| {
        b.iter(|| {
            black_box(design.response_normalized(-3.0, 3.0, 2001).unwrap())
        })
    });
}

criterion_group!(benches, bench_synthesis, bench_topology, bench_response);
criterion_main!(benches);
