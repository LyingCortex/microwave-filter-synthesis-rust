//! Native Rust (Criterion) synthesis benchmark.
//!
//! Aligned with `fuqtgui/src/scripts/bench_synthesis.py` so the native Rust
//! numbers are directly comparable to the Python (numpy) and Python-via-FFI
//! (`mfs` module) numbers. The same three stages are benchmarked per case:
//!
//!   - synthesis : generalized Chebyshev polynomial + transversal matrix build
//!   - folding   : reduce to folded canonical form (Givens rotations)
//!   - response  : S-parameter sweep over a normalized prototype grid
//!
//! The stage names and their order come from the standard's `gates.json`
//! (`data::stages("twoport")`), and the case set comes from the standard's
//! `cases/twoport.json` (`data::cases("twoport")`). No case is hardcoded here:
//! adding a case to the standard is enough to benchmark it (see
//! `cases/README`-style provenance in the JSON). Timing policy (warm-up,
//! sample count) is read from `harness.json` per case.
//!
//! See `fuqtgui/docs/analysis/benchmark-python-vs-rust-2026-09.md` for the
//! methodology; `tests/bench_support/` holds the loader shared with the
//! `benchmark_standard` conformance test.
//!
//! Run: `cargo bench --bench synthesis`

#[path = "../tests/bench_support/mod.rs"]
mod bench_support;

use std::path::Path;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::design::FilterDesign;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::response::{ResponseSolver, SParameterResponse};

use bench_support::data::{self, Case};
use bench_support::gates::{
    self, apply_harness_policy, TWOPORT_RESPONSE_BOUND, TWOPORT_RESPONSE_POINTS,
};

/// A program-generated order/transmission-zero scan.
///
/// **Deliberately not part of the standard.** The standard's cases are either
/// publication-anchored (`Cameron 2018` figures, with printed matrices in
/// `golden/`) or synthetic cases carrying declared physical invariants; a
/// parameter sweep has neither provenance nor invariants, so it cannot be
/// expressed as a `cases/<family>.json` entry. It is kept here as a scaling
/// probe for mfs's own regression interest, and is benchmarked with the
/// standard's *default* harness policy.
struct SweepCase {
    name: &'static str,
    order: usize,
    rl: f64,
    zeros: &'static [f64],
}

const SWEEP_CASES: &[SweepCase] = &[
    SweepCase {
        name: "order4 tz0",
        order: 4,
        rl: 22.0,
        zeros: &[],
    },
    SweepCase {
        name: "order4 tz2",
        order: 4,
        rl: 22.0,
        zeros: &[1.8, -1.8],
    },
    SweepCase {
        name: "order6 tz2",
        order: 6,
        rl: 22.0,
        zeros: &[1.8, -1.8],
    },
    SweepCase {
        name: "order8 tz4",
        order: 8,
        rl: 23.0,
        zeros: &[1.5, -1.5, 2.2, -2.2],
    },
    SweepCase {
        name: "order10 tz4",
        order: 10,
        rl: 22.0,
        zeros: &[1.5, -1.5, 2.2, -2.2],
    },
    SweepCase {
        name: "order12 tz6",
        order: 12,
        rl: 22.0,
        zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5],
    },
    SweepCase {
        name: "order16 tz8",
        order: 16,
        rl: 22.0,
        zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5, 3.0, -3.0],
    },
    SweepCase {
        name: "order20 tz8",
        order: 20,
        rl: 22.0,
        zeros: &[1.5, -1.5, 2.0, -2.0, 2.5, -2.5, 3.0, -3.0],
    },
];

/// Standard stage names, asserted against the contract so a rename upstream
/// fails loudly instead of silently relabelling the benchmark groups.
///
/// `data::stages` returns the family's stages as a set (the contract stores them
/// as a JSON object, so the iteration order is alphabetical, not pipeline order);
/// the pipeline order below is the one `gates.json`'s per-stage descriptions
/// document, and each name is required to be present.
fn twoport_stages() -> [String; 3] {
    let stages = data::stages("twoport").expect("standard defines twoport stages");
    let found: std::collections::BTreeSet<&str> = stages.iter().map(String::as_str).collect();
    assert_eq!(
        found,
        ["synthesis", "folding", "response"].into_iter().collect(),
        "standard twoport stage set changed"
    );
    ["synthesis", "folding", "response"].map(str::to_string)
}

fn build_design(order: usize, rl: f64, zeros: &[f64]) -> mfs::design::FilterDesign {
    FilterDesign::prototype(order, rl)
        .zeros(zeros.iter().copied())
        .synthesize()
        .unwrap()
}

fn bench_standard_cases(c: &mut Criterion) {
    let [synthesis_stage, folding_stage, response_stage] = twoport_stages();
    let cases = data::cases("twoport").expect("standard twoport cases");
    let grid = FrequencyGrid::linspace(
        -TWOPORT_RESPONSE_BOUND,
        TWOPORT_RESPONSE_BOUND,
        TWOPORT_RESPONSE_POINTS,
    )
    .expect("contract response grid");

    let mut report = data::Report::new("rust-native", "rust", &data::platform());

    for case in &cases {
        let (order, rl, zeros) = case_inputs(case);
        let row = report.add_case(&case.qualid);

        // Criterion timing (kept for regression visibility).
        {
            let mut group = c.benchmark_group(format!("{synthesis_stage}/{}", case.qualid));
            apply_harness_policy(&mut group, &case.qualid);
            group.bench_function(&case.qualid, |b| b.iter(|| build_design(order, rl, &zeros)));
            group.finish();
        }

        // Report wall-clock + gates.
        let harness = data::harness(&case.qualid);
        let synth_timing = data::time_call(
            || {
                let _ = build_design(order, rl, &zeros);
            },
            harness.timing.warmup,
            harness.timing.repeats,
        );
        row.set_stage(&synthesis_stage, synth_timing);

        let design = build_design(order, rl, &zeros);
        let transversal = design.matrix().clone();
        let folded = design.to_folded().expect("folded form");

        let fold_timing = data::time_call(
            || {
                let _ = build_design(order, rl, &zeros).to_folded();
            },
            harness.timing.warmup,
            harness.timing.repeats,
        );
        row.set_stage(&folding_stage, fold_timing);

        let resp_timing = data::time_call(
            || {
                let _ = ResponseSolver.evaluate_normalized(design.matrix(), &grid);
            },
            harness.timing.warmup,
            harness.timing.repeats,
        );
        row.set_stage(&response_stage, resp_timing);

        set_twoport_gates(row, case, order, rl, &transversal, &folded, &grid);
    }

    report
        .write(Path::new("target/bench-synthesis-report.json"))
        .expect("write synthesis bench report");
}

/// The contract's four twoport gates on mfs's own synthesis + folding,
/// evaluated on the normalized prototype grid. Mirrors
/// `fuqtgui/src/scripts/bench_synthesis.py::run_twoport_gates` so the Python and
/// Rust reports carry the same gate set for a per-case diff.
fn set_twoport_gates(
    row: &mut data::CaseResult,
    case: &Case,
    order: usize,
    rl: f64,
    transversal: &CouplingMatrix,
    folded: &CouplingMatrix,
    grid: &FrequencyGrid,
) {
    let tols = data::gate_tols();
    let n_points = grid.as_slice().len();
    let eval = |m: &CouplingMatrix| -> SParameterResponse {
        ResponseSolver
            .evaluate_normalized(m, grid)
            .expect("response")
    };
    let r_t = eval(transversal);
    let r_f = eval(folded);

    // 1. response_invariance (hard, 1e-10)
    let inv = gates::response_invariance_metric(&r_t, &r_f).expect("equal-length");
    let tol_inv = tols.get("response_invariance").copied().unwrap_or(1e-10);
    row.set_gate(
        "response_invariance",
        data::GateResult {
            status: if inv < tol_inv { "PASS" } else { "FAIL" }.to_string(),
            metric: Some(inv),
            tol: Some(tol_inv),
            detail: Some(format!(
                "transversal vs folded S11/S21 max|Δ| at {n_points} points"
            )),
        },
    );

    // 2. folded_mask (hard, 0.0)
    let worst = gates::folded_mask_violations(folded);
    let tol_mask = tols.get("folded_mask").copied().unwrap_or(0.0);
    row.set_gate(
        "folded_mask",
        data::GateResult {
            status: if worst <= tol_mask.max(1e-6) {
                "PASS"
            } else {
                "FAIL"
            }
            .to_string(),
            metric: Some(worst),
            tol: Some(tol_mask),
            detail: Some(format!("max|off-mask|={worst:e} (either convention)")),
        },
    );

    // 3. book_matrix (soft, 0.02) — only when the anchor is itself a folded form.
    let tol_coup = tols.get("coupling_abs").copied().unwrap_or(0.02);
    match case.anchor_id.as_deref() {
        None => row.set_gate(
            "book_matrix",
            data::GateResult::skip("case has no anchor; book matrix not transcribed"),
        ),
        Some(id) => {
            let anchors = data::anchors("twoport").expect("twoport anchors");
            let entry = anchors.get(id);
            match entry {
                None => row.set_gate(
                    "book_matrix",
                    data::GateResult::skip(format!(
                        "anchor {id} not in golden/twoport_anchors.json"
                    )),
                ),
                Some(e) if e["folded"].is_null() => row.set_gate(
                    "book_matrix",
                    data::GateResult::skip(format!(
                        "anchor {id} has no transcribed folded matrix (folded = null)"
                    )),
                ),
                Some(e) if e["topology"].as_str() != Some("folded") => row.set_gate(
                    "book_matrix",
                    data::GateResult::skip(format!(
                        "anchor {id} topology {:?} is not a folded form",
                        e["topology"].as_str().unwrap_or("?")
                    )),
                ),
                Some(e) => {
                    let flat: Vec<f64> = e["folded"]
                        .as_array()
                        .expect("folded anchor is an array")
                        .iter()
                        .flat_map(|r| {
                            r.as_array()
                                .expect("row is an array")
                                .iter()
                                .map(|v| v.as_f64().expect("numeric"))
                        })
                        .collect();
                    let book = CouplingMatrix::new(order, flat).expect("book folded matrix");
                    match gates::coupling_set_deviation(folded, &book) {
                        Ok(dev) => row.set_gate(
                            "book_matrix",
                            data::GateResult {
                                status: if dev < tol_coup { "PASS" } else { "FAIL" }.to_string(),
                                metric: Some(dev),
                                tol: Some(tol_coup),
                                detail: Some(format!(
                                    "sorted |main-chain|/|cross| deviation from book folded ({id})"
                                )),
                            },
                        ),
                        Err(e) => row.set_gate(
                            "book_matrix",
                            data::GateResult::skip(format!("book matrix comparison: {e}")),
                        ),
                    }
                }
            }
        }
    };

    // 4. return_loss (soft, 1.0) — worst passband RL vs spec.rl
    let tol_rl = tols.get("return_loss_db").copied().unwrap_or(1.0);
    let worst_rl = gates::worst_passband_return_loss_db(&r_t);
    let rl_dev = (worst_rl - rl).abs();
    row.set_gate(
        "return_loss",
        data::GateResult {
            status: if rl_dev < tol_rl { "PASS" } else { "FAIL" }.to_string(),
            metric: Some(rl_dev),
            tol: Some(tol_rl),
            detail: Some(format!(
                "worst passband RL {worst_rl:.2}dB vs spec {rl:.2}dB"
            )),
        },
    );
}

/// Case spec -> synthesis inputs (`order`, `rl`, normalized transmission zeros).
fn case_inputs(case: &Case) -> (usize, f64, Vec<f64>) {
    (
        gates::spec_usize(&case.spec, "order").expect("spec.order"),
        gates::spec_f64(&case.spec, "rl").expect("spec.rl"),
        gates::spec_f64_vec(&case.spec, "zeros").expect("spec.zeros"),
    )
}

fn bench_sweep_cases(c: &mut Criterion) {
    let [synthesis_stage, folding_stage, response_stage] = twoport_stages();
    let grid = FrequencyGrid::linspace(
        -TWOPORT_RESPONSE_BOUND,
        TWOPORT_RESPONSE_BOUND,
        TWOPORT_RESPONSE_POINTS,
    )
    .expect("contract response grid");

    for case in SWEEP_CASES {
        let label = format!("sweep/{}", case.name);

        let mut group = c.benchmark_group(format!("{synthesis_stage}/{label}"));
        apply_harness_policy(&mut group, "twoport_sweep");
        group.bench_function(&label, |b| {
            b.iter(|| build_design(case.order, case.rl, case.zeros))
        });
        group.finish();

        let design = build_design(case.order, case.rl, case.zeros);

        let mut group = c.benchmark_group(format!("{folding_stage}/{label}"));
        apply_harness_policy(&mut group, "twoport_sweep");
        group.bench_function(&label, |b| {
            b.iter(|| black_box(design.to_folded().unwrap()))
        });
        group.finish();

        let mut group = c.benchmark_group(format!("{response_stage}/{label}"));
        apply_harness_policy(&mut group, "twoport_sweep");
        group.bench_function(&label, |b| {
            b.iter(|| {
                black_box(
                    ResponseSolver
                        .evaluate_normalized(design.matrix(), &grid)
                        .unwrap(),
                )
            })
        });
        group.finish();
    }
}

criterion_group!(benches, bench_standard_cases, bench_sweep_cases);
criterion_main!(benches);
