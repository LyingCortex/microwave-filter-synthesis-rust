//! Native Rust (Criterion) multiplexer benchmark.
//!
//! Case set and policy come from the standard (`data::cases("multiplexer")`,
//! `data::harness(qualid)`, `data::stages("multiplexer")`); no case constant is
//! hardcoded here.
//!
//! ## Capability and honest stage reporting
//!
//! The contract's `multiplexer` stages are `poly` → `cm` → `star` → `response`
//! (`gates.json`), describing the model-paper flow implemented by fuqtgui
//! (Zhao 2026 star-junction synthesis: rational-function iteration, residue-matrix
//! CM, similarity transform to the star topology). mfs implements a *different*
//! algorithm — Macchiarella–Tamiazzo per-channel iterative perturbation
//! (`src/multiplexer.rs`) followed by a direct per-channel response evaluation.
//!
//! Two of the four stages have an mfs counterpart and are reported:
//!
//! * `poly` — the per-channel iterative synthesis (`MultiplexerSynthesis::
//!   synthesize`, which produces each channel's S-polynomial roots and coupling
//!   matrix). It is the same *phase* of the flow as the contract's `poly`, but a
//!   different algorithm, so its numbers are only comparable across mfs
//!   backends, not against fuqtgui.
//! * `response` — the multiplexer response sweep.
//!
//! `cm` and `star` (residue-matrix analytic CM; similarity transform to the
//! star junction) are **not produced** by mfs and are omitted from the report
//! entirely. Fabricating timings for them from the per-channel work would be a
//! false claim about the algorithm; every gate that requires one of them is
//! recorded as `SKIP`, never `PASS`.
//!
//! A case is *capability-matched* when its tags contain `type1` or `synthetic`
//! (non-resonant junction — what `JunctionType::TypeI` implements). Cases tagged
//! `type2` (resonant junction), `contiguous`, or `highorder` are reported as
//! `SKIP` with the deciding tag in `detail`; they are never silently dropped.
//!
//! Run: `cargo bench --bench multiplexer`

#[path = "../tests/bench_support/mod.rs"]
mod bench_support;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mfs::freq::BandPassMapping;
use mfs::multiplexer::{Multiplexer, MultiplexerSynthesis};
use mfs::prelude::FrequencyMapping;
use mfs::{ChannelSpec, JunctionType};

use bench_support::data::{self, Case};
use bench_support::gates::{
    self, apply_harness_policy, channel_inputs, spec_f64, MULTIPLEXER_RESPONSE_POINTS,
};

/// Stage mfs can actually time for this family (see module docs).
const MFS_MULTIPLEXER_STAGE: &str = "response";
/// Contract stages mfs has no counterpart for at all.
const UNAVAILABLE_STAGES: [&str; 2] = ["cm", "star"];

/// Why a case cannot be run by mfs's multiplexer, or `None` when it can.
fn capability_blocker(case: &Case) -> Option<String> {
    if case.has_tag("type2") {
        return Some(
            "type2: resonant junction (JunctionType::TypeII star node) is not what mfs's \
             multiplexer implements; fuqtgui's star-junction algorithm covers it"
                .to_string(),
        );
    }
    if case.has_tag("contiguous") && !case.has_tag("synthetic") {
        return Some(
            "contiguous: guard band = 0 breaks the per-channel loading approximation".to_string(),
        );
    }
    if case.has_tag("highorder") {
        return Some(
            "highorder: runs at --stages poly with warmup 0 / repeats 1 (see harness.json note); \
             a Criterion bench is the wrong harness for it"
                .to_string(),
        );
    }
    if case.has_tag("type1") || case.has_tag("synthetic") {
        // A matched case must be fully expressible as mfs channels; an extra
        // resonator-zero channel would be silently dropped otherwise.
        let extra = case.spec.get("extra_rz_channel").filter(|v| !v.is_null());
        if let Some(extra) = extra {
            return Some(format!(
                "spec declares extra_rz_channel={extra}, which mfs's ChannelSpec cannot express"
            ));
        }
        return None;
    }
    Some(format!("no capability tag matched (tags: {:?})", case.tags))
}

/// Case spec -> mfs channel specs.
///
/// `center`/`bandwidth` follow the standard's edge convention
/// (`center = √(f_low·f_high)`, `bw = f_high - f_low` — the same convention as
/// mfs `freq::edges_to_center_bandwidth`, under which Ω(f_low) = -1 and
/// Ω(f_high) = +1). Frequencies in the spec are MHz.
///
/// `tz` (MHz) are mapped into the channel's own normalized prototype domain via
/// `BandPassMapping`, which is what `ChannelSpec::zeros` expects; being
/// out-of-band, they land at |Ω| > 1.
fn to_channel_specs(case: &Case) -> Vec<ChannelSpec> {
    let channels = case
        .spec
        .get("channels")
        .and_then(serde_json::Value::as_array)
        .expect("multiplexer spec.channels");

    let mut specs = Vec::with_capacity(channels.len());
    // `channel_inputs` yields (order, rl_db, tz) straight from the standard.
    for (raw, (order, return_loss_db, tz)) in channels
        .iter()
        .zip(channel_inputs(&case.spec).expect("channel inputs"))
    {
        let f_low = spec_f64(raw, "f_low").expect("f_low") * 1e6;
        let f_high = spec_f64(raw, "f_high").expect("f_high") * 1e6;
        let center_hz = (f_low * f_high).sqrt();
        let bandwidth_hz = f_high - f_low;
        let mapping = BandPassMapping::new(center_hz, bandwidth_hz).expect("channel mapping");
        let zeros = tz
            .iter()
            .map(|hz| {
                mapping
                    .map_hz_to_normalized(*hz * 1e6)
                    .expect("transmission zero mapping")
                    .omega
            })
            .collect();
        specs.push(ChannelSpec {
            order,
            return_loss_db,
            zeros,
            center_hz,
            bandwidth_hz,
        });
    }
    specs
}

/// Junction type declared by the case spec, restricted to what mfs implements.
fn junction(case: &Case) -> JunctionType {
    match case
        .spec
        .get("junction")
        .and_then(serde_json::Value::as_str)
    {
        Some("nonresonant") | None => JunctionType::TypeI { n: 1.0, b0: 0.0 },
        Some(other) => panic!(
            "{}: junction {other:?} has no mfs mapping (capability filter should have skipped it)",
            case.qualid
        ),
    }
}

/// Synthesizes a multiplexer and returns the channel designs.
fn synthesize_mux(case: &Case) -> Vec<mfs::design::FilterDesign> {
    let specs = to_channel_specs(case);
    let harness = data::harness(&case.qualid);
    let mut synthesis = MultiplexerSynthesis::new()
        .channels(specs)
        .junction(junction(case));
    if let Some(max_iter) = harness.max_iter {
        synthesis = synthesis.max_iterations(max_iter);
    }
    synthesis
        .synthesize()
        .unwrap_or_else(|e| panic!("{}: multiplexer synthesis failed: {e}", case.qualid))
        .channels
}

/// Full span of the case's bands, in Hz, used for the response sweep.
fn response_span_hz(case: &Case) -> (f64, f64) {
    let channels = case
        .spec
        .get("channels")
        .and_then(serde_json::Value::as_array)
        .expect("multiplexer spec.channels");
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for raw in channels {
        let f_low = spec_f64(raw, "f_low").expect("f_low") * 1e6;
        let f_high = spec_f64(raw, "f_high").expect("f_high") * 1e6;
        low = low.min(f_low);
        high = high.max(f_high);
    }
    // 10 % margin so the response is not clipped at the band edges.
    let margin = 0.1 * (high - low);
    (low - margin, high + margin)
}

fn bench_multiplexer(c: &mut Criterion) {
    let mut report = data::Report::new("mfs", "rust", &data::platform());
    let stages = data::stages("multiplexer").expect("standard defines multiplexer stages");
    for stage in &UNAVAILABLE_STAGES {
        assert!(
            stages.contains(&stage.to_string()),
            "standard multiplexer stages changed: {stages:?}"
        );
    }
    assert!(
        stages.contains(&MFS_MULTIPLEXER_STAGE.to_string()),
        "standard multiplexer stages changed: {stages:?}"
    );
    assert!(
        stages.contains(&"poly".to_string()),
        "standard multiplexer stages changed: {stages:?}"
    );

    let cases = data::cases("multiplexer").expect("standard multiplexer cases");

    for case in &cases {
        let row = report.add_case(&case.qualid);

        if let Some(blocker) = capability_blocker(case) {
            println!("SKIP {}: {blocker}", case.qualid);
            row.set_gate(
                "convergence",
                data::GateResult::skip(format!("{blocker}; no poly stage")),
            );
            row.set_gate(
                "unitarity",
                data::GateResult::skip(format!(
                    "{blocker}; mfs produces none of {}",
                    UNAVAILABLE_STAGES.join("/")
                )),
            );
            continue;
        }

        // mfs has a `poly` analogue (the iterative per-channel synthesis), so
        // `convergence` is evaluated on it; `rz_distribution` needs the contract's
        // `poly` outputs (per-channel rational-function poles/zeros), which mfs's
        // algorithm does not expose in that form. Gates needing `cm`/`star` stay
        // SKIP - mfs produces neither.
        row.set_gate(
            "rz_distribution",
            data::GateResult::skip(
                "mfs's `poly` analogue is iterative per-channel synthesis and exposes \
                 coupling matrices, not the contract's rational-function pole/zero set",
            ),
        );
        for gate in ["matrix_size", "unitarity", "equiripple", "s11_match"] {
            row.set_gate(
                gate,
                data::GateResult::skip(format!(
                    "mfs has no {} stage",
                    UNAVAILABLE_STAGES.join("/")
                )),
            );
        }

        let harness = data::harness(&case.qualid);
        let (start_hz, stop_hz) = response_span_hz(case);

        // `poly` analogue: mfs's per-channel iterative synthesis, run through
        // the FULL Macchiarella–Tamiazzo spectral path (B2/B3). The contract's
        // `poly` is fuqtgui's rational-function iteration, so this is the same
        // phase of the flow reached by a different algorithm — the timing is
        // comparable only across mfs backends, and the report says so. The
        // full path is what actually reaches the contract's convergence gate
        // (the simplified perturbation stalls at a bounded residual by design).
        let specs = to_channel_specs(case);
        let mut synthesis = MultiplexerSynthesis::new()
            .channels(specs.clone())
            .junction(junction(case))
            .full_iteration(true);
        if let Some(max_iter) = harness.max_iter {
            synthesis = synthesis.max_iterations(max_iter);
        }
        let mut group = c.benchmark_group(format!("poly/{}", case.qualid));
        apply_harness_policy(&mut group, &case.qualid);
        group.bench_function(&case.qualid, |b| {
            b.iter(|| black_box(synthesis.clone().synthesize().unwrap()))
        });
        group.finish();
        // The contract declares no `convergence` tolerance, so the gate compares
        // against mfs's own criterion; the budget is set from the harness.
        let converged = synthesis.clone().synthesize().unwrap();
        row.set_gate("convergence", gates::convergence_gate(&converged));
        let poly_timing = data::time_call(
            || {
                let _ = synthesis.clone().synthesize();
            },
            harness.timing.warmup,
            harness.timing.repeats,
        );
        println!(
            "{}: poly median {:.3} ms (n={}) [mfs-specific algorithm]",
            case.qualid,
            poly_timing.median_ms,
            poly_timing.n.unwrap_or(0)
        );
        row.set_stage("poly", poly_timing);

        let channels = synthesize_mux(case);
        let mux = Multiplexer::new(channels).expect("multiplexer from channels");

        let mut group = c.benchmark_group(format!("{MFS_MULTIPLEXER_STAGE}/{}", case.qualid));
        apply_harness_policy(&mut group, &case.qualid);
        group.bench_function(&case.qualid, |b| {
            b.iter(|| {
                black_box(
                    mux.response(start_hz, stop_hz, MULTIPLEXER_RESPONSE_POINTS)
                        .unwrap(),
                )
            })
        });
        group.finish();

        let timing = data::time_call(
            || {
                let _ = mux.response(start_hz, stop_hz, MULTIPLEXER_RESPONSE_POINTS);
            },
            harness.timing.warmup,
            harness.timing.repeats,
        );
        println!(
            "{}: {} median {:.3} ms (n={})",
            case.qualid,
            MFS_MULTIPLEXER_STAGE,
            timing.median_ms,
            timing.n.unwrap_or(0)
        );
        row.set_stage(MFS_MULTIPLEXER_STAGE, timing);
        // `max_iter` from the harness is honoured by `synthesize_mux`; surface it
        // so a reader can tell which policy produced the timing.
        if let Some(max_iter) = harness.max_iter {
            println!("{}: max_iter={max_iter} (harness.json)", case.qualid);
        }
    }

    report
        .write(std::path::Path::new("target/bench-multiplexer-report.json"))
        .expect("write multiplexer bench report");
}

criterion_group!(benches, bench_multiplexer);
criterion_main!(benches);
