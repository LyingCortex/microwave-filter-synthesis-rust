//! mfs-side evaluation of the standard's gates, on mfs domain objects.
//!
//! Each helper implements exactly one gate from `gates.json` (`metric` field)
//! so the benches and the conformance test measure the same quantity the
//! contract defines. Nothing here compensates for missing capability: epic
//! stages a backend does not produce are reported as `SKIP` by the callers.

use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::response::SParameterResponse;

use super::data::{self, Result};

/// Grid used by the standard's `twoport` response stage (201 points, [-1.3, 1.3]).
pub const TWOPORT_RESPONSE_POINTS: usize = 201;
/// Grid bounds used by the standard's `twoport` response stage.
pub const TWOPORT_RESPONSE_BOUND: f64 = 1.3;
/// Grid used by the standard's `multiplexer` response stage (401 points, [-1.3, 1.3]).
pub const MULTIPLEXER_RESPONSE_POINTS: usize = 401;
/// Grid bounds used by the standard's `multiplexer` response stage.
pub const MULTIPLEXER_RESPONSE_BOUND: f64 = 1.3;

/// Normalized evaluation grid of the contract (`[-bound, bound]`, `points` samples).
pub fn normalized_grid(bound: f64, points: usize) -> Result<FrequencyGrid> {
    Ok(FrequencyGrid::linspace(-bound, bound, points)?)
}

// ─── spec access ─────────────────────────────────────────────────────────────

/// Reads a `f64` field of a case spec.
pub fn spec_f64(spec: &serde_json::Value, key: &str) -> Result<f64> {
    spec.get(key)
        .and_then(serde_json::Value::as_f64)
        .ok_or_else(|| format!("case spec is missing numeric field {key:?}").into())
}

/// Reads a `usize` field of a case spec.
pub fn spec_usize(spec: &serde_json::Value, key: &str) -> Result<usize> {
    spec.get(key)
        .and_then(serde_json::Value::as_u64)
        .map(|v| v as usize)
        .ok_or_else(|| format!("case spec is missing integer field {key:?}").into())
}

/// Reads a numeric array field of a case spec.
pub fn spec_f64_vec(spec: &serde_json::Value, key: &str) -> Result<Vec<f64>> {
    spec.get(key)
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| format!("case spec is missing array field {key:?}"))?
        .iter()
        .map(|v| {
            v.as_f64().ok_or_else(|| {
                format!("case spec field {key:?} contains a non-numeric entry").into()
            })
        })
        .collect()
}

/// Reads a nested `spec.channels` entry as `(order, rl_db, tz)` — the mfs
/// `ChannelSpec` inputs that come straight from the standard.
pub fn channel_inputs(spec: &serde_json::Value) -> Result<Vec<(usize, f64, Vec<f64>)>> {
    spec.get("channels")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| "multiplexer case spec has no `channels` array".to_string())?
        .iter()
        .map(|ch| {
            Ok((
                spec_usize(ch, "n_rz")?,
                spec_f64(ch, "rl_db")?,
                spec_f64_vec(ch, "tz")?,
            ))
        })
        .collect()
}

// ─── gates ───────────────────────────────────────────────────────────────────

/// `response_invariance` — max S-parameter deviation between two responses,
/// `max|S11_a - S11_b|` or the sign-symmetric `min(max|S21_a-S21_b|, max|S21_a+S21_b|)`.
///
/// S21 is compared up to a global sign: an odd number of resonator-order
/// reversals in a similarity rotation flips it, which is physical.
pub fn response_invariance_metric(a: &SParameterResponse, b: &SParameterResponse) -> Result<f64> {
    if a.samples.len() != b.samples.len() {
        return Err(format!(
            "response_invariance needs equal-length responses: {} vs {}",
            a.samples.len(),
            b.samples.len()
        )
        .into());
    }
    let mut max_s11 = 0.0_f64;
    let mut max_s21 = 0.0_f64;
    let mut max_s21_neg = 0.0_f64;
    for (x, y) in a.samples.iter().zip(b.samples.iter()) {
        max_s11 = max_s11.max((x.s11_mag() - y.s11_mag()).abs());
        max_s21 = max_s21.max((x.s21_mag() - y.s21_mag()).abs());
        max_s21_neg = max_s21_neg.max((x.s21_mag() + y.s21_mag()).abs());
    }
    Ok(max_s11.max(max_s21.min(max_s21_neg)))
}

/// Criterion refuses sample sizes below 10; the standard's `repeats` is a
/// contract figure for the report harness, so it is honoured as a floor.
pub const CRITERION_MIN_SAMPLES: usize = 10;

/// Applies the standard's per-case timing policy (`harness.json`) to a Criterion
/// group: warm-up time and sample count.
pub fn apply_harness_policy(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    qualid: &str,
) {
    let harness = data::harness(qualid);
    group.warm_up_time(std::time::Duration::from_secs(
        harness.timing.warmup.max(1) as u64
    ));
    group.sample_size(harness.timing.repeats.max(CRITERION_MIN_SAMPLES));
}

/// `gate("response_invariance")` tolerance from the standard (`1e-10`).
pub fn response_invariance_tol() -> f64 {
    data::gate_tols()
        .get("response_invariance")
        .copied()
        .unwrap_or(1e-10)
}

/// `folded_mask` (hard, tol 0.0) — largest |M[i][j]| over the upper triangle
/// outside the folded mask, using the contract's default `"either"` convention:
/// the union of the source-end (`canonical`) and load-end (`mirrored`) masks.
///
/// Normative source: `data/gates.json:gates.folded_mask` (semantics + evidence)
/// and `CONTRACT.md`. Cameron's folded form is not unique (§8.3), so the gate
/// accepts both mirror images; mfs's `to_folded` emits the mirrored one for even
/// orders, the book's printed folded matrices the canonical one, fuqtgui's
/// computed folded matrices the canonical one.
pub fn folded_mask_violations(matrix: &CouplingMatrix) -> f64 {
    let side = matrix.side();
    let order = matrix.order();
    let mut worst = 0.0_f64;
    for row in 0..side {
        for col in (row + 1)..side {
            let mirrored = (mirror_index(col, order), mirror_index(row, order));
            if folded_mask_keeps(row, col, order)
                || folded_mask_keeps(mirrored.0, mirrored.1, order)
            {
                continue;
            }
            worst = worst.max(matrix.at(row, col).unwrap_or_default().abs());
        }
    }
    worst
}

/// Centre reflection of a matrix index: `i -> N+1-i`, with source (0) and load
/// (N+1) swapped. This is the mirror relating the two legal folded conventions.
fn mirror_index(index: usize, order: usize) -> usize {
    order + 1 - index
}

/// The folded matrix reflected through the filter's centre, i.e. the same
/// structure expressed in the other legal folded convention.
pub fn reversal_mirror(matrix: &CouplingMatrix) -> CouplingMatrix {
    let side = matrix.side();
    let order = matrix.order();
    let mut flat = vec![0.0_f64; side * side];
    for row in 0..side {
        for col in 0..side {
            flat[mirror_index(row, order) * side + mirror_index(col, order)] =
                matrix.at(row, col).unwrap_or_default();
        }
    }
    CouplingMatrix::new(order, flat).expect("mirrored matrix keeps the same order")
}

/// One convention's positional predicate of the folded mask (0 = source,
/// `order+1` = load), transcribed from the contract's `folded_mask` semantics:
/// the main chain `(i, i+1)`, the symmetric cross-diagonal `i+j = N+1+shift`,
/// the asymmetric cross-diagonal `i+j = N+shift` with `shift = N%2`, and the
/// multiple I/O couplings `(0, j)` / `(i, N+1)` that those diagonals contain.
/// Self-couplings (the main diagonal) are never checked: the gate looks at the
/// upper triangle with `j > i` only.
fn folded_mask_keeps(row: usize, col: usize, order: usize) -> bool {
    if col == row + 1 {
        return col <= order + 1;
    }
    let shift = order % 2;
    let load = order + 1;
    col > row + 1 && col <= load && (row + col == order + 1 + shift || row + col == order + shift)
}

/// `convergence` gate (hard) for the `multiplexer` family: the outer iteration
/// must stop on its own criterion rather than run out of iterations.
///
/// `gates.json` names the metric ("converged flag / final root displacement")
/// but declares no `tol_key`, so there is no contract tolerance to compare
/// against. The gate reads the synthesis's own `converged` flag: `true` when
/// the loop stops on its convergence criterion — the strict root-displacement
/// tolerance, or, for the full Macchiarella–Tamiazzo spectral path, the
/// high-degree stagnation accept (the displacement noise floor exceeds the
/// strict tolerance while the physical fixed point has converged) — and
/// `false` when the iteration budget is exhausted. An exhausted budget is a
/// non-convergence whatever the residual. The final root displacement is
/// reported as the metric for transparency.
pub fn convergence_gate(result: &mfs::multiplexer::MultiplexerSynthesisResult) -> data::GateResult {
    let converged = result.converged;
    data::GateResult {
        status: if converged { "PASS" } else { "FAIL" }.to_string(),
        metric: Some(result.convergence_error),
        tol: None,
        detail: Some(format!(
            "mfs multiplexer synthesis: {} of {} iteration(s), converged on its \
             own criterion: {}, final root displacement {:.3e}",
            result.iterations,
            result.max_iterations,
            converged,
            result.convergence_error
        )),
    }
}

/// `return_loss` — worst passband return loss over the normalized passband
/// (`|Ω| <= 1`), reported as a positive dB figure.
pub fn worst_passband_return_loss_db(response: &SParameterResponse) -> f64 {
    response
        .samples
        .iter()
        .filter(|s| s.normalized_omega.abs() <= 1.0)
        .map(|s| -s.s11_db())
        .fold(f64::INFINITY, f64::min)
}

/// `book_matrix` — deviation between two matrices' sorted main-chain and sorted
/// cross-coupling magnitude sets (the folded form is defined up to a plane
/// rotation, so element-wise equality does not hold).
pub fn coupling_set_deviation(a: &CouplingMatrix, b: &CouplingMatrix) -> Result<f64> {
    if a.order() != b.order() {
        return Err(format!(
            "coupling_set_deviation needs equal orders: {} vs {}",
            a.order(),
            b.order()
        )
        .into());
    }
    let (main_a, cross_a) = coupling_sets(a);
    let (main_b, cross_b) = coupling_sets(b);
    Ok(sorted_set_deviation(&main_a, &main_b).max(sorted_set_deviation(&cross_a, &cross_b)))
}

fn coupling_sets(matrix: &CouplingMatrix) -> (Vec<f64>, Vec<f64>) {
    let side = matrix.side();
    let mut main = Vec::new();
    let mut cross = Vec::new();
    for row in 0..side {
        for col in (row + 1)..side {
            let value = matrix.at(row, col).unwrap_or_default().abs();
            if col == row + 1 || (row == 0 && col == side - 1) {
                main.push(value);
            } else {
                cross.push(value);
            }
        }
    }
    (main, cross)
}

fn sorted_set_deviation(a: &[f64], b: &[f64]) -> f64 {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    a.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    b.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let len = a.len().max(b.len());
    let mut worst = 0.0_f64;
    for index in 0..len {
        let x = a.get(index).copied().unwrap_or(0.0);
        let y = b.get(index).copied().unwrap_or(0.0);
        worst = worst.max((x - y).abs());
    }
    worst
}
