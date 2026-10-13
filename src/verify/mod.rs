use crate::error::Result;
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::response::{ResponseSample, SParameterResponse};
use crate::spec::FilterSpec;

/// Tolerance profile used when comparing two sampled responses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseTolerance {
    /// Maximum allowed difference between physical frequency samples in Hz.
    pub frequency_hz: f64,
    /// Maximum allowed difference between normalized frequency samples.
    pub normalized_omega: f64,
    /// Maximum allowed difference in `|S11|`.
    pub s11_magnitude: f64,
    /// Maximum allowed difference in `|S21|`.
    pub s21_magnitude: f64,
    /// Maximum allowed difference in group delay.
    pub group_delay: f64,
}

impl ResponseTolerance {
    /// A strict but practical default for regression-style response checks.
    pub fn strict() -> Self {
        Self {
            frequency_hz: 1e-9,
            normalized_omega: 1e-12,
            s11_magnitude: 1e-8,
            s21_magnitude: 1e-8,
            group_delay: 1e-7,
        }
    }
}

impl Default for ResponseTolerance {
    fn default() -> Self {
        Self::strict()
    }
}

/// Tolerance profile used when interpreting whether a matrix entry is effectively zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatrixPatternTolerance {
    /// Maximum absolute value still treated as structurally zero.
    pub structural_zero: f64,
}

impl MatrixPatternTolerance {
    /// A practical default for topology-pattern checks on transformed matrices.
    pub fn strict() -> Self {
        Self {
            structural_zero: 1e-6,
        }
    }
}

impl Default for MatrixPatternTolerance {
    fn default() -> Self {
        Self::strict()
    }
}

/// Summary of how closely two responses match.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseComparison {
    /// Number of samples that were compared.
    pub samples_compared: usize,
    /// Largest observed frequency-axis mismatch in Hz.
    pub max_frequency_hz_deviation: f64,
    /// Largest observed normalized-frequency mismatch.
    pub max_normalized_omega_deviation: f64,
    /// Largest observed deviation in `|S11|`.
    pub max_s11_magnitude_deviation: f64,
    /// Largest observed deviation in `|S21|`.
    pub max_s21_magnitude_deviation: f64,
    /// Largest observed group-delay deviation.
    pub max_group_delay_deviation: f64,
    /// Whether group delay was actually measured at **any** sample pair. When
    /// both responses leave `group_delay = None` everywhere (the pole-expansion
    /// fast path), the group-delay deviation is a silent no-op (0.0) and this
    /// stays `false`, so a tolerance that demands a finite group delay fails
    /// rather than passing on an unmeasured value (review P2-5). Group delay is
    /// legitimately `None` at individual samples (transmission zeros, where it
    /// diverges); a partial measurement still counts as compared.
    pub group_delay_compared: bool,
}

/// The result of a response-invariance check. Skipped is kept distinct from
/// Pass so a workflow that never ran the electrical check is not reported as a
/// passing one (review P3 #2: `passes()` used to conflate "not checked" with
/// "passed").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckOutcome {
    /// The response check was not performed.
    Skipped,
    /// The response check ran and met the requested tolerance.
    Pass,
    /// The response check ran and exceeded the requested tolerance.
    Fail,
}

impl std::fmt::Display for CheckOutcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            CheckOutcome::Skipped => "skipped",
            CheckOutcome::Pass => "pass",
            CheckOutcome::Fail => "fail",
        };
        f.write_str(s)
    }
}

/// Optional response-check summary that can be embedded in transform reports.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseCheckReport {
    /// Optional comparison details when a response sweep was performed.
    pub comparison: Option<ResponseComparison>,
    /// Whether the optional comparison passed the requested tolerance.
    pub invariant: Option<bool>,
}

impl ResponseCheckReport {
    /// Creates an empty response-check report for workflows that skipped electrical checking.
    pub fn skipped() -> Self {
        Self {
            comparison: None,
            invariant: None,
        }
    }

    /// Creates a populated response-check report from a comparison result and tolerance.
    pub fn from_comparison(comparison: ResponseComparison, tolerance: ResponseTolerance) -> Self {
        let invariant = comparison.passes(tolerance);
        Self {
            comparison: Some(comparison),
            invariant: Some(invariant),
        }
    }

    /// Returns whether the response check passed, or `true` when no response check was requested.
    pub fn passes(&self) -> bool {
        self.invariant.unwrap_or(true)
    }

    /// Returns the three-state outcome of the response check, keeping a skipped
    /// check distinct from a passing one (review P3 #2).
    pub fn outcome(&self) -> CheckOutcome {
        match self.invariant {
            Some(true) => CheckOutcome::Pass,
            Some(false) => CheckOutcome::Fail,
            None => CheckOutcome::Skipped,
        }
    }
}

/// Summary of structural checks performed for an extracted section workflow.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SectionVerificationReport {
    /// Number of entries expected to be effectively zero.
    pub zero_checks: usize,
    /// Number of entries expected to stay meaningfully non-zero.
    pub nonzero_checks: usize,
    /// Number of zero checks that failed.
    pub failed_zero_checks: usize,
    /// Number of non-zero checks that failed.
    pub failed_nonzero_checks: usize,
}

impl SectionVerificationReport {
    /// Returns whether every expected section-pattern check passed.
    pub fn passes(&self) -> bool {
        self.failed_zero_checks == 0 && self.failed_nonzero_checks == 0
    }
}

impl ResponseComparison {
    /// Returns whether all tracked deviations fall within the supplied tolerance.
    pub fn passes(&self, tolerance: ResponseTolerance) -> bool {
        // P2-5: when the tolerance constrains group delay (finite bound), the
        // comparison must have actually measured it. A response path that left
        // `group_delay = None` makes the deviation a silent no-op (0.0), which
        // must not satisfy a finite group-delay bound.
        if tolerance.group_delay.is_finite() && !self.group_delay_compared {
            return false;
        }
        self.max_frequency_hz_deviation <= tolerance.frequency_hz
            && self.max_normalized_omega_deviation <= tolerance.normalized_omega
            && self.max_s11_magnitude_deviation <= tolerance.s11_magnitude
            && self.max_s21_magnitude_deviation <= tolerance.s21_magnitude
            && self.max_group_delay_deviation <= tolerance.group_delay
    }
}

/// Compares two sampled responses and reports their maximum deviations.
pub fn compare_responses(
    lhs: &SParameterResponse,
    rhs: &SParameterResponse,
) -> Result<ResponseComparison> {
    let expected = lhs.samples.len();
    let actual = rhs.samples.len();
    if expected != actual {
        return Err(crate::error::MfsError::DimensionMismatch { expected, actual });
    }

    let mut comparison = ResponseComparison {
        samples_compared: expected,
        max_frequency_hz_deviation: 0.0,
        max_normalized_omega_deviation: 0.0,
        max_s11_magnitude_deviation: 0.0,
        max_s21_magnitude_deviation: 0.0,
        max_group_delay_deviation: 0.0,
        // P2-5: set to true once at least one sample pair has a measurable
        // group delay on both sides. Stays `false` when group delay is
        // absent everywhere (the pole-expansion fast path leaves every
        // sample `None`), which makes the group-delay deviation a silent
        // no-op that `passes()` must not treat as 0.0 (review P2-5).
        group_delay_compared: false,
    };

    for (left, right) in lhs.samples.iter().zip(rhs.samples.iter()) {
        comparison.max_frequency_hz_deviation = comparison
            .max_frequency_hz_deviation
            .max((left.frequency_hz - right.frequency_hz).abs());
        comparison.max_normalized_omega_deviation = comparison
            .max_normalized_omega_deviation
            .max((left.normalized_omega - right.normalized_omega).abs());
        comparison.max_s11_magnitude_deviation = comparison.max_s11_magnitude_deviation.max(
            (magnitude(left.s11_re, left.s11_im) - magnitude(right.s11_re, right.s11_im)).abs(),
        );
        comparison.max_s21_magnitude_deviation = comparison.max_s21_magnitude_deviation.max(
            (magnitude(left.s21_re, left.s21_im) - magnitude(right.s21_re, right.s21_im)).abs(),
        );
        // Group delay is only compared where both samples have it. It is
        // legitimately `None` at transmission zeros (where it diverges) and on
        // the pole-expansion fast path (which leaves it `None` everywhere).
        // Compare at every available pair and remember that at least one was
        // measured, so an entirely-unmeasured comparison is surfaced by
        // `passes()` rather than passing on a silent 0.0 (review P2-5).
        if let (Some(l), Some(r)) = (left.group_delay, right.group_delay) {
            comparison.group_delay_compared = true;
            comparison.max_group_delay_deviation =
                comparison.max_group_delay_deviation.max((l - r).abs());
        }
    }

    Ok(comparison)
}


// ─── Absolute spec checks (review §D.3) ──────────────────────────────────────
/// Margin (dB) applied to the spec return-loss gate. The worst passband RL
/// sits at the band edges (the ripple), where the solver matches the spec to
/// ~1e-4 dB, so a strict `>=` comparison is a coin flip at the floating-point
/// boundary; a small margin avoids failing a design that meets its edge RL.
pub const RETURN_LOSS_GATE_MARGIN_DB: f64 = 0.1;
/// Maximum allowed `|S21|` (dB) inside the normalized passband. A lossless
/// matched design sits at 0 dB modulo ripple; +0.5 dB absorbs sampling
/// placement near the band edges and light `unloaded_q` loss. Heuristic gate
/// threshold, not physical law.
pub const PASSBAND_MAX_IL_DB: f64 = 0.5;
/// Minimum allowed `|S21|` (dB) inside the normalized passband. −1 dB is
/// below the ripple of any realistic Chebyshev passband and far above the
/// depth of a genuine in-band null (which would indicate a misplaced zero).
/// Heuristic gate threshold, not physical law.
pub const PASSBAND_MIN_IL_DB: f64 = -1.0;
/// `|S21|` (dB) a prescribed transmission zero must dip to (i.e. ≤ this)
/// within its search window. −40 dB is far below any passband ripple and
/// far above the solver's numerical floor, so a missing or misplaced zero
/// fails by many tens of dB while a true zero passes with margin. Heuristic
/// gate threshold, not physical law.
pub const TRANSMISSION_ZERO_DEPTH_DB: f64 = -40.0;
/// Half-width of the search window around a prescribed transmission zero, in
/// grid steps (the zero's dip must occur within ± this many steps of it).
pub const TRANSMISSION_ZERO_WINDOW_STEPS: usize = 2;

/// The result of one absolute-spec gate (e.g. passband return loss).
#[derive(Debug, Clone, PartialEq)]
pub struct SpecGateResult {
    /// Stable gate name (e.g. `passband_return_loss`, `transmission_zero_depth_1`).
    pub gate: String,
    /// Three-state outcome; reuses [`CheckOutcome`] so "data unavailable"
    /// stays `Skipped`, never a silent `Pass` (review P3 #2 / §D.3).
    pub outcome: CheckOutcome,
    /// Key measured scalar in dB that the gate judged; `None` when skipped.
    /// For `passband_insertion_loss` this is the violated side, or the floor
    /// value when the gate passes.
    pub measured: Option<f64>,
    /// Threshold applied to `measured` (direction depends on the gate);
    /// `None` when skipped.
    pub threshold: Option<f64>,
    /// Human-readable context: both measured values where a gate checks two,
    /// or the reason the gate was skipped.
    pub detail: String,
}

/// Report of the absolute spec checks performed by [`verify_against_spec`].
#[derive(Debug, Clone, PartialEq)]
pub struct SpecCheckReport {
    /// One entry per gate, in the order the gates run.
    pub gates: Vec<SpecGateResult>,
}

impl SpecCheckReport {
    /// Aggregates the gate outcomes: any `Fail` → `Fail`; every gate
    /// `Skipped` → `Skipped` (no data was available to check against);
    /// otherwise `Pass`.
    pub fn outcome(&self) -> CheckOutcome {
        if self.gates.iter().any(|g| g.outcome == CheckOutcome::Fail) {
            return CheckOutcome::Fail;
        }
        if self.gates.iter().all(|g| g.outcome == CheckOutcome::Skipped) {
            return CheckOutcome::Skipped;
        }
        CheckOutcome::Pass
    }

    /// Returns `true` when no gate failed (a fully-skipped report is `true`,
    /// mirroring [`ResponseCheckReport::passes`]). Callers that must keep
    /// "never checked" distinct from "passed" should prefer
    /// [`Self::outcome`] (review P3 #2).
    pub fn passes(&self) -> bool {
        !self.gates.iter().any(|g| g.outcome == CheckOutcome::Fail)
    }

    /// Looks up the result of a gate by name.
    pub fn gate(&self, name: &str) -> Option<&SpecGateResult> {
        self.gates.iter().find(|g| g.gate == name)
    }

    /// One reportable line per gate — `name = outcome (detail)` — the same
    /// row shape the terminal/markdown renderers use for
    /// [`CheckOutcome`]-backed checks, so callers can print the new check
    /// type with the existing rendering conventions.
    pub fn format_rows(&self) -> String {
        self.gates
            .iter()
            .map(|g| format!("{:<28} = {:<7}  {}", g.gate, g.outcome, g.detail))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

/// Checks a sampled response **against a [`FilterSpec`]** — the absolute
/// gates that [`compare_responses`] (two responses agree) never covers
/// (review §D.3): passband return loss, passband insertion loss,
/// transmission-zero depth, and (by convention, skipped) stopband rejection.
///
/// The response is expected in **normalized prototype coordinates**: the
/// passband is `|ω| ≤ 1` (the convention `FilterSpec` stores; it carries no
/// explicit edges) and the spec's transmission zeros are already normalized.
/// Physical-grid responses work as long as `normalized_omega` is populated,
/// which the response solver always does.
///
/// Each gate returns a three-state [`CheckOutcome`]: a gate whose data is
/// unavailable (no passband samples, zero at infinity, empty response, ...)
/// is `Skipped`, never a silent `Pass`. The stopband gate is always
/// `Skipped`: `FilterSpec` carries no stopband requirement in normalized
/// coordinates (its optional out-of-band windows are physical Hz and would
/// need a band mapping that is unavailable here), so there is nothing to
/// gate on without inventing a requirement.
///
/// ```rust
/// # use mfs::prelude::*;
/// let design = FilterDesign::prototype(4, 20.0)
///     .zeros([-2.0, 2.0])
///     .synthesize()?;
/// let response = design.response_normalized(-3.0, 3.0, 601)?;
/// let report = verify_against_spec(design.spec(), &response);
/// assert!(report.passes());
/// # Ok::<(), MfsError>(())
/// ```
///
/// A report whose gates `Fail` is a *result*, not an error: it carries the
/// measured value and the threshold so the caller can decide (tighten the
/// design, relax the spec, or reject it). `Ok(())`/no-error does not mean
/// the design met the spec — check [`SpecCheckReport::passes`].
///
/// ```rust
/// # use mfs::prelude::*;
/// let design = FilterDesign::prototype(4, 20.0)
///     .zeros([-2.0, 2.0])
///     .synthesize()?;
/// let response = design.response_normalized(-3.0, 3.0, 601)?;
/// // A tightened 33 dB RL spec the ~20 dB design cannot meet.
/// let tightened = FilterSpec::new(4, 33.0)?.with_normalized_transmission_zeros([-2.0, 2.0]);
/// let report = verify_against_spec(&tightened, &response);
/// assert!(!report.passes());
/// let rl = report.gate("passband_return_loss").expect("RL gate present");
/// assert_eq!(rl.outcome, CheckOutcome::Fail);
/// // The failure is actionable: measured ~20 dB vs the 33 dB requirement.
/// assert_eq!(rl.threshold, Some(33.0));
/// # Ok::<(), MfsError>(())
/// ```
pub fn verify_against_spec(spec: &FilterSpec, response: &SParameterResponse) -> SpecCheckReport {
    let samples = &response.samples;
    let passband: Vec<&ResponseSample> = samples
        .iter()
        .filter(|s| s.normalized_omega.abs() <= 1.0)
        .collect();

    let mut gates = Vec::new();

    // Gate 1: passband return loss against the spec's own `return_loss_db`.
    if passband.is_empty() {
        gates.push(skipped_gate(
            "passband_return_loss",
            "no samples inside the normalized passband |ω| ≤ 1",
        ));
    } else {
        // RL = −20·log10(|S11|); the spec constrains the *worst* passband
        // point (|S11| is largest at the band edges, where it equals the
        // ripple), so the gate compares the min passband RL against the spec
        // with a boundary margin (see RETURN_LOSS_GATE_MARGIN_DB). The center
        // passband can dip far deeper (reflection nulls); a max would pass a
        // design that misses its RL at the band edges.
        let min_rl_db = passband
            .iter()
            .map(|s| -s.s11_db())
            .fold(f64::INFINITY, f64::min);
        let required = spec.return_loss_db() - RETURN_LOSS_GATE_MARGIN_DB;
        gates.push(measured_gate(
            "passband_return_loss",
            min_rl_db,
            Some(spec.return_loss_db()),
            min_rl_db >= required,
            format!(
                "min passband RL {:.1} dB, spec requires {:.1} dB",
                min_rl_db,
                spec.return_loss_db()
            ),
        ));
    }

    // Gate 2: passband insertion loss — |S21| must stay near 0 dB.
    if passband.is_empty() {
        gates.push(skipped_gate(
            "passband_insertion_loss",
            "no samples inside the normalized passband |ω| ≤ 1",
        ));
    } else {
        let s21_dbs: Vec<f64> = passband
            .iter()
            .map(|s| s.s21_db())
            .collect();
        let max_db = s21_dbs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_db = s21_dbs.iter().cloned().fold(f64::INFINITY, f64::min);
        let too_hot = max_db > PASSBAND_MAX_IL_DB;
        let too_cold = min_db < PASSBAND_MIN_IL_DB;
        // `measured` reports the violated side when one failed, otherwise the
        // floor (the side a misplaced in-band zero would violate).
        let (measured, threshold) = if too_hot {
            (max_db, PASSBAND_MAX_IL_DB)
        } else {
            (min_db, PASSBAND_MIN_IL_DB)
        };
        gates.push(measured_gate(
            "passband_insertion_loss",
            measured,
            Some(threshold),
            !too_hot && !too_cold,
            format!(
                "passband |S21| {min_db:.1} .. {max_db:.1} dB (allowed {PASSBAND_MIN_IL_DB:.1} .. {PASSBAND_MAX_IL_DB:.1} dB)"
            ),
        ));
    }

    // Gate 3: transmission-zero depth, one sub-gate per prescribed zero.
    if spec.transmission_zeros.is_empty() {
        gates.push(skipped_gate(
            "transmission_zero_depth",
            "spec prescribes no transmission zeros",
        ));
    } else {
        for (index, zero) in spec.transmission_zeros.iter().enumerate() {
            let name = format!("transmission_zero_depth_{index}");
            if !zero.value.is_finite() {
                // Spec convention: f64::INFINITY is a zero at infinity — it
                // imposes no in-band dip requirement, so the sub-gate is
                // Skipped (not Pass: no dip was actually verified).
                gates.push(skipped_gate(&name, "zero at infinity — no dip requirement"));
            } else if samples.is_empty() {
                gates.push(skipped_gate(&name, "response has no samples"));
            } else {
                match zero_window_indices(samples, zero.value) {
                    Some(indices) if indices.is_empty() => {
                        gates.push(skipped_gate(
                            &name,
                            format!(
                                "no samples within ±{TRANSMISSION_ZERO_WINDOW_STEPS} grid steps of ω = {:.3}",
                                zero.value
                            ),
                        ));
                    }
                    Some(indices) => {
                        let depth_db = indices
                            .iter()
                            .map(|&i| samples[i].s21_db())
                            .fold(f64::INFINITY, f64::min);
                        gates.push(measured_gate(
                            &name,
                            depth_db,
                            Some(TRANSMISSION_ZERO_DEPTH_DB),
                            depth_db <= TRANSMISSION_ZERO_DEPTH_DB,
                            format!(
                                "min |S21| {depth_db:.1} dB within ±{TRANSMISSION_ZERO_WINDOW_STEPS} steps of ω = {:.3} (floor {TRANSMISSION_ZERO_DEPTH_DB:.1} dB)",
                                zero.value
                            ),
                        ));
                    }
                    None => unreachable!("zero_window_indices returns None only for an empty sample set"),
                }
            }
        }
    }

    // Gate 4: stopband rejection — always Skipped by design. FilterSpec has
    // no normalized-domain stopband requirement (see function docs), so
    // inventing a floor here would gate on a requirement the spec never made.
    gates.push(skipped_gate(
        "stopband_rejection",
        "spec carries no stopband requirement in normalized coordinates \
         (out-of-band windows are physical Hz and need a band mapping)",
    ));

    SpecCheckReport { gates }
}

/// Builds a Pass/Fail gate result (Skipped is built by [`skipped_gate`]).
fn measured_gate(
    gate: &str,
    measured: f64,
    threshold: Option<f64>,
    ok: bool,
    detail: String,
) -> SpecGateResult {
    SpecGateResult {
        gate: gate.to_string(),
        outcome: if ok { CheckOutcome::Pass } else { CheckOutcome::Fail },
        measured: Some(measured),
        threshold,
        detail,
    }
}

/// Builds a Skipped gate result: data unavailable, so the gate is reported
/// as `Skipped`, never as `Pass` (review P3 #2 / §D.3).
fn skipped_gate(gate: &str, detail: impl Into<String>) -> SpecGateResult {
    SpecGateResult {
        gate: gate.to_string(),
        outcome: CheckOutcome::Skipped,
        measured: None,
        threshold: None,
        detail: detail.into(),
    }
}

/// Indices of the samples within ±[`TRANSMISSION_ZERO_WINDOW_STEPS`] grid
/// steps of `zero_omega`, assuming a monotonic grid (the solver's
/// linspace grids are). For a degenerate/non-monotonic grid the window
/// falls back to the nearest sample and its two neighbors each side.
fn zero_window_indices(samples: &[ResponseSample], zero_omega: f64) -> Option<Vec<usize>> {
    let n = samples.len();
    if n == 0 {
        return None;
    }

    let spacing = if n >= 2 {
        (samples[n - 1].normalized_omega - samples[0].normalized_omega).abs() / (n - 1) as f64
    } else {
        0.0
    };
    if spacing > 0.0 && spacing.is_finite() {
        let radius = spacing * TRANSMISSION_ZERO_WINDOW_STEPS as f64;
        return Some((0..n).filter(|&i| {
            (samples[i].normalized_omega - zero_omega).abs() <= radius
        }).collect());
    }

    let mut nearest = 0usize;
    for i in 1..n {
        if (samples[i].normalized_omega - zero_omega).abs()
            < (samples[nearest].normalized_omega - zero_omega).abs()
        {
            nearest = i;
        }
    }
    let lo = nearest.saturating_sub(TRANSMISSION_ZERO_WINDOW_STEPS);
    let hi = (nearest + TRANSMISSION_ZERO_WINDOW_STEPS + 1).min(n);
    Some((lo..hi).collect())
}

/// Returns whether a matrix matches the expected folded source-side sparsity pattern.
pub fn matches_folded_pattern(matrix: &CouplingMatrix, tolerance: MatrixPatternTolerance) -> bool {
    for col in 2..=matrix.order() {
        if !is_structural_zero(matrix.at(0, col).unwrap_or_default(), tolerance) {
            return false;
        }
    }

    for row in 0..matrix.side() {
        for col in 0..matrix.side() {
            let lhs = matrix.at(row, col).unwrap_or_default();
            let rhs = matrix.at(col, row).unwrap_or_default();
            if (lhs - rhs).abs() > tolerance.structural_zero {
                return false;
            }
        }
    }

    true
}

/// Returns whether a matrix matches the current arrow-reduction sparsity pattern.
pub fn matches_arrow_pattern(matrix: &CouplingMatrix, tolerance: MatrixPatternTolerance) -> bool {
    for row in 0..matrix.order().saturating_sub(1) {
        for col in (row + 2)..=matrix.order() {
            if !is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
                return false;
            }
        }
    }

    true
}

/// Returns whether a matrix matches the expected sparsity pattern for a topology.
pub fn matches_topology_pattern(
    matrix: &CouplingMatrix,
    topology: MatrixTopology,
    tolerance: MatrixPatternTolerance,
) -> bool {
    match topology {
        MatrixTopology::Transversal => true,
        MatrixTopology::Folded => matches_folded_pattern(matrix, tolerance),
        // The current wheel reduction follows the same sparsity target as arrow.
        #[allow(deprecated)]
        MatrixTopology::Arrow | MatrixTopology::Wheel => matches_arrow_pattern(matrix, tolerance),
    }
}

/// Verifies the current structural expectations for a triplet extraction.
pub fn verify_triplet_extraction(
    matrix: &CouplingMatrix,
    center_resonator: usize,
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    if center_resonator < 2 || center_resonator + 2 > matrix.order() {
        return Err(crate::error::MfsError::PreconditionViolation(format!(
            "triplet verification requires center in [2, {}), got {center_resonator}",
            matrix.order()
        )));
    }

    // The triplet extraction should zero out the coupling that bypasses the
    // extracted section. The cross-coupling at (center-1, center+1) is not
    // universally nonzero — it depends on the polynomial set — so we only
    // verify the structural zero.
    verify_required_entries(
        matrix,
        &[(center_resonator + 1, center_resonator + 3)],
        &[],
        tolerance,
    )
}

/// Verifies the current structural expectations for a quadruplet extraction.
pub fn verify_quadruplet_extraction(
    matrix: &CouplingMatrix,
    position: usize,
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    if position < 2 || position + 2 > matrix.order() {
        return Err(crate::error::MfsError::PreconditionViolation(format!(
            "quadruplet verification requires position in [2, {}), got {position}",
            matrix.order().saturating_sub(1)
        )));
    }

    verify_required_entries(
        matrix,
        &[(position + 2, position)],
        &[(position + 1, position - 1)],
        tolerance,
    )
}

/// Verifies the current structural expectations for a trisection extraction.
pub fn verify_trisection_extraction(
    matrix: &CouplingMatrix,
    zero_positions: (usize, usize),
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    let (start, end) = zero_positions;
    if start < 2 || end > matrix.order() || end != start + 2 {
        return Err(crate::error::MfsError::PreconditionViolation(format!(
            "trisection verification requires ordered positions spanning one center resonator, got ({start}, {end})"
        )));
    }

    let center = (start + end) / 2;
    verify_required_entries(
        matrix,
        &[(start, end + 1), (center, end + 1)],
        &[],
        tolerance,
    )
}

fn magnitude(re: f64, im: f64) -> f64 {
    (re * re + im * im).sqrt()
}

fn is_structural_zero(value: f64, tolerance: MatrixPatternTolerance) -> bool {
    value.abs() <= tolerance.structural_zero
}

fn verify_required_entries(
    matrix: &CouplingMatrix,
    required_zeros: &[(usize, usize)],
    required_nonzeros: &[(usize, usize)],
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    let side = matrix.side();
    for &(row, col) in required_zeros.iter().chain(required_nonzeros.iter()) {
        if row >= side || col >= side {
            return Err(crate::error::MfsError::DimensionMismatch {
                expected: side,
                actual: row.max(col) + 1,
            });
        }
    }

    let mut report = SectionVerificationReport {
        zero_checks: required_zeros.len(),
        nonzero_checks: required_nonzeros.len(),
        failed_zero_checks: 0,
        failed_nonzero_checks: 0,
    };

    for &(row, col) in required_zeros {
        if !is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
            report.failed_zero_checks += 1;
        }
    }

    for &(row, col) in required_nonzeros {
        if is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
            report.failed_nonzero_checks += 1;
        }
    }

    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::CouplingMatrixBuilder;
    use crate::response::{ResponseSample, SParameterResponse};

    fn sample(
        frequency_hz: f64,
        normalized_omega: f64,
        s11_re: f64,
        s21_re: f64,
    ) -> ResponseSample {
        ResponseSample {
            frequency_hz,
            normalized_omega,
            group_delay: Some(0.25),
            s11_re,
            s11_im: 0.0,
            s21_re,
            s21_im: 0.0,
        }
    }

    #[test]
    fn compare_responses_reports_zero_deviation_for_identical_data() -> Result<()> {
        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9), sample(2.0, 0.0, 0.1, 0.95)],
        };
        let rhs = lhs.clone();

        let comparison = compare_responses(&lhs, &rhs)?;
        assert_eq!(comparison.samples_compared, 2);
        assert_eq!(comparison.max_s11_magnitude_deviation, 0.0);
        assert!(comparison.passes(ResponseTolerance::default()));
        Ok(())
    }

    #[test]
    fn response_check_report_tracks_skipped_and_checked_states() -> Result<()> {
        // A skipped check is `true` for the legacy boolean `passes()` but must
        // report as `Skipped`, not `Pass` (review P3 #2).
        let skipped = ResponseCheckReport::skipped();
        assert!(skipped.passes());
        assert_eq!(skipped.outcome(), CheckOutcome::Skipped);

        // A measured comparison that meets the tolerance is a genuine Pass.
        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9)],
        };
        let rhs = lhs.clone();
        let comparison = compare_responses(&lhs, &rhs)?;
        let checked =
            ResponseCheckReport::from_comparison(comparison, ResponseTolerance::default());
        assert_eq!(checked.invariant, Some(true));
        assert!(checked.passes());
        assert_eq!(checked.outcome(), CheckOutcome::Pass);

        // The same comparison against an impossible tolerance is a Fail.
        let failing = ResponseCheckReport::from_comparison(
            compare_responses(&lhs, &rhs)?,
            ResponseTolerance {
                frequency_hz: -1.0,
                ..ResponseTolerance::default()
            },
        );
        assert_eq!(failing.invariant, Some(false));
        assert!(!failing.passes());
        assert_eq!(failing.outcome(), CheckOutcome::Fail);
        Ok(())
    }

    #[test]
    fn compare_responses_detects_s_parameter_magnitude_change() -> Result<()> {
        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9)],
        };
        let rhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.23, 0.87)],
        };

        let comparison = compare_responses(&lhs, &rhs)?;
        assert!(comparison.max_s11_magnitude_deviation > 0.0);
        assert!(!comparison.passes(ResponseTolerance::default()));
        Ok(())
    }

    /// Regression for review P2-5: the pole-expansion fast path leaves
    /// `group_delay = None`, so without the `group_delay_compared` flag a
    /// transform that preserves magnitudes but shifts group delay would pass a
    /// strict tolerance on an unmeasured (0.0) deviation. A finite group-delay
    /// bound must fail when the group delay was never actually compared.
    #[test]
    fn compare_responses_uncompared_group_delay_fails_finite_tolerance() -> Result<()> {
        fn sample_none_group_delay() -> ResponseSample {
            ResponseSample {
                frequency_hz: 1.0,
                normalized_omega: -1.0,
                group_delay: None,
                s11_re: 0.2,
                s11_im: 0.0,
                s21_re: 0.9,
                s21_im: 0.0,
            }
        }
        let lhs = SParameterResponse {
            samples: vec![sample_none_group_delay()],
        };
        let rhs = lhs.clone();

        let comparison = compare_responses(&lhs, &rhs)?;
        assert!(!comparison.group_delay_compared);
        // A finite group-delay bound must NOT be satisfied by an unmeasured
        // group delay.
        assert!(
            !comparison.passes(ResponseTolerance::default()),
            "finite group-delay tolerance must fail when group delay was not compared"
        );

        // When the tolerance does not constrain group delay (infinite bound),
        // the unmeasured comparison is tolerated and the (identical)
        // magnitudes pass.
        let lax = ResponseTolerance {
            group_delay: f64::INFINITY,
            ..ResponseTolerance::default()
        };
        assert!(comparison.passes(lax));
        Ok(())
    }

    #[test]
    fn folded_pattern_check_accepts_folded_source_sparsity() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(matches_folded_pattern(
            &matrix,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    fn arrow_pattern_check_rejects_unreduced_source_coupling() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.4)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(!matches_arrow_pattern(
            &matrix,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    #[allow(deprecated)]
    fn topology_pattern_dispatch_accepts_wheel_using_current_arrow_rules() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(matches_topology_pattern(
            &matrix,
            MatrixTopology::Wheel,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    fn triplet_verification_report_passes_for_current_expected_pattern() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?
            .extract_triplet(-1.3, 2)?;

        let report = verify_triplet_extraction(&matrix, 2, MatrixPatternTolerance::default())?;
        assert!(report.passes());
        Ok(())
    }

    #[test]
    fn trisection_verification_report_passes_for_current_expected_pattern() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .topology(MatrixTopology::Arrow)
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.86)?
            .set_symmetric(1, 5, 0.25)?
            .set_symmetric(2, 3, 0.78)?
            .set_symmetric(3, 4, 0.69)?
            .set_symmetric(4, 5, 0.58)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.18)?
            .build()?
            .extract_trisection(-1.25, (2, 4))?;

        let report =
            verify_trisection_extraction(&matrix, (2, 4), MatrixPatternTolerance::default())?;
        assert!(report.passes());
        Ok(())
    }
    // ─── Absolute spec checks (review §D.3) ───────────────────────────────

    use crate::design::FilterDesign;
    use crate::spec::TransmissionZero;

    fn good_sixth_order_design() -> Result<FilterDesign> {
        FilterDesign::prototype(6, 23.0)
            .zeros([-2.0, 2.0])
            .synthesize()
    }

    /// (a) A good synthesized 6th-order design passes every measurable gate.
    #[test]
    fn verify_against_spec_good_design_passes_all_gates() -> Result<()> {
        let design = good_sixth_order_design()?;
        let response = design.response_normalized(-3.0, 3.0, 601)?;
        let report = verify_against_spec(design.spec(), &response);

        assert_eq!(report.outcome(), CheckOutcome::Pass);
        assert!(report.passes());

        let rl = report.gate("passband_return_loss").expect("RL gate present");
        assert_eq!(rl.outcome, CheckOutcome::Pass, "RL detail: {}", rl.detail);
        let measured_rl = rl.measured.expect("pass carries the measured RL");
        assert!(
            (22.5..=28.0).contains(&measured_rl),
            "measured passband RL {measured_rl} dB out of expected band"
        );

        let il = report.gate("passband_insertion_loss").expect("IL gate present");
        assert_eq!(il.outcome, CheckOutcome::Pass, "IL detail: {}", il.detail);
        let measured_il = il.measured.expect("pass carries the measured IL");
        assert!(
            (PASSBAND_MIN_IL_DB..=PASSBAND_MAX_IL_DB).contains(&measured_il),
            "measured passband IL {measured_il} dB out of allowed band"
        );

        for name in ["transmission_zero_depth_0", "transmission_zero_depth_1"] {
            let tz = report.gate(name).expect("TZ gate present");
            assert_eq!(tz.outcome, CheckOutcome::Pass, "{name}: {}", tz.detail);
            assert!(
                tz.measured.unwrap() <= TRANSMISSION_ZERO_DEPTH_DB,
                "{name}: measured depth above floor"
            );
        }

        // The stopband gate is documented as always Skipped: FilterSpec
        // carries no stopband requirement in normalized coordinates.
        let sb = report.gate("stopband_rejection").expect("stopband gate present");
        assert_eq!(sb.outcome, CheckOutcome::Skipped);
        Ok(())
    }

    /// (b) A spec 10 dB tighter than the design can achieve fails the RL gate
    /// with the measured (achievable) value, not the unachievable spec.
    #[test]
    fn verify_against_spec_tightened_rl_spec_fails_return_loss_gate() -> Result<()> {
        let design = good_sixth_order_design()?;
        let response = design.response_normalized(-3.0, 3.0, 601)?;
        let tightened = FilterSpec::new(6, 33.0)?.with_normalized_transmission_zeros([-2.0, 2.0]);

        let report = verify_against_spec(&tightened, &response);
        assert_eq!(report.outcome(), CheckOutcome::Fail);
        assert!(!report.passes());

        let rl = report.gate("passband_return_loss").expect("RL gate present");
        assert_eq!(rl.outcome, CheckOutcome::Fail, "RL detail: {}", rl.detail);
        let measured = rl.measured.expect("fail carries the measured value");
        assert!(
            (22.5..=28.0).contains(&measured),
            "measured RL {measured} dB should be the achievable ~23, not ~33"
        );
        assert_eq!(rl.threshold, Some(33.0));

        // The other gates are unaffected by the tightened RL spec.
        assert_eq!(
            report.gate("passband_insertion_loss").unwrap().outcome,
            CheckOutcome::Pass
        );
        Ok(())
    }

    /// (b′) A prescribed zero the response does not have fails the depth gate
    /// with the shallow (off-grid) measured depth.
    #[test]
    fn verify_against_spec_prescribed_zero_off_grid_fails_depth_gate() -> Result<()> {
        // A zero-free 4th-order design: its passband has no in-band dip, so a
        // prescribed in-band zero the response does not realize must fail.
        let design = FilterDesign::prototype(4, 20.0).synthesize()?;
        let response = design.response_normalized(-3.0, 3.0, 601)?;
        // The spec prescribes a zero at ω = 0.5 (inside the passband) that the
        // synthesized response does not have.
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros([0.5]);

        let report = verify_against_spec(&spec, &response);
        assert_eq!(report.outcome(), CheckOutcome::Fail);

        let missing = report.gate("transmission_zero_depth_0").unwrap();
        assert_eq!(missing.outcome, CheckOutcome::Fail, "detail: {}", missing.detail);
        let depth = missing.measured.expect("fail carries the measured depth");
        assert!(
            depth > TRANSMISSION_ZERO_DEPTH_DB,
            "depth {depth} dB should not reach the -40 dB floor"
        );
        Ok(())
    }

    /// (c) Missing data → Skipped, never a silent Pass (review P3 #2 / §D.3).
    #[test]
    fn verify_against_spec_missing_data_is_skipped_not_pass() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros([-2.0]);

        // Empty response: every gate Skipped, aggregate Skipped.
        let empty = SParameterResponse { samples: vec![] };
        let report = verify_against_spec(&spec, &empty);
        assert_eq!(report.outcome(), CheckOutcome::Skipped);
        assert_ne!(report.outcome(), CheckOutcome::Pass);
        for gate in &report.gates {
            assert_eq!(gate.outcome, CheckOutcome::Skipped, "{}", gate.gate);
            assert!(gate.measured.is_none(), "{}", gate.gate);
        }
        // The reportable rows still render the three-state outcome.
        assert!(report.format_rows().contains("skipped"));

        // Samples exist but none in the normalized passband: the passband
        // gates are Skipped, and the zero at ω = -2 is outside the window of
        // every sample (grid spacing 0.5 → radius 1.0, nearest ω = 1.5), so
        // the zero gate is Skipped too. Nothing passes silently.
        let out_of_band = SParameterResponse {
            samples: vec![sample(10.0, 1.5, 0.01, 0.9), sample(20.0, 2.0, 0.01, 0.2)],
        };
        let report = verify_against_spec(&spec, &out_of_band);
        assert_eq!(
            report.gate("passband_return_loss").unwrap().outcome,
            CheckOutcome::Skipped
        );
        assert_eq!(
            report.gate("passband_insertion_loss").unwrap().outcome,
            CheckOutcome::Skipped
        );
        assert_eq!(
            report.gate("transmission_zero_depth_0").unwrap().outcome,
            CheckOutcome::Skipped
        );
        Ok(())
    }

    /// A zero at infinity (spec convention: `f64::INFINITY`) imposes no dip
    /// requirement → its sub-gate is Skipped, while a real zero still Passes.
    #[test]
    fn verify_against_spec_zero_at_infinity_is_skipped() -> Result<()> {
        let design = FilterDesign::prototype(4, 20.0)
            .zeros([-2.0])
            .synthesize()?;
        let response = design.response_normalized(-3.0, 3.0, 601)?;
        let spec = FilterSpec::new(4, 20.0)?.with_transmission_zeros([
            TransmissionZero::normalized(-2.0),
            TransmissionZero::infinite(),
        ]);

        let report = verify_against_spec(&spec, &response);
        assert_eq!(
            report.gate("transmission_zero_depth_0").unwrap().outcome,
            CheckOutcome::Pass
        );
        let infinite = report.gate("transmission_zero_depth_1").unwrap();
        assert_eq!(infinite.outcome, CheckOutcome::Skipped);
        assert!(infinite.detail.contains("infinity"));
        Ok(())
    }

}
