// Feature: synthesis-numerical-fix
// Regression tests for the 11 previously-failing filter configurations.
// **Validates: Requirements 5.1, 5.2, 5.3**
//
// These configurations previously failed with:
//   NumericalFailure("y12 residue is unexpectedly complex in the current real-valued synthesis path")
//
// The fix extends the synthesis path to handle complex-conjugate residue pairs.
//
// Note on order-2 all-pole: These configurations fail with DimensionMismatch because
// the admittance polynomial denominator degree drops below N due to leading coefficient
// cancellation in even-order all-pole filters. This is a separate known limitation
// (not the NumericalFailure that this spec fixes). The tests document this behavior.
//
// Note on triplet/trisection: These section extractions require order >= 5 (triplet
// needs center_resonator in [2, order-2) and trisection needs positions spanning one
// center resonator). For order-3 configs with TZs, we verify folded topology
// transformation succeeds instead, which is the prerequisite for section extraction.

use mfs::approx::generalized_chebyshev_polynomials;
use mfs::error::{MfsError, Result};
use mfs::freq::BandPassMapping;
use mfs::spec::FilterSpec;
use mfs::synthesis::MatrixSynthesisEngine;
use mfs::transform::{transform_matrix, TopologyKind};

// ---------------------------------------------------------------------------
// Helper: build a FilterSpec from order, return loss, and optional TZ list
// ---------------------------------------------------------------------------

fn build_spec(order: usize, return_loss_db: f64, tzs: &[f64]) -> FilterSpec {
    FilterSpec::new(order, return_loss_db)
        .expect("valid order and return loss")
        .with_normalized_transmission_zeros(tzs.to_vec())
}

// ---------------------------------------------------------------------------
// Helper: run the full synthesis pipeline and return the coupling matrix
// ---------------------------------------------------------------------------

fn synthesize_matrix(spec: &FilterSpec) -> mfs::matrix::CouplingMatrix {
    let polynomials =
        generalized_chebyshev_polynomials(spec).expect("polynomial generation should succeed");
    MatrixSynthesisEngine
        .synthesize(&polynomials)
        .expect("synthesis should succeed for previously-failing configuration")
}

// ===========================================================================
// Requirement 5.1 / Core: Synthesis succeeds for previously-failing configs
// ===========================================================================

// --- Order 2, all-pole (known limitation: DimensionMismatch) ---

/// Order 2 all-pole filters previously failed with DimensionMismatch due to
/// leading coefficient cancellation. This is now fixed by normalizing E(s)
/// when the leading coefficients would cancel.
#[test]
fn regression_order2_allpole_15db_known_dimension_mismatch() {
    let spec = build_spec(2, 15.0, &[]);
    let polynomials = generalized_chebyshev_polynomials(&spec).unwrap();
    let result = MatrixSynthesisEngine.synthesize(&polynomials);
    assert!(
        result.is_ok(),
        "order 2 all-pole should now synthesize successfully, got: {:?}",
        result.unwrap_err()
    );
}

#[test]
fn regression_order2_allpole_20db_known_dimension_mismatch() {
    let spec = build_spec(2, 20.0, &[]);
    let polynomials = generalized_chebyshev_polynomials(&spec).unwrap();
    let result = MatrixSynthesisEngine.synthesize(&polynomials);
    assert!(
        result.is_ok(),
        "order 2 all-pole should now synthesize successfully, got: {:?}",
        result.unwrap_err()
    );
}

#[test]
fn regression_order2_allpole_25db_known_dimension_mismatch() {
    let spec = build_spec(2, 25.0, &[]);
    let polynomials = generalized_chebyshev_polynomials(&spec).unwrap();
    let result = MatrixSynthesisEngine.synthesize(&polynomials);
    assert!(
        result.is_ok(),
        "order 2 all-pole should now synthesize successfully, got: {:?}",
        result.unwrap_err()
    );
}

// --- Order 3, all-pole (previously failed with NumericalFailure, now fixed) ---

#[test]
fn regression_order3_allpole_15db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 15.0, &[]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

#[test]
fn regression_order3_allpole_20db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 20.0, &[]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

#[test]
fn regression_order3_allpole_25db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 25.0, &[]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

// --- Order 3, 2 TZs at ±1.5 (previously failed with NumericalFailure, now fixed) ---

#[test]
fn regression_order3_2tz_pm1_5_15db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 15.0, &[-1.5, 1.5]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

#[test]
fn regression_order3_2tz_pm1_5_20db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 20.0, &[-1.5, 1.5]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

#[test]
fn regression_order3_2tz_pm1_5_25db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 25.0, &[-1.5, 1.5]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

// --- Order 3, 2 TZs at ±2.0 (previously failed with NumericalFailure, now fixed) ---

#[test]
fn regression_order3_2tz_pm2_0_20db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 20.0, &[-2.0, 2.0]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

// --- Order 3, 1 TZ at 1.5 (previously failed with NumericalFailure, now fixed) ---

#[test]
fn regression_order3_1tz_1_5_20db_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(3, 20.0, &[1.5]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 3);
    assert_eq!(matrix.side(), 5);
    Ok(())
}

// ===========================================================================
// Requirement 5.1: Bandpass scaling succeeds on previously-failing configs
// ===========================================================================

fn assert_bandpass_scaling_succeeds(spec: &FilterSpec) {
    let matrix = synthesize_matrix(spec);
    let mapping = BandPassMapping::new(6.75e9, 300.0e6).expect("valid bandpass mapping");
    let scaled = matrix
        .denormalize_bandpass(&mapping)
        .expect("denormalize_bandpass should succeed on synthesized matrix");
    assert_eq!(scaled.order(), matrix.order());
    assert_eq!(scaled.side(), matrix.side());
}

#[test]
fn regression_bandpass_scaling_order3_allpole_15db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 15.0, &[]));
}

#[test]
fn regression_bandpass_scaling_order3_allpole_20db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 20.0, &[]));
}

#[test]
fn regression_bandpass_scaling_order3_allpole_25db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 25.0, &[]));
}

#[test]
fn regression_bandpass_scaling_order3_2tz_pm1_5_15db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 15.0, &[-1.5, 1.5]));
}

#[test]
fn regression_bandpass_scaling_order3_2tz_pm1_5_20db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 20.0, &[-1.5, 1.5]));
}

#[test]
fn regression_bandpass_scaling_order3_2tz_pm1_5_25db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 25.0, &[-1.5, 1.5]));
}

#[test]
fn regression_bandpass_scaling_order3_2tz_pm2_0_20db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 20.0, &[-2.0, 2.0]));
}

#[test]
fn regression_bandpass_scaling_order3_1tz_1_5_20db() {
    assert_bandpass_scaling_succeeds(&build_spec(3, 20.0, &[1.5]));
}

// ===========================================================================
// Requirement 5.2: External Q extraction produces positive values
// ===========================================================================

fn assert_external_q_positive(spec: &FilterSpec) {
    let matrix = synthesize_matrix(spec);
    let mapping = BandPassMapping::new(6.75e9, 300.0e6).expect("valid bandpass mapping");
    let scaled = matrix
        .denormalize_bandpass_with_external_q(&mapping)
        .expect("denormalize_bandpass_with_external_q should succeed");
    assert!(
        scaled.source_external_q() > 0.0,
        "source external Q should be positive, got {}",
        scaled.source_external_q()
    );
    assert!(
        scaled.load_external_q() > 0.0,
        "load external Q should be positive, got {}",
        scaled.load_external_q()
    );
}

#[test]
fn regression_external_q_order3_allpole_15db() {
    assert_external_q_positive(&build_spec(3, 15.0, &[]));
}

#[test]
fn regression_external_q_order3_allpole_20db() {
    assert_external_q_positive(&build_spec(3, 20.0, &[]));
}

#[test]
fn regression_external_q_order3_allpole_25db() {
    assert_external_q_positive(&build_spec(3, 25.0, &[]));
}

#[test]
fn regression_external_q_order3_2tz_pm1_5_15db() {
    assert_external_q_positive(&build_spec(3, 15.0, &[-1.5, 1.5]));
}

#[test]
fn regression_external_q_order3_2tz_pm1_5_20db() {
    assert_external_q_positive(&build_spec(3, 20.0, &[-1.5, 1.5]));
}

#[test]
fn regression_external_q_order3_2tz_pm1_5_25db() {
    assert_external_q_positive(&build_spec(3, 25.0, &[-1.5, 1.5]));
}

#[test]
fn regression_external_q_order3_2tz_pm2_0_20db() {
    assert_external_q_positive(&build_spec(3, 20.0, &[-2.0, 2.0]));
}

#[test]
fn regression_external_q_order3_1tz_1_5_20db() {
    assert_external_q_positive(&build_spec(3, 20.0, &[1.5]));
}

// ===========================================================================
// Requirement 5.3: Topology transformation succeeds for configs with TZs
// ===========================================================================
//
// Triplet/trisection section extraction requires order >= 5 (triplet needs
// center_resonator in [2, order-2), trisection needs positions (start, end)
// where start >= 2, end <= order, end == start + 2).
//
// For order-3 configs with TZs, we verify:
// 1. Folded topology transformation succeeds (prerequisite for section extraction)
// 2. For a higher-order filter (order 5) with the same TZ pattern, triplet and
//    trisection extraction succeed — proving the synthesis output is compatible
//    with downstream section operations.

fn assert_folded_transform_succeeds(spec: &FilterSpec) {
    let matrix = synthesize_matrix(spec);
    let outcome = transform_matrix(&matrix, TopologyKind::Folded)
        .expect("folded topology transformation should succeed");
    assert_eq!(outcome.topology, TopologyKind::Folded);
    assert!(
        outcome.report.passes(),
        "folded transform report should pass verification"
    );
}

#[test]
fn regression_folded_transform_order3_2tz_pm1_5_15db() {
    assert_folded_transform_succeeds(&build_spec(3, 15.0, &[-1.5, 1.5]));
}

#[test]
fn regression_folded_transform_order3_2tz_pm1_5_20db() {
    assert_folded_transform_succeeds(&build_spec(3, 20.0, &[-1.5, 1.5]));
}

#[test]
fn regression_folded_transform_order3_2tz_pm1_5_25db() {
    assert_folded_transform_succeeds(&build_spec(3, 25.0, &[-1.5, 1.5]));
}

#[test]
fn regression_folded_transform_order3_2tz_pm2_0_20db() {
    assert_folded_transform_succeeds(&build_spec(3, 20.0, &[-2.0, 2.0]));
}

#[test]
fn regression_folded_transform_order3_1tz_1_5_20db() {
    assert_folded_transform_succeeds(&build_spec(3, 20.0, &[1.5]));
}

// ---------------------------------------------------------------------------
// Higher-order triplet/trisection extraction using the same TZ patterns
// ---------------------------------------------------------------------------
//
// These tests verify that the synthesis fix produces matrices compatible with
// triplet and trisection extraction by using order-5 filters with the same
// TZ magnitudes as the previously-failing order-3 configs.

use mfs::synthesis::SectionSynthesis;

#[test]
fn regression_triplet_extraction_order5_tz_1_5() -> Result<()> {
    let spec = build_spec(5, 20.0, &[-1.5]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let sections = SectionSynthesis {};
    // Triplet extraction succeeds (returns Ok) — the synthesized matrix is
    // compatible with the section extraction pipeline.
    let outcome = sections.synthesize_triplet(&polynomials, -1.5, 2)?;
    assert_eq!(outcome.matrix.order(), 5);
    Ok(())
}

#[test]
fn regression_triplet_extraction_order5_tz_2_0() -> Result<()> {
    let spec = build_spec(5, 20.0, &[-2.0]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let sections = SectionSynthesis {};
    // Triplet extraction succeeds (returns Ok) — the synthesized matrix is
    // compatible with the section extraction pipeline.
    let outcome = sections.synthesize_triplet(&polynomials, -2.0, 2)?;
    assert_eq!(outcome.matrix.order(), 5);
    Ok(())
}

#[test]
fn regression_trisection_extraction_order5_tz_1_5() -> Result<()> {
    let spec = build_spec(5, 20.0, &[-1.5]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let sections = SectionSynthesis {};
    let outcome = sections.synthesize_trisection(&polynomials, -1.5, (2, 4))?;
    assert!(
        outcome.verification.passes(),
        "trisection extraction should pass verification"
    );
    Ok(())
}

#[test]
fn regression_trisection_extraction_order5_tz_2_0() -> Result<()> {
    let spec = build_spec(5, 20.0, &[-2.0]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    let sections = SectionSynthesis {};
    let outcome = sections.synthesize_trisection(&polynomials, -2.0, (2, 4))?;
    assert!(
        outcome.verification.passes(),
        "trisection extraction should pass verification"
    );
    Ok(())
}

// ===========================================================================
// All-finite zero-set feasibility domain (review P1-1)
// ===========================================================================
//
// The all-finite formulation uses eps_r = eps / sqrt(eps^2 - 1), which is
// real-valued only for eps > 1. When every transmission zero is finite and
// lies near the band edge (|z| ~ 1) or the return loss is high, eps <= 1 and
// the spec is infeasible for the all-finite path. Before the P1-1 fix this
// propagated a NaN 1/eps_r into the E-polynomial stage and surfaced later as
// the misleading generic "coefficients must be finite" error. `find_eps` now
// rejects with a domain-specific PreconditionViolation naming the feasible
// alternative (a mixed zero set, which forces eps_r = 1).

/// A 4th-order all-finite spec with well-separated zeros but RL = 34 dB has
/// eps <= 1 and must be rejected with the domain-specific message (the review
/// measured RL <= 33 feasible, RL >= 34 infeasible for these zeros) — not the
/// generic "coefficients must be finite" finiteness error.
#[test]
fn all_finite_infeasible_high_rl_returns_typed_error() {
    let spec = build_spec(4, 34.0, &[-1.5, 1.5, -2.0, 2.0]);
    let err = generalized_chebyshev_polynomials(&spec).unwrap_err();
    assert!(
        matches!(err, MfsError::PreconditionViolation(_)),
        "expected PreconditionViolation for infeasible all-finite spec, got: {err:?}"
    );
    let msg = format!("{err}");
    assert!(
        msg.contains("infeasible"),
        "error must name the all-finite domain, got: {msg}"
    );
    assert!(
        !msg.contains("must be finite"),
        "must not fall through to the generic finiteness error: {msg}"
    );
}

/// Control: the same 4th-order all-finite zeros at a lower return loss (RL = 20
/// dB) are feasible (eps > 1) and the polynomial stage succeeds with a finite,
/// positive eps_r — proving the guard above is specific to the infeasible
/// domain and does not reject working all-finite designs.
#[test]
fn all_finite_feasible_low_rl_succeeds() -> Result<()> {
    let spec = build_spec(4, 20.0, &[-1.5, 1.5, -2.0, 2.0]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    // Feasible all-finite: eps_r = eps/sqrt(eps^2-1) is a finite positive real,
    // and the E polynomial is finite. This is the boundary the guard must NOT
    // trip.
    assert!(
        polynomials.eps_r.is_finite() && polynomials.eps_r > 0.0,
        "feasible all-finite spec must yield a finite positive eps_r, got {}",
        polynomials.eps_r
    );
    assert!(
        polynomials.e.coefficients.iter().all(|c| c.re.is_finite() && c.im.is_finite()),
        "E coefficients must all be finite for a feasible spec"
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Pole-coalescence regression (7p3z, RL = 23 dB)
// ---------------------------------------------------------------------------
// A transmission zero can drive two denominator poles of the Chebyshev
// prototype into numerical near-coalescence. For order 7 with
// z = [1.1362, 2.3415, -1.4132] the two poles at s ≈ j·1.1362 are distinct
// but 4.95e-6 apart — below the former 1e-5 cluster tolerance, which merged
// them into a single double pole and failed the transversal order check
// ("expected 7, got 6"). The tolerance is now anchored to the root
// solver's resolution floor (3e-7), so these distinct poles stay separate
// and the design synthesizes. See docs/internal/pole-count-defect-7p3z-
// 2026-10-26.md for the full root-cause record.

/// The recorded defect input: order 7, RL = 23 dB, z0 = 1.1362.
#[test]
fn regression_order7_7p3z_pole_coalescence_synthesis_succeeds() -> Result<()> {
    let spec = build_spec(7, 23.0, &[1.1362, 2.3415, -1.4132]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 7);
    assert!(
        matrix.as_slice().iter().all(|v| v.is_finite()),
        "coalescence-regression matrix must be fully finite"
    );
    // The near-coincident pair must appear as two distinct resonator
    // diagonals at ≈ -1.1362, not a single merged pole.
    let near = (0..matrix.order())
        .filter_map(|i| matrix.resonator_diagonal(i))
        .filter(|d| (d - -1.1362).abs() < 5e-3)
        .count();
    assert_eq!(
        near, 2,
        "expected two distinct resonators near -1.1362, got {near}"
    );
    Ok(())
}

/// Inside the former failure window: z0 = 1.13621 puts the pole pair at a
/// 5.85e-7 gap — still above the 3e-7 resolution floor, so it must synthesize
/// even though the old 1e-5 tolerance rejected it.
#[test]
fn regression_order7_7p3z_coalescence_window_edge_succeeds() -> Result<()> {
    let spec = build_spec(7, 23.0, &[1.13621, 2.3415, -1.4132]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 7);
    Ok(())
}

/// Control just outside the coalescence region: z0 = 1.13624 (gap 1.25e-5)
/// synthesized under the old tolerance and must continue to.
#[test]
fn regression_order7_7p3z_control_outside_window_succeeds() -> Result<()> {
    let spec = build_spec(7, 23.0, &[1.13624, 2.3415, -1.4132]);
    let matrix = synthesize_matrix(&spec);
    assert_eq!(matrix.order(), 7);
    Ok(())
}
