//! C3 verification surface for true multi-passband (multi-band) synthesis —
//! the RED→GREEN anchor for section C of `docs/continuation-plan.md`.
//!
//! The acceptance fixture is Cameron 2011 Fig. 25
//! (`Cameron2011_10_5_dual_band_symmetric`, `tests/cameron2011_fixtures.json`):
//! a dual-band symmetric design, order 10, 22 dB return loss, passbands
//! `[-1,-0.35]` and `[0.35,1]` in the shared normalized prototype axis, with
//! 5 transmission zeros (2 outer at ±1.1573, 3 inter-band at
//! `{0, ±0.2542}`).
//!
//! Route-α (C0, `docs/internal/multi-band-synthesis-design.md`) synthesizes a
//! **transversal** $N+2$ matrix via the existing `MatrixSynthesisEngine`,
//! while the fixture's printed matrix is a `trisection_cascade_with_extracted_
//! zero` topology. The two matrices are *different topologies of the same
//! network*, so "reproduces the fixture" is asserted at the **response**
//! level (topologically invariant), plus absolute physical checks that pin
//! the design to the spec independent of the fixture.
//!
//! Gates (all measured on a 4001-point grid over `-2..2`):
//! 1. **Response invariance** — max|ΔS11|, max|ΔS21| < 1e-3 vs the fixture.
//! 2. **Per-band equiripple** — worst in-band |S11| = −22 dB (the RL spec) ±1.
//! 3. **Unitarity** — max‖|S11|²+|S21|²−1‖ < 1e-10 (lossless, matched 2-port).
//! 4. **Transmission zeros** — |S21| < −40 dB in a ±0.005 window at each TZ.
//! 5. **Outer stopband** — worst |S21| < −15 dB over |w|∈[1.10, 2.0] (away
//!    from the band edge, where the Chebyshev lobe is weakest).

use mfs::approx::multi_band_chebyshev_polynomials;
use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::response::ResponseSolver;
use mfs::synthesis::MatrixSynthesisEngine;

/// Cameron 2011 Fig. 25 dual-band symmetric spec (normalized prototype domain).
const BANDS: &[(f64, f64)] = &[(-1.0, -0.35), (0.35, 1.0)];
const ORDER: usize = 10;
const RL: f64 = 22.0;
const TZ: &[f64] = &[0.0, 0.25420871, -0.25420871, 1.15728568, -1.15728568];

/// Loads the printed $N+2$ coupling matrix for the Fig. 25 fixture.
fn fixture_matrix() -> Result<CouplingMatrix> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/cameron2011_fixtures.json");
    let raw = std::fs::read_to_string(&path).expect("read fixtures JSON");
    let val: serde_json::Value = serde_json::from_str(&raw).expect("parse fixtures JSON");
    let case = val["case"]
        .as_array()
        .ok_or_else(|| {
            mfs::MfsError::PreconditionViolation("fixtures: missing `case` array".into())
        })?
        .iter()
        .find(|c| c["case_id"] == "Cameron2011_10_5_dual_band_symmetric")
        .ok_or_else(|| {
            mfs::MfsError::PreconditionViolation("fixtures: Fig.25 case not found".into())
        })?;
    let rows = case["coupling_matrix"]["data"]
        .as_array()
        .ok_or_else(|| {
            mfs::MfsError::PreconditionViolation("fixtures: missing matrix `data`".into())
        })?;
    let side = rows.len();
    let mut flat = Vec::with_capacity(side * side);
    for r in rows {
        let cols = r
            .as_array()
            .ok_or_else(|| mfs::MfsError::PreconditionViolation("fixtures: data row".into()))?;
        if cols.len() != side {
            return Err(mfs::MfsError::PreconditionViolation(
                "fixtures: matrix row/column length mismatch".into(),
            ));
        }
        for v in cols {
            flat.push(v.as_f64().ok_or_else(|| {
                mfs::MfsError::PreconditionViolation("fixtures: non-numeric entry".into())
            })?);
        }
    }
    // The printed matrix is (N+2)×(N+2); the API's `order` is the resonator
    // count (side − 2 = 10).
    CouplingMatrix::new(side - 2, flat)
}

/// The C3 anchor: the synthesized multi-band response reproduces the
/// fixture's response (topologically invariant) and meets every absolute
/// spec-level quantity on the Fig. 25 design.
#[test]
fn dual_band_response_reproduces_cameron_fig25_fixture() -> Result<()> {
    let set = multi_band_chebyshev_polynomials(BANDS, ORDER, RL, TZ)?;
    let synthesized = MatrixSynthesisEngine {}.synthesize(&set)?;
    let fixture = fixture_matrix()?;

    let grid = FrequencyGrid::linspace(-2.0, 2.0, 4001)?;
    let omegas = grid.as_slice().to_vec();
    let solver = ResponseSolver;
    let rs = solver.evaluate_normalized(&synthesized, &grid)?;
    let rf = solver.evaluate_normalized(&fixture, &grid)?;

    // Gate 1 — response invariance vs the fixture (topologically invariant).
    let mut d11 = 0.0_f64;
    let mut d21 = 0.0_f64;
    for (a, b) in rs.samples.iter().zip(rf.samples.iter()) {
        d11 = d11.max((a.s11_mag() - b.s11_mag()).abs());
        d21 = d21.max((a.s21_mag() - b.s21_mag()).abs());
    }
    assert!(
        d11 < 1e-3,
        "response invariance: max|ΔS11| = {d11:.3e} ≥ 1e-3 (synthesized vs fixture)"
    );
    assert!(
        d21 < 1e-3,
        "response invariance: max|ΔS21| = {d21:.3e} ≥ 1e-3 (synthesized vs fixture)"
    );

    // Gate 2 — per-band equiripple: worst in-band |S11| is the RL floor.
    for (lo, hi) in [(0.35_f64, 1.0_f64), (-1.0_f64, -0.35_f64)] {
        let worst = omegas
            .iter()
            .enumerate()
            .filter(|(_, w)| **w >= lo && **w <= hi)
            .map(|(i, _)| rs.samples[i].s11_db())
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (worst + RL).abs() <= 1.0,
            "band [{lo},{hi}]: worst in-band |S11| = {worst:.2} dB, expected -{RL} ± 1 dB"
        );
    }

    // Gate 3 — unitarity: lossless matched 2-port (|S11|² + |S21|² = 1).
    let mut dunit = 0.0_f64;
    for s in &rs.samples {
        let a = s.s11_mag();
        let b = s.s21_mag();
        dunit = dunit.max((a * a + b * b - 1.0).abs());
    }
    assert!(
        dunit < 1e-10,
        "unitarity: max‖|S11|²+|S21|²−1‖ = {dunit:.3e} ≥ 1e-10"
    );

    // Gate 4 — transmission zeros: a deep |S21| null at each spec TZ.
    for tz in [0.0_f64, 0.25420871, 1.15728568] {
        let min_db = omegas
            .iter()
            .enumerate()
            .filter(|(_, w)| (**w - tz).abs() <= 5e-3)
            .map(|(i, _)| rs.samples[i].s21_db())
            .fold(f64::INFINITY, f64::min);
        assert!(
            min_db < -40.0,
            "TZ at w={tz}: min |S21| = {min_db:.1} dB, expected < -40 dB"
        );
    }

    // Gate 5 — outer stopband: rejection away from the band edge. The
    // Chebyshev stopband lobe is weakest right at the band edge (|w|≈1.0,
    // ≈ −3.7 dB) and deepens outward (worst −16.6 dB over |w|∈[1.10, 2.0],
    // with the outer TZ null at ±1.1573). The window excludes the edge so the
    // gate is not dominated by the passband-edge transition.
    let outer = omegas
        .iter()
        .enumerate()
        .filter(|(_, w)| w.abs() >= 1.10 && w.abs() <= 2.0)
        .map(|(i, _)| rs.samples[i].s21_db())
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        outer < -15.0,
        "outer stopband: worst |S21| over |w|∈[1.10, 2.0] = {outer:.2} dB, expected < -15 dB"
    );

    Ok(())
}
