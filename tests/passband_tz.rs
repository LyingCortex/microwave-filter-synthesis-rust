//! Integration coverage for the E1 investigation of
//! `docs/continuation-plan.md` section E item 1 ("passband transmission
//! zeros, `|zero| < 1`").
//!
//! The investigation reached a **negative** result, and these tests pin it:
//!
//! 1. `generalized_chebyshev_polynomials` rejects a finite zero with
//!    `|zero| < 1` with a **typed** `MfsError::InvalidTransmissionZero` (not a
//!    helper-level `NotImplemented`). A literal relaxation of the Cameron
//!    `sqrt(1 - 1/zero^2)` term to an imaginary value produces a complex
//!    non-Hermitian `F(s)` and breaks the Feldtkeller identity
//!    `|E|^2 = |F|^2 + |P|^2/eps^2` at O(1); a valid in-band null would need
//!    the *multifunction* characteristic, not this helper.
//!
//! 2. Multi-band route-α's inter-band gap zeros are **ordinary stopband
//!    zeros**, not passband zeros. A `BandPassMapping`'s `|Omega| = 1`
//!    boundary *is* that band's own passband edge, so every inter-band gap
//!    frequency maps to `|Omega| >= 1` on **both** adjacent bands. This
//!    corrects the C0 design doc's original claim and means route-α does not
//!    depend on E1.

use mfs::approx::generalized_chebyshev_polynomials;
use mfs::error::{MfsError, Result};
use mfs::freq::{BandPassMapping, FrequencyMapping};
use mfs::spec::{FilterSpec, TransmissionZero};

/// A zero strictly inside the passband (`|zero| < 1`) must be rejected at the
/// approximation entry point with the typed error, for any sign.
#[test]
fn passband_transmission_zero_is_rejected_with_typed_error() {
    for zero in [0.75_f64, 0.5, 0.1, -0.9, -0.25] {
        let spec = FilterSpec::new(4, 20.0)
            .expect("valid order and return loss")
            .with_transmission_zeros(vec![TransmissionZero::normalized(zero)]);
        match generalized_chebyshev_polynomials(&spec) {
            Err(MfsError::InvalidTransmissionZero(message)) => {
                assert!(
                    message.contains("passband"),
                    "zero {zero}: unexpected message {message:?}"
                );
            }
            other => panic!("zero {zero}: expected InvalidTransmissionZero, got {other:?}"),
        }
    }
}

/// A zero **on** the passband edge (`|zero| = 1`) is also rejected: it is the
/// degenerate case where the ripple zero meets the transmission zero.
#[test]
fn unit_magnitude_transmission_zero_is_rejected() {
    let spec = FilterSpec::new(4, 20.0)
        .expect("valid order and return loss")
        .with_transmission_zeros(vec![TransmissionZero::normalized(1.0)]);
    match generalized_chebyshev_polynomials(&spec) {
        Err(MfsError::InvalidTransmissionZero(_)) => {}
        other => panic!("expected InvalidTransmissionZero for |zero| = 1, got {other:?}"),
    }
}

/// `|zero| > 1` (stopband) remains fully supported — the E1 rejection must not
/// regress the shipped path.
#[test]
fn stopband_transmission_zero_still_synthesizes() -> Result<()> {
    let spec = FilterSpec::new(4, 20.0)?
        .with_normalized_transmission_zeros(vec![1.5, -2.5]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;
    assert!(polynomials.generalized.is_some());
    Ok(())
}

/// The corrected TZ law: for two disjoint passbands, every frequency in the
/// inter-band gap maps to `|Omega| >= 1` on **both** adjacent band prototypes.
/// Worked example from the investigation: bands 10-11 GHz and 12-13 GHz, gap
/// 11-12 GHz.
#[test]
fn inter_band_gap_maps_to_stopband_on_both_bands() -> Result<()> {
    use mfs::freq::edges_to_center_bandwidth;
    let (c_lo, bw_lo) = edges_to_center_bandwidth(10.0e9, 11.0e9);
    let (c_hi, bw_hi) = edges_to_center_bandwidth(12.0e9, 13.0e9);
    let lower = BandPassMapping::new(c_lo, bw_lo)?;
    let upper = BandPassMapping::new(c_hi, bw_hi)?;

    // Passband edges map exactly to |Omega| = 1 for their own band (that is the
    // definition of the transform's edges).
    assert!((lower.map_hz_to_normalized(10.0e9)?.omega + 1.0).abs() < 1e-9);
    assert!((lower.map_hz_to_normalized(11.0e9)?.omega - 1.0).abs() < 1e-9);
    assert!((upper.map_hz_to_normalized(12.0e9)?.omega + 1.0).abs() < 1e-9);
    assert!((upper.map_hz_to_normalized(13.0e9)?.omega - 1.0).abs() < 1e-9);

    // Every gap sample is >= 1 in magnitude on both bands.
    for f in [11.25e9, 11.5e9, 11.75e9] {
        let ol = lower.map_hz_to_normalized(f)?.omega;
        let ou = upper.map_hz_to_normalized(f)?.omega;
        assert!(
            ol.abs() >= 1.0 - 1e-9,
            "lower band: gap {f} mapped to Omega={ol} (expected |Omega| >= 1)"
        );
        assert!(
            ou.abs() >= 1.0 - 1e-9,
            "upper band: gap {f} mapped to Omega={ou} (expected |Omega| >= 1)"
        );
        // The two sides have opposite sign (below the upper band's edge).
        assert!(ol > 0.0 && ou < 0.0, "gap {f}: Omega_l={ol}, Omega_u={ou}");
    }
    Ok(())
}
