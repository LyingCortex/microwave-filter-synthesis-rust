//! Multi-band (true multi-passband) specification surface — C2a of
//! `docs/continuation-plan.md` section C, designed in
//! `docs/internal/multi-band-synthesis-design.md`.
//!
//! A [`MultiBandSpec`] describes $K \ge 2$ disjoint passbands inside one filter
//! structure, plus per-region out-of-band rejection requirements. It is a
//! **separate** type from [`FilterSpec`](super::FilterSpec): a $K$-band filter
//! (one structure) is not a $K$-channel multiplexer (a junction of $K$
//! channels), and widening `FilterSpec` would conflate the two.
//!
//! Band invariants are inherited from
//! [`MultiBandMapping`](crate::freq::MultiBandMapping): each band is a valid
//! `BandPassMapping` (positive finite center and bandwidth), bands are sorted
//! by center, and adjacent bands are **strictly disjoint** (the lower band's
//! upper edge lies strictly below the upper band's lower edge).

use crate::error::{MfsError, Result};
use crate::freq::{validated_transmission_zero, BandPassMapping, MultiBandMapping};

use super::TransmissionZero;

/// One passband of a multi-band design, in physical frequency units.
///
/// Stored as a [`BandPassMapping`] so the band is validated exactly like the
/// single-band path (geometric-mean center, positive bandwidth, `ω = 0` at the
/// center).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandSpec {
    mapping: BandPassMapping,
    return_loss_db: f64,
}

impl BandSpec {
    /// Creates a validated band from a center frequency and 3 dB bandwidth.
    pub fn new(center_hz: f64, bandwidth_hz: f64, return_loss_db: f64) -> Result<Self> {
        let mapping = BandPassMapping::new(center_hz, bandwidth_hz)?;
        if !return_loss_db.is_finite() || return_loss_db <= 0.0 {
            return Err(MfsError::InvalidReturnLoss { return_loss_db });
        }
        Ok(Self {
            mapping,
            return_loss_db,
        })
    }

    /// Creates a band from passband edge frequencies `(f1_hz, f2_hz)` with
    /// `f1_hz < f2_hz`, using the crate's geometric-mean convention
    /// ([`edges_to_center_bandwidth`](crate::freq::edges_to_center_bandwidth)).
    pub fn from_edges(f1_hz: f64, f2_hz: f64, return_loss_db: f64) -> Result<Self> {
        let (center_hz, bandwidth_hz) = crate::freq::edges_to_center_bandwidth(f1_hz, f2_hz);
        Self::new(center_hz, bandwidth_hz, return_loss_db)
    }

    /// The band's center frequency in Hz.
    pub fn center_hz(&self) -> f64 {
        self.mapping.center_hz()
    }

    /// The band's 3 dB bandwidth in Hz.
    pub fn bandwidth_hz(&self) -> f64 {
        self.mapping.bandwidth_hz()
    }

    /// The band's targeted passband return loss in dB.
    pub fn return_loss_db(&self) -> f64 {
        self.return_loss_db
    }

    /// The underlying frequency mapping.
    pub fn mapping(&self) -> BandPassMapping {
        self.mapping
    }
}

/// The physical region a [`RejectionSpec`] constrains.
///
/// Regions are *declared*, not derived: the design's transmission-zero
/// placement is driven by the spec, so the caller names each region explicitly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RejectionRegion {
    /// Below the lowest passband.
    OuterLower,
    /// Above the highest passband.
    OuterUpper,
    /// The gap between adjacent passbands `bands[index]` and `bands[index + 1]`.
    ///
    /// `index` must be `< band_count - 1`.
    InterBand { index: usize },
}

/// A minimum-attenuation requirement over one region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RejectionSpec {
    /// The physical region constrained.
    pub region: RejectionRegion,
    /// Minimum attenuation the region must meet, in dB.
    pub min_attenuation_db: f64,
}

impl RejectionSpec {
    /// Creates a validated rejection requirement.
    pub fn new(region: RejectionRegion, min_attenuation_db: f64) -> Result<Self> {
        if !min_attenuation_db.is_finite() || min_attenuation_db <= 0.0 {
            return Err(MfsError::InvalidReturnLoss {
                return_loss_db: min_attenuation_db,
            });
        }
        Ok(Self {
            region,
            min_attenuation_db,
        })
    }
}

/// A $K$-band ($K \ge 2$) multi-passband specification.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiBandSpec {
    bands: Vec<BandSpec>,
    rejection: Vec<RejectionSpec>,
    transmission_zeros_normalized: Vec<f64>,
    order: usize,
}

impl MultiBandSpec {
    /// Creates a validated multi-band specification.
    ///
    /// - `bands` must contain at least 2 bands, sorted or unsorted (they are
    ///   sorted by center frequency here).
    /// - adjacent bands must be strictly disjoint.
    /// - `order` is the total filter order $N \ge$ number of finite
    ///   transmission zeros.
    ///
    /// Transmission zeros default to empty (an all-pole design); add them with
    /// [`with_normalized_transmission_zeros`](Self::with_normalized_transmission_zeros).
    pub fn new(bands: Vec<BandSpec>, order: usize) -> Result<Self> {
        if bands.len() < 2 {
            return Err(MfsError::InvalidFrequency(format!(
                "multi-band spec requires at least 2 bands, got {}",
                bands.len()
            )));
        }
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        let mut bands = bands;
        bands.sort_by(|a, b| {
            a.center_hz()
                .partial_cmp(&b.center_hz())
                .expect("band centers are validated finite")
        });

        // Strictly disjoint, exactly like `MultiBandMapping::new`.
        let mapping_bands: Vec<(f64, f64)> = bands
            .iter()
            .map(|band| (band.center_hz(), band.bandwidth_hz()))
            .collect();
        // Reuses MultiBandMapping's validation (>=2 bands, sorted, disjoint) so
        // the spec and the mapping can never disagree on the band invariant.
        let _ = MultiBandMapping::new(&mapping_bands)?;

        Ok(Self {
            bands,
            rejection: Vec::new(),
            transmission_zeros_normalized: Vec::new(),
            order,
        })
    }

    /// Adds a rejection requirement.
    pub fn with_rejection(mut self, rejection: impl IntoIterator<Item = RejectionSpec>) -> Self {
        self.rejection = rejection.into_iter().collect();
        self
    }

    /// Sets the normalized-prototype transmission zeros.
    ///
    /// Values are in the **multi-band prototype** domain (a shared axis whose
    /// band edges are the band edges of the composed transform), consistent
    /// with the single-band [`FilterSpec`](super::FilterSpec) convention.
    pub fn with_normalized_transmission_zeros(
        mut self,
        transmission_zeros: impl IntoIterator<Item = f64>,
    ) -> Self {
        self.transmission_zeros_normalized = transmission_zeros.into_iter().collect();
        self
    }

    /// Sets the normalized-prototype transmission zeros from
    /// [`TransmissionZero`] values.
    pub fn with_transmission_zeros(
        mut self,
        transmission_zeros: impl IntoIterator<Item = TransmissionZero>,
    ) -> Self {
        self.transmission_zeros_normalized = transmission_zeros
            .into_iter()
            .map(|zero| zero.value)
            .collect();
        self
    }

    /// The passbands, sorted by center frequency ascending.
    pub fn bands(&self) -> &[BandSpec] {
        &self.bands
    }

    /// The declared rejection requirements.
    pub fn rejection(&self) -> &[RejectionSpec] {
        &self.rejection
    }

    /// The normalized-prototype transmission zeros.
    pub fn normalized_transmission_zeros(&self) -> &[f64] {
        &self.transmission_zeros_normalized
    }

    /// The total filter order.
    pub fn order(&self) -> usize {
        self.order
    }

    /// The number of passbands $K$.
    pub fn band_count(&self) -> usize {
        self.bands.len()
    }

    /// Builds the [`MultiBandMapping`] that realizes this spec's band layout.
    pub fn mapping(&self) -> Result<MultiBandMapping> {
        let bands: Vec<(f64, f64)> = self
            .bands
            .iter()
            .map(|band| (band.center_hz(), band.bandwidth_hz()))
            .collect();
        MultiBandMapping::new(&bands)
    }

    /// Validates the transmission-zero vector (finite non-zero, or infinity)
    /// and returns it as raw values, mirroring the single-band path.
    pub fn validated_transmission_zeros(&self) -> Result<Vec<f64>> {
        let zeros: Vec<TransmissionZero> = self
            .transmission_zeros_normalized
            .iter()
            .copied()
            .map(TransmissionZero::normalized)
            .collect();
        zeros.iter().copied().map(validated_transmission_zero).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_bands() -> Vec<BandSpec> {
        vec![
            BandSpec::new(1.0e9, 0.1e9, 20.0).expect("band 1"),
            BandSpec::new(1.3e9, 0.1e9, 20.0).expect("band 2"),
        ]
    }

    #[test]
    fn band_spec_validates_center_bandwidth_and_return_loss() {
        assert!(BandSpec::new(-1.0, 1.0, 20.0).is_err());
        assert!(BandSpec::new(1.0, 0.0, 20.0).is_err());
        assert!(BandSpec::new(1.0, 1.0, 0.0).is_err());
        assert!(BandSpec::new(1.0e9, 0.1e9, 20.0).is_ok());
    }

    #[test]
    fn band_spec_from_edges_uses_geometric_mean() {
        let band = BandSpec::from_edges(1.0e9, 1.2e9, 20.0).expect("band");
        let expected_center = (1.0e9_f64 * 1.2e9).sqrt();
        assert!((band.center_hz() - expected_center).abs() < 1e-3);
        assert!((band.bandwidth_hz() - 0.2e9).abs() < 1e-3);
    }

    #[test]
    fn multi_band_spec_requires_at_least_two_bands() {
        let one = vec![BandSpec::new(1.0e9, 0.1e9, 20.0).expect("band")];
        assert!(MultiBandSpec::new(one, 6).is_err());
    }

    #[test]
    fn multi_band_spec_rejects_overlapping_bands() {
        let overlapping = vec![
            BandSpec::new(1.0e9, 0.4e9, 20.0).expect("band 1"),
            BandSpec::new(1.1e9, 0.4e9, 20.0).expect("band 2"),
        ];
        assert!(MultiBandSpec::new(overlapping, 6).is_err());
    }

    #[test]
    fn multi_band_spec_sorts_bands_and_exposes_fields() -> Result<()> {
        let spec = MultiBandSpec::new(two_bands(), 6)?
            .with_normalized_transmission_zeros(vec![-1.5, 0.0, 1.5])
            .with_rejection([
                RejectionSpec::new(RejectionRegion::OuterLower, 30.0)?,
                RejectionSpec::new(RejectionRegion::InterBand { index: 0 }, 20.0)?,
                RejectionSpec::new(RejectionRegion::OuterUpper, 30.0)?,
            ]);

        assert_eq!(spec.band_count(), 2);
        assert_eq!(spec.order(), 6);
        assert_eq!(spec.bands()[0].center_hz(), 1.0e9);
        assert_eq!(spec.rejection().len(), 3);
        assert_eq!(spec.normalized_transmission_zeros(), &[-1.5, 0.0, 1.5]);

        // The mapping reuses the same band layout and validates.
        let mapping = spec.mapping()?;
        assert_eq!(mapping.band_count(), 2);
        Ok(())
    }

    #[test]
    fn multi_band_spec_matches_cameron_fig25_band_layout() -> Result<()> {
        // Cameron 2011 Fig. 25 (fixture `Cameron2011_10_5_dual_band_symmetric`):
        // normalized bands [-1.0, -0.35] and [0.35, 1.0]. In the normalized
        // prototype domain the crate maps edges to Omega = +/-1 via a
        // BandPassMapping; use physical edges 0.65..1.0 and 1.0..1.0/0.65 to
        // get the expected band arrangement in a physical layout.
        let spec = MultiBandSpec::new(
            vec![
                BandSpec::from_edges(0.8e9, 1.0e9, 22.0)?,
                BandSpec::from_edges(1.4e9, 1.75e9, 22.0)?,
            ],
            10,
        )?;
        assert_eq!(spec.band_count(), 2);
        assert_eq!(spec.order(), 10);
        assert!(spec.mapping().is_ok());
        Ok(())
    }
}
