//! Frequency-domain helpers for prototype and physical-frequency workflows.
//!
//! Convention:
//! - transmission zeros stored in [`crate::spec::FilterSpec`] are already
//!   normalized prototype values
//! - [`FrequencyMapping`] is used to relate physical-frequency grids to the
//!   normalized prototype axis
//! - if your transmission zeros start in Hz, convert them before building a
//!   spec instead of expecting synthesis to normalize them implicitly
//! - prototype sign conventions: low-pass maps `f` to `f/fc`; high-pass maps
//!   `f` to `-fc/f`; band-stop maps `f` to `-1 / Ω_BPF(f)`. Forward mappings
//!   clamp at their singular physical points (`f = 0` for high-pass,
//!   `f = f0` for band-stop) to `-SINGULAR_CLAMP` instead of failing,
//!   matching the convention used by the cross-language benchmark standard.

use crate::error::{MfsError, Result};
use crate::spec::TransmissionZero;

/// One sample expressed in the normalized prototype frequency domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalizedSample {
    /// Normalized angular frequency used by approximation and response code.
    pub omega: f64,
}

/// Maps between physical frequency units and normalized prototype coordinates.
pub trait FrequencyMapping {
    /// Converts one physical-frequency sample into prototype coordinates.
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample>;

    /// Converts one normalized prototype sample back into physical frequency.
    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64>;

    /// Converts every sample in a grid through the mapping in one pass.
    fn map_grid_hz_to_normalized(&self, grid: &FrequencyGrid) -> Result<Vec<NormalizedSample>> {
        grid.as_slice()
            .iter()
            .copied()
            .map(|frequency_hz| self.map_hz_to_normalized(frequency_hz))
            .collect()
    }
}

/// Validates one normalized transmission zero and returns its stored coordinate.
///
/// Explicit zeros at infinity are accepted. `NaN` remains invalid.
pub fn validated_transmission_zero(zero: TransmissionZero) -> Result<f64> {
    if zero.value.is_nan() {
        return Err(MfsError::InvalidTransmissionZero(
            "transmission zero must be a real value or infinity".to_string(),
        ));
    }

    Ok(zero.value)
}

/// Collects and validates the normalized transmission zeros already stored in a spec.
pub fn validated_transmission_zeros(zeros: &[TransmissionZero]) -> Result<Vec<f64>> {
    zeros
        .iter()
        .copied()
        .map(validated_transmission_zero)
        .collect()
}

/// Low-pass frequency mapping parameterized by a single cutoff frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowPassMapping {
    cutoff_hz: f64,
}

impl LowPassMapping {
    /// Creates a validated low-pass mapping.
    pub fn new(cutoff_hz: f64) -> Result<Self> {
        if !cutoff_hz.is_finite() || cutoff_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "cutoff must be > 0, got {cutoff_hz}"
            )));
        }
        Ok(Self { cutoff_hz })
    }

    /// Returns the cutoff frequency used for normalization.
    pub fn cutoff_hz(&self) -> f64 {
        self.cutoff_hz
    }
}

impl FrequencyMapping for LowPassMapping {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz < 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "sample must be >= 0, got {frequency_hz}"
            )));
        }

        Ok(NormalizedSample {
            // A low-pass prototype normalizes directly by the cutoff frequency.
            omega: frequency_hz / self.cutoff_hz,
        })
    }

    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64> {
        if !sample.omega.is_finite() || sample.omega < 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "normalized low-pass frequency must be >= 0, got {}",
                sample.omega
            )));
        }

        Ok(sample.omega * self.cutoff_hz)
    }
}

/// Band-pass frequency mapping parameterized by center frequency and bandwidth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandPassMapping {
    center_hz: f64,
    bandwidth_hz: f64,
}

impl BandPassMapping {
    /// Creates a validated band-pass mapping.
    pub fn new(center_hz: f64, bandwidth_hz: f64) -> Result<Self> {
        if !center_hz.is_finite() || center_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "center frequency must be > 0, got {center_hz}"
            )));
        }
        if !bandwidth_hz.is_finite() || bandwidth_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "bandwidth must be > 0, got {bandwidth_hz}"
            )));
        }

        Ok(Self {
            center_hz,
            bandwidth_hz,
        })
    }

    /// Returns the center frequency of the mapped passband.
    pub fn center_hz(&self) -> f64 {
        self.center_hz
    }

    /// Returns the absolute 3 dB bandwidth used by the mapping.
    pub fn bandwidth_hz(&self) -> f64 {
        self.bandwidth_hz
    }
}

impl FrequencyMapping for BandPassMapping {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "band-pass sample must be > 0, got {frequency_hz}"
            )));
        }

        let fractional_bw = self.bandwidth_hz / self.center_hz;
        // This is the standard low-pass to band-pass frequency transformation.
        let omega =
            (1.0 / fractional_bw) * (frequency_hz / self.center_hz - self.center_hz / frequency_hz);
        Ok(NormalizedSample { omega })
    }

    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64> {
        if !sample.omega.is_finite() {
            return Err(MfsError::InvalidFrequency(
                "normalized band-pass frequency must be finite".to_string(),
            ));
        }

        // Use the positive-frequency branch of the inverse quadratic mapping.
        let discriminant =
            (sample.omega * self.bandwidth_hz).powi(2) + 4.0 * self.center_hz.powi(2);
        let frequency_hz = (sample.omega * self.bandwidth_hz + discriminant.sqrt()) / 2.0;

        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "band-pass inverse mapping produced invalid value {frequency_hz}"
            )));
        }

        Ok(frequency_hz)
    }
}

/// Normalized-magnitude clamp used by the high-pass and band-stop forward
/// mappings at their singular physical points (`f = 0`, `f = f0`), where the
/// exact image has magnitude infinity. The clamp takes the sign of the side
/// the forward map approaches from (`-SINGULAR_CLAMP` for both singular
/// points). It keeps response evaluation finite without changing the filter
/// shape in the neighborhood of the singularity.
pub const SINGULAR_CLAMP: f64 = 1.0e6;

/// High-pass frequency mapping parameterized by a single cutoff frequency.
///
/// The prototype sign convention is the frequency inversion
/// `ω = -fc / f`: passband frequencies (`f >= fc`) map to `ω <= -1` and
/// `ω(f = fc) = -1` exactly, so the low-pass prototype's stopband
/// (`|ω| > 1`) becomes the high-pass passband. `f = 0` is singular and
/// clamps to `-SINGULAR_CLAMP` (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HighPassMapping {
    cutoff_hz: f64,
}

impl HighPassMapping {
    /// Creates a validated high-pass mapping.
    pub fn new(cutoff_hz: f64) -> Result<Self> {
        if !cutoff_hz.is_finite() || cutoff_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "cutoff must be > 0, got {cutoff_hz}"
            )));
        }
        Ok(Self { cutoff_hz })
    }

    /// Returns the cutoff frequency used for normalization.
    pub fn cutoff_hz(&self) -> f64 {
        self.cutoff_hz
    }
}

impl FrequencyMapping for HighPassMapping {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz < 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "sample must be >= 0, got {frequency_hz}"
            )));
        }
        // f = 0 maps to -inf under -fc/f; clamp instead of failing.
        let omega = if frequency_hz == 0.0 {
            -SINGULAR_CLAMP
        } else {
            -self.cutoff_hz / frequency_hz
        };
        Ok(NormalizedSample { omega })
    }

    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64> {
        // The image of the forward mapping is (-inf, 0); anything else has
        // no physical high-pass preimage.
        if !sample.omega.is_finite() || sample.omega >= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "normalized high-pass frequency must be < 0 and finite, got {}",
                sample.omega
            )));
        }
        // Inputs beyond the clamp bound (e.g. a grid endpoint exactly at
        // -SINGULAR_CLAMP) invert to 0 or a negative frequency; clamp to the
        // nearest physical value instead of failing.
        let omega = sample.omega.max(-SINGULAR_CLAMP);
        let frequency_hz = -self.cutoff_hz / omega;

        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "high-pass inverse mapping produced invalid value {frequency_hz}"
            )));
        }

        Ok(frequency_hz)
    }
}

/// Band-stop frequency mapping parameterized by stopband center frequency and
/// stopband width.
///
/// The mapping is the inverse of the band-pass transformation,
/// `ω = -1 / Ω_BPF(f)`, with the same geometric-mean center convention: the
/// stopband edges map to `+1` (lower edge) and `-1` (upper edge) exactly,
/// `f = f0` is singular (`|ω| -> ∞` from both sides, clamped to `-SINGULAR_CLAMP`), and the
/// two passbands map to `|ω| < 1`. The forward image covers every `ω ≠ 0`,
/// so the inverse is defined for both signs of `ω`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandStopMapping {
    center_hz: f64,
    stopband_width_hz: f64,
}

impl BandStopMapping {
    /// Creates a validated band-stop mapping.
    ///
    /// `center_hz` is the geometric mean of the stopband edges and
    /// `stopband_width_hz` their difference, mirroring
    /// [`BandPassMapping::new`] (see [`edges_to_center_bandwidth`]).
    pub fn new(center_hz: f64, stopband_width_hz: f64) -> Result<Self> {
        if !center_hz.is_finite() || center_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "stopband center frequency must be > 0, got {center_hz}"
            )));
        }
        if !stopband_width_hz.is_finite() || stopband_width_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "stopband width must be > 0, got {stopband_width_hz}"
            )));
        }

        Ok(Self {
            center_hz,
            stopband_width_hz,
        })
    }

    /// Returns the stopband center frequency of the mapped stopband.
    pub fn center_hz(&self) -> f64 {
        self.center_hz
    }

    /// Returns the stopband width used by the mapping.
    pub fn stopband_width_hz(&self) -> f64 {
        self.stopband_width_hz
    }
}

impl FrequencyMapping for BandStopMapping {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz < 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "sample must be >= 0, got {frequency_hz}"
            )));
        }

        let fractional_bw = self.stopband_width_hz / self.center_hz;
        // The band-stop transformation is the reciprocal (with a sign flip)
        // of the band-pass transformation. f = f0 is its singular point;
        // frequencies arbitrarily close to it map to arbitrarily large
        // |ω|, so clamp to the singular bound instead of overflowing.
        let bp_omega =
            (1.0 / fractional_bw) * (frequency_hz / self.center_hz - self.center_hz / frequency_hz);
        let omega = if bp_omega == 0.0 {
            -SINGULAR_CLAMP
        } else {
            (-1.0 / bp_omega).clamp(-SINGULAR_CLAMP, SINGULAR_CLAMP)
        };
        Ok(NormalizedSample { omega })
    }

    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64> {
        // ω = 0 is the limit of the forward mapping as f -> 0+ (ω -> 0+)
        // or f -> +inf (ω -> 0-); it has no finite preimage.
        if !sample.omega.is_finite() || sample.omega == 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "normalized band-stop frequency must be nonzero and finite, got {}",
                sample.omega
            )));
        }
        // Inputs beyond the clamp bound (e.g. a grid endpoint exactly at
        // ±SINGULAR_CLAMP) invert past the stopband center, where the
        // mapping is singular; clamp to the nearest physical value instead
        // of failing.
        let omega = sample.omega.clamp(-SINGULAR_CLAMP, SINGULAR_CLAMP);

        // Both branches solve f² + (W/ω)·f − c² = 0; the root forms below
        // are algebraically identical, chosen to avoid cancellation.
        let c = self.center_hz;
        let w = self.stopband_width_hz;
        let discriminant = (w * w + 4.0 * omega * omega * c * c).sqrt();
        let frequency_hz = if omega > 0.0 {
            // f < c (lower passband / lower stopband half); rationalized
            // form stays accurate as ω -> 0+.
            2.0 * omega * c * c / (w + discriminant)
        } else {
            // f > c (upper stopband half / upper passband); no cancellation.
            (w + discriminant) / (-2.0 * omega)
        };

        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "band-stop inverse mapping produced invalid value {frequency_hz}"
            )));
        }

        Ok(frequency_hz)
    }
}

/// Multi-band frequency mapping: one band-pass prototype per passband.
///
/// `MultiBandMapping` holds a disjoint, center-ascending list of
/// [`BandPassMapping`]s — one per passband — and is the mapping-layer
/// counterpart of a multi-passband spec. Each band's physical edges map to
/// `-1` (lower) and `+1` (upper) exactly in that band's own normalized
/// domain, under the same geometric-mean convention as
/// [`BandPassMapping::new`].
///
/// The forward map is piecewise: a physical frequency is assigned to the
/// band whose prototype-transform distance `|Ω_b(f)|` is smallest (ties go
/// to the lower band), so every `f > 0` — passbands, inter-band gaps, and
/// the out-of-band regions — has a well-defined owner. The inverse is
/// band-local: it inverts one band's transform, so it takes a band index.
/// Because the band index is required, this type deliberately does not
/// implement [`FrequencyMapping`]; the trait's inverse signature has no
/// place for it.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiBandMapping {
    /// Sorted by center frequency ascending; adjacent bands are disjoint.
    bands: Vec<BandPassMapping>,
}

impl MultiBandMapping {
    /// Creates a mapping from `(center_hz, bandwidth_hz)` pairs.
    ///
    /// At least two bands are required. Bands are validated like
    /// [`BandPassMapping::new`], sorted by center frequency, and must be
    /// disjoint: each band's upper edge must lie strictly below the next
    /// band's lower edge, because a shared edge would make the piecewise
    /// forward map ambiguous at that point.
    pub fn new(bands: &[(f64, f64)]) -> Result<Self> {
        if bands.len() < 2 {
            return Err(MfsError::InvalidFrequency(format!(
                "multi-band mapping requires at least 2 bands, got {}",
                bands.len()
            )));
        }

        let mut parsed = Vec::with_capacity(bands.len());
        for &(center_hz, bandwidth_hz) in bands {
            parsed.push(BandPassMapping::new(center_hz, bandwidth_hz)?);
        }
        parsed.sort_by(|a, b| {
            a.center_hz()
                .partial_cmp(&b.center_hz())
                .expect("centers are validated finite")
        });

        for window in parsed.windows(2) {
            let (lo, hi) = (window[0], window[1]);
            let (_, lo_upper) = band_edge_pair(lo.center_hz(), lo.bandwidth_hz());
            let (hi_lower, _) = band_edge_pair(hi.center_hz(), hi.bandwidth_hz());
            // Both edges are validated finite, so the total comparison is
            // well-defined.
            if lo_upper >= hi_lower {
                return Err(MfsError::InvalidFrequency(format!(
                    "multi-band bands must be disjoint: upper edge {} Hz of the band centered at {} Hz meets lower edge {} Hz of the band centered at {} Hz",
                    lo_upper,
                    lo.center_hz(),
                    hi_lower,
                    hi.center_hz()
                )));
            }
        }

        Ok(Self { bands: parsed })
    }

    /// Creates a mapping from `(f1_hz, f2_hz)` passband edge pairs, using
    /// [`edges_to_center_bandwidth`]'s geometric-mean convention (S4).
    pub fn from_edges(edges: &[(f64, f64)]) -> Result<Self> {
        if edges.len() < 2 {
            return Err(MfsError::InvalidFrequency(format!(
                "multi-band mapping requires at least 2 bands, got {}",
                edges.len()
            )));
        }

        let center_bandwidth: Vec<(f64, f64)> = edges
            .iter()
            .map(|&(f1_hz, f2_hz)| {
                if !f1_hz.is_finite() || !f2_hz.is_finite() || f1_hz <= 0.0 || f2_hz <= f1_hz {
                    return Err(MfsError::InvalidFrequency(format!(
                        "band edges must satisfy 0 < f1 < f2 (strictly ordered), got f1={f1_hz}, f2={f2_hz}"
                    )));
                }
                Ok(edges_to_center_bandwidth(f1_hz, f2_hz))
            })
            .collect::<Result<Vec<_>>>()?;

        Self::new(&center_bandwidth)
    }

    /// Returns the number of bands.
    pub fn band_count(&self) -> usize {
        self.bands.len()
    }

    /// Returns the bands in center-ascending order.
    pub fn bands(&self) -> &[BandPassMapping] {
        &self.bands
    }

    /// Returns the `(lower, upper)` physical edge frequencies of one band.
    pub fn band_edges(&self, band: usize) -> Result<(f64, f64)> {
        let band = self
            .bands
            .get(band)
            .ok_or_else(|| MfsError::InvalidFrequency(format!(
                "band index {band} out of range ({} bands)",
                self.bands.len()
            )))?;
        Ok(band_edge_pair(band.center_hz(), band.bandwidth_hz()))
    }

    /// Piecewise forward map: assigns `frequency_hz` to the band with the
    /// smallest prototype-transform distance `|Ω_b(f)|` and returns that
    /// band's index together with the band-local sample.
    pub fn map_hz_to_normalized(
        &self,
        frequency_hz: f64,
    ) -> Result<(usize, NormalizedSample)> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "multi-band sample must be > 0, got {frequency_hz}"
            )));
        }

        // One pass over the bands; a strict comparison keeps the lower band
        // (bands are center-ascending) on exact ties.
        let mut best: Option<(usize, NormalizedSample)> = None;
        for (index, band) in self.bands.iter().enumerate() {
            let sample = band.map_hz_to_normalized(frequency_hz)?;
            let closer = best
                .as_ref()
                .is_none_or(|(_, owner)| sample.omega.abs() < owner.omega.abs());
            if closer {
                best = Some((index, sample));
            }
        }

        best.ok_or_else(|| {
            MfsError::InvalidFrequency("multi-band mapping has no bands".to_string())
        })
    }

    /// Forward map with an explicit band (no nearest-band assignment).
    pub fn map_band_hz_to_normalized(
        &self,
        band: usize,
        frequency_hz: f64,
    ) -> Result<NormalizedSample> {
        let band = self
            .bands
            .get(band)
            .ok_or_else(|| MfsError::InvalidFrequency(format!(
                "band index {band} out of range ({} bands)",
                self.bands.len()
            )))?;
        band.map_hz_to_normalized(frequency_hz)
    }

    /// Band-local inverse: inverts one band's transform back to physical
    /// frequency.
    pub fn map_normalized_to_hz(
        &self,
        band: usize,
        sample: NormalizedSample,
    ) -> Result<f64> {
        let band = self
            .bands
            .get(band)
            .ok_or_else(|| MfsError::InvalidFrequency(format!(
                "band index {band} out of range ({} bands)",
                self.bands.len()
            )))?;
        band.map_normalized_to_hz(sample)
    }

    /// Converts every grid sample through the piecewise forward map.
    pub fn map_grid_hz_to_normalized(
        &self,
        grid: &FrequencyGrid,
    ) -> Result<Vec<(usize, NormalizedSample)>> {
        grid.as_slice()
            .iter()
            .copied()
            .map(|frequency_hz| self.map_hz_to_normalized(frequency_hz))
            .collect()
    }
}

/// Returns the `(lower, upper)` physical edges of a band-pass mapping with
/// geometric-mean center `center_hz` and bandwidth `bandwidth_hz`.
///
/// The single source of truth for passband edges (geometric-mean convention);
/// consumers that need the Ω=±1 edges derive them from this rather than the
/// arithmetic `center ∓ bandwidth/2`.
pub fn band_edge_pair(center_hz: f64, bandwidth_hz: f64) -> (f64, f64) {
    let spread = (bandwidth_hz * bandwidth_hz + 4.0 * center_hz * center_hz).sqrt();
    let lower = (-bandwidth_hz + spread) / 2.0;
    (lower, lower + bandwidth_hz)
}

/// Upper bound on the response-evaluation grid point count (review P2-11).
/// A single `points` value is a linear allocation; 100k samples keeps a
/// hostile request bounded without touching legitimate dense grids.
pub const MAX_GRID_POINTS: usize = 100_000;

/// Concrete frequency samples used when evaluating a synthesized response.
#[derive(Debug, Clone, PartialEq)]
pub struct FrequencyGrid {
    samples_hz: Vec<f64>,
}

impl FrequencyGrid {
    /// Builds an evenly spaced grid from `start_hz` to `stop_hz`, inclusive.
    pub fn linspace(start_hz: f64, stop_hz: f64, points: usize) -> Result<Self> {
        if points < 2 {
            return Err(MfsError::InvalidGridSize { points });
        }
        if points > MAX_GRID_POINTS {
            return Err(MfsError::InvalidGridSize { points });
        }
        if !start_hz.is_finite() || !stop_hz.is_finite() || stop_hz <= start_hz {
            return Err(MfsError::InvalidFrequency(format!(
                "invalid grid range: start={start_hz}, stop={stop_hz}"
            )));
        }

        let step = (stop_hz - start_hz) / (points.saturating_sub(1) as f64);
        let samples_hz = (0..points)
            .map(|index| start_hz + step * index as f64)
            .collect();

        Ok(Self { samples_hz })
    }

    /// Returns the underlying grid samples as a borrowed slice.
    pub fn as_slice(&self) -> &[f64] {
        &self.samples_hz
    }
}

/// Converts passband edge frequencies into (center, bandwidth).
///
/// The center is the **geometric mean** `√(f1·f2)`, the convention under which
/// `Ω(f1) = -1` and `Ω(f2) = +1` hold exactly in the band-pass mapping. This is
/// the single definition of the edge→center/bandwidth convention; all
/// edge-based entry points (`FilterDesign::bandpass_edges` and the CLI `design`
/// subcommand) go through it so the convention cannot drift. Input validity
/// (positive, ordered edges) is the caller's responsibility.
pub fn edges_to_center_bandwidth(f1_hz: f64, f2_hz: f64) -> (f64, f64) {
    let center_hz = (f1_hz * f2_hz).sqrt();
    let bandwidth_hz = f2_hz - f1_hz;
    (center_hz, bandwidth_hz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::FilterDesign;
    use crate::spec::TransmissionZero;

    #[test]
    fn edges_to_center_bandwidth_uses_geometric_mean_convention() -> Result<()> {
        let (f1, f2) = (6.6e9, 6.9e9);
        let (center, bw) = edges_to_center_bandwidth(f1, f2);
        approx_eq(center, (f1 * f2).sqrt(), 1e-6);
        approx_eq(bw, f2 - f1, 1e-6);
        // The high-level entry point must use the identical convention (S4):
        // this cross-check is what prevents the CLI / bandpass_edges from
        // drifting apart on the band-center definition.
        let d = FilterDesign::bandpass_edges(4, 20.0, f1, f2).synthesize()?;
        assert_eq!(d.center_hz(), Some(center));
        assert_eq!(d.bandwidth_hz(), Some(bw));
        Ok(())
    }

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn low_pass_forward_and_inverse_mapping_round_trip() -> Result<()> {
        let mapping = LowPassMapping::new(2.0e9)?;
        let normalized = mapping.map_hz_to_normalized(1.0e9)?;
        approx_eq(normalized.omega, 0.5, 1e-12);

        let restored = mapping.map_normalized_to_hz(normalized)?;
        approx_eq(restored, 1.0e9, 1e-3);
        Ok(())
    }

    #[test]
    fn band_pass_center_frequency_maps_to_zero() -> Result<()> {
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let normalized = mapping.map_hz_to_normalized(mapping.center_hz())?;
        approx_eq(normalized.omega, 0.0, 1e-12);
        Ok(())
    }

    #[test]
    fn band_pass_inverse_matches_python_positive_branch() -> Result<()> {
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;
        let sample = NormalizedSample { omega: 1.25 };
        let mapped = mapping.map_normalized_to_hz(sample)?;

        let expected = (sample.omega * mapping.bandwidth_hz()
            + ((sample.omega * mapping.bandwidth_hz()).powi(2)
                + 4.0 * mapping.center_hz().powi(2))
            .sqrt())
            / 2.0;
        approx_eq(mapped, expected, 1e-3);
        Ok(())
    }

    #[test]
    fn grid_mapping_preserves_length() -> Result<()> {
        let mapping = LowPassMapping::new(1.0e9)?;
        let grid = FrequencyGrid::linspace(0.5e9, 1.5e9, 5)?;
        let normalized = mapping.map_grid_hz_to_normalized(&grid)?;

        assert_eq!(normalized.len(), 5);
        approx_eq(normalized[0].omega, 0.5, 1e-12);
        approx_eq(normalized[4].omega, 1.5, 1e-12);
        Ok(())
    }

    #[test]
    fn validated_transmission_zero_returns_stored_value() -> Result<()> {
        let normalized =
            validated_transmission_zero(TransmissionZero::normalized(0.9891304347826066))?;
        approx_eq(normalized, 0.9891304347826066, 1e-12);
        Ok(())
    }

    #[test]
    fn validated_transmission_zero_accepts_infinity() -> Result<()> {
        let normalized = validated_transmission_zero(TransmissionZero::infinite())?;
        assert!(normalized.is_infinite());
        Ok(())
    }

    #[test]
    fn validated_transmission_zero_rejects_nan() {
        let error = validated_transmission_zero(TransmissionZero::normalized(f64::NAN))
            .expect_err("non-finite transmission zero must fail");
        assert!(matches!(error, MfsError::InvalidTransmissionZero(_)));
    }

    #[test]
    fn high_pass_cutoff_maps_to_minus_one() -> Result<()> {
        let mapping = HighPassMapping::new(2.0e9)?;
        let normalized = mapping.map_hz_to_normalized(mapping.cutoff_hz())?;
        approx_eq(normalized.omega, -1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn high_pass_dc_clamps_to_singular() -> Result<()> {
        let mapping = HighPassMapping::new(2.0e9)?;
        let normalized = mapping.map_hz_to_normalized(0.0)?;
        // The forward image is (-inf, 0]; f = 0 sits at the -inf end, so the
        // clamp is negative (deep in the prototype stopband), not positive.
        approx_eq(normalized.omega, -SINGULAR_CLAMP, 0.0);

        // The inverse of the clamped endpoint returns to the nearest
        // physical frequency instead of failing.
        let restored = mapping.map_normalized_to_hz(normalized)?;
        approx_eq(restored, 2.0e9 / SINGULAR_CLAMP, 1e-21);
        Ok(())
    }

    #[test]
    fn high_pass_forward_inverse_round_trip_both_branches() -> Result<()> {
        let mapping = HighPassMapping::new(1.0e9)?;
        // Stopband side (positive prototype frequency) and passband side
        // (negative prototype frequency) must both invert.
        for frequency_hz in [0.25e9, 0.75e9, 1.2e9, 3.0e9] {
            let normalized = mapping.map_hz_to_normalized(frequency_hz)?;
            let restored = mapping.map_normalized_to_hz(normalized)?;
            approx_eq(restored, frequency_hz, frequency_hz * 1e-9);
        }
        // ω = -2.5 inverts to fc / 2.5 exactly (frequency inversion).
        let restored = mapping.map_normalized_to_hz(NormalizedSample { omega: -2.5 })?;
        approx_eq(restored, 0.4e9, 1e-12);
        Ok(())
    }

    #[test]
    fn high_pass_inverse_rejects_zero_and_nonfinite() {
        let mapping = HighPassMapping::new(1.0e9).expect("valid mapping");
        assert!(mapping
            .map_normalized_to_hz(NormalizedSample { omega: 0.0 })
            .is_err());
        assert!(mapping
            .map_normalized_to_hz(NormalizedSample { omega: f64::NAN })
            .is_err());
    }

    #[test]
    fn band_stop_center_maps_to_minus_clamp() -> Result<()> {
        let mapping = BandStopMapping::new(12.0e9, 200.0e6)?;
        let normalized = mapping.map_hz_to_normalized(mapping.center_hz())?;
        approx_eq(normalized.omega, -SINGULAR_CLAMP, 0.0);

        // The inverse of the clamped endpoint returns to the stopband center.
        let restored = mapping.map_normalized_to_hz(normalized)?;
        approx_eq(restored, 12.0e9, 12.0e9 * 1e-6);
        Ok(())
    }

    #[test]
    fn band_stop_edges_map_to_plus_minus_one() -> Result<()> {
        let (f1, f2) = (11.9e9, 12.1e9);
        let (center, width) = edges_to_center_bandwidth(f1, f2);
        let mapping = BandStopMapping::new(center, width)?;

        let lower = mapping.map_hz_to_normalized(f1)?;
        let upper = mapping.map_hz_to_normalized(f2)?;
        // Geometric-mean convention (S4): the lower stopband edge (f1) hits
        // +1 and the upper edge (f2) hits -1 exactly, because
        // ω = -1/Ω_BPF flips the band-pass edge signs.
        approx_eq(lower.omega, 1.0, 1e-9);
        approx_eq(upper.omega, -1.0, 1e-9);
        Ok(())
    }

    #[test]
    fn band_stop_forward_inverse_round_trip() -> Result<()> {
        let mapping = BandStopMapping::new(12.0e9, 200.0e6)?;
        // Samples across the stopband and both passbands.
        for frequency_hz in [9.0e9, 11.5e9, 11.95e9, 12.05e9, 12.5e9, 15.0e9] {
            let normalized = mapping.map_hz_to_normalized(frequency_hz)?;
            let restored = mapping.map_normalized_to_hz(normalized)?;
            approx_eq(restored, frequency_hz, frequency_hz * 1e-9);
        }
        Ok(())
    }

    #[test]
    fn band_stop_inverse_rejects_zero_and_nonfinite() {
        let mapping = BandStopMapping::new(12.0e9, 200.0e6).expect("valid mapping");
        // ω = 0 is the limit of the forward map as f -> 0+ or f -> +inf; it
        // has no finite preimage. Both signs of nonzero ω are physical.
        assert!(mapping
            .map_normalized_to_hz(NormalizedSample { omega: 0.0 })
            .is_err());
        assert!(mapping
            .map_normalized_to_hz(NormalizedSample { omega: f64::NAN })
            .is_err());
        // Positive ω (lower passband / lower stopband half) must invert.
        let lower_passband = mapping
            .map_normalized_to_hz(NormalizedSample { omega: 0.5 })
            .expect("positive omega is a valid band-stop preimage");
        assert!(lower_passband < 12.0e9);
    }

    #[test]
    fn multi_band_requires_at_least_two_bands() {
        assert!(MultiBandMapping::new(&[]).is_err());
        assert!(MultiBandMapping::new(&[(6.75e9, 300.0e6)]).is_err());
        assert!(MultiBandMapping::from_edges(&[(6.6e9, 6.9e9)]).is_err());
    }

    #[test]
    fn multi_band_from_edges_uses_geometric_mean_convention() -> Result<()> {
        let edges = [(6.6e9, 6.9e9), (11.9e9, 12.1e9)];
        let mapping = MultiBandMapping::from_edges(&edges)?;

        for (index, (band, (f1, f2))) in mapping.bands().iter().zip(edges).enumerate() {
            let (center, bw) = edges_to_center_bandwidth(f1, f2);
            // S4 invariant: the edge-based entry point must produce exactly
            // the geometric-mean center/bandwidth convention.
            assert_eq!(band.center_hz(), center);
            assert_eq!(band.bandwidth_hz(), bw);

            // Edges recover to the physical frequencies within relative 1e-6.
            let (lo, hi) = mapping.band_edges(index)?;
            approx_eq(lo, f1, f1 * 1e-6);
            approx_eq(hi, f2, f2 * 1e-6);
        }
        Ok(())
    }

    #[test]
    fn multi_band_band_edges_map_to_minus_plus_one() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;

        for (index, (f1, f2)) in
            [(6.6e9, 6.9e9), (11.9e9, 12.1e9)].into_iter().enumerate()
        {
            // Geometric-mean convention (S4): each band's edges hit -1/+1
            // exactly in its own normalized domain.
            let lower = mapping.map_band_hz_to_normalized(index, f1)?;
            let upper = mapping.map_band_hz_to_normalized(index, f2)?;
            approx_eq(lower.omega, -1.0, 1e-9);
            approx_eq(upper.omega, 1.0, 1e-9);
        }
        Ok(())
    }

    #[test]
    fn multi_band_per_band_round_trip() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;

        // In-band samples round-trip through the owning band.
        for (owner, samples) in [
            (0, [6.55e9, 6.7e9, 6.75e9, 6.8e9, 6.95e9]),
            (1, [11.85e9, 11.95e9, 12.0e9, 12.05e9, 12.15e9]),
        ] {
            for frequency_hz in samples {
                let (assigned, normalized) = mapping.map_hz_to_normalized(frequency_hz)?;
                assert_eq!(assigned, owner);
                let restored = mapping.map_normalized_to_hz(owner, normalized)?;
                approx_eq(restored, frequency_hz, frequency_hz * 1e-9);
            }
        }
        Ok(())
    }

    #[test]
    fn multi_band_nearest_band_assignment() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;

        // (frequency, expected owner): the owner is the band with the
        // smallest prototype-transform distance |Ω_b(f)|.
        for (frequency_hz, expected) in [
            (1.0e9, 0),   // far below: band 0 (wider) dominates
            (9.0e9, 0),   // inter-band gap, closer to band 0
            (13.0e9, 1),  // gap above band 1
            (25.0e9, 0),  // far above: |Ω_1| grows faster (narrower band),
                          // so the wider band 0 reclaims ownership
        ] {
            let (owner, _) = mapping.map_hz_to_normalized(frequency_hz)?;
            assert_eq!(owner, expected, "owner of {frequency_hz} Hz");
        }
        Ok(())
    }

    #[test]
    fn multi_band_grid_round_trip() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;
        let grid = FrequencyGrid::linspace(6.0e9, 13.0e9, 7)?;
        let samples = mapping.map_grid_hz_to_normalized(&grid)?;

        assert_eq!(samples.len(), 7);
        for ((owner, normalized), &frequency_hz) in samples.iter().zip(grid.as_slice()) {
            let restored = mapping.map_normalized_to_hz(*owner, *normalized)?;
            approx_eq(restored, frequency_hz, frequency_hz * 1e-9);
        }
        Ok(())
    }

    #[test]
    fn multi_band_rejects_overlapping_bands() {
        // Band A spans 6.6516–6.9516 GHz, band B 6.7516–7.0516 GHz: the
        // shared-edge region would make the piecewise map ambiguous.
        assert!(MultiBandMapping::new(&[(6.8e9, 300.0e6), (6.9e9, 300.0e6)]).is_err());
    }

    #[test]
    fn multi_band_sorts_unsorted_input() -> Result<()> {
        let mapping =
            MultiBandMapping::new(&[(11.9995833e9, 200.0e6), (6.7483335e9, 300.0e6)])?;
        assert_eq!(mapping.bands().len(), 2);
        assert!(mapping.bands()[0].center_hz() < mapping.bands()[1].center_hz());
        Ok(())
    }

    #[test]
    fn multi_band_forward_rejects_nonpositive_and_nonfinite() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;
        assert!(mapping.map_hz_to_normalized(0.0).is_err());
        assert!(mapping.map_hz_to_normalized(-1.0e9).is_err());
        assert!(mapping.map_hz_to_normalized(f64::NAN).is_err());
        Ok(())
    }

    #[test]
    fn multi_band_band_index_out_of_range_errors() -> Result<()> {
        let mapping = MultiBandMapping::from_edges(&[(6.6e9, 6.9e9), (11.9e9, 12.1e9)])?;
        assert!(mapping.band_edges(2).is_err());
        assert!(mapping
            .map_normalized_to_hz(
                2,
                NormalizedSample { omega: 0.5 }
            )
            .is_err());
        assert!(mapping.map_band_hz_to_normalized(1, -1.0e9).is_err());
        Ok(())
    }

    #[test]
    fn linspace_rejects_points_above_max() {
        // P2-11: the choke point must reject an OOM-hostile point count for
        // every entry (not just the JSON validator).
        let error = FrequencyGrid::linspace(1.0, 2.0, MAX_GRID_POINTS + 1)
            .expect_err("points above MAX_GRID_POINTS must be rejected");
        assert!(matches!(error, MfsError::InvalidGridSize { .. }));
        // The boundary itself is still accepted.
        assert!(FrequencyGrid::linspace(1.0, 2.0, MAX_GRID_POINTS).is_ok());
    }

}
