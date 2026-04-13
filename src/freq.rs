//! Frequency-domain helpers for prototype and physical-frequency workflows.
//!
//! Convention:
//! - transmission zeros stored in [`crate::spec::FilterSpec`] are already
//!   normalized prototype values
//! - [`FrequencyMapping`] is used to relate physical-frequency grids to the
//!   normalized prototype axis
//! - if your transmission zeros start in Hz, convert them before building a
//!   spec instead of expecting synthesis to normalize them implicitly

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
pub fn validated_transmission_zero(zero: TransmissionZero) -> Result<f64> {
    if !zero.value.is_finite() {
        return Err(MfsError::InvalidTransmissionZero(
            "transmission zero must be finite".to_string(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::TransmissionZero;

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
            + ((sample.omega * mapping.bandwidth_hz()).powi(2) + 4.0 * mapping.center_hz().powi(2))
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
        let normalized = validated_transmission_zero(TransmissionZero::normalized(0.9891304347826066))?;
        approx_eq(normalized, 0.9891304347826066, 1e-12);
        Ok(())
    }

    #[test]
    fn validated_transmission_zero_rejects_non_finite_values() {
        let error = validated_transmission_zero(TransmissionZero::normalized(f64::NAN))
            .expect_err("non-finite transmission zero must fail");

        assert!(matches!(error, MfsError::InvalidTransmissionZero(_)));
    }
}
