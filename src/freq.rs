use crate::error::{MfsError, Result};
use crate::spec::{TransmissionZero, TransmissionZeroDomain};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalizedSample {
    pub omega: f64,
}

pub trait FrequencyPlan {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample>;

    fn map_normalized_to_hz(&self, sample: NormalizedSample) -> Result<f64>;

    fn map_grid_hz_to_normalized(&self, grid: &FrequencyGrid) -> Result<Vec<NormalizedSample>> {
        grid.as_slice()
            .iter()
            .copied()
            .map(|frequency_hz| self.map_hz_to_normalized(frequency_hz))
            .collect()
    }
}

pub fn normalize_transmission_zero(
    zero: TransmissionZero,
    plan: &impl FrequencyPlan,
) -> Result<f64> {
    if !zero.value.is_finite() {
        return Err(MfsError::InvalidTransmissionZero(
            "transmission zero must be finite".to_string(),
        ));
    }

    match zero.domain {
        TransmissionZeroDomain::Normalized => Ok(zero.value),
        TransmissionZeroDomain::PhysicalHz => Ok(plan.map_hz_to_normalized(zero.value)?.omega),
    }
}

pub fn normalize_transmission_zeros(
    zeros: &[TransmissionZero],
    plan: &impl FrequencyPlan,
) -> Result<Vec<f64>> {
    zeros
        .iter()
        .copied()
        .map(|zero| normalize_transmission_zero(zero, plan))
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowPassPlan {
    cutoff_hz: f64,
}

impl LowPassPlan {
    pub fn new(cutoff_hz: f64) -> Result<Self> {
        if !cutoff_hz.is_finite() || cutoff_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "cutoff must be > 0, got {cutoff_hz}"
            )));
        }
        Ok(Self { cutoff_hz })
    }

    pub fn cutoff_hz(&self) -> f64 {
        self.cutoff_hz
    }
}

impl FrequencyPlan for LowPassPlan {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz < 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "sample must be >= 0, got {frequency_hz}"
            )));
        }

        Ok(NormalizedSample {
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandPassPlan {
    center_hz: f64,
    bandwidth_hz: f64,
}

impl BandPassPlan {
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

    pub fn center_hz(&self) -> f64 {
        self.center_hz
    }

    pub fn bandwidth_hz(&self) -> f64 {
        self.bandwidth_hz
    }
}

impl FrequencyPlan for BandPassPlan {
    fn map_hz_to_normalized(&self, frequency_hz: f64) -> Result<NormalizedSample> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(MfsError::InvalidFrequency(format!(
                "band-pass sample must be > 0, got {frequency_hz}"
            )));
        }

        let fractional_bw = self.bandwidth_hz / self.center_hz;
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

#[derive(Debug, Clone, PartialEq)]
pub struct FrequencyGrid {
    samples_hz: Vec<f64>,
}

impl FrequencyGrid {
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
        let plan = LowPassPlan::new(2.0e9)?;
        let normalized = plan.map_hz_to_normalized(1.0e9)?;
        approx_eq(normalized.omega, 0.5, 1e-12);

        let restored = plan.map_normalized_to_hz(normalized)?;
        approx_eq(restored, 1.0e9, 1e-3);
        Ok(())
    }

    #[test]
    fn band_pass_center_frequency_maps_to_zero() -> Result<()> {
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
        let normalized = plan.map_hz_to_normalized(plan.center_hz())?;
        approx_eq(normalized.omega, 0.0, 1e-12);
        Ok(())
    }

    #[test]
    fn band_pass_inverse_matches_python_positive_branch() -> Result<()> {
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
        let sample = NormalizedSample { omega: 1.25 };
        let mapped = plan.map_normalized_to_hz(sample)?;

        let expected = (sample.omega * plan.bandwidth_hz()
            + ((sample.omega * plan.bandwidth_hz()).powi(2) + 4.0 * plan.center_hz().powi(2))
                .sqrt())
            / 2.0;
        approx_eq(mapped, expected, 1e-3);
        Ok(())
    }

    #[test]
    fn grid_mapping_preserves_length() -> Result<()> {
        let plan = LowPassPlan::new(1.0e9)?;
        let grid = FrequencyGrid::linspace(0.5e9, 1.5e9, 5)?;
        let normalized = plan.map_grid_hz_to_normalized(&grid)?;

        assert_eq!(normalized.len(), 5);
        approx_eq(normalized[0].omega, 0.5, 1e-12);
        approx_eq(normalized[4].omega, 1.5, 1e-12);
        Ok(())
    }

    #[test]
    fn transmission_zero_normalization_uses_frequency_plan() -> Result<()> {
        let plan = BandPassPlan::new(6.75e9, 300.0e6)?;
        let normalized = normalize_transmission_zero(TransmissionZero::physical_hz(6.9e9), &plan)?;

        approx_eq(normalized, 0.9891304347826066, 1e-12);
        Ok(())
    }

    #[test]
    fn transmission_zero_normalization_rejects_non_finite_values() {
        let plan = LowPassPlan::new(1.0).expect("valid plan");
        let error = normalize_transmission_zero(TransmissionZero::normalized(f64::NAN), &plan)
            .expect_err("non-finite transmission zero must fail");

        assert!(matches!(error, MfsError::InvalidTransmissionZero(_)));
    }
}
