use std::fs;
use std::path::Path;

use crate::error::{MfsError, Result};

use super::FilterSpecBuilder;

/// Transmission-zero and high-level specification types.
///
/// `FilterSpec` follows a normalized-prototype convention: transmission zeros
/// are stored as normalized prototype coordinates, not physical Hz values.
/// Convert physical zeros before constructing the spec.

/// Defines one finite transmission zero in normalized low-pass prototype coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionZero {
    /// Numeric value of the zero in normalized low-pass prototype coordinates.
    pub value: f64,
}

impl TransmissionZero {
    /// Creates a finite zero in normalized low-pass coordinates.
    pub fn finite(normalized_position: f64) -> Self {
        Self::normalized(normalized_position)
    }

    /// Creates a zero that is already normalized.
    pub fn normalized(value: f64) -> Self {
        Self { value }
    }
}

impl From<f64> for TransmissionZero {
    fn from(value: f64) -> Self {
        TransmissionZero::normalized(value)
    }
}

/// One out-of-band attenuation window in physical frequency units.
#[derive(Debug, Clone, PartialEq)]
pub struct OutOfBandAttenuationWindow {
    /// Inclusive lower edge of the attenuation window in Hz.
    pub start_freq_hz: f64,
    /// Inclusive upper edge of the attenuation window in Hz.
    pub stop_freq_hz: f64,
    /// Minimum required attenuation across this window in dB.
    pub attenuation_db: f64,
}

impl OutOfBandAttenuationWindow {
    /// Creates a validated out-of-band attenuation window.
    pub fn new(start_freq_hz: f64, stop_freq_hz: f64, attenuation_db: f64) -> Result<Self> {
        if !start_freq_hz.is_finite() || !stop_freq_hz.is_finite() || stop_freq_hz <= start_freq_hz
        {
            return Err(MfsError::InvalidFrequency(format!(
                "invalid out-of-band attenuation window: start={start_freq_hz}, stop={stop_freq_hz}"
            )));
        }
        if !attenuation_db.is_finite() || attenuation_db <= 0.0 {
            return Err(MfsError::Unsupported(format!(
                "invalid out-of-band attenuation: {attenuation_db}"
            )));
        }

        Ok(Self {
            start_freq_hz,
            stop_freq_hz,
            attenuation_db,
        })
    }
}

/// Collection of out-of-band attenuation requirements.
///
/// File format (three columns):
/// `start_freq_hz  stop_freq_hz  attenuation_db`
/// - whitespace or commas are accepted as separators
/// - blank lines and lines starting with `#` or `//` are ignored
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OutOfBandAttenuationSpec {
    /// One or more attenuation windows supplied by the user or a fixture.
    pub windows: Vec<OutOfBandAttenuationWindow>,
}

impl OutOfBandAttenuationSpec {
    /// Creates an out-of-band attenuation specification from validated windows.
    pub fn new(windows: impl IntoIterator<Item = OutOfBandAttenuationWindow>) -> Self {
        Self {
            windows: windows.into_iter().collect(),
        }
    }

    /// Loads an out-of-band attenuation specification from a three-column file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path).map_err(|error| {
            MfsError::Unsupported(format!(
                "failed to read out-of-band attenuation file {}: {error}",
                path.display()
            ))
        })?;

        let mut windows = Vec::new();
        for (line_index, line) in content.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("//") {
                continue;
            }

            let cleaned = trimmed.replace(',', " ");
            let parts: Vec<&str> = cleaned.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(MfsError::Unsupported(format!(
                    "invalid out-of-band attenuation line {} in {} (expected 3 columns)",
                    line_index + 1,
                    path.display()
                )));
            }

            let start_freq_hz = parts[0].parse::<f64>().map_err(|error| {
                MfsError::Unsupported(format!(
                    "invalid out-of-band attenuation start frequency on line {} in {}: {error}",
                    line_index + 1,
                    path.display()
                ))
            })?;
            let stop_freq_hz = parts[1].parse::<f64>().map_err(|error| {
                MfsError::Unsupported(format!(
                    "invalid out-of-band attenuation stop frequency on line {} in {}: {error}",
                    line_index + 1,
                    path.display()
                ))
            })?;
            let attenuation_db = parts[2].parse::<f64>().map_err(|error| {
                MfsError::Unsupported(format!(
                    "invalid out-of-band attenuation value on line {} in {}: {error}",
                    line_index + 1,
                    path.display()
                ))
            })?;

            windows.push(OutOfBandAttenuationWindow::new(
                start_freq_hz,
                stop_freq_hz,
                attenuation_db,
            )?);
        }

        Ok(Self { windows })
    }
}

/// Minimal normalized internal specification used by the synthesis pipeline.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedFilterSpec {
    /// Number of resonators in the synthesized network.
    pub order: usize,
    /// Minimum passband return loss in dB.
    pub return_loss_db: f64,
    /// Finite transmission zeros expressed in normalized prototype coordinates.
    pub transmission_zeros: Vec<TransmissionZero>,
    /// Optional unloaded Q shared by the design at this stage.
    pub unloaded_q: Option<f64>,
    /// Optional out-of-band attenuation requirements for the current band-pass flow.
    pub out_of_band_attenuation: Option<OutOfBandAttenuationSpec>,
}

/// Backward-compatible public name for the current normalized internal spec.
pub type FilterSpec = NormalizedFilterSpec;

impl NormalizedFilterSpec {
    /// Creates a validated normalized internal specification.
    pub fn new(order: usize, return_loss_db: f64) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }
        if !return_loss_db.is_finite() || return_loss_db <= 0.0 {
            return Err(MfsError::InvalidReturnLoss { return_loss_db });
        }

        Ok(Self {
            order,
            return_loss_db,
            transmission_zeros: Vec::new(),
            unloaded_q: None,
            out_of_band_attenuation: None,
        })
    }

    /// Starts a builder-based filter specification workflow.
    pub fn builder() -> FilterSpecBuilder {
        FilterSpecBuilder::default()
    }

    /// Returns the requested passband return loss in dB.
    pub fn return_loss_db(&self) -> f64 {
        self.return_loss_db
    }

    /// Returns a copy of the spec with the provided transmission zeros.
    pub fn with_transmission_zeros(
        mut self,
        transmission_zeros: impl IntoIterator<Item = TransmissionZero>,
    ) -> Self {
        self.transmission_zeros = transmission_zeros.into_iter().collect();
        self
    }

    /// Returns a copy of the spec with normalized transmission zeros.
    pub fn with_normalized_transmission_zeros(
        self,
        transmission_zeros: impl IntoIterator<Item = f64>,
    ) -> Self {
        self.with_transmission_zeros(
            transmission_zeros
                .into_iter()
                .map(TransmissionZero::normalized),
        )
    }

    /// Returns a copy of the spec with a single unloaded Q value.
    pub fn with_unloaded_q(mut self, unloaded_q: f64) -> Self {
        self.unloaded_q = Some(unloaded_q);
        self
    }

    /// Returns a copy of the spec with an out-of-band attenuation specification.
    pub fn with_out_of_band_attenuation(
        mut self,
        out_of_band_attenuation: OutOfBandAttenuationSpec,
    ) -> Self {
        self.out_of_band_attenuation = Some(out_of_band_attenuation);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spec_stores_normalized_transmission_zeros() -> Result<()> {
        let spec = NormalizedFilterSpec::new(4, 20.0)?.with_transmission_zeros([
            TransmissionZero::normalized(-1.5),
            TransmissionZero::normalized(0.75),
        ]);

        assert_eq!(spec.order, 4);
        assert_eq!(spec.return_loss_db, 20.0);
        assert_eq!(spec.transmission_zeros.len(), 2);
        assert_eq!(spec.transmission_zeros[0].value, -1.5);
        assert_eq!(spec.transmission_zeros[1].value, 0.75);
        Ok(())
    }

    #[test]
    fn builder_constructs_spec_with_transmission_zeros() -> Result<()> {
        let spec = NormalizedFilterSpec::builder()
            .order(4)
            .return_loss_db(20.0)
            .normalized_transmission_zeros(vec![-1.5])
            .build()?;

        assert_eq!(spec.order, 4);
        assert_eq!(spec.transmission_zeros.len(), 1);
        Ok(())
    }

    #[test]
    fn spec_validates_return_loss() {
        let error = NormalizedFilterSpec::new(4, 0.0)
            .expect_err("return loss must be positive");
        assert!(matches!(error, MfsError::InvalidReturnLoss { .. }));
    }

    #[test]
    fn spec_can_store_out_of_band_attenuation_and_unloaded_q() -> Result<()> {
        let out_of_band_attenuation = OutOfBandAttenuationSpec::new([
            OutOfBandAttenuationWindow::new(7.0e9, 7.2e9, 40.0)?,
        ]);
        let spec = NormalizedFilterSpec::new(5, 22.0)?
            .with_unloaded_q(3500.0)
            .with_out_of_band_attenuation(out_of_band_attenuation.clone());

        assert_eq!(spec.unloaded_q, Some(3500.0));
        assert_eq!(
            spec.out_of_band_attenuation,
            Some(out_of_band_attenuation)
        );
        Ok(())
    }
}
