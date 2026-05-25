use std::path::Path;

use crate::error::Result;

use super::{FilterSpec, OutOfBandAttenuationSpec, TransmissionZero};

/// Builder-style constructor for filter specifications.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FilterSpecBuilder {
    order: Option<usize>,
    return_loss_db: Option<f64>,
    transmission_zeros: Vec<TransmissionZero>,
    unloaded_q: Option<f64>,
    out_of_band_attenuation: Option<OutOfBandAttenuationSpec>,
}

impl FilterSpecBuilder {
    /// Sets the filter order.
    pub fn order(mut self, order: usize) -> Self {
        self.order = Some(order);
        self
    }

    /// Sets the passband return loss in dB.
    pub fn return_loss_db(mut self, return_loss_db: f64) -> Self {
        self.return_loss_db = Some(return_loss_db);
        self
    }

    /// Sets normalized transmission zeros.
    pub fn transmission_zeros(
        mut self,
        transmission_zeros: impl IntoIterator<Item = TransmissionZero>,
    ) -> Self {
        self.transmission_zeros = transmission_zeros.into_iter().collect();
        self
    }

    /// Sets normalized transmission zeros.
    pub fn normalized_transmission_zeros(
        self,
        transmission_zeros: impl IntoIterator<Item = f64>,
    ) -> Self {
        self.transmission_zeros(
            transmission_zeros
                .into_iter()
                .map(TransmissionZero::normalized),
        )
    }

    /// Sets one unloaded Q value for the synthesized design.
    pub fn unloaded_q(mut self, unloaded_q: f64) -> Self {
        self.unloaded_q = Some(unloaded_q);
        self
    }

    /// Sets an out-of-band attenuation specification.
    pub fn out_of_band_attenuation(
        mut self,
        out_of_band_attenuation: OutOfBandAttenuationSpec,
    ) -> Self {
        self.out_of_band_attenuation = Some(out_of_band_attenuation);
        self
    }

    /// Loads an out-of-band attenuation specification from a three-column file.
    pub fn out_of_band_attenuation_file(
        mut self,
        path: impl AsRef<Path>,
    ) -> Result<Self> {
        self.out_of_band_attenuation = Some(OutOfBandAttenuationSpec::from_file(path)?);
        Ok(self)
    }

    /// Builds a validated filter specification.
    pub fn build(self) -> Result<FilterSpec> {
        let order = self.order.unwrap_or(0);
        let return_loss_db = self.return_loss_db.unwrap_or(0.0);
        let mut spec =
            FilterSpec::new(order, return_loss_db)?.with_transmission_zeros(self.transmission_zeros);
        if let Some(unloaded_q) = self.unloaded_q {
            spec = spec.with_unloaded_q(unloaded_q);
        }
        if let Some(out_of_band_attenuation) = self.out_of_band_attenuation {
            spec = spec.with_out_of_band_attenuation(out_of_band_attenuation);
        }
        Ok(spec)
    }
}
