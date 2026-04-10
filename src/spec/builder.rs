use crate::error::Result;

use super::{ApproximationFamily, FilterClass, FilterSpec, ReturnLossSpec, TransmissionZero};

/// Builder-style constructor for filter specifications.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterSpecBuilder {
    order: Option<usize>,
    filter_class: FilterClass,
    approximation_family: ApproximationFamily,
    return_loss_db: Option<f64>,
    transmission_zeros: Vec<TransmissionZero>,
}

impl Default for FilterSpecBuilder {
    fn default() -> Self {
        Self {
            order: None,
            filter_class: FilterClass::BandPass,
            approximation_family: ApproximationFamily::Chebyshev,
            return_loss_db: None,
            transmission_zeros: Vec::new(),
        }
    }
}

impl FilterSpecBuilder {
    /// Sets the filter order.
    pub fn order(mut self, order: usize) -> Self {
        self.order = Some(order);
        self
    }

    /// Sets the requested physical filter class.
    pub fn filter_class(mut self, filter_class: FilterClass) -> Self {
        self.filter_class = filter_class;
        self
    }

    /// Sets the approximation family.
    pub fn approximation_family(mut self, approximation_family: ApproximationFamily) -> Self {
        self.approximation_family = approximation_family;
        self
    }

    /// Sets the passband return loss in dB.
    pub fn return_loss_db(mut self, return_loss_db: f64) -> Self {
        self.return_loss_db = Some(return_loss_db);
        self
    }

    /// Sets the requested transmission zeros.
    pub fn transmission_zeros(mut self, transmission_zeros: Vec<TransmissionZero>) -> Self {
        self.transmission_zeros = transmission_zeros;
        self
    }

    /// Builds a validated filter specification.
    pub fn build(self) -> Result<FilterSpec> {
        let order = self.order.unwrap_or(0);
        let return_loss = ReturnLossSpec::new(self.return_loss_db.unwrap_or(0.0))?;
        let mut spec = FilterSpec::new(
            order,
            self.filter_class,
            self.approximation_family,
            return_loss,
        )?;
        spec.transmission_zeros = self.transmission_zeros;
        Ok(spec)
    }
}
