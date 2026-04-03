use crate::error::{MfsError, Result};

/// Describes the physical topology of the filter being synthesized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterClass {
    /// Low-pass response with a single upper cutoff.
    LowPass,
    /// High-pass response with a single lower cutoff.
    HighPass,
    /// Band-pass response around a center frequency.
    BandPass,
    /// Band-stop response that rejects a finite band.
    BandStop,
    /// Multi-band pass response with multiple passbands.
    MultiBandPass,
    /// Duplexer-style network with multiple channels.
    Duplexer,
}

/// Backward-compatible alias kept for earlier API naming.

/// Selects which approximation family defines the prototype polynomials.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproximationFamily {
    /// Equal-ripple Chebyshev response.
    Chebyshev,
}

/// Captures the passband return-loss requirement for the design.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReturnLossSpec {
    /// Minimum return loss in dB across the passband.
    pub return_loss_db: f64,
}

impl ReturnLossSpec {
    /// Creates a validated return-loss specification.
    pub fn new(return_loss_db: f64) -> Result<Self> {
        if !return_loss_db.is_finite() || return_loss_db <= 0.0 {
            return Err(MfsError::InvalidReturnLoss { return_loss_db });
        }

        Ok(Self { return_loss_db })
    }
}

/// Indicates whether a transmission zero is already normalized or still in Hz.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransmissionZeroDomain {
    /// Zero is expressed in the normalized low-pass prototype domain.
    Normalized,
    /// Zero is expressed in physical frequency units.
    PhysicalHz,
}

/// Defines one finite transmission zero for the target response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionZero {
    /// Numeric value of the zero in the domain selected below.
    pub value: f64,
    /// Domain that tells the synthesizer how to interpret `value`.
    pub domain: TransmissionZeroDomain,
}

impl TransmissionZero {
    /// Creates a finite zero in normalized low-pass coordinates.
    pub fn finite(normalized_position: f64) -> Self {
        Self::normalized(normalized_position)
    }

    /// Creates a zero that is already normalized.
    pub fn normalized(value: f64) -> Self {
        Self {
            value,
            domain: TransmissionZeroDomain::Normalized,
        }
    }

    /// Creates a zero specified in physical frequency units.
    pub fn physical_hz(value: f64) -> Self {
        Self {
            value,
            domain: TransmissionZeroDomain::PhysicalHz,
        }
    }
}

/// Full parameter bundle required to define a synthesis prototype.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterParameter {
    /// Number of resonators in the synthesized network.
    pub order: usize,
    /// Desired physical filter topology.
    pub filter_class: FilterClass,
    /// Approximation family used to build the prototype.
    pub approximation_family: ApproximationFamily,
    /// Passband return-loss requirement.
    pub return_loss: ReturnLossSpec,
    /// Optional finite transmission zeros for generalized responses.
    pub transmission_zeros: Vec<TransmissionZero>,
}

impl FilterParameter {
    /// Creates a validated parameter bundle with explicit semantic axes.
    pub fn new(
        order: usize,
        filter_class: FilterClass,
        approximation_family: ApproximationFamily,
        return_loss: ReturnLossSpec,
    ) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        Ok(Self {
            order,
            filter_class,
            approximation_family,
            return_loss,
            transmission_zeros: Vec::new(),
        })
    }

    /// Convenience constructor for the common Chebyshev band-pass case.
    pub fn chebyshev(order: usize, return_loss_db: f64) -> Result<Self> {
        Self::new(
            order,
            FilterClass::BandPass,
            ApproximationFamily::Chebyshev,
            ReturnLossSpec::new(return_loss_db)?,
        )
    }

    /// Returns a copy of the spec with a different filter class.
    pub fn with_filter_class(mut self, filter_class: FilterClass) -> Self {
        self.filter_class = filter_class;
        self
    }


    /// Returns the requested passband return loss in dB.
    pub fn return_loss_db(&self) -> f64 {
        self.return_loss.return_loss_db
    }


    /// Returns the requested filter topology.
    pub fn filter_class(&self) -> FilterClass {
        self.filter_class
    }

    /// Returns a copy of the spec with the provided transmission zeros.
    pub fn with_transmission_zeros(mut self, transmission_zeros: Vec<TransmissionZero>) -> Self {
        self.transmission_zeros = transmission_zeros;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chebyshev_spec_sets_independent_semantic_axes() -> Result<()> {
        let spec = FilterParameter::chebyshev(4, 20.0)?;

        assert_eq!(spec.order, 4);
        assert_eq!(spec.filter_class, FilterClass::BandPass);
        assert_eq!(spec.approximation_family, ApproximationFamily::Chebyshev);
        assert_eq!(spec.return_loss.return_loss_db, 20.0);
        Ok(())
    }

    #[test]
    fn generic_constructor_supports_explicit_semantic_axes() -> Result<()> {
        let spec = FilterParameter::new(
            5,
            FilterClass::LowPass,
            ApproximationFamily::Chebyshev,
            ReturnLossSpec::new(19.5)?,
        )?;

        assert_eq!(spec.order, 5);
        assert_eq!(spec.filter_class(), FilterClass::LowPass);
        assert_eq!(spec.approximation_family, ApproximationFamily::Chebyshev);
        assert_eq!(spec.return_loss_db(), 19.5);
        Ok(())
    }

    #[test]
    fn return_loss_spec_validates_return_loss() {
        let error = ReturnLossSpec::new(0.0).expect_err("return loss must be positive");
        assert!(matches!(error, MfsError::InvalidReturnLoss { .. }));
    }
}
