use crate::error::{MfsError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterClass {
    LowPass,
    HighPass,
    BandPass,
    BandStop,
    MultiBandPass,
    Duplexer,
}

pub type FilterType = FilterClass;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproximationFamily {
    Chebyshev,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PerformanceSpec {
    pub return_loss_db: f64,
}

impl PerformanceSpec {
    pub fn new(return_loss_db: f64) -> Result<Self> {
        if !return_loss_db.is_finite() || return_loss_db <= 0.0 {
            return Err(MfsError::InvalidReturnLoss { return_loss_db });
        }

        Ok(Self { return_loss_db })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransmissionZeroDomain {
    Normalized,
    PhysicalHz,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionZero {
    pub value: f64,
    pub domain: TransmissionZeroDomain,
}

impl TransmissionZero {
    pub fn finite(normalized_position: f64) -> Self {
        Self::normalized(normalized_position)
    }

    pub fn normalized(value: f64) -> Self {
        Self {
            value,
            domain: TransmissionZeroDomain::Normalized,
        }
    }

    pub fn physical_hz(value: f64) -> Self {
        Self {
            value,
            domain: TransmissionZeroDomain::PhysicalHz,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FilterSpec {
    pub order: usize,
    pub filter_class: FilterClass,
    pub approximation_family: ApproximationFamily,
    pub performance: PerformanceSpec,
    pub transmission_zeros: Vec<TransmissionZero>,
}

impl FilterSpec {
    pub fn new(
        order: usize,
        filter_class: FilterClass,
        approximation_family: ApproximationFamily,
        performance: PerformanceSpec,
    ) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        Ok(Self {
            order,
            filter_class,
            approximation_family,
            performance,
            transmission_zeros: Vec::new(),
        })
    }

    pub fn chebyshev(order: usize, return_loss_db: f64) -> Result<Self> {
        Self::new(
            order,
            FilterClass::BandPass,
            ApproximationFamily::Chebyshev,
            PerformanceSpec::new(return_loss_db)?,
        )
    }

    pub fn with_filter_class(mut self, filter_class: FilterClass) -> Self {
        self.filter_class = filter_class;
        self
    }

    pub fn with_filter_type(self, filter_type: FilterType) -> Self {
        self.with_filter_class(filter_type)
    }

    pub fn return_loss_db(&self) -> f64 {
        self.performance.return_loss_db
    }

    pub fn filter_type(&self) -> FilterType {
        self.filter_class
    }

    pub fn filter_class(&self) -> FilterClass {
        self.filter_class
    }

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
        let spec = FilterSpec::chebyshev(4, 20.0)?;

        assert_eq!(spec.order, 4);
        assert_eq!(spec.filter_class, FilterClass::BandPass);
        assert_eq!(spec.approximation_family, ApproximationFamily::Chebyshev);
        assert_eq!(spec.performance.return_loss_db, 20.0);
        Ok(())
    }

    #[test]
    fn generic_constructor_supports_explicit_semantic_axes() -> Result<()> {
        let spec = FilterSpec::new(
            5,
            FilterClass::LowPass,
            ApproximationFamily::Chebyshev,
            PerformanceSpec::new(19.5)?,
        )?;

        assert_eq!(spec.order, 5);
        assert_eq!(spec.filter_class(), FilterClass::LowPass);
        assert_eq!(spec.approximation_family, ApproximationFamily::Chebyshev);
        assert_eq!(spec.return_loss_db(), 19.5);
        Ok(())
    }

    #[test]
    fn filter_type_compatibility_helper_still_updates_filter_class() -> Result<()> {
        let spec = FilterSpec::chebyshev(3, 18.0)?.with_filter_type(FilterType::LowPass);

        assert_eq!(spec.filter_class(), FilterClass::LowPass);
        assert_eq!(spec.filter_type(), FilterType::LowPass);
        Ok(())
    }

    #[test]
    fn performance_spec_validates_return_loss() {
        let error = PerformanceSpec::new(0.0).expect_err("return loss must be positive");
        assert!(matches!(error, MfsError::InvalidReturnLoss { .. }));
    }
}
