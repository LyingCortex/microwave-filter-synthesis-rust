mod chebyshev;
mod generalized_chebyshev;
mod polynomial;

use crate::error::Result;
use crate::freq::FrequencyPlan;
use crate::spec::FilterSpec;

pub use chebyshev::ChebyshevApproximation;
pub use generalized_chebyshev::{
    CameronRecurrence, ComplexCoefficient, ComplexPolynomial, GeneralizedChebyshevData,
    PaddedTransmissionZeros, cameron_recursive, find_a_polynomial, find_e_polynomial, find_eps,
    find_p_polynomial, pad_transmission_zeros, synthesize_generalized_chebyshev_data,
};
pub use polynomial::{PolynomialSet, chebyshev_ripple_factor, monic_polynomial_from_real_roots};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrototypePoint {
    pub x: f64,
    pub y: f64,
}

pub trait ApproximationEngine {
    fn synthesize(&self, spec: &FilterSpec, plan: &impl FrequencyPlan) -> Result<PolynomialSet>;
}
