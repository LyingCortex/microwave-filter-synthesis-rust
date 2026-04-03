mod chebyshev;
mod generalized_chebyshev;
mod polynomial;

use crate::error::Result;
use crate::freq::FrequencyMapping;
use crate::spec::FilterParameter;

pub use chebyshev::ChebyshevApproximation;
pub use generalized_chebyshev::{
    CameronRecurrence, ComplexCoefficient, ComplexPolynomial, GeneralizedChebyshevData,
    PaddedTransmissionZeros, cameron_recursive, find_a_polynomial, find_e_polynomial, find_eps,
    find_p_polynomial, pad_transmission_zeros, synthesize_generalized_chebyshev_data,
};
pub use polynomial::{PolynomialSet, chebyshev_ripple_factor, monic_polynomial_from_real_roots};

/// Generic point on a prototype response curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrototypePoint {
    /// Horizontal coordinate, typically normalized frequency.
    pub x: f64,
    /// Vertical coordinate, typically amplitude or attenuation.
    pub y: f64,
}

/// Trait implemented by approximation engines that generate prototype polynomials.
pub trait ApproximationEngine {
    /// Synthesizes a polynomial set from a validated filter specification.
    fn synthesize(
        &self,
        spec: &FilterParameter,
        mapping: &impl FrequencyMapping,
    ) -> Result<PolynomialSet>;
}
