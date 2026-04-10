//! Approximation-layer building blocks for prototype polynomials.
//!
//! This layer is intentionally split into:
//! - `chebyshev`: high-level approximation engines and main-flow orchestration
//! - `complex_poly`: reusable complex-polynomial storage and root-solving
//! - `generalized_chebyshev`: Cameron/generalized helper recurrences and
//!   domain-specific polynomial transforms
//! - `generalized_ops`: shared `w <-> s` transforms and recurrence helpers for
//!   the generalized path
//! - `polynomial`: real-valued approximation output containers and projections

mod chebyshev;
mod complex_poly;
mod generalized_chebyshev;
mod generalized_ops;
mod polynomial;

use crate::error::Result;
use crate::freq::FrequencyMapping;
use crate::spec::FilterSpec;

pub use complex_poly::{
    ComplexCoefficient, ComplexPolynomial, ComplexRootSolver, DurandKernerRootSolver,
};
pub use chebyshev::{ChebyshevApproximation, GeneralizedChebyshevApproximation};
pub use generalized_chebyshev::{
    APolynomialStage, CameronRecurrence, EPolynomialStage, GeneralizedChebyshevData,
    PaddedTransmissionZeros, build_a_polynomial_stage, build_e_polynomial_stage,
    cameron_recursive, find_a_polynomial, find_e_polynomial, find_eps, find_p_polynomial,
    pad_transmission_zeros, synthesize_generalized_chebyshev_data,
};

// Complex-polynomial primitives live in `complex_poly`, while generalized
// Chebyshev helpers reuse them to build domain-specific E/F/P/A polynomials.
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
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
    ) -> Result<PolynomialSet>;
}
