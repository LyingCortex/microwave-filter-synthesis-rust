//! Approximation-layer building blocks for prototype polynomials.
//!
//! This layer is intentionally split into:
//! - `generalized_chebyshev`: spec-facing generalized Chebyshev polynomial assembly
//! - `complex_poly`: reusable complex-polynomial storage and root-solving
//! - `helpers`: Cameron/generalized helper recurrences and domain-specific
//!   polynomial transforms
//! - `generalized_ops`: shared `w <-> s` transforms and recurrence helpers for
//!   the generalized path
//! - `polynomial`: real-valued approximation output containers and projections

mod complex_poly;
mod generalized_chebyshev;
mod generalized_chebyshev_helpers;
mod generalized_ops;
mod polynomial;

pub use complex_poly::{
    AberthRootSolver, AdaptiveRootSolver, CompanionMatrixRootSolver,
    ComplexCoefficient, ComplexPolynomial, ComplexRootSolver, DurandKernerRootSolver,
};
pub use generalized_chebyshev::generalized_chebyshev_polynomials;
pub use polynomial::{PolynomialSet, chebyshev_ripple_factor, monic_polynomial_from_real_roots};

/// Advanced Cameron/generalized helper primitives.
///
/// Most users should start with [`generalized_chebyshev_polynomials`]. This
/// namespace keeps the literature-oriented helper stages explicit without
/// flattening them into the main `approx` surface.
pub mod helpers {
    pub use super::generalized_chebyshev_helpers::{
        APolynomialStage, CameronRecurrence, EPolynomialStage, GeneralizedChebyshevData,
        PaddedTransmissionZeros, build_a_polynomial_stage, build_e_polynomial_stage,
        cameron_recursive, find_a_polynomial, find_e_polynomial, find_eps, find_p_polynomial,
        pad_transmission_zeros, synthesize_generalized_chebyshev_data,
    };
}

// Complex-polynomial primitives live in `complex_poly`, while generalized
// Chebyshev helpers live under `approx::helpers` to keep the common entry
// points smaller.
