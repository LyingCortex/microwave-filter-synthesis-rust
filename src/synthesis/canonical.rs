use crate::approx::PolynomialSet;
use crate::error::Result;
use crate::matrix::{CouplingMatrix, MatrixTopology};

pub use super::residues::{AdmittancePolynomials, ResidueExpansion, ResiduePole};
pub use super::engine::MatrixSynthesisOutcome;
use super::MatrixSynthesisEngine;

/// Top-level facade for coupling-matrix synthesis operations.
#[derive(Debug, Default, Clone, Copy)]
pub struct CanonicalMatrixSynthesis {
    engine: MatrixSynthesisEngine,
}

impl CanonicalMatrixSynthesis {
    /// Synthesizes the canonical matrix implied by the approximation output.
    pub fn synthesize(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        self.engine.synthesize(polynomials)
    }

    /// Synthesizes the canonical matrix and reports which internal path was used.
    pub fn synthesize_with_details(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<MatrixSynthesisOutcome> {
        self.engine.synthesize_with_details(polynomials)
    }

    /// Synthesizes a matrix and converts it into the requested topology.
    pub fn synthesize_with_topology(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
    ) -> Result<CouplingMatrix> {
        self.engine.synthesize_with_topology(polynomials, topology)
    }
}
