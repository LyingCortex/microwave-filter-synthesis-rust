//! Residue-based coupling-matrix synthesis.
//!
//! This module is split into three seams:
//! - [`admittance`]: input-parsing — Y-parameter polynomials from the approximation output.
//! - [`expansion`]: residue-computation — partial-fraction expansions and classification.
//! - [`transversal`]: matrix-assembly — building the transversal coupling matrix.

mod admittance;
mod expansion;
mod transversal;

pub use admittance::{AdmittancePolynomials, synthesize_admittance_polynomials};
pub use expansion::{
    classify_residues, residue_real_component, synthesize_residue_expansions, ClassifiedResidues,
    ResidueClassification, ResidueExpansion, ResiduePole,
};
pub(crate) use transversal::build_transversal_from_residues;

#[allow(unused_imports)]
pub(crate) use transversal::combine_conjugate_pair;
