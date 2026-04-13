mod canonical;
mod engine;
mod orchestration;
mod placeholder;
mod residues;
mod sections;

pub use canonical::{
    AdmittancePolynomials, CanonicalMatrixSynthesis, ResidueExpansion, ResiduePole,
};
pub use engine::{MatrixSynthesisEngine, MatrixSynthesisMethod, MatrixSynthesisOutcome};
pub use orchestration::{
    synthesize_and_evaluate_with_mapping, synthesize_generalized_chebyshev, EvaluationOutcome,
    SynthesisOutcome,
};
pub use sections::SectionSynthesis;
pub use sections::VerifiedSectionSynthesis;

pub(crate) use placeholder::synthesize_placeholder_matrix;
pub(crate) use residues::{
    build_transversal_from_residues, synthesize_admittance_polynomials,
    synthesize_residue_expansions,
};
