mod engine;
mod orchestration;
mod residues;
mod sections;

pub use residues::{AdmittancePolynomials, ResidueExpansion, ResiduePole};
pub use residues::{ClassifiedResidues, ResidueClassification, classify_residues, residue_real_component};
pub use residues::synthesize_residue_expansions;
pub use engine::{MatrixSynthesisEngine, MatrixSynthesisMethod, MatrixSynthesisOutcome};
pub use orchestration::{
    build_filter_spec, synthesize_and_evaluate_with_mapping, synthesize_generalized_chebyshev,
    ApproximationKind, EvaluationOutcome, SynthesisOutcome,
};
pub use sections::SectionSynthesis;

pub use residues::synthesize_admittance_polynomials;

pub(crate) use residues::build_transversal_from_residues;
