mod engine;
mod machiarella;
mod orchestration;
mod recovery;
mod residues;
mod sections;

pub use engine::{MatrixSynthesisEngine, MatrixSynthesisMethod, MatrixSynthesisOutcome};
pub use machiarella::{
    junction_response, ChannelState, FullIterationState, JunctionCommon, JunctionResponse,
    Prepared, prepare, run_full_iteration, synthesize_full,
};
pub use orchestration::{
    build_filter_spec, synthesize_and_evaluate_with_mapping, synthesize_generalized_chebyshev,
    ApproximationKind, EvaluationOutcome, SynthesisOutcome,
};
pub use residues::synthesize_residue_expansions;
pub use residues::{
    classify_residues, residue_real_component, ClassifiedResidues, ResidueClassification,
};
pub use residues::{AdmittancePolynomials, ResidueExpansion, ResiduePole};
pub use sections::SectionSynthesis;

pub use residues::synthesize_admittance_polynomials;

pub(crate) use recovery::recover_channel_design;
pub(crate) use residues::build_transversal_from_residues;
