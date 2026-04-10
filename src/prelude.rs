//! Convenient re-exports for common end-to-end library workflows.

pub use crate::approx::{
    ApproximationEngine, ChebyshevApproximation, GeneralizedChebyshevApproximation, PolynomialSet,
    PrototypePoint,
};
pub use crate::error::{MfsError, Result};
pub use crate::freq::{
    BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping, NormalizedSample,
};
pub use crate::matrix::{
    BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder, MatrixShape, MatrixTopology,
};
pub use crate::response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use crate::spec::{
    ApproximationFamily, FilterClass, FilterSpec, FilterSpecBuilder, ReturnLossSpec, TransmissionZero,
    TransmissionZeroDomain,
};
pub use crate::synthesis::{
    ApproximationStageKind, CanonicalMatrixSynthesis, ChebyshevSynthesis, EvaluationOutcome,
    SectionSynthesis, SynthesisOutcome, VerifiedSectionSynthesis,
};
pub use crate::transform::{
    extract_quadruplet_section, extract_quadruplet_section_with_report,
    extract_quadruplet_section_with_response_check, extract_triplet_section,
    extract_triplet_section_with_report, extract_triplet_section_with_response_check,
    extract_trisection_section, extract_trisection_section_with_report,
    extract_trisection_section_with_response_check, to_arrow, to_folded, to_wheel,
    transform_matrix, transform_matrix_with_response_check, SectionTransformEngine,
    SectionTransformOutcome, TopologyKind, TransformEngine, TransformReport,
};
pub use crate::verify::{
    compare_responses, matches_arrow_pattern, matches_folded_pattern, matches_topology_pattern,
    verify_quadruplet_extraction, verify_triplet_extraction, verify_trisection_extraction,
    MatrixPatternTolerance, ResponseCheckReport, ResponseComparison, ResponseTolerance,
    SectionVerificationReport,
};
