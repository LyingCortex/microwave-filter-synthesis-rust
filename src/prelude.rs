//! Convenient re-exports for common end-to-end library workflows.

pub use crate::approx::{
    ApproximationEngine, GeneralizedChebyshevApproximation, PolynomialSet, PrototypePoint,
};
pub use crate::error::{MfsError, Result};
pub use crate::freq::{
    BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping, NormalizedSample,
};
pub use crate::matrix::{
    BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder, MatrixShape, MatrixTopology,
};
pub use crate::output::{
    Table, format_aligned_summary, format_aligned_summary_with_width, format_box_table, format_complex_scalar_parts, format_decimal_scalar, format_key_value_table_data,
    format_markdown_table, format_matrix_table_data, format_polynomial_table_data,
    format_out_of_band_attenuation_table_data,
    print_terminal_synthesis_report,
    format_reference_actual_polynomial_table_data, format_response_samples_table_data,
    format_root_comparison_table_data, format_section_title, format_singularity_table_data,
    render_markdown_synthesis_report, render_terminal_filter_database_report,
    render_terminal_synthesis_report,
};
pub use crate::response::{ResponseSample, ResponseSolver, SParameterResponse};
pub use crate::spec::{
    FilterSpec, FilterSpecBuilder, OutOfBandAttenuationSpec, OutOfBandAttenuationWindow,
    TransmissionZero,
};
pub use crate::synthesis::{
    synthesize_and_evaluate_with_mapping, synthesize_generalized_chebyshev,
    CanonicalMatrixSynthesis,
    EvaluationOutcome, SectionSynthesis, SynthesisOutcome, VerifiedSectionSynthesis,
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
pub use crate::{
    bandpass, filter_spec, generalized_chebyshev, generalized_chebyshev_polynomials,
    generalized_chebyshev_with_response, lowpass, normalize_transmission_zeros_hz,
    synthesize_and_evaluate_generalized_chebyshev,
};
