//! Convenient re-exports for common library workflows.
//!
//! ```rust
//! use mfs::prelude::*;
//!
//! let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
//!     .zeros_hz([6.5e9, 7.0e9])
//!     .synthesize()?;
//!
//! let folded = design.to_folded()?;
//! let response = design.response(6.5e9, 7.0e9, 201)?;
//! # Ok::<(), MfsError>(())
//! ```

// Primary API
pub use crate::design::FilterDesign;
pub use crate::error::{MfsError, Result};
pub use crate::matrix::{CouplingMatrix, MatrixTopology};
pub use crate::response::{ResponseSample, SParameterResponse};

// Supporting types (commonly needed)
pub use crate::freq::{BandPassMapping, FrequencyGrid, FrequencyMapping, LowPassMapping};
pub use crate::spec::{FilterSpec, FilterSpecBuilder, TransmissionZero};
pub use crate::response::ResponseSolver;
pub use crate::verify::ResponseTolerance;
pub use crate::matrix::MatrixTopology as TopologyKind;

// Legacy functions (for existing code)
pub use crate::{
    bandpass, filter_spec, generalized_chebyshev, generalized_chebyshev_polynomials,
    generalized_chebyshev_with_response, lowpass, normalize_transmission_zeros_hz,
};
pub use crate::synthesis::{
    ApproximationKind, EvaluationOutcome, SectionSynthesis, SynthesisOutcome,
};
pub use crate::transform::{
    transform_matrix, transform_matrix_with_response_check, SectionTransformOutcome,
};
pub use crate::output::{
    print_terminal_synthesis_report, render_markdown_synthesis_report,
    render_terminal_filter_database_report, render_terminal_synthesis_report,
};
