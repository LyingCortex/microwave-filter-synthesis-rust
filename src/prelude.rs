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

pub use crate::design::FilterDesign;
pub use crate::error::{MfsError, Result};
pub use crate::matrix::{CouplingMatrix, MatrixTopology};
pub use crate::response::{ResponseSample, SParameterResponse};
