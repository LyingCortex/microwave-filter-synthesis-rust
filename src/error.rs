use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Error type shared across synthesis, mapping, and response evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum MfsError {
    /// The requested filter order is zero or otherwise invalid.
    InvalidOrder { order: usize },
    /// The requested return loss is non-positive or non-finite.
    InvalidReturnLoss { return_loss_db: f64 },
    /// A physical or normalized frequency input failed validation.
    InvalidFrequency(String),
    /// The requested frequency grid is too small to be useful.
    InvalidGridSize { points: usize },
    /// A transmission zero value or its placement is invalid.
    InvalidTransmissionZero(String),
    /// Two related vectors or matrices do not share the expected size.
    DimensionMismatch { expected: usize, actual: usize },
    /// A numerical algorithm failed to converge or produced degenerate results.
    NumericalFailure(String),
    /// The requested feature is recognized but not yet implemented.
    NotImplemented(String),
    /// A required precondition was not met (wrong input state, missing data, invalid parameters).
    PreconditionViolation(String),
}

impl Display for MfsError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidOrder { order } => write!(f, "filter order must be positive, got {order}"),
            Self::InvalidReturnLoss { return_loss_db } => {
                write!(f, "return loss must be positive, got {return_loss_db}")
            }
            Self::InvalidFrequency(message) => write!(f, "invalid frequency: {message}"),
            Self::InvalidGridSize { points } => {
                write!(f, "frequency grid requires at least 2 points, got {points}")
            }
            Self::InvalidTransmissionZero(message) => {
                write!(f, "invalid transmission zero: {message}")
            }
            Self::DimensionMismatch { expected, actual } => {
                write!(f, "dimension mismatch: expected {expected}, got {actual}")
            }
            Self::NumericalFailure(message) => write!(f, "numerical failure: {message}"),
            Self::NotImplemented(message) => write!(f, "not implemented: {message}"),
            Self::PreconditionViolation(message) => {
                write!(f, "precondition violation: {message}")
            }
        }
    }
}

impl Error for MfsError {}

/// Crate-wide result alias used by public APIs and internal helpers.
pub type Result<T> = std::result::Result<T, MfsError>;
