use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, PartialEq)]
pub enum MfsError {
    InvalidOrder { order: usize },
    InvalidReturnLoss { return_loss_db: f64 },
    InvalidFrequency(String),
    InvalidGridSize { points: usize },
    InvalidTransmissionZero(String),
    DimensionMismatch { expected: usize, actual: usize },
    Unsupported(String),
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
            Self::Unsupported(message) => write!(f, "unsupported operation: {message}"),
        }
    }
}

impl Error for MfsError {}

pub type Result<T> = std::result::Result<T, MfsError>;
