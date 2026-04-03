use crate::error::{MfsError, Result};

use super::CouplingMatrix;

/// Builder for assembling a dense coupling matrix entry by entry.
#[derive(Debug, Default)]
pub struct CouplingMatrixBuilder {
    order: usize,
    data: Vec<f64>,
}

impl CouplingMatrixBuilder {
    /// Allocates a zero-initialized matrix for the requested filter order.
    pub fn new(order: usize) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        let side = order + 2;
        Ok(Self {
            order,
            data: vec![0.0; side * side],
        })
    }

    /// Sets one matrix entry.
    pub fn set(mut self, row: usize, col: usize, value: f64) -> Result<Self> {
        let side = self.order + 2;
        if row >= side || col >= side {
            return Err(MfsError::DimensionMismatch {
                expected: side,
                actual: row.max(col) + 1,
            });
        }

        self.data[row * side + col] = value;
        Ok(self)
    }

    /// Writes the same value into symmetric off-diagonal positions.
    pub fn set_symmetric(self, row: usize, col: usize, value: f64) -> Result<Self> {
        self.set(row, col, value)?.set(col, row, value)
    }

    /// Finalizes the builder and validates the resulting matrix dimensions.
    pub fn build(self) -> Result<CouplingMatrix> {
        CouplingMatrix::new(self.order, self.data)
    }
}
