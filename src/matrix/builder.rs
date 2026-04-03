use crate::error::{MfsError, Result};

use super::CouplingMatrix;

#[derive(Debug, Default)]
pub struct CouplingMatrixBuilder {
    order: usize,
    data: Vec<f64>,
}

impl CouplingMatrixBuilder {
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

    pub fn set_symmetric(self, row: usize, col: usize, value: f64) -> Result<Self> {
        self.set(row, col, value)?.set(col, row, value)
    }

    pub fn build(self) -> Result<CouplingMatrix> {
        CouplingMatrix::new(self.order, self.data)
    }
}
