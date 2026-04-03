use crate::error::{MfsError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatrixShape {
    pub rows: usize,
    pub cols: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CouplingMatrix {
    shape: MatrixShape,
    data: Vec<f64>,
}

impl CouplingMatrix {
    pub fn new(order: usize, data: Vec<f64>) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        let side = order + 2;
        let expected = side * side;
        if data.len() != expected {
            return Err(MfsError::DimensionMismatch {
                expected,
                actual: data.len(),
            });
        }

        Ok(Self {
            shape: MatrixShape {
                rows: side,
                cols: side,
            },
            data,
        })
    }

    pub fn identity(order: usize) -> Result<Self> {
        if order == 0 {
            return Err(MfsError::InvalidOrder { order });
        }

        let side = order + 2;
        let mut data = vec![0.0; side * side];
        for index in 0..side {
            data[index * side + index] = 1.0;
        }

        Self::new(order, data)
    }

    pub fn order(&self) -> usize {
        self.shape.rows - 2
    }

    pub fn shape(&self) -> MatrixShape {
        self.shape
    }

    pub fn at(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.shape.rows || col >= self.shape.cols {
            return None;
        }

        Some(self.data[row * self.shape.cols + col])
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub fn source_coupling(&self) -> f64 {
        self.at(0, 1).unwrap_or_default().abs()
    }

    pub fn load_coupling(&self) -> f64 {
        self.at(self.order(), self.order() + 1)
            .unwrap_or_default()
            .abs()
    }

    pub fn resonator_diagonal(&self, resonator_index: usize) -> Option<f64> {
        if resonator_index >= self.order() {
            return None;
        }

        self.at(resonator_index + 1, resonator_index + 1)
    }

    pub fn chain_couplings(&self) -> Vec<f64> {
        (0..self.order().saturating_sub(1))
            .filter_map(|step| self.at(step + 1, step + 2))
            .map(f64::abs)
            .collect()
    }
}
