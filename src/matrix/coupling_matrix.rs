use crate::error::{MfsError, Result};
use nalgebra::DMatrix;
use num_complex::Complex64;

/// Simple shape metadata for a dense coupling matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatrixShape {
    /// Number of matrix rows.
    pub rows: usize,
    /// Number of matrix columns.
    pub cols: usize,
}

/// Dense coupling matrix including source and load rows/columns.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingMatrix {
    order: usize,
    data: Vec<f64>,
}

impl CouplingMatrix {
    /// Creates a coupling matrix from flattened row-major data.
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
            order,
            data,
        })
    }

    /// Creates an identity matrix of the correct source-load augmented size.
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

    /// Returns the resonator count represented by this matrix.
    pub fn order(&self) -> usize {
        self.order
    }

    /// Returns the physical matrix side length including source and load nodes.
    pub fn side(&self) -> usize {
        self.order + 2
    }

    /// Returns the matrix dimensions.
    pub fn shape(&self) -> MatrixShape {
        let side = self.side();
        MatrixShape {
            rows: side,
            cols: side,
        }
    }

    /// Returns one matrix entry if the indices are in range.
    pub fn at(&self, row: usize, col: usize) -> Option<f64> {
        let side = self.side();
        if row >= side || col >= side {
            return None;
        }

        Some(self.data[row * side + col])
    }

    /// Returns the underlying row-major storage.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Returns the matrix as a dense complex matrix for numerical solver backends.
    pub(crate) fn to_complex_dense(&self) -> DMatrix<Complex64> {
        let side = self.side();
        DMatrix::from_row_slice(
            side,
            side,
            &self
                .data
                .iter()
                .copied()
                .map(Complex64::from)
                .collect::<Vec<_>>(),
        )
    }

    /// Returns the source-to-first-resonator coupling magnitude.
    pub fn source_coupling(&self) -> f64 {
        self.at(0, 1).unwrap_or_default().abs()
    }

    /// Returns the last-resonator-to-load coupling magnitude.
    pub fn load_coupling(&self) -> f64 {
        self.at(self.order(), self.side() - 1).unwrap_or_default().abs()
    }

    /// Returns the diagonal detuning term for one resonator.
    pub fn resonator_diagonal(&self, resonator_index: usize) -> Option<f64> {
        if resonator_index >= self.order() {
            return None;
        }

        self.at(resonator_index + 1, resonator_index + 1)
    }

    /// Returns the nearest-neighbor coupling magnitudes along the resonator chain.
    pub fn chain_couplings(&self) -> Vec<f64> {
        (0..self.order().saturating_sub(1))
            .filter_map(|step| self.at(step + 1, step + 2))
            .map(f64::abs)
            .collect()
    }
}
