use crate::error::{MfsError, Result};
use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Supported coupling-matrix topologies exposed by the library.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum MatrixTopology {
    /// Standard transversal or otherwise untransformed matrix form.
    #[default]
    Transversal,
    /// Folded topology obtained by similarity rotations.
    Folded,
    /// Arrow topology obtained by similarity rotations.
    Arrow,
    /// Wheel topology: currently implemented as Arrow with a relabeled topology tag.
    /// A true Wheel reduction (with distinct sparsity pattern) is not yet available.
    #[deprecated(note = "Use Arrow topology instead. Real Wheel implementation is not yet available.")]
    Wheel,
}

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
    pub(crate) order: usize,
    pub(crate) topology: MatrixTopology,
    pub(crate) data: Vec<f64>,
}

/// Physical-frequency view of a normalized coupling matrix after band-pass scaling.
#[derive(Debug, Clone, PartialEq)]
pub struct BandPassScaledCouplingMatrix {
    pub(crate) matrix_hz: CouplingMatrix,
    pub(crate) source_external_q: f64,
    pub(crate) load_external_q: f64,
}

impl BandPassScaledCouplingMatrix {
    /// Returns the scaled dense matrix with resonator couplings in Hz.
    pub fn matrix_hz(&self) -> &CouplingMatrix {
        &self.matrix_hz
    }

    /// Returns the source external quality factor implied by the normalized matrix.
    pub fn source_external_q(&self) -> f64 {
        self.source_external_q
    }

    /// Returns the load external quality factor implied by the normalized matrix.
    pub fn load_external_q(&self) -> f64 {
        self.load_external_q
    }
}

impl CouplingMatrix {
    /// Creates a coupling matrix from flattened row-major data.
    pub fn new(order: usize, data: Vec<f64>) -> Result<Self> {
        Self::new_with_topology(order, MatrixTopology::Transversal, data)
    }

    /// Creates a coupling matrix with an explicit topology label.
    pub fn new_with_topology(
        order: usize,
        topology: MatrixTopology,
        data: Vec<f64>,
    ) -> Result<Self> {
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
            topology,
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

    /// Returns the topology label currently attached to the matrix.
    pub fn topology(&self) -> MatrixTopology {
        self.topology
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

    /// Debug-panicking accessor for internal use. Panics on OOB in debug, unchecked in release.
    #[inline]
    pub(crate) fn get(&self, row: usize, col: usize) -> f64 {
        debug_assert!(
            row < self.side() && col < self.side(),
            "matrix access out of bounds: ({row}, {col}) for side {}",
            self.side()
        );
        unsafe { *self.data.get_unchecked(row * self.side() + col) }
    }

    /// Returns the underlying row-major storage.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Returns mutable row-major storage for internal synthesis helpers.
    pub(crate) fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Returns a new matrix transformed into the requested topology when supported.
    pub(crate) fn transform_topology(&self, topology: MatrixTopology) -> Result<Self> {
        match topology {
            MatrixTopology::Transversal => Ok(self.clone()),
            MatrixTopology::Folded => Ok(self.to_folded()),
            MatrixTopology::Arrow => Ok(self.to_arrow()),
            #[allow(deprecated)]
            MatrixTopology::Wheel => {
                let mut matrix = self.to_arrow();
                matrix.topology = MatrixTopology::Wheel;
                Ok(matrix)
            }
        }
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
        self.get(0, 1).abs()
    }

    /// Returns the last-resonator-to-load coupling magnitude.
    pub fn load_coupling(&self) -> f64 {
        self.get(self.order(), self.side() - 1).abs()
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

    pub(crate) fn set_entry(&mut self, row: usize, col: usize, value: f64) {
        let side = self.side();
        self.data[row * side + col] = value;
    }

    pub(crate) fn clean_small_values(&mut self) {
        for value in &mut self.data {
            if value.abs() <= 1e-10 {
                *value = 0.0;
            }
        }
    }

    #[allow(dead_code)]
    pub(crate) fn multiply(&self, rhs: &Self) -> Self {
        let side = self.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for k in 0..side {
                let lhs_val = self.data[row * side + k];
                if lhs_val == 0.0 {
                    continue;
                }
                for col in 0..side {
                    data[row * side + col] += lhs_val * rhs.data[k * side + col];
                }
            }
        }
        Self {
            order: self.order,
            topology: self.topology,
            data,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn transpose(&self) -> Self {
        let side = self.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for col in 0..side {
                data[col * side + row] = self.data[row * side + col];
            }
        }
        Self {
            order: self.order,
            topology: self.topology,
            data,
        }
    }

    /// Applies a Givens similarity rotation in-place: M' = R * M * R^T
    /// where R differs from identity only at entries (p,p), (p,q), (q,p), (q,q).
    /// This is O(n) per rotation instead of O(n³) for full matrix multiply.
    pub(crate) fn apply_givens_similarity_inplace(&mut self, pivot_a: usize, pivot_b: usize, cosine: f64, sine: f64) {
        let side = self.side();

        // Step 1: M' = R * M (affects rows pivot_a and pivot_b)
        // row_a' = cos * row_a - sin * row_b
        // row_b' = sin * row_a + cos * row_b
        for col in 0..side {
            let a = self.data[pivot_a * side + col];
            let b = self.data[pivot_b * side + col];
            self.data[pivot_a * side + col] = cosine * a - sine * b;
            self.data[pivot_b * side + col] = sine * a + cosine * b;
        }

        // Step 2: M'' = M' * R^T (affects columns pivot_a and pivot_b)
        // col_a' = cos * col_a - sin * col_b  (of M')
        // col_b' = sin * col_a + cos * col_b  (of M')
        // Note: R^T has cos on diagonal, but transposed off-diag: R^T[a,b] = sin, R^T[b,a] = -sin
        // So col_a'' = cos * col_a' + sin * col_b'
        //    col_b'' = -sin * col_a' + cos * col_b'
        for row in 0..side {
            let a = self.data[row * side + pivot_a];
            let b = self.data[row * side + pivot_b];
            self.data[row * side + pivot_a] = cosine * a - sine * b;
            self.data[row * side + pivot_b] = sine * a + cosine * b;
        }
    }

    fn flip_sign(&self, diagonal_index: usize) -> Self {
        let side = self.side();
        let mut result = self.clone();
        // Negate row `diagonal_index`
        for col in 0..side {
            result.data[diagonal_index * side + col] = -result.data[diagonal_index * side + col];
        }
        // Negate column `diagonal_index`
        for row in 0..side {
            result.data[row * side + diagonal_index] = -result.data[row * side + diagonal_index];
        }
        // The diagonal element was negated twice, restore it
        result.data[diagonal_index * side + diagonal_index] =
            -result.data[diagonal_index * side + diagonal_index];
        result
    }

    fn to_folded(&self) -> Self {
        use super::rotations::safe_angle;

        let mut matrix = self.clone();
        let side = matrix.side();
        let filter_order = matrix.order();
        let end_oper = if filter_order % 2 == 0 {
            filter_order / 2
        } else {
            (filter_order - 1) / 2
        };

        for row_oper_num in 0..end_oper {
            for col_in_row_oper in ((row_oper_num + 2)..=filter_order).rev() {
                if matrix
                    .get(row_oper_num, col_in_row_oper)
                    .abs()
                    > 1e-7
                {
                    // Column operation: zero out (row_oper_num, col_in_row_oper)
                    let pivot_a = col_in_row_oper - 1;
                    let pivot_b = col_in_row_oper;
                    let numerator = -matrix.get(row_oper_num, col_in_row_oper);
                    let denominator = matrix.get(row_oper_num, col_in_row_oper - 1);
                    let theta = safe_angle(numerator, denominator);
                    matrix.apply_givens_similarity_inplace(pivot_a, pivot_b, theta.cos(), theta.sin());
                } else {
                    matrix.set_entry(row_oper_num, col_in_row_oper, 0.0);
                    matrix.set_entry(col_in_row_oper, row_oper_num, 0.0);
                }
            }

            let col_oper_num = side - 1 - row_oper_num;
            for row_in_col_oper in (row_oper_num + 2)..=(col_oper_num.saturating_sub(2)) {
                if matrix
                    .get(row_in_col_oper, col_oper_num)
                    .abs()
                    > 1e-7
                {
                    // Row operation: zero out (row_in_col_oper, col_oper_num)
                    let pivot_a = row_in_col_oper;
                    let pivot_b = row_in_col_oper + 1;
                    let numerator = matrix.get(row_in_col_oper, col_oper_num);
                    let denominator = matrix.get(row_in_col_oper + 1, col_oper_num);
                    let theta = safe_angle(numerator, denominator);
                    matrix.apply_givens_similarity_inplace(pivot_a, pivot_b, theta.cos(), theta.sin());
                } else {
                    matrix.set_entry(row_in_col_oper, col_oper_num, 0.0);
                    matrix.set_entry(col_oper_num, row_in_col_oper, 0.0);
                }
            }
        }

        // Keep nearest-neighbor couplings positive after the orthogonal rotations.
        for index in 0..(side - 1) {
            if matrix.get(index, index + 1) < 0.0 {
                matrix = matrix.flip_sign(index + 1);
            }
        }

        matrix.clean_small_values();
        matrix.topology = MatrixTopology::Folded;
        matrix
    }

    fn to_arrow(&self) -> Self {
        use super::rotations::{diagonal_rotation_angle, safe_angle};

        let mut matrix = self.clone();
        let order = matrix.order();

        for resonator in 1..order {
            for target in (resonator + 1)..=order {
                let target_index = resonator - 1;
                let pivot_a = resonator;
                let pivot_b = target;
                let sign = -1.0_f64;

                let theta = if target_index != pivot_a {
                    let numerator = matrix.get(target_index, pivot_b);
                    let denominator = matrix.get(target_index, pivot_a);
                    safe_angle(sign * numerator, denominator)
                } else {
                    diagonal_rotation_angle(&matrix, pivot_a, pivot_b)
                };

                matrix.apply_givens_similarity_inplace(pivot_a, pivot_b, theta.cos(), theta.sin());
            }
        }

        matrix.clean_small_values();
        matrix.topology = MatrixTopology::Arrow;
        matrix
    }

    #[allow(dead_code)]
    pub(crate) fn rotate_matrix(&self, row: usize, col: usize, column_operation: bool) -> Self {
        use super::rotations::safe_angle;

        let (pivot_a, pivot_b, theta) = if column_operation {
            let pivot_a = col - 1;
            let pivot_b = col;
            let numerator = -self.get(row, col);
            let denominator = self.get(row, col - 1);
            (pivot_a, pivot_b, safe_angle(numerator, denominator))
        } else {
            let pivot_a = row;
            let pivot_b = row + 1;
            let numerator = self.get(row, col);
            let denominator = self.get(row + 1, col);
            (pivot_a, pivot_b, safe_angle(numerator, denominator))
        };

        let cosine = theta.cos();
        let sine = theta.sin();

        let mut result = self.clone();
        result.apply_givens_similarity_inplace(pivot_a, pivot_b, cosine, sine);
        result
    }

    pub(crate) fn rotate_matrix_with_indices(
        &self,
        target_index: usize,
        pivot_a: usize,
        pivot_b: usize,
        sign: f64,
        axis: super::rotations::RotationAxis,
    ) -> Self {
        use super::rotations::{diagonal_rotation_angle, safe_angle, RotationAxis};

        let theta = match axis {
            RotationAxis::Row => {
                if target_index != pivot_b {
                    let numerator = self.get(pivot_a, target_index);
                    let denominator = self.get(pivot_b, target_index);
                    safe_angle(sign * numerator, denominator)
                } else {
                    diagonal_rotation_angle(self, pivot_a, pivot_b)
                }
            }
            RotationAxis::Column => {
                if target_index != pivot_a {
                    let numerator = self.get(target_index, pivot_b);
                    let denominator = self.get(target_index, pivot_a);
                    safe_angle(sign * numerator, denominator)
                } else {
                    diagonal_rotation_angle(self, pivot_a, pivot_b)
                }
            }
        };

        let cosine = theta.cos();
        let sine = theta.sin();

        let mut result = self.clone();
        result.apply_givens_similarity_inplace(pivot_a, pivot_b, cosine, sine);
        result
    }
}
