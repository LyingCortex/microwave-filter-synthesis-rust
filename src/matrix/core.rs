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
    #[deprecated(
        note = "Use Arrow topology instead. Real Wheel implementation is not yet available."
    )]
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
    ///
    /// Rejects a one-sided off-diagonal entry (`M[i][j]` set with `M[j][i]`
    /// not equal): the response solver and every other consumer assume a
    /// symmetric real matrix — `backend.rs` solves `A x = e_0` once and reads
    /// the last *row* of A⁻¹ as the transpose of the last *column*, which
    /// holds only for a symmetric A. An asymmetric matrix would produce a
    /// silently wrong S21/S12, not a loud failure (the same silent-corruption
    /// class the `flip_sign` diagonal bug was). Use
    /// [`Self::new_with_topology_asymmetric`] for the deliberately asymmetric
    /// matrices the folded-mask gate fixtures build.
    pub fn new_with_topology(
        order: usize,
        topology: MatrixTopology,
        data: Vec<f64>,
    ) -> Result<Self> {
        let m = Self::new_with_topology_asymmetric(order, topology, data)?;
        m.assert_symmetric()?;
        Ok(m)
    }

    /// Rejects a one-sided off-diagonal entry (`M[i][j]` not equal to
    /// `M[j][i]`).
    ///
    /// The response solver, multiplication, and section extraction all assume
    /// a symmetric real matrix: `backend.rs` solves `A x = e_0` once and reads
    /// the last *row* of A⁻¹ as the transpose of the last *column*, which
    /// holds only for a symmetric A. An asymmetric matrix would produce
    /// silently wrong S21/S12 and group delay, not a loud failure. This is the
    /// consumer-side half of that invariant (the constructors enforce it via
    /// [`Self::new_with_topology`]).
    ///
    /// Tolerance-based and scale-aware: a folded (Givens-rotated) matrix is
    /// symmetric only to within rounding (~1e-15), so exact `a != b` would
    /// false-reject legitimate rotated matrices. The threshold rejects a gross
    /// one-sided entry (the silent-corruption bug) while permitting the
    /// rotation/rounding asymmetry every rotated or perturbed matrix carries.
    pub(crate) fn assert_symmetric(&self) -> Result<()> {
        let side = self.side();
        for i in 0..side {
            for j in (i + 1)..side {
                let a = self.data[i * side + j];
                let b = self.data[j * side + i];
                let tol = 1e-12 * 1.0_f64.max(a.abs().max(b.abs()));
                if (a - b).abs() > tol {
                    return Err(MfsError::PreconditionViolation(format!(
                        "asymmetric coupling matrix: M[{i}][{j}] = {a} != \
                         M[{j}][{i}] = {b} (coupling matrices are symmetric; \
                         write both entries, or use new_with_topology_asymmetric \
                         for a deliberately asymmetric matrix)"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Creates a coupling matrix with an explicit topology label, **without**
    /// the symmetry check. For the deliberately asymmetric matrices the
    /// folded-mask gate fixtures build to exercise `folded_mask_violations`.
    /// Never use this on a matrix that feeds the response solver, multiplication,
    /// or section extraction — those assume symmetry.
    pub fn new_with_topology_asymmetric(
        order: usize,
        topology: MatrixTopology,
        data: Vec<f64>,
    ) -> Result<Self> {
        let _side = validate_and_side(order, &data)?;
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
    pub(crate) fn apply_givens_similarity_inplace(
        &mut self,
        pivot_a: usize,
        pivot_b: usize,
        cosine: f64,
        sine: f64,
    ) {
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
        // With R^T[a,a] = cos, R^T[a,b] = sin, R^T[b,a] = -sin, R^T[b,b] = cos:
        // col_a'' = cos * col_a' - sin * col_b'  (of M')
        //    col_b'' = sin * col_a' + cos * col_b'  (of M')
        for row in 0..side {
            let a = self.data[row * side + pivot_a];
            let b = self.data[row * side + pivot_b];
            self.data[row * side + pivot_a] = cosine * a - sine * b;
            self.data[row * side + pivot_b] = sine * a + cosine * b;
        }
    }

    // TODO: This mirrors `CouplingMatrix::flip_sign_inplace` in
    // `src/matrix/sections.rs` (same row/column negation, same diagonal
    // preservation), differing only in returning a new matrix instead of
    // mutating in place. Keep the two in sync: a similarity sign flip must
    // leave the diagonal unchanged (R*M*R^T with R[i][i] = -1). If one is
    // edited, the other must be checked too. Consider deduplicating into a
    // single shared helper.
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
        // The diagonal element has now been negated twice (once by the row
        // negation, once by the column negation), so it is back to its original
        // value. A proper similarity flip R*M*R^T (R[i][i] = -1) leaves the
        // diagonal unchanged, so we must NOT negate it a third time. Negating
        // it again here was a bug that flipped the self-coupling sign of any
        // flipped node, corrupting the response (lost transmission zeros) for
        // even-order asymmetric filters.
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
            for col_in_row_oper in ((row_oper_num + 2)..=(filter_order - row_oper_num)).rev() {
                if matrix.get(row_oper_num, col_in_row_oper).abs() > 1e-7 {
                    // Column operation: zero out (row_oper_num, col_in_row_oper)
                    let pivot_a = col_in_row_oper - 1;
                    let pivot_b = col_in_row_oper;
                    let numerator = -matrix.get(row_oper_num, col_in_row_oper);
                    let denominator = matrix.get(row_oper_num, col_in_row_oper - 1);
                    let theta = safe_angle(numerator, denominator);
                    matrix.apply_givens_similarity_inplace(
                        pivot_a,
                        pivot_b,
                        theta.cos(),
                        theta.sin(),
                    );
                } else {
                    matrix.set_entry(row_oper_num, col_in_row_oper, 0.0);
                    matrix.set_entry(col_in_row_oper, row_oper_num, 0.0);
                }
            }

            let col_oper_num = side - 1 - row_oper_num;
            for row_in_col_oper in (row_oper_num + 2)..=(col_oper_num.saturating_sub(2)) {
                if matrix.get(row_in_col_oper, col_oper_num).abs() > 1e-7 {
                    // Row operation: zero out (row_in_col_oper, col_oper_num)
                    let pivot_a = row_in_col_oper;
                    let pivot_b = row_in_col_oper + 1;
                    let numerator = matrix.get(row_in_col_oper, col_oper_num);
                    let denominator = matrix.get(row_in_col_oper + 1, col_oper_num);
                    let theta = safe_angle(numerator, denominator);
                    matrix.apply_givens_similarity_inplace(
                        pivot_a,
                        pivot_b,
                        theta.cos(),
                        theta.sin(),
                    );
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

/// Validates order (non-zero) and data length (side²); returns the side.
fn validate_and_side(order: usize, data: &[f64]) -> Result<usize> {
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
    Ok(side)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A symmetric order-2 (side 4) matrix: source couplings 0.8 / 0.5, chain
    /// 0.6, load couplings 0.7 / 0.4. All off-diagonal pairs equal.
    fn symmetric_data() -> Vec<f64> {
        // Order-2, side-4. Indices are row*4+col; every off-diagonal pair
        // below is written twice with the same value (symmetric).
        let mut d = vec![0.0_f64; 16];
        // source couplings (row 0): M[0][k] == M[k][0]
        d[1] = 0.8; d[4] = 0.8; // M[0][1]/M[1][0]
        d[2] = 0.5; d[8] = 0.5; // M[0][2]/M[2][0]
        // chain: M[1][2] == M[2][1]
        d[6] = 0.6; d[9] = 0.6;
        // load couplings (col 3): M[k][3] == M[3][k]
        d[7] = 0.7; d[13] = 0.7; // M[1][3]/M[3][1]
        d[11] = 0.4; d[14] = 0.4; // M[2][3]/M[3][2]
        d
    }

    #[test]
    fn new_with_topology_accepts_symmetric() -> Result<()> {
        let m = CouplingMatrix::new_with_topology(2, MatrixTopology::Transversal, symmetric_data())?;
        assert_eq!(m.order(), 2);
        assert_eq!(m.topology(), MatrixTopology::Transversal);
        assert!((m.at(1, 2).unwrap() - 0.6).abs() < 1e-15);
        Ok(())
    }

    #[test]
    fn new_with_topology_rejects_one_sided_off_diagonal() {
        let mut d = symmetric_data();
        d[9] = 0.9; // M[2][1] = 0.9 — one-sided: M[2][1] != M[1][2]=0.6
        let err = CouplingMatrix::new_with_topology(2, MatrixTopology::Transversal, d).unwrap_err();
        assert!(
            matches!(err, MfsError::PreconditionViolation(_)),
            "a one-sided off-diagonal entry must be rejected, got {err:?}"
        );
        // The error must name the offending pair so a user can find it.
        assert!(
            err.to_string().contains("M[1][2]") && err.to_string().contains("M[2][1]"),
            "error must name both entries of the asymmetric pair: {err}"
        );
    }

    #[test]
    fn new_with_topology_asymmetric_permits_one_sided_off_diagonal() -> Result<()> {
        let mut d = symmetric_data();
        d[9] = 0.9; // M[2][1] = 0.9 (M[1][2] stays 0.6 — deliberately one-sided)
        let m = CouplingMatrix::new_with_topology_asymmetric(2, MatrixTopology::Transversal, d)?;
        assert!((m.at(1, 2).unwrap() - 0.6).abs() < 1e-15);
        assert!((m.at(2, 1).unwrap() - 0.9).abs() < 1e-15);
        Ok(())
    }

    #[test]
    fn new_with_topology_accepts_diagonal_only() -> Result<()> {
        // A diagonal-only (all-off-diagonal-zero) matrix is trivially symmetric.
        let d = vec![0.0_f64; 9]; // order 1, side 3
        let m = CouplingMatrix::new_with_topology(1, MatrixTopology::Transversal, d).unwrap();
        assert_eq!(m.order(), 1);
        Ok(())
    }

    /// Pins the *lower* side of the symmetry tolerance: a relative difference
    /// of ~1e-13 — the level a Givens-rotated (folded) matrix legitimately
    /// carries after rounding — must be **accepted**. If this test starts
    /// failing, the tolerance has been over-tightened (e.g. back to exact
    /// `a != b`), which would false-reject every rotated/perturbed matrix.
    #[test]
    fn new_with_topology_permits_rotation_level_asymmetry() -> Result<()> {
        let mut d = symmetric_data();
        // M[2][1] differs from M[1][2]=0.6 by a relative 1e-13.
        d[9] = 0.6 * (1.0 + 1e-13);
        let m = CouplingMatrix::new_with_topology(2, MatrixTopology::Transversal, d)?;
        assert_eq!(m.order(), 2);
        Ok(())
    }

    /// Pins the *upper* side of the symmetry tolerance: a relative difference
    /// of ~1e-6 — far above rotation/rounding noise, but far below a gross
    /// one-sided entry — must be **rejected**. The threshold
    /// (`1e-12 * max(1,|a|,|b|)`) sits between the two boundary tests, so a
    /// loosened tolerance that would let real asymmetry through is caught here.
    #[test]
    fn new_with_topology_rejects_relative_asymmetry_above_noise() {
        let mut d = symmetric_data();
        // M[2][1] differs from M[1][2]=0.6 by a relative 1e-6.
        d[9] = 0.6 * (1.0 + 1e-6);
        let err = CouplingMatrix::new_with_topology(2, MatrixTopology::Transversal, d).unwrap_err();
        assert!(
            matches!(err, MfsError::PreconditionViolation(_)),
            "a relative-1e-6 asymmetric pair must be rejected, got {err:?}"
        );
    }
}
