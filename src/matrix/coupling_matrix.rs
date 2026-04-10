use crate::error::{MfsError, Result};
use crate::freq::BandPassMapping;
use nalgebra::DMatrix;
use num_complex::Complex64;

/// Supported coupling-matrix topologies exposed by the library.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MatrixTopology {
    /// Standard transversal or otherwise untransformed matrix form.
    #[default]
    Transversal,
    /// Folded topology obtained by similarity rotations.
    Folded,
    /// Arrow topology obtained by similarity rotations.
    Arrow,
    /// Wheel topology obtained by similarity rotations.
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
    order: usize,
    topology: MatrixTopology,
    data: Vec<f64>,
}

/// Physical-frequency view of a normalized coupling matrix after band-pass scaling.
#[derive(Debug, Clone, PartialEq)]
pub struct BandPassScaledCouplingMatrix {
    matrix_hz: CouplingMatrix,
    source_external_q: f64,
    load_external_q: f64,
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

    /// Returns the underlying row-major storage.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Returns mutable row-major storage for internal synthesis helpers.
    pub(crate) fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Returns a new matrix transformed into the requested topology when supported.
    pub fn transform_topology(&self, topology: MatrixTopology) -> Result<Self> {
        match topology {
            MatrixTopology::Transversal => Ok(self.clone()),
            MatrixTopology::Folded => Ok(self.to_folded()),
            MatrixTopology::Arrow => Ok(self.to_arrow()),
            MatrixTopology::Wheel => {
                let mut matrix = self.to_arrow();
                matrix.topology = MatrixTopology::Wheel;
                Ok(matrix)
            }
        }
    }

    /// Extracts one trisection from the tail of the matrix and moves it to the requested center.
    ///
    /// `center_resonator` uses 1-based resonator numbering, excluding source and load.
    /// For example, `center_resonator = 2` targets the resonator triplet `(1, 2, 3)`.
    pub(crate) fn extract_triplet(
        &self,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<Self> {
        validate_triplet_center(self.order, center_resonator)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "triplet transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.at(tail, tail).unwrap_or_default();
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.at(tail - 1, tail).unwrap_or_default(), denominator)
        };

        let rotation = rotation_matrix_basic(matrix.order(), theta, tail - 1, tail)?;
        matrix = rotation.multiply(&matrix).multiply(&rotation.transpose());

        let move_steps = order - center_resonator - 1;
        for step in 0..move_steps {
            let pivot_a = order - step - 2;
            let pivot_b = order - step - 1;
            matrix = matrix.rotate_matrix_with_indices(
                order - step,
                pivot_a,
                pivot_b,
                1.0,
                RotationAxis::Row,
            );
        }

        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Extracts two neighboring trisections and merges them into a quadruplet.
    ///
    /// `position` matches the first triplet center in 1-based resonator numbering.
    /// `common_resonator` must be `1` or `4`, matching the two elimination formulas
    /// used by the Python prototype.
    pub(crate) fn extract_quadruplet(
        &self,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<Self> {
        validate_quadruplet_position(self.order, position)?;
        if common_resonator != 1 && common_resonator != 4 {
            return Err(MfsError::Unsupported(
                "common resonator for quadruplet extraction must be 1 or 4".to_string(),
            ));
        }

        let mut matrix = if swap_zero_order {
            self.extract_triplet(second_zero, position)?
                .extract_triplet(first_zero, position + 1)?
        } else {
            self.extract_triplet(first_zero, position)?
                .extract_triplet(second_zero, position + 1)?
        };

        let tail = position + 2;
        let theta = if common_resonator == 4 {
            -safe_angle(
                matrix.at(tail - 1, tail - 3).unwrap_or_default(),
                matrix.at(tail - 3, tail - 2).unwrap_or_default(),
            )
        } else {
            safe_angle(
                matrix.at(tail, tail - 2).unwrap_or_default(),
                matrix.at(tail - 1, tail).unwrap_or_default(),
            )
        };

        let rotation = rotation_matrix_basic(matrix.order(), theta, tail - 2, tail - 1)?;
        matrix = rotation.multiply(&matrix).multiply(&rotation.transpose());
        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Converts an arrow-style matrix into a trisection-centered topology.
    ///
    /// `zero_positions` uses 1-based resonator numbering and must span exactly
    /// one center resonator, for example `(2, 4)` to target a trisection centered
    /// on resonator `3`.
    pub(crate) fn extract_trisection(
        &self,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<Self> {
        if self.topology != MatrixTopology::Arrow {
            return Err(MfsError::Unsupported(format!(
                "trisection extraction requires Arrow input, got {:?}",
                self.topology
            )));
        }
        validate_trisection_positions(self.order, zero_positions)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "trisection transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.at(tail, tail).unwrap_or_default();
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.at(tail - 1, tail).unwrap_or_default(), denominator)
        };

        let rotation = rotation_matrix_basic(matrix.order(), theta, tail - 1, tail)?;
        matrix = rotation.multiply(&matrix).multiply(&rotation.transpose());

        let center_resonator = (zero_positions.0 + zero_positions.1) / 2;
        let pull_steps = order - 1 - center_resonator;
        for step in 0..pull_steps {
            matrix = matrix.rotate_matrix_with_indices(
                order - step,
                order - step - 2,
                order - step - 1,
                1.0,
                RotationAxis::Row,
            );
        }

        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Scales a normalized band-pass matrix into physical-frequency units.
    ///
    /// Internal resonator couplings and diagonal terms are converted into Hz.
    /// Source/load couplings remain as couplings, matching the Python helper's
    /// bandwidth-scaled mode rather than the external-Q representation.
    pub fn denormalize_bandpass(&self, mapping: &BandPassMapping) -> Result<Self> {
        let center_hz = mapping.center_hz();
        let bandwidth_hz = mapping.bandwidth_hz();
        let side = self.side();
        let mut data = vec![0.0; side * side];

        for row in 0..side {
            for col in 0..side {
                let value = self.at(row, col).unwrap_or_default();
                let scaled = if row == col {
                    if row == 0 || row == side - 1 {
                        0.0
                    } else {
                        denormalize_resonator_frequency(value, center_hz, bandwidth_hz)
                    }
                } else {
                    value * bandwidth_hz
                };
                data[row * side + col] = scaled;
            }
        }

        Self::new_with_topology(self.order, self.topology, data)
    }

    /// Converts a normalized matrix into a physical band-pass matrix plus port Q values.
    pub fn denormalize_bandpass_with_external_q(
        &self,
        mapping: &BandPassMapping,
    ) -> Result<BandPassScaledCouplingMatrix> {
        let mut matrix_hz = self.denormalize_bandpass(mapping)?;
        let fractional_bw = mapping.bandwidth_hz() / mapping.center_hz();
        let source_coupling = self.at(0, 1).unwrap_or_default();
        let load_coupling = self.at(self.order(), self.side() - 1).unwrap_or_default();
        let source_external_q = external_q_from_normalized_coupling(source_coupling, fractional_bw)?;
        let load_external_q = external_q_from_normalized_coupling(load_coupling, fractional_bw)?;

        matrix_hz.set_entry(0, 1, source_external_q);
        matrix_hz.set_entry(1, 0, source_external_q);
        matrix_hz.set_entry(self.order(), self.side() - 1, load_external_q);
        matrix_hz.set_entry(self.side() - 1, self.order(), load_external_q);

        Ok(BandPassScaledCouplingMatrix {
            matrix_hz,
            source_external_q,
            load_external_q,
        })
    }

    /// Converts a physical band-pass matrix back into normalized units.
    ///
    /// This expects source/load entries to still be couplings, not external Q values.
    pub fn normalize_bandpass(&self, mapping: &BandPassMapping) -> Result<Self> {
        let center_hz = mapping.center_hz();
        let bandwidth_hz = mapping.bandwidth_hz();
        let side = self.side();
        let mut data = vec![0.0; side * side];

        for row in 0..side {
            for col in 0..side {
                let value = self.at(row, col).unwrap_or_default();
                let normalized = if row == col {
                    if row == 0 || row == side - 1 {
                        0.0
                    } else {
                        normalize_resonator_frequency(value, center_hz, bandwidth_hz)?
                    }
                } else {
                    value / bandwidth_hz
                };
                data[row * side + col] = normalized;
            }
        }

        Self::new_with_topology(self.order, self.topology, data)
    }

    /// Converts a physical band-pass matrix that stores external Q values back into normalized form.
    pub fn normalize_bandpass_with_external_q(&self, mapping: &BandPassMapping) -> Result<Self> {
        let mut normalized = self.normalize_bandpass(mapping)?;
        let fractional_bw = mapping.bandwidth_hz() / mapping.center_hz();
        let source_q = self.at(0, 1).unwrap_or_default();
        let load_q = self.at(self.order(), self.side() - 1).unwrap_or_default();

        let source_coupling = normalized_coupling_from_external_q(source_q, fractional_bw)?;
        let load_coupling = normalized_coupling_from_external_q(load_q, fractional_bw)?;
        normalized.set_entry(0, 1, source_coupling);
        normalized.set_entry(1, 0, source_coupling);
        normalized.set_entry(self.order(), self.side() - 1, load_coupling);
        normalized.set_entry(self.side() - 1, self.order(), load_coupling);

        Ok(normalized)
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

    fn to_folded(&self) -> Self {
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
                    .at(row_oper_num, col_in_row_oper)
                    .unwrap_or_default()
                    .abs()
                    > 1e-7
                {
                    matrix = matrix.rotate_matrix(row_oper_num, col_in_row_oper, true);
                } else {
                    matrix.set_entry(row_oper_num, col_in_row_oper, 0.0);
                    matrix.set_entry(col_in_row_oper, row_oper_num, 0.0);
                }
            }

            let col_oper_num = side - 1 - row_oper_num;
            for row_in_col_oper in (row_oper_num + 2)..=(col_oper_num.saturating_sub(2)) {
                if matrix
                    .at(row_in_col_oper, col_oper_num)
                    .unwrap_or_default()
                    .abs()
                    > 1e-7
                {
                    matrix = matrix.rotate_matrix(row_in_col_oper, col_oper_num, false);
                } else {
                    matrix.set_entry(row_in_col_oper, col_oper_num, 0.0);
                    matrix.set_entry(col_oper_num, row_in_col_oper, 0.0);
                }
            }
        }

        // Keep nearest-neighbor couplings positive after the orthogonal rotations.
        for index in 0..(side - 1) {
            if matrix.at(index, index + 1).unwrap_or_default() < 0.0 {
                matrix = matrix.flip_sign(index + 1);
            }
        }

        matrix.clean_small_values();
        matrix.topology = MatrixTopology::Folded;
        matrix
    }

    fn to_arrow(&self) -> Self {
        let mut matrix = self.clone();
        let order = matrix.order();

        for resonator in 1..order {
            for target in (resonator + 1)..=order {
                matrix = matrix.rotate_matrix_with_indices(
                    resonator - 1,
                    resonator,
                    target,
                    -1.0,
                    RotationAxis::Column,
                );
            }
        }

        matrix.clean_small_values();
        matrix.topology = MatrixTopology::Arrow;
        matrix
    }

    fn rotate_matrix(&self, row: usize, col: usize, column_operation: bool) -> Self {
        let mut rotation = Self::identity(self.order()).expect("matrix order already validated");
        let (pivot_a, pivot_b, theta) = if column_operation {
            let pivot_a = col - 1;
            let pivot_b = col;
            let numerator = -self.at(row, col).unwrap_or_default();
            let denominator = self.at(row, col - 1).unwrap_or_default();
            (pivot_a, pivot_b, safe_angle(numerator, denominator))
        } else {
            let pivot_a = row;
            let pivot_b = row + 1;
            let numerator = self.at(row, col).unwrap_or_default();
            let denominator = self.at(row + 1, col).unwrap_or_default();
            (pivot_a, pivot_b, safe_angle(numerator, denominator))
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        rotation.set_entry(pivot_a, pivot_a, cosine);
        rotation.set_entry(pivot_b, pivot_b, cosine);
        rotation.set_entry(pivot_b, pivot_a, sine);
        rotation.set_entry(pivot_a, pivot_b, -sine);

        rotation
            .multiply(self)
            .multiply(&rotation.transpose())
    }

    fn rotate_matrix_with_indices(
        &self,
        target_index: usize,
        pivot_a: usize,
        pivot_b: usize,
        sign: f64,
        axis: RotationAxis,
    ) -> Self {
        let mut rotation = Self::identity(self.order()).expect("matrix order already validated");
        let theta = match axis {
            RotationAxis::Row => {
                if target_index != pivot_b {
                    let numerator = self.at(pivot_a, target_index).unwrap_or_default();
                    let denominator = self.at(pivot_b, target_index).unwrap_or_default();
                    safe_angle(sign * numerator, denominator)
                } else {
                    diagonal_rotation_angle(self, pivot_a, pivot_b)
                }
            }
            RotationAxis::Column => {
                if target_index != pivot_a {
                    let numerator = self.at(target_index, pivot_b).unwrap_or_default();
                    let denominator = self.at(target_index, pivot_a).unwrap_or_default();
                    safe_angle(sign * numerator, denominator)
                } else {
                    diagonal_rotation_angle(self, pivot_a, pivot_b)
                }
            }
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        rotation.set_entry(pivot_a, pivot_a, cosine);
        rotation.set_entry(pivot_b, pivot_b, cosine);
        rotation.set_entry(pivot_b, pivot_a, sine);
        rotation.set_entry(pivot_a, pivot_b, -sine);

        rotation
            .multiply(self)
            .multiply(&rotation.transpose())
    }

    fn multiply(&self, rhs: &Self) -> Self {
        let side = self.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for col in 0..side {
                let mut acc = 0.0;
                for inner in 0..side {
                    acc += self.at(row, inner).unwrap_or_default() * rhs.at(inner, col).unwrap_or_default();
                }
                data[row * side + col] = acc;
            }
        }
        Self {
            order: self.order,
            topology: self.topology,
            data,
        }
    }

    fn transpose(&self) -> Self {
        let side = self.side();
        let mut data = vec![0.0; side * side];
        for row in 0..side {
            for col in 0..side {
                data[col * side + row] = self.at(row, col).unwrap_or_default();
            }
        }
        Self {
            order: self.order,
            topology: self.topology,
            data,
        }
    }

    fn flip_sign(&self, diagonal_index: usize) -> Self {
        let mut reflection = Self::identity(self.order()).expect("matrix order already validated");
        reflection.set_entry(diagonal_index, diagonal_index, -1.0);
        reflection
            .multiply(self)
            .multiply(&reflection.transpose())
    }

    fn set_entry(&mut self, row: usize, col: usize, value: f64) {
        let side = self.side();
        self.data[row * side + col] = value;
    }

    fn clean_small_values(&mut self) {
        for value in &mut self.data {
            if value.abs() <= 1e-10 {
                *value = 0.0;
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum RotationAxis {
    Row,
    Column,
}

fn rotation_matrix_basic(order: usize, theta: f64, row: usize, col: usize) -> Result<CouplingMatrix> {
    let mut rotation = CouplingMatrix::identity(order)?;
    let cosine = theta.cos();
    let sine = theta.sin();
    rotation.set_entry(row, row, cosine);
    rotation.set_entry(col, col, cosine);
    rotation.set_entry(row, col, -sine);
    rotation.set_entry(col, row, sine);
    Ok(rotation)
}

fn safe_angle(y: f64, x: f64) -> f64 {
    if x.abs() < 1e-10 {
        if y.is_sign_positive() {
            -std::f64::consts::FRAC_PI_2
        } else {
            std::f64::consts::FRAC_PI_2
        }
    } else {
        (y / x).atan()
    }
}

fn diagonal_rotation_angle(matrix: &CouplingMatrix, index_a: usize, index_b: usize) -> f64 {
    let diagonal_delta = matrix
        .at(index_b, index_b)
        .unwrap_or_default()
        - matrix.at(index_a, index_a).unwrap_or_default();
    if diagonal_delta.abs() < 1e-10 {
        0.0
    } else {
        let ratio = (2.0 * matrix.at(index_a, index_b).unwrap_or_default()) / diagonal_delta;
        0.5 * safe_angle(ratio, 1.0)
    }
}

fn denormalize_resonator_frequency(normalized: f64, center_hz: f64, bandwidth_hz: f64) -> f64 {
    let fractional_bw = bandwidth_hz / center_hz;
    center_hz
        * ((1.0 + (normalized * fractional_bw / 2.0).powi(2)).sqrt()
            - normalized * fractional_bw / 2.0)
}

fn normalize_resonator_frequency(physical_hz: f64, center_hz: f64, bandwidth_hz: f64) -> Result<f64> {
    if !physical_hz.is_finite() || physical_hz <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "physical resonator frequency must be > 0, got {physical_hz}"
        )));
    }

    Ok((center_hz / bandwidth_hz) * (center_hz / physical_hz - physical_hz / center_hz))
}

fn external_q_from_normalized_coupling(coupling: f64, fractional_bw: f64) -> Result<f64> {
    if !coupling.is_finite() || coupling.abs() <= 1e-12 {
        return Err(MfsError::InvalidFrequency(
            "normalized source/load coupling must be non-zero when converting to external Q"
                .to_string(),
        ));
    }

    Ok(1.0 / (coupling * coupling * fractional_bw))
}

fn normalized_coupling_from_external_q(external_q: f64, fractional_bw: f64) -> Result<f64> {
    if !external_q.is_finite() || external_q <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "external Q must be > 0, got {external_q}"
        )));
    }

    Ok((1.0 / (external_q * fractional_bw)).sqrt())
}

fn validate_triplet_center(order: usize, center_resonator: usize) -> Result<()> {
    if center_resonator < 2 || center_resonator >= order {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "triplet center must be in [2, {}), got {center_resonator}",
            order
        )));
    }

    Ok(())
}

fn validate_quadruplet_position(order: usize, position: usize) -> Result<()> {
    if position < 2 || position + 1 >= order {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "quadruplet position must leave room for two adjacent triplets, got {position} for order {order}"
        )));
    }

    Ok(())
}

fn validate_trisection_positions(order: usize, zero_positions: (usize, usize)) -> Result<()> {
    let (start, end) = zero_positions;
    if start < 1 || end > order || start >= end {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "trisection zero positions must be ordered resonator indices within 1..={order}, got ({start}, {end})"
        )));
    }
    if end - start != 2 {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "trisection zero positions must differ by exactly 2, got ({start}, {end})"
        )));
    }

    Ok(())
}
