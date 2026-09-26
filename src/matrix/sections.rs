use crate::error::{MfsError, Result};

use super::core::{CouplingMatrix, MatrixTopology};
use super::rotations::{safe_angle, RotationAxis};

impl CouplingMatrix {
    /// Extracts one trisection from the tail of the matrix and moves it to the requested center.
    ///
    /// `center_resonator` uses 1-based resonator numbering, excluding source and load.
    /// For example, `center_resonator = 2` targets the resonator triplet `(1, 2, 3)`.
    pub(crate) fn extract_triplet(
        &self,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<Self> {
        validate_triplet_center(self.order(), center_resonator)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "triplet transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.get(tail, tail);
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.get(tail - 1, tail), denominator)
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 1, tail, cosine, sine);

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
        validate_quadruplet_position(self.order(), position)?;
        if common_resonator != 1 && common_resonator != 4 {
            return Err(MfsError::PreconditionViolation(
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
                matrix.get(tail - 1, tail - 3),
                matrix.get(tail - 3, tail - 2),
            )
        } else {
            safe_angle(
                matrix.get(tail, tail - 2),
                matrix.get(tail - 1, tail),
            )
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 2, tail - 1, cosine, sine);
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
        if self.topology() != MatrixTopology::Arrow {
            return Err(MfsError::PreconditionViolation(
                "extract_trisection requires an Arrow-topology matrix (trisection patterns only exist in the arrow reduction); convert with transform_topology(MatrixTopology::Arrow) first".to_string(),
            ));
        }
        validate_trisection_positions(self.order(), zero_positions)?;
        if !transmission_zero.is_finite() {
            return Err(MfsError::InvalidTransmissionZero(
                "trisection transmission zero must be finite".to_string(),
            ));
        }

        let mut matrix = self.clone();
        let order = matrix.order();
        let tail = order;
        let denominator = transmission_zero + matrix.get(tail, tail);
        let theta = if denominator.abs() < 1e-10 {
            std::f64::consts::FRAC_PI_2
        } else {
            safe_angle(matrix.get(tail - 1, tail), denominator)
        };

        let cosine = theta.cos();
        let sine = theta.sin();
        matrix.apply_givens_similarity_inplace(tail - 1, tail, cosine, sine);

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
}

fn validate_triplet_center(order: usize, center_resonator: usize) -> Result<()> {
    if center_resonator < 1 || center_resonator >= order {
        return Err(MfsError::InvalidTransmissionZero(format!(
            "triplet center must be in [1, {}), got {center_resonator}",
            order
        )));
    }

    Ok(())
}

fn validate_quadruplet_position(order: usize, position: usize) -> Result<()> {
    if position < 1 || position + 1 >= order {
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

impl CouplingMatrix {
    /// Extract cascade trisections (CT topology) for multiple zeros.
    ///
    /// Applies trisection extraction sequentially for each zero.
    /// Input should be in Arrow form.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    pub(crate) fn extract_cascade_trisections(
        &self,
        zeros: &[f64],
        positions: &[(usize, usize)],
    ) -> Result<Self> {
        if zeros.len() != positions.len() {
            return Err(MfsError::DimensionMismatch {
                expected: zeros.len(),
                actual: positions.len(),
            });
        }
        let mut matrix = self.clone();
        for (tz, pos) in zeros.iter().zip(positions.iter()) {
            matrix = matrix.extract_trisection(*tz, *pos)?;
        }
        matrix.ensure_positive_main_couplings();
        Ok(matrix)
    }

    /// Extract cascade quadruplets (CQ topology) for paired zeros.
    ///
    /// Processes zeros in pairs as quadruplets. An odd remaining zero
    /// is extracted as a standalone trisection.
    /// Input should be in Arrow form.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    pub(crate) fn extract_cascade_quadruplets(
        &self,
        zeros: &[f64],
        z_locs: &[[usize; 2]],
    ) -> Result<Self> {
        if zeros.len() != z_locs.len() {
            return Err(MfsError::DimensionMismatch {
                expected: zeros.len(),
                actual: z_locs.len(),
            });
        }

        let mut matrix = self.clone();
        let _num_zeros = zeros.len();

        // Classify z_locs into CQ pairs and CT singles
        // span=3 → quadruplet diagonal, span=2 → trisection
        let mut quad_indices: Vec<(usize, usize, usize)> = Vec::new(); // (zero_idx, start, end)
        let mut tri_indices: Vec<(usize, usize, usize)> = Vec::new();  // (zero_idx, start, end)

        for (i, z) in z_locs.iter().enumerate() {
            let span = z[1].abs_diff(z[0]);
            let start = z[0].min(z[1]);
            let end = z[0].max(z[1]);
            match span {
                3 => quad_indices.push((i, start, end)),
                2 => tri_indices.push((i, start, end)),
                _ => return Err(MfsError::InvalidTransmissionZero(
                    format!("z_loc[{}]=[{},{}] has invalid span {}. Must be 2 (CT) or 3 (CQ).", i, z[0], z[1], span)
                )),
            }
        }

        // Match trisections to quadruplets
        let mut used_tri = vec![false; tri_indices.len()];

        // Process quadruplets (from high position to low = tail first)
        let mut cq_ops: Vec<(f64, f64, usize, usize)> = Vec::new(); // (z1, z2, position, common_res)
        for &(q_zero_idx, q_start, _q_end) in quad_indices.iter() {
            // Find matching trisection within this quadruplet range
            let mut partner = None;
            for (ti, &(t_zero_idx, t_start, t_end)) in tri_indices.iter().enumerate() {
                if used_tri[ti] { continue; }
                if t_start >= q_start && t_end <= _q_end {
                    partner = Some((ti, t_zero_idx, t_start));
                    break;
                }
            }
            if let Some((ti, t_zero_idx, t_start)) = partner {
                used_tri[ti] = true;
                let position = q_start + 1; // triplet center
                let common_res = if t_start == q_start { 1 } else { 4 };
                cq_ops.push((zeros[t_zero_idx], zeros[q_zero_idx], position, common_res));
            } else {
                return Err(MfsError::InvalidTransmissionZero(
                    format!("Quadruplet ({},{}) has no matching trisection partner", q_start, _q_end)
                ));
            }
        }

        // Collect remaining CTs
        let mut ct_ops: Vec<(f64, usize)> = Vec::new(); // (zero, center)
        for (ti, &(t_zero_idx, t_start, t_end)) in tri_indices.iter().enumerate() {
            if !used_tri[ti] {
                let center = (t_start + t_end) / 2;
                ct_ops.push((zeros[t_zero_idx], center));
            }
        }

        // Sort all operations by position descending (extract from tail first)
        let mut all_ops: Vec<(usize, bool)> = Vec::new(); // (position, is_cq)
        for &(_, _, pos, _) in cq_ops.iter() {
            all_ops.push((pos, true));
        }
        for &(_, center) in ct_ops.iter() {
            all_ops.push((center, false));
        }
        // We need indices after sorting
        let mut cq_order: Vec<usize> = (0..cq_ops.len()).collect();
        cq_order.sort_by(|a, b| cq_ops[*b].2.cmp(&cq_ops[*a].2));
        let mut ct_order: Vec<usize> = (0..ct_ops.len()).collect();
        ct_order.sort_by(|a, b| ct_ops[*b].1.cmp(&ct_ops[*a].1));

        // Interleave by position (highest first)
        let mut cq_iter = cq_order.iter().peekable();
        let mut ct_iter = ct_order.iter().peekable();

        loop {
            let cq_pos = cq_iter.peek().map(|&&i| cq_ops[i].2);
            let ct_pos = ct_iter.peek().map(|&&i| ct_ops[i].1);
            match (cq_pos, ct_pos) {
                (Some(cp), Some(tp)) if cp >= tp => {
                    let &i = cq_iter.next().unwrap();
                    let (z1, z2, pos, cr) = cq_ops[i];
                    matrix = matrix.extract_quadruplet(z1, z2, pos, cr, false)?;
                }
                (_, Some(_)) => {
                    let &i = ct_iter.next().unwrap();
                    let (z, center) = ct_ops[i];
                    matrix = matrix.extract_triplet(z, center)?;
                }
                (Some(_), None) => {
                    let &i = cq_iter.next().unwrap();
                    let (z1, z2, pos, cr) = cq_ops[i];
                    matrix = matrix.extract_quadruplet(z1, z2, pos, cr, false)?;
                }
                (None, None) => break,
            }
        }

        matrix.ensure_positive_main_couplings();
        Ok(matrix)
    }

    /// Ensure all main (nearest-neighbor) couplings are positive.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    fn ensure_positive_main_couplings(&mut self) {
        let side = self.side();
        for i in 0..(side - 1) {
            if self.get(i, i + 1) < 0.0 {
                self.flip_sign_inplace(i + 1);
            }
        }
    }

    // TODO: This mirrors `CouplingMatrix::flip_sign` in `src/matrix/core.rs`
    // (same row/column negation, same diagonal preservation), differing only
    // in mutating in place instead of returning a new matrix. Keep the two in
    // sync: a similarity sign flip must leave the diagonal unchanged
    // (R*M*R^T with R[i][i] = -1). If one is edited, the other must be checked
    // too. Consider deduplicating into a single shared helper.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    fn flip_sign_inplace(&mut self, index: usize) {
        let side = self.side();
        for col in 0..side {
            let val = self.data[index * side + col];
            self.data[index * side + col] = -val;
        }
        for row in 0..side {
            let val = self.data[row * side + index];
            self.data[row * side + index] = -val;
        }
    }
}

impl CouplingMatrix {
    /// Extract BOX topology from Arrow matrix.
    ///
    /// BOX is a 4-resonator structure with cross-couplings M(k,k+2) and M(k+1,k+3),
    /// and eliminated inline M(k+1,k+2) = 0.
    ///
    /// Algorithm: extract_triplet → Jacobi rotation to eliminate M[center, center+1]
    ///
    /// `span` is (start, end) with end-start=3, representing the 4-resonator block.
    /// Internally: center = start+1, used for triplet extraction.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    pub(crate) fn extract_box_section(
        &self,
        transmission_zero: f64,
        span_start: usize,
    ) -> Result<Self> {
        let center = span_start + 1;

        // Step 1: extract trisection
        let mut matrix = self.extract_triplet(transmission_zero, center)?;

        // Step 2: Jacobi rotation in plane (center, center+1) to eliminate M[center, center+1]
        let i = center;
        let j = center + 1;
        let mij = matrix.get(i, j);
        let mii = matrix.get(i, i);
        let mjj = matrix.get(j, j);

        let denom = mii - mjj;
        let theta = if denom.abs() > 1e-12 {
            -0.5 * (2.0 * mij / denom).atan()
        } else {
            if mij > 0.0 { -std::f64::consts::FRAC_PI_4 } else { std::f64::consts::FRAC_PI_4 }
        };

        matrix.apply_givens_similarity_inplace(i, j, theta.cos(), theta.sin());
        matrix.clean_small_values();
        Ok(matrix)
    }

    /// Extract multiple BOX sections.
    /// Each span is [start, end] with end-start=3.
    // kept: real algorithm, not yet wired to a topology path (see P7 dead-code policy)
    #[allow(dead_code)]
    pub(crate) fn extract_box_sections(
        &self,
        zeros: &[f64],
        spans: &[[usize; 2]],
    ) -> Result<Self> {
        if zeros.len() != spans.len() {
            return Err(MfsError::DimensionMismatch {
                expected: zeros.len(),
                actual: spans.len(),
            });
        }
        let mut matrix = self.clone();
        // Process from high position to low (tail first)
        let mut ops: Vec<(f64, usize)> = zeros.iter()
            .zip(spans.iter())
            .map(|(&z, s)| (z, s[0].min(s[1])))
            .collect();
        ops.sort_by_key(|op| std::cmp::Reverse(op.1));

        for (tz, start) in ops {
            matrix = matrix.extract_box_section(tz, start)?;
        }
        matrix.ensure_positive_main_couplings();
        Ok(matrix)
    }
}
