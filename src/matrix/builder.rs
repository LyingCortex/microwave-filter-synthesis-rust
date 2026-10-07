use crate::error::{MfsError, Result};

use super::{CouplingMatrix, MatrixTopology};

/// Builder for assembling a dense coupling matrix entry by entry.
#[derive(Debug, Default)]
pub struct CouplingMatrixBuilder {
    order: usize,
    topology: MatrixTopology,
    allow_asymmetric: bool,
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
            topology: MatrixTopology::Transversal,
            allow_asymmetric: false,
            data: vec![0.0; side * side],
        })
    }

    /// Overrides the topology metadata attached to the built matrix.
    pub fn topology(mut self, topology: MatrixTopology) -> Self {
        self.topology = topology;
        self
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

    /// Permits `build()` to return a matrix whose off-diagonal entries are
    /// not symmetric. The default (`false`) rejects a one-sided off-diagonal
    /// `set`, which would otherwise silently produce a matrix every downstream
    /// consumer (response solver, multiplication, section extraction) assumes
    /// is symmetric. Intended for the folded-mask gate fixtures, which build
    /// a deliberately asymmetric matrix. The symmetry check itself lives in
    /// `CouplingMatrix::new_with_topology`; this flag routes the build to the
    /// unchecked `new_with_topology_asymmetric` instead.
    pub fn allow_asymmetric(mut self, allowed: bool) -> Self {
        self.allow_asymmetric = allowed;
        self
    }

    /// Finalizes the builder. By default it rejects a one-sided off-diagonal
    /// entry (an asymmetric coupling matrix); opt out with `allow_asymmetric`.
    pub fn build(self) -> Result<CouplingMatrix> {
        if self.allow_asymmetric {
            CouplingMatrix::new_with_topology_asymmetric(self.order, self.topology, self.data)
        } else {
            CouplingMatrix::new_with_topology(self.order, self.topology, self.data)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_symmetric_builds_symmetric_matrix() -> Result<()> {
        let m = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.8)?
            .set(1, 1, -0.3)?
            .build()?;
        assert!((m.get(0, 1) - 0.8).abs() < 1e-15);
        assert!((m.get(1, 0) - 0.8).abs() < 1e-15);
        assert!((m.get(1, 1) + 0.3).abs() < 1e-15);
        Ok(())
    }

    #[test]
    fn one_sided_off_diagonal_set_is_rejected_by_default() {
        let err = CouplingMatrixBuilder::new(2)
            .unwrap()
            .set(0, 1, 0.8)
            .unwrap()
            .build();
        assert!(
            matches!(err, Err(MfsError::PreconditionViolation(_))),
            "a one-sided off-diagonal entry must be rejected by default, got {err:?}"
        );
    }

    #[test]
    fn allow_asymmetric_permits_one_sided_off_diagonal_set() -> Result<()> {
        let m = CouplingMatrixBuilder::new(2)?
            .set(0, 1, 0.8)?
            .allow_asymmetric(true)
            .build()?;
        assert!((m.get(0, 1) - 0.8).abs() < 1e-15);
        assert!((m.get(1, 0)).abs() < 1e-15, "the other entry must stay zero");
        Ok(())
    }
}
