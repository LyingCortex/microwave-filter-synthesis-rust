mod builder;
mod coupling_matrix;
mod synthesis;

pub use builder::CouplingMatrixBuilder;
pub use coupling_matrix::{CouplingMatrix, MatrixShape};
pub use synthesis::CouplingMatrixSynthesizer;

#[cfg(test)]
mod tests {
    use crate::approx::PolynomialSet;
    use crate::error::Result;

    use super::{CouplingMatrixBuilder, CouplingMatrixSynthesizer, MatrixShape};

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn builder_can_set_symmetric_entries() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(2)?
            .set_symmetric(0, 1, 0.75)?
            .build()?;

        approx_eq(matrix.at(0, 1).unwrap_or_default(), 0.75, 1e-12);
        approx_eq(matrix.at(1, 0).unwrap_or_default(), 0.75, 1e-12);
        Ok(())
    }

    #[test]
    fn matrix_synthesizer_builds_matrix_from_polynomials() -> Result<()> {
        let polynomials = PolynomialSet::new(
            3,
            0.1,
            vec![-1.5, 1.5],
            vec![1.0, 0.2, 0.3, 0.4],
            vec![0.8, 0.6, 0.4, 0.2],
            vec![1.0, 0.5, -2.25],
        )?;

        let matrix = CouplingMatrixSynthesizer.synthesize(&polynomials)?;
        assert_eq!(matrix.shape(), MatrixShape { rows: 5, cols: 5 });
        approx_eq(matrix.at(0, 1).unwrap_or_default(), 0.8, 1e-12);
        approx_eq(matrix.at(1, 2).unwrap_or_default(), 0.4, 1e-12);
        approx_eq(matrix.at(2, 2).unwrap_or_default(), 0.5, 1e-12);
        approx_eq(matrix.at(3, 3).unwrap_or_default(), -2.25, 1e-12);
        approx_eq(matrix.at(3, 4).unwrap_or_default(), 1.0, 1e-12);
        Ok(())
    }
}
