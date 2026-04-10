use crate::approx::PolynomialSet;
use crate::error::Result;
use crate::matrix::{CouplingMatrix, CouplingMatrixBuilder};

/// Builds the current placeholder chain-style coupling matrix from polynomial metadata.
pub(crate) fn synthesize_placeholder_matrix(polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
    let order = polynomials.order;
    let mut builder = CouplingMatrixBuilder::new(order)?;
    let projected_f = polynomials.f_real_projection();
    let projected_e = polynomials.e_real_projection();
    let projected_p = polynomials.p_real_projection();

    let source_coupling = projected_f
        .get(0)
        .copied()
        .unwrap_or(1.0)
        .abs()
        .max(1e-12);
    // Keep the source/load couplings non-zero so the response matrix remains invertible.
    let load_coupling = projected_e
        .get(0)
        .copied()
        .unwrap_or(1.0)
        .abs()
        .max(1e-12);

    builder = builder.set_symmetric(0, 1, source_coupling)?;
    builder = builder.set_symmetric(order, order + 1, load_coupling)?;

    for resonator in 0..order {
        let diagonal = *projected_p.get(resonator).unwrap_or(&0.0);
        builder = builder.set(resonator + 1, resonator + 1, diagonal)?;
    }

    for step in 0..order.saturating_sub(1) {
        let e_coeff = projected_e
            .get(step + 1)
            .copied()
            .unwrap_or_default()
            .abs();
        let f_coeff = projected_f
            .get(step + 1)
            .copied()
            .unwrap_or_default()
            .abs();
        // Average the neighboring E/F coefficients to obtain a simple chain coupling.
        let coupling = ((e_coeff + f_coeff) / 2.0).max(1e-12);
        builder = builder.set_symmetric(step + 1, step + 2, coupling)?;
    }

    builder.build()
}
