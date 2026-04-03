use crate::approx::PolynomialSet;
use crate::error::Result;

use super::{CouplingMatrix, CouplingMatrixBuilder};

#[derive(Debug, Default, Clone, Copy)]
pub struct CouplingMatrixSynthesizer;

impl CouplingMatrixSynthesizer {
    pub fn synthesize(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        let order = polynomials.order;
        let mut builder = CouplingMatrixBuilder::new(order)?;

        let source_coupling = polynomials
            .f
            .get(0)
            .copied()
            .unwrap_or(1.0)
            .abs()
            .max(1e-12);
        let load_coupling = polynomials
            .e
            .get(0)
            .copied()
            .unwrap_or(1.0)
            .abs()
            .max(1e-12);

        builder = builder.set_symmetric(0, 1, source_coupling)?;
        builder = builder.set_symmetric(order, order + 1, load_coupling)?;

        for resonator in 0..order {
            let diagonal = *polynomials.p.get(resonator).unwrap_or(&0.0);
            builder = builder.set(resonator + 1, resonator + 1, diagonal)?;
        }

        for step in 0..order.saturating_sub(1) {
            let e_coeff = polynomials
                .e
                .get(step + 1)
                .copied()
                .unwrap_or_default()
                .abs();
            let f_coeff = polynomials
                .f
                .get(step + 1)
                .copied()
                .unwrap_or_default()
                .abs();
            let coupling = ((e_coeff + f_coeff) / 2.0).max(1e-12);
            builder = builder.set_symmetric(step + 1, step + 2, coupling)?;
        }

        builder.build()
    }
}
