use crate::approx::{ComplexCoefficient, ComplexPolynomial, PolynomialSet};
use crate::error::{MfsError, Result};
use crate::freq::BandPassMapping;

use super::{BandPassScaledCouplingMatrix, CouplingMatrix, CouplingMatrixBuilder, MatrixTopology};

/// Indicates which matrix-construction path produced the final matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixSynthesisMethod {
    /// Matrix was reconstructed from Y-parameter poles and residues.
    ResidueExpansion,
    /// Matrix fell back to the placeholder chain builder.
    PlaceholderFallback,
}

/// Detailed result of one coupling-matrix synthesis attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixSynthesisOutcome {
    /// Final synthesized matrix.
    pub matrix: CouplingMatrix,
    /// Construction path that produced the matrix.
    pub method: MatrixSynthesisMethod,
}

/// Polynomial-form Y parameters derived from the approximation output.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittancePolynomials {
    /// Common denominator polynomial used by all Y parameters.
    pub denominator: ComplexPolynomial,
    /// Numerator of normalized input admittance.
    pub y11: ComplexPolynomial,
    /// Numerator of transfer admittance.
    pub y12: ComplexPolynomial,
    /// Numerator of output admittance.
    pub y22: ComplexPolynomial,
}

/// One simple pole with its associated residue.
#[derive(Debug, Clone, PartialEq)]
pub struct ResiduePole {
    /// Pole location in the `s` plane.
    pub pole: ComplexCoefficient,
    /// Residue associated with that pole.
    pub residue: ComplexCoefficient,
}

/// Partial-fraction data extracted from one rational Y parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueExpansion {
    /// Simple-pole residues sorted by pole imaginary part.
    pub residues: Vec<ResiduePole>,
    /// Constant term left over after removing the simple-pole part.
    pub constant_term: Option<ComplexCoefficient>,
}

/// Maps polynomial data into the placeholder coupling-matrix representation.
#[derive(Debug, Default, Clone, Copy)]
pub struct CouplingMatrixSynthesizer;

impl CouplingMatrixSynthesizer {
    /// Builds polynomial-form Y parameters from generalized Chebyshev helper data.
    pub fn synthesize_admittance_polynomials(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<AdmittancePolynomials> {
        let generalized = polynomials.generalized.as_ref().ok_or_else(|| {
            MfsError::Unsupported(
                "admittance synthesis requires generalized Chebyshev helper data".to_string(),
            )
        })?;
        let e_s = generalized.e_s.as_ref().ok_or_else(|| {
            MfsError::Unsupported("generalized helper data is missing E(s)".to_string())
        })?;
        let f_over_eps_r = generalized
            .f_s
            .scale(ComplexCoefficient::from_real(1.0 / polynomials.eps_r))?;
        let p_transfer = generalized
            .p_s
            .scale(ComplexCoefficient::from_real(-2.0 / polynomials.eps))?;
        let conjugated_e = e_s.alternating_conjugate()?;
        let conjugated_f = f_over_eps_r.alternating_conjugate()?;

        let denominator = e_s
            .add(&f_over_eps_r)?
            .add(&conjugated_f)?
            .add(&conjugated_e)?;
        let y11 = e_s
            .sub(&f_over_eps_r)?
            .add(&conjugated_f)?
            .sub(&conjugated_e)?;
        let y22 = e_s
            .add(&f_over_eps_r)?
            .sub(&conjugated_f)?
            .sub(&conjugated_e)?;

        Ok(AdmittancePolynomials {
            denominator,
            y11,
            y12: p_transfer,
            y22,
        })
    }

    /// Splits the Y-parameter numerators into simple residues over the shared denominator.
    pub fn synthesize_residue_expansions(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<(ResidueExpansion, ResidueExpansion, ResidueExpansion)> {
        let admittance = self.synthesize_admittance_polynomials(polynomials)?;
        let y11 = residue_expansion(&admittance.y11, &admittance.denominator)?;
        let y12 = residue_expansion(&admittance.y12, &admittance.denominator)?;
        let y22 = residue_expansion(&admittance.y22, &admittance.denominator)?;
        Ok((y11, y12, y22))
    }

    /// Produces a coupling matrix and reports which synthesis path succeeded.
    pub fn synthesize_with_details(
        &self,
        polynomials: &PolynomialSet,
    ) -> Result<MatrixSynthesisOutcome> {
        if polynomials.generalized.is_some() {
            match self
                .synthesize_residue_expansions(polynomials)
                .and_then(|(y11, y12, y22)| {
                    build_transversal_from_residues(polynomials, &y11, &y12, &y22)
                }) {
                Ok(matrix) => {
                    return Ok(MatrixSynthesisOutcome {
                        matrix,
                        method: MatrixSynthesisMethod::ResidueExpansion,
                    })
                }
                Err(MfsError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }

        Ok(MatrixSynthesisOutcome {
            matrix: synthesize_placeholder_matrix(polynomials)?,
            method: MatrixSynthesisMethod::PlaceholderFallback,
        })
    }

    /// Produces a coupling matrix that is structurally compatible with the current pipeline.
    pub fn synthesize(&self, polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        Ok(self.synthesize_with_details(polynomials)?.matrix)
    }

    /// Synthesizes a matrix and immediately applies one topology transformation.
    pub fn synthesize_with_topology(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
    ) -> Result<CouplingMatrix> {
        self.synthesize(polynomials)?.transform_topology(topology)
    }

    /// Synthesizes a matrix, applies the requested topology, and scales it into band-pass units.
    pub fn synthesize_bandpass(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
        mapping: &BandPassMapping,
    ) -> Result<CouplingMatrix> {
        self.synthesize_with_topology(polynomials, topology)?
            .denormalize_bandpass(mapping)
    }

    /// Synthesizes a matrix, applies topology, and returns the band-pass/Qe representation.
    pub fn synthesize_bandpass_with_external_q(
        &self,
        polynomials: &PolynomialSet,
        topology: MatrixTopology,
        mapping: &BandPassMapping,
    ) -> Result<BandPassScaledCouplingMatrix> {
        self.synthesize_with_topology(polynomials, topology)?
            .denormalize_bandpass_with_external_q(mapping)
    }

    /// Synthesizes a matrix and extracts one triplet section at the requested center.
    pub fn synthesize_triplet(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        center_resonator: usize,
    ) -> Result<CouplingMatrix> {
        self.synthesize(polynomials)?
            .extract_triplet(transmission_zero, center_resonator)
    }

    /// Synthesizes a matrix and extracts a quadruplet section from two adjacent triplets.
    pub fn synthesize_quadruplet(
        &self,
        polynomials: &PolynomialSet,
        first_zero: f64,
        second_zero: f64,
        position: usize,
        common_resonator: usize,
        swap_zero_order: bool,
    ) -> Result<CouplingMatrix> {
        self.synthesize(polynomials)?.extract_quadruplet(
            first_zero,
            second_zero,
            position,
            common_resonator,
            swap_zero_order,
        )
    }

    /// Synthesizes a matrix and pulls one trisection into the requested resonator window.
    pub fn synthesize_trisection(
        &self,
        polynomials: &PolynomialSet,
        transmission_zero: f64,
        zero_positions: (usize, usize),
    ) -> Result<CouplingMatrix> {
        self.synthesize(polynomials)?
            .extract_trisection(transmission_zero, zero_positions)
    }
}

fn synthesize_placeholder_matrix(polynomials: &PolynomialSet) -> Result<CouplingMatrix> {
        let order = polynomials.order;
        let mut builder = CouplingMatrixBuilder::new(order)?;

        let source_coupling = polynomials
            .f
            .get(0)
            .copied()
            .unwrap_or(1.0)
            .abs()
            .max(1e-12);
        // Keep the source/load couplings non-zero so the response matrix remains invertible.
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
            // Average the neighboring E/F coefficients to obtain a simple chain coupling.
            let coupling = ((e_coeff + f_coeff) / 2.0).max(1e-12);
            builder = builder.set_symmetric(step + 1, step + 2, coupling)?;
        }

        builder.build()
}

fn residue_expansion(
    numerator: &ComplexPolynomial,
    denominator: &ComplexPolynomial,
) -> Result<ResidueExpansion> {
    if numerator.degree() > denominator.degree() {
        return Err(MfsError::Unsupported(
            "residue expansion currently requires a proper or constant-offset rational function"
                .to_string(),
        ));
    }

    let constant_term = if numerator.degree() == denominator.degree() {
        Some(numerator.leading_coefficient() / denominator.leading_coefficient())
    } else {
        None
    };

    let adjusted_numerator = if let Some(constant) = constant_term {
        numerator.sub(&denominator.scale(constant)?)?
    } else {
        numerator.clone()
    };

    let derivative = denominator.derivative()?;
    let mut residues = denominator
        .roots()?
        .into_iter()
        .map(|pole| {
            let slope = derivative.evaluate(pole);
            if slope.norm_sqr() <= 1e-20 {
                return Err(MfsError::Unsupported(
                    "repeated poles are not yet supported in residue expansion".to_string(),
                ));
            }

            Ok(ResiduePole {
                pole,
                residue: adjusted_numerator.evaluate(pole) / slope,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    residues.sort_by(|left, right| {
        left.pole
            .im
            .partial_cmp(&right.pole.im)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(ResidueExpansion {
        residues,
        constant_term,
    })
}

fn build_transversal_from_residues(
    polynomials: &PolynomialSet,
    y11: &ResidueExpansion,
    y12: &ResidueExpansion,
    y22: &ResidueExpansion,
) -> Result<CouplingMatrix> {
    let order = polynomials.order;
    if y11.residues.len() != order || y12.residues.len() != order || y22.residues.len() != order {
        return Err(MfsError::DimensionMismatch {
            expected: order,
            actual: y11
                .residues
                .len()
                .min(y12.residues.len())
                .min(y22.residues.len()),
        });
    }

    let mut builder = CouplingMatrixBuilder::new(order)?;
    for index in 0..order {
        let pole = y11.residues[index].pole;
        if pole.re.abs() > 1e-6 {
            return Err(MfsError::Unsupported(
                "transversal synthesis currently expects poles close to the imaginary axis"
                    .to_string(),
            ));
        }

        let pole_imag = pole.im;
        let residue_11 = real_part_if_almost_real(y11.residues[index].residue, "y11 residue")?;
        let residue_12 = real_part_if_almost_real(y12.residues[index].residue, "y12 residue")?;
        let residue_22 = real_part_if_almost_real(y22.residues[index].residue, "y22 residue")?;
        let use_r11 = residue_11.abs() >= residue_22.abs();

        let source = if use_r11 {
            nonzero_sqrt_abs(residue_11, "y11 residue")?
        } else {
            residue_12 / nonzero_sqrt_abs(residue_22, "y22 residue")?
        };
        let load = if use_r11 {
            residue_12 / nonzero_sqrt_abs(residue_11, "y11 residue")?
        } else {
            nonzero_sqrt_abs(residue_22, "y22 residue")?
        };

        builder = builder.set(index + 1, index + 1, -pole_imag)?;
        builder = builder.set_symmetric(0, index + 1, source)?;
        builder = builder.set_symmetric(index + 1, order + 1, load)?;
    }

    let mut matrix = builder.build()?;
    if let Some(constant) = y12.constant_term {
        if constant.re.abs() > 1e-6 {
            return Err(MfsError::Unsupported(
                "direct source-load term must be purely imaginary in current synthesis path"
                    .to_string(),
            ));
        }
        let direct = constant.im;
        if direct.abs() > 1e-12 {
            let side = matrix.side();
            set_matrix_entry(&mut matrix, 0, side - 1, direct);
            set_matrix_entry(&mut matrix, side - 1, 0, direct);
        }
    } else if polynomials.transmission_zeros_normalized.len() == order {
        let direct = polynomials.eps * (polynomials.eps_r - 1.0) / polynomials.eps_r;
        if direct.abs() > 1e-12 {
            let side = matrix.side();
            set_matrix_entry(&mut matrix, 0, side - 1, direct);
            set_matrix_entry(&mut matrix, side - 1, 0, direct);
        }
    }

    Ok(matrix)
}

fn nonzero_sqrt_abs(value: f64, label: &str) -> Result<f64> {
    if value.abs() <= 1e-12 {
        return Err(MfsError::Unsupported(format!(
            "{label} is too small to derive a stable coupling value"
        )));
    }

    Ok(value.abs().sqrt())
}

fn real_part_if_almost_real(value: ComplexCoefficient, label: &str) -> Result<f64> {
    if value.im.abs() > 1e-6 {
        return Err(MfsError::Unsupported(format!(
            "{label} is unexpectedly complex in the current real-valued synthesis path"
        )));
    }

    Ok(value.re)
}

fn set_matrix_entry(matrix: &mut CouplingMatrix, row: usize, col: usize, value: f64) {
    let side = matrix.side();
    matrix.as_mut_slice()[row * side + col] = value;
}
