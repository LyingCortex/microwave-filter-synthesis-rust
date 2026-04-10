use crate::approx::{
    ComplexCoefficient, ComplexPolynomial, DurandKernerRootSolver, PolynomialSet,
};
use crate::error::{MfsError, Result};
use crate::matrix::{CouplingMatrix, CouplingMatrixBuilder};

fn complex_from_real(value: f64) -> ComplexCoefficient {
    ComplexCoefficient::new(value, 0.0)
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

/// Builds polynomial-form Y parameters from generalized Chebyshev helper data.
pub(crate) fn synthesize_admittance_polynomials(
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
        .scale(complex_from_real(1.0 / polynomials.eps_r))?;
    let p_transfer = generalized
        .p_s
        .scale(complex_from_real(-2.0 / polynomials.eps))?;
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
pub(crate) fn synthesize_residue_expansions(
    polynomials: &PolynomialSet,
) -> Result<(ResidueExpansion, ResidueExpansion, ResidueExpansion)> {
    let admittance = synthesize_admittance_polynomials(polynomials)?;
    let y11 = residue_expansion(&admittance.y11, &admittance.denominator)?;
    let y12 = residue_expansion(&admittance.y12, &admittance.denominator)?;
    let y22 = residue_expansion(&admittance.y22, &admittance.denominator)?;
    Ok((y11, y12, y22))
}

pub(crate) fn build_transversal_from_residues(
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

fn residue_expansion(
    numerator: &ComplexPolynomial,
    denominator: &ComplexPolynomial,
) -> Result<ResidueExpansion> {
    let solver = DurandKernerRootSolver;
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
        .roots_with(&solver)?
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
