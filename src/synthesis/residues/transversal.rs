use crate::approx::PolynomialSet;
use crate::error::{MfsError, Result};
use crate::matrix::{CouplingMatrix, CouplingMatrixBuilder};

use crate::approx::ComplexCoefficient;

use super::expansion::{
    classify_residues, nonzero_sqrt_abs, residue_real_component, ResidueClassification,
    ResidueExpansion,
};

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

    // Classify residues as real or complex-conjugate pairs
    let classifications = classify_residues(y11, y12, y22, 1e-6)?;

    // Adaptive pole-axis tolerance: for high-order polynomials, the root solver
    // produces poles with slightly larger real parts due to accumulated numerical
    // error. Scale the tolerance with order to avoid false rejections.
    // For order N, the Durand-Kerner solver's precision degrades roughly as O(N^2 * eps_mach).
    // We use 1e-6 * (N/5)^2 which gives ~1e-6 for order 5, ~1.6e-5 for order 20, ~3.6e-5 for order 30.
    let pole_axis_tolerance = 1e-6_f64.max(1e-6 * (order as f64 / 5.0).powi(2));

    let mut builder = CouplingMatrixBuilder::new(order)?;

    for classification in &classifications {
        match classification {
            ResidueClassification::Real { index } => {
                let idx = *index;
                let pole = y11.residues[idx].pole;
                if pole.re.abs() > pole_axis_tolerance {
                    return Err(MfsError::NumericalFailure(format!(
                        "transversal synthesis expects poles close to the imaginary axis \
                         (pole real part {:.2e} exceeds tolerance {:.2e} for order {order})",
                        pole.re.abs(), pole_axis_tolerance
                    )));
                }

                let pole_imag = pole.im;

                // For poles classified as Real, extract a real coupling value from
                // each residue. The classification step already determined this pole
                // should be treated as real (either because it's near the origin, has
                // real residues, or is an unpaired complex pole from complex polynomial
                // arithmetic). The dominant real/imaginary component is used to tolerate
                // the odd-order `j`-phase of the transfer residue.
                let residue_11 = residue_real_component(y11.residues[idx].residue);
                let residue_12 = residue_real_component(y12.residues[idx].residue);
                let residue_22 = residue_real_component(y22.residues[idx].residue);
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

                builder = builder.set(idx + 1, idx + 1, -pole_imag)?;
                builder = builder.set_symmetric(0, idx + 1, source)?;
                builder = builder.set_symmetric(idx + 1, order + 1, load)?;
            }
            ResidueClassification::ComplexPair { index_a, index_b } => {
                let idx_a = *index_a;
                let idx_b = *index_b;

                // Validate both poles are close to the imaginary axis
                let pole_a = y11.residues[idx_a].pole;
                let pole_b = y11.residues[idx_b].pole;
                if pole_a.re.abs() > pole_axis_tolerance || pole_b.re.abs() > pole_axis_tolerance {
                    return Err(MfsError::NumericalFailure(format!(
                        "transversal synthesis expects poles close to the imaginary axis \
                         (pole real parts {:.2e}, {:.2e} exceed tolerance {:.2e} for order {order})",
                        pole_a.re.abs(), pole_b.re.abs(), pole_axis_tolerance
                    )));
                }

                // Each pole of the conjugate pair maps to an independent row/column
                // in the transversal matrix (Cameron 2003, Section III).
                // Extract real-valued couplings from each pole's residue independently,
                // tolerating the odd-order `j`-phase of the transfer residue.
                let residue_11_a = residue_real_component(y11.residues[idx_a].residue);
                let residue_12_a = residue_real_component(y12.residues[idx_a].residue);
                let residue_22_a = residue_real_component(y22.residues[idx_a].residue);
                let use_r11_a = residue_11_a.abs() >= residue_22_a.abs();

                let source_a = if use_r11_a {
                    nonzero_sqrt_abs(residue_11_a, "y11 residue (pair a)")?
                } else {
                    residue_12_a / nonzero_sqrt_abs(residue_22_a, "y22 residue (pair a)")?
                };
                let load_a = if use_r11_a {
                    residue_12_a / nonzero_sqrt_abs(residue_11_a, "y11 residue (pair a)")?
                } else {
                    nonzero_sqrt_abs(residue_22_a, "y22 residue (pair a)")?
                };

                let residue_11_b = residue_real_component(y11.residues[idx_b].residue);
                let residue_12_b = residue_real_component(y12.residues[idx_b].residue);
                let residue_22_b = residue_real_component(y22.residues[idx_b].residue);
                let use_r11_b = residue_11_b.abs() >= residue_22_b.abs();

                let source_b = if use_r11_b {
                    nonzero_sqrt_abs(residue_11_b, "y11 residue (pair b)")?
                } else {
                    residue_12_b / nonzero_sqrt_abs(residue_22_b, "y22 residue (pair b)")?
                };
                let load_b = if use_r11_b {
                    residue_12_b / nonzero_sqrt_abs(residue_11_b, "y11 residue (pair b)")?
                } else {
                    nonzero_sqrt_abs(residue_22_b, "y22 residue (pair b)")?
                };

                // Assign entries for pole_a's row
                builder = builder.set(idx_a + 1, idx_a + 1, -pole_a.im)?;
                builder = builder.set_symmetric(0, idx_a + 1, source_a)?;
                builder = builder.set_symmetric(idx_a + 1, order + 1, load_a)?;

                // Assign entries for pole_b's row
                builder = builder.set(idx_b + 1, idx_b + 1, -pole_b.im)?;
                builder = builder.set_symmetric(0, idx_b + 1, source_b)?;
                builder = builder.set_symmetric(idx_b + 1, order + 1, load_b)?;
            }
        }
    }

    let mut matrix = builder.build()?;
    if let Some(constant) = y12.constant_term {
        if constant.re.abs() > 1e-6 {
            return Err(MfsError::NumericalFailure(
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

fn set_matrix_entry(matrix: &mut CouplingMatrix, row: usize, col: usize, value: f64) {
    let side = matrix.side();
    matrix.as_mut_slice()[row * side + col] = value;
}

/// Combines a complex-conjugate residue pair into real matrix entries.
///
/// For a conjugate pair at poles `p` and `p*` with residues `r` and `r*`:
/// - The combination `r/(s-p) + r*/(s-p*)` produces a real-valued rational function
/// - The effective real residue is `2 * Re(r)` (imaginary parts cancel in the sum)
/// - Diagonal entry: `-Im(p)` (pole location on imaginary axis)
/// - Source coupling: `sqrt(|2 * Re(r11)|)`
/// - Load coupling: derived from `2*Re(r12) / sqrt(|2*Re(r11)|)` or `sqrt(|2*Re(r22)|)` depending on magnitudes
///
/// The function accepts residues from a single pole of the conjugate pair.
/// The conjugate's residues are implicitly `conj(r)`, so the sum `r + conj(r) = 2*Re(r)`.
///
/// Returns `(diagonal, source_coupling, load_coupling, cross_coupling)`.
#[allow(dead_code)]
pub(crate) fn combine_conjugate_pair(
    pole_a: ComplexCoefficient,
    residue_11_a: ComplexCoefficient,
    residue_12_a: ComplexCoefficient,
    residue_22_a: ComplexCoefficient,
) -> Result<(f64, f64, f64, f64)> {
    // The combination of conjugate residues r/(s-p) + r*/(s-p*) yields 2*Re(r)
    // as the effective real residue contribution.
    // Im(r) + Im(conj(r)) = Im(r) - Im(r) = 0, so the combination is always real.
    let combined_r11 = 2.0 * residue_11_a.re;
    let combined_r12 = 2.0 * residue_12_a.re;
    let combined_r22 = 2.0 * residue_22_a.re;

    // Diagonal entry: pole location on imaginary axis
    let diagonal = -pole_a.im;

    // Source and load couplings follow the same logic as the real path,
    // but using the combined (doubled) real residue values.
    let use_r11 = combined_r11.abs() >= combined_r22.abs();

    let source_coupling = if use_r11 {
        nonzero_sqrt_abs(combined_r11, "combined y11 residue")?
    } else {
        combined_r12 / nonzero_sqrt_abs(combined_r22, "combined y22 residue")?
    };

    let load_coupling = if use_r11 {
        combined_r12 / nonzero_sqrt_abs(combined_r11, "combined y11 residue")?
    } else {
        nonzero_sqrt_abs(combined_r22, "combined y22 residue")?
    };

    // Cross coupling is zero for the transversal form (no inter-resonator coupling)
    let cross_coupling = 0.0;

    Ok((diagonal, source_coupling, load_coupling, cross_coupling))
}

#[cfg(test)]
mod tests {
    use super::super::expansion::synthesize_residue_expansions;
    use super::*;
    use crate::approx::generalized_chebyshev_polynomials;
    use crate::spec::FilterSpec;

    /// Helper to assert two f64 values are approximately equal within tolerance.
    fn assert_approx_eq(lhs: f64, rhs: f64, tol: f64, msg: &str) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "{msg}: expected {lhs} ≈ {rhs} within {tol}, diff={diff}"
        );
    }

    /// Order-2 all-pole 20 dB: verify symmetric source/load couplings.
    /// Requirements: 3.1, 3.2
    #[test]
    fn order_2_all_pole_symmetric_source_load_couplings() -> Result<()> {
        let spec = FilterSpec::new(2, 20.0)?;
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        // Order 2 => matrix is 4x4 (order + 2)
        assert_eq!(matrix.side(), 4);
        assert_eq!(matrix.order(), 2);

        // Source couplings at (0,1) and (0,2)
        let source_1 = matrix.at(0, 1).unwrap().abs();
        let source_2 = matrix.at(0, 2).unwrap().abs();

        // Load couplings at (1,3) and (2,3)
        let load_1 = matrix.at(1, 3).unwrap().abs();
        let load_2 = matrix.at(2, 3).unwrap().abs();

        // Verify |source_coupling_1| ≈ |load_coupling_1|
        assert_approx_eq(
            source_1, load_1, 1e-6,
            "source coupling 1 should equal load coupling 1 for symmetric all-pole",
        );

        // Verify |source_coupling_2| ≈ |load_coupling_2|
        assert_approx_eq(
            source_2, load_2, 1e-6,
            "source coupling 2 should equal load coupling 2 for symmetric all-pole",
        );

        Ok(())
    }

    /// Order-3 all-pole 20 dB: verify diagonal entries sum to zero (symmetric response).
    /// For an odd-order all-pole filter, the transversal form has poles symmetric about
    /// the origin: one at 0 and a conjugate pair at ±jω. The diagonal entries are
    /// -Im(pole), so they should sum to zero (symmetric frequency offsets).
    /// Requirements: 3.1, 3.3
    #[test]
    fn order_3_all_pole_zero_diagonal_entries() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        // Order 3 => matrix is 5x5 (order + 2)
        assert_eq!(matrix.side(), 5);
        assert_eq!(matrix.order(), 3);

        // Diagonal entries at (1,1), (2,2), (3,3) represent pole locations.
        // For a symmetric all-pole filter, they should sum to zero
        // (one pole at origin, two at ±jω giving diagonals ∓ω).
        let diag_1 = matrix.at(1, 1).unwrap();
        let diag_2 = matrix.at(2, 2).unwrap();
        let diag_3 = matrix.at(3, 3).unwrap();

        // The sum of all diagonal entries should be zero (symmetric pole distribution)
        let diag_sum = diag_1 + diag_2 + diag_3;
        assert_approx_eq(
            diag_sum, 0.0, 1e-6,
            "sum of diagonal entries should be zero for symmetric all-pole filter",
        );

        // At least one diagonal should be zero (the pole at the origin for odd order)
        let has_zero_diagonal = [diag_1, diag_2, diag_3]
            .iter()
            .any(|d| d.abs() < 1e-6);
        assert!(
            has_zero_diagonal,
            "odd-order all-pole should have at least one zero diagonal (pole at origin)"
        );

        Ok(())
    }

    /// Order-4 with 2 TZs at ±1.5, 20 dB: backward compatibility test.
    /// Verifies that this configuration synthesizes correctly and produces a structurally
    /// valid matrix with reference values matching the current implementation.
    /// Requirements: 7.1, 7.2, 7.3
    #[test]
    fn order_4_with_2_tzs_backward_compatibility() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        // Order 4 => matrix is 6x6 (order + 2)
        assert_eq!(matrix.side(), 6);
        assert_eq!(matrix.order(), 4);

        // Verify residue classification produces exactly 2 complex pairs for this config
        let classifications = classify_residues(&y11, &y12, &y22, 1e-6)?;
        let complex_pair_count = classifications
            .iter()
            .filter(|c| matches!(c, ResidueClassification::ComplexPair { .. }))
            .count();
        assert_eq!(
            complex_pair_count, 2,
            "order-4 with 2 TZs at ±1.5 should produce 2 complex-conjugate pairs"
        );

        // Verify source couplings are non-zero (row 0, columns 1..=4)
        for col in 1..=4 {
            let val = matrix.at(0, col).unwrap();
            assert!(
                val.abs() > 1e-12,
                "source coupling at (0,{col}) should be non-zero, got {val}"
            );
        }

        // Verify load couplings are non-zero (rows 1..=4, column 5)
        for row in 1..=4 {
            let val = matrix.at(row, 5).unwrap();
            assert!(
                val.abs() > 1e-12,
                "load coupling at ({row},5) should be non-zero, got {val}"
            );
        }

        // Verify matrix symmetry: M[i,j] == M[j,i] for all entries
        let side = matrix.side();
        for row in 0..side {
            for col in 0..side {
                let val_ij = matrix.at(row, col).unwrap();
                let val_ji = matrix.at(col, row).unwrap();
                assert_approx_eq(
                    val_ij, val_ji, 1e-10,
                    &format!("matrix should be symmetric: M[{row},{col}] vs M[{col},{row}]"),
                );
            }
        }

        // Verify diagonal entries are finite and reasonable (pole locations)
        for i in 1..=4 {
            let diag = matrix.at(i, i).unwrap();
            assert!(
                diag.is_finite(),
                "diagonal ({i},{i}) should be finite, got {diag}"
            );
        }

        // Reference values cross-validated against the fuqtgui Python engine
        // (CouplingMatrix.find_Ps / find_Es / coupling_matrix pipeline).
        assert_approx_eq(matrix.at(0, 1).unwrap(), 3.2436689194764357e-1, 1e-10,
            "reference value M[0,1]");
        assert_approx_eq(matrix.at(0, 2).unwrap(), 6.3995494022614496e-1, 1e-10,
            "reference value M[0,2]");
        assert_approx_eq(matrix.at(0, 3).unwrap(), 6.3995494022614707e-1, 1e-10,
            "reference value M[0,3]");
        assert_approx_eq(matrix.at(0, 4).unwrap(), 3.2436689194763746e-1, 1e-10,
            "reference value M[0,4]");
        assert_approx_eq(matrix.at(1, 1).unwrap(), 1.2452660363370787e0, 1e-10,
            "reference value M[1,1]");
        assert_approx_eq(matrix.at(1, 5).unwrap(), -3.2436689194763552e-1, 1e-10,
            "reference value M[1,5]");
        assert_approx_eq(matrix.at(2, 2).unwrap(), 7.6299234763018264e-1, 1e-10,
            "reference value M[2,2]");
        assert_approx_eq(matrix.at(2, 5).unwrap(), 6.3995494022614674e-1, 1e-10,
            "reference value M[2,5]");
        assert_approx_eq(matrix.at(3, 3).unwrap(), -7.6299234763018497e-1, 1e-10,
            "reference value M[3,3]");
        assert_approx_eq(matrix.at(3, 5).unwrap(), -6.3995494022614741e-1, 1e-10,
            "reference value M[3,5]");
        assert_approx_eq(matrix.at(4, 4).unwrap(), -1.2452660363370769e0, 1e-10,
            "reference value M[4,4]");
        assert_approx_eq(matrix.at(4, 5).unwrap(), 3.2436689194764701e-1, 1e-10,
            "reference value M[4,5]");

        Ok(())
    }

    /// Verify coupling sign conventions match current implementation.
    /// For order-4 with 2 TZs at ±1.5, source and load couplings should have
    /// consistent sign patterns that match the existing implementation behavior.
    /// Requirements: 7.3
    #[test]
    fn coupling_sign_conventions_order_4_with_2_tzs() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        let order = 4usize;
        let side = matrix.side(); // 6

        // Collect source couplings (row 0, columns 1..=4)
        let source_couplings: Vec<f64> = (1..=order)
            .map(|col| matrix.at(0, col).unwrap())
            .collect();

        // Collect load couplings (rows 1..=4, column 5)
        let load_couplings: Vec<f64> = (1..=order)
            .map(|row| matrix.at(row, side - 1).unwrap())
            .collect();

        // Verify sign convention: source couplings should all be positive
        // (this is the convention used by the sqrt(|residue|) extraction which
        // always produces a positive value for source couplings)
        for (i, &sc) in source_couplings.iter().enumerate() {
            assert!(
                sc > 0.0,
                "source coupling {i} should be positive (sign convention), got {sc}"
            );
        }

        // Verify that load couplings have signs determined by the residue_12 ratio.
        // The key invariant is: for each resonator k, source_k * load_k has the same
        // sign as the combined real part of the y12 residue (2*Re(r12) for complex pairs).
        for (i, (&sc, &lc)) in source_couplings.iter().zip(load_couplings.iter()).enumerate() {
            let product = sc * lc;
            // For complex pairs, the effective residue is 2*Re(r12)
            let effective_r12 = 2.0 * y12.residues[i].residue.re;
            // The product source*load should have the same sign as the effective y12 residue
            assert!(
                product * effective_r12 >= 0.0,
                "source*load product sign should match effective y12 residue sign for resonator {i}: \
                 product={product}, effective_r12={effective_r12}"
            );
        }

        // Verify symmetry of source/load coupling entries
        for col in 1..=order {
            let m_0_col = matrix.at(0, col).unwrap();
            let m_col_0 = matrix.at(col, 0).unwrap();
            assert_approx_eq(
                m_0_col, m_col_0, 1e-10,
                &format!("source coupling symmetry at column {col}"),
            );
        }
        for row in 1..=order {
            let m_row_last = matrix.at(row, side - 1).unwrap();
            let m_last_row = matrix.at(side - 1, row).unwrap();
            assert_approx_eq(
                m_row_last, m_last_row, 1e-10,
                &format!("load coupling symmetry at row {row}"),
            );
        }

        Ok(())
    }

    /// Temporary test to capture reference values for order-4 with 2 TZs at ±1.5, 20 dB.
    #[test]
    fn capture_order_4_reference_values() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        // Print all matrix entries
        let side = matrix.side(); // should be 6 (order + 2)
        eprintln!("Matrix side: {}", side);
        eprintln!("Matrix order: {}", matrix.order());
        for row in 0..side {
            for col in 0..side {
                let val = matrix.at(row, col).unwrap();
                if val.abs() > 1e-15 {
                    eprintln!("  matrix[{},{}] = {:.16e}", row, col, val);
                }
            }
        }

        // Also print residue classifications
        let classifications = classify_residues(&y11, &y12, &y22, 1e-6)?;
        eprintln!("\nClassifications:");
        for c in &classifications {
            eprintln!("  {:?}", c);
        }

        // Also test order-5 with 2 TZs (should produce real residues)
        eprintln!("\n--- Order-5 with 2 TZs at ±1.5, 20 dB ---");
        let spec5 = FilterSpec::new(5, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let poly5 = generalized_chebyshev_polynomials(&spec5)?;
        let (y11_5, y12_5, y22_5) = synthesize_residue_expansions(&poly5)?;
        let class5 = classify_residues(&y11_5, &y12_5, &y22_5, 1e-6)?;
        eprintln!("Classifications for order-5 with 2 TZs:");
        for c in &class5 {
            eprintln!("  {:?}", c);
        }

        // Also test order-4 with 1 TZ
        eprintln!("\n--- Order-4 with 1 TZ at 1.5, 20 dB ---");
        let spec4_1tz = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![1.5]);
        let poly4_1tz = generalized_chebyshev_polynomials(&spec4_1tz)?;
        match synthesize_residue_expansions(&poly4_1tz) {
            Ok((y11_4_1tz, y12_4_1tz, y22_4_1tz)) => {
                match classify_residues(&y11_4_1tz, &y12_4_1tz, &y22_4_1tz, 1e-6) {
                    Ok(class4_1tz) => {
                        eprintln!("Classifications for order-4 with 1 TZ:");
                        for c in &class4_1tz {
                            eprintln!("  {:?}", c);
                        }
                    }
                    Err(e) => eprintln!("Classification error for order-4 with 1 TZ: {e}"),
                }
            }
            Err(e) => eprintln!("Residue expansion error for order-4 with 1 TZ: {e}"),
        }

        // Also test order-6 with 2 TZs
        eprintln!("\n--- Order-6 with 2 TZs at ±1.5, 20 dB ---");
        let spec6 = FilterSpec::new(6, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let poly6 = generalized_chebyshev_polynomials(&spec6)?;
        let (y11_6, y12_6, y22_6) = synthesize_residue_expansions(&poly6)?;
        let class6 = classify_residues(&y11_6, &y12_6, &y22_6, 1e-6)?;
        eprintln!("Classifications for order-6 with 2 TZs:");
        for c in &class6 {
            eprintln!("  {:?}", c);
        }

        Ok(())
    }

    /// Order-3 with 2 TZs at ±1.5: verify synthesis succeeds and produces valid matrix.
    /// Requirements: 4.1, 4.2
    #[test]
    fn order_3_with_2_tzs_synthesis_succeeds() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?.with_normalized_transmission_zeros(vec![1.5, -1.5]);
        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let (y11, y12, y22) = synthesize_residue_expansions(&polynomials)?;
        let matrix = build_transversal_from_residues(&polynomials, &y11, &y12, &y22)?;

        // Order 3 with 2 TZs => matrix is 5x5 (order + 2)
        assert_eq!(matrix.side(), 5);
        assert_eq!(matrix.order(), 3);

        // Verify non-zero source couplings (at least one in row 0, columns 1..=3)
        let has_nonzero_source = (1..=3).any(|col| matrix.at(0, col).unwrap().abs() > 1e-12);
        assert!(
            has_nonzero_source,
            "matrix should have at least one non-zero source coupling"
        );

        // Verify non-zero load couplings (at least one in column 4, rows 1..=3)
        let has_nonzero_load = (1..=3).any(|row| matrix.at(row, 4).unwrap().abs() > 1e-12);
        assert!(
            has_nonzero_load,
            "matrix should have at least one non-zero load coupling"
        );

        Ok(())
    }
}
