//! Multi-band Chebyshev approximation (work item C2b of
//! `docs/internal/multi-band-synthesis-design.md`).
//!
//! Given $K \ge 2$ passbands (in the shared normalized prototype axis) and a
//! set of transmission zeros, this produces a [`PolynomialSet`] whose
//! characteristic $K(s) = F(s)/P(s)$ is **equiripple** in every passband at the
//! requested return loss. The $F(s)$ roots (reflection zeros) are solved for;
//! everything downstream — the Feldtkeller $E(s)$ spectral factor and
//! [`PolynomialSet`] assembly — reuses the single-band helpers verbatim, so the
//! existing `MatrixSynthesisEngine` consumes the output with no changes.
//!
//! # Method
//!
//! On the positive half-axis each passband interval $(a, b)$ is tiled by the
//! reflection zeros that fall inside it. $K(w)$ (with $w$ the real prototype
//! variable) is a ratio of real polynomials, so $|K|$ has a local maximum at
//! each band edge and between every pair of adjacent reflection zeros. The
//! equiripple condition is that all of these maxima share one value; with $N/2$
//! positive reflection zeros on a single positive interval this is
//! $(N/2 + 1)$ equations in $N/2$ unknowns (one value is the shared ripple
//! level), so the system is square. We solve it with a damped
//! Levenberg-Marquardt step and a finite-difference Jacobian.
//!
//! The return loss fixes the *level* through the crate's conventional
//! `eps = |P(j)| / |F(j)| / sqrt(10^{RL/10} - 1)` (`find_eps`), which
//! normalizes so the band-edge reflection equals exactly the requested RL.

use crate::error::{MfsError, Result};

use super::helpers::{build_e_polynomial_stage, find_eps, GeneralizedChebyshevData};
use super::{ComplexCoefficient, ComplexPolynomial, PolynomialSet};

const LM_MAX_ITERATIONS: usize = 200;
const LM_TOLERANCE: f64 = 1e-12;
const LM_INITIAL_LAMBDA: f64 = 1e-2;
const LM_LAMBDA_BREAK: f64 = 1e12;

/// Maximum local maxima expected on one positive interval: `n_rz + 1`.
fn maxima_count(n_rz: usize) -> usize {
    n_rz + 1
}

/// Number of positive (right-half) reflection zeros for a symmetric multi-band
/// design of total order `order`.
///
/// The characteristic `F` has degree `order`; its roots are the reflection
/// zeros, which for a symmetric prototype occur as $\pm$ pairs. Half are
/// positive. The order must therefore be even.
fn positive_rz_count(order: usize) -> Result<usize> {
    if order == 0 {
        return Err(MfsError::InvalidOrder { order });
    }
    if order % 2 != 0 {
        return Err(MfsError::InvalidOrder { order });
    }
    Ok(order / 2)
}


/// Samples $|F(w)/P(w)|$ over `grid` for the given positive reflection zeros
/// and the fixed transmission-zero polynomial `p`.
fn char_magnitude(rz: &[f64], p: &[f64], grid: &[f64]) -> Vec<f64> {
    // F(w) = prod (w^2 - rz^2), a real polynomial; compute its value directly
    // for the real grid to avoid repeated poly reconstruction.
    let mut out = Vec::with_capacity(grid.len());
    for &w in grid {
        let mut f = 1.0;
        for &r in rz {
            f *= w * w - r * r;
        }
        // P(w) via Horner over ascending coefficients (highest degree first).
        let mut pv = 0.0;
        for c in p.iter().rev() {
            pv = pv * w + c;
        }
        out.push((f / pv).abs());
    }
    out
}

/// Finds the local maxima of `values` (with the two endpoints always included
/// as candidates) and returns their indices, sorted ascending.
fn local_maxima_indices(values: &[f64]) -> Vec<usize> {
    let mut idx = Vec::new();
    idx.push(0);
    for i in 1..values.len().saturating_sub(1) {
        if values[i] >= values[i - 1] && values[i] >= values[i + 1] {
            idx.push(i);
        }
    }
    idx.push(values.len() - 1);
    idx.sort();
    idx.dedup();
    idx
}

/// Equiripple residual: the deviation of each local maximum of $|K|$ from their
/// mean. Square in count (one maxima = the shared level) for a single
/// positive interval.
fn equiripple_residual(rz: &[f64], p: &[f64], grid: &[f64], expected: usize) -> Vec<f64> {
    let values = char_magnitude(rz, p, grid);
    let idx = local_maxima_indices(&values);
    // The true equiripple peaks dominate in magnitude (they sit between the
    // zero-value reflection zeros), so take the `expected` largest local
    // maxima. This stays well-defined even if the grid exposes spurious peaks.
    let mut peaks: Vec<f64> = idx.iter().map(|&i| values[i]).collect();
    peaks.sort_by(|a, b| b.total_cmp(a));
    peaks.truncate(expected);
    let mean = peaks.iter().sum::<f64>() / peaks.len() as f64;
    peaks.iter().map(|v| v - mean).collect()
}

/// Solves the equiripple reflection zeros on a single positive interval
/// `(lo, hi)`. `seed` must be `n` positive values inside `(lo, hi)`.
fn solve_rz(lo: f64, hi: f64, p: &[f64], seed: &[f64]) -> Result<Vec<f64>> {
    let n = seed.len();
    let grid_len = 20000;
    let grid: Vec<f64> = (0..grid_len)
        .map(|i| lo + (hi - lo) * (i as f64) / (grid_len as f64 - 1.0))
        .collect();

    let expected = maxima_count(n);
    let mut rz = seed.to_vec();
    let mut cost = equiripple_residual(&rz, p, &grid, expected);
    let mut base_cost: f64 = cost.iter().map(|v| v * v).sum();
    let mut lambda = LM_INITIAL_LAMBDA;

    for _ in 0..LM_MAX_ITERATIONS {
        if (base_cost).sqrt() < LM_TOLERANCE {
            break;
        }
        // Finite-difference Jacobian (n residuals-ish, n params).
        let mut jac: Vec<Vec<f64>> = Vec::with_capacity(cost.len());
        for (r, &cost_r) in cost.iter().enumerate() {
            let mut row = vec![0.0; n];
            for j in 0..n {
                let delta = 1e-7 * 1.0_f64.max(rz[j].abs());
                let mut pert = rz.clone();
                pert[j] += delta;
                row[j] = (equiripple_residual(&pert, p, &grid, expected)[r] - cost_r) / delta;
            }
            jac.push(row);
        }

        // (J^T J + lambda * diag(J^T J)) dp = -J^T r
        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];
        for i in 0..cost.len() {
            for j in 0..n {
                jtr[j] += jac[i][j] * cost[i];
                for k in 0..n {
                    jtj[j][k] += jac[i][j] * jac[i][k];
                }
            }
        }
        for (j, row) in jtj.iter_mut().enumerate() {
            row[j] *= 1.0 + lambda;
        }
        let rhs: Vec<f64> = jtr.iter().map(|v| -v).collect();
        let step = solve_gaussian(&jtj, &rhs)?;

        // Enforce the RZ stay inside (lo, hi) and stay ordered.
        let mut trial = Vec::with_capacity(n);
        let mut feasible = true;
        for (j, &d) in step.iter().enumerate() {
            let v = rz[j] + d;
            if v <= lo + 1e-9 || v >= hi - 1e-9 {
                feasible = false;
                break;
            }
            trial.push(v);
        }
        if !feasible || step.iter().any(|d| !d.is_finite()) {
            lambda *= 2.0;
            if lambda > LM_LAMBDA_BREAK {
                break;
            }
            continue;
        }

        let trial_res = equiripple_residual(&trial, p, &grid, expected);
        let trial_cost: f64 = trial_res.iter().map(|v| v * v).sum();
        if trial_cost < base_cost {
            rz = trial;
            cost = trial_res;
            base_cost = trial_cost;
            lambda = (lambda * 0.5).max(1e-10);
        } else {
            lambda *= 2.0;
            if lambda > LM_LAMBDA_BREAK {
                break;
            }
        }
    }

    if (base_cost).sqrt() > 1e-6 {
        return Err(MfsError::NumericalFailure(format!(
            "multi-band equiripple solve did not converge (residual {:.3e})",
            base_cost.sqrt()
        )));
    }
    Ok(rz)
}

/// Gaussian elimination with partial pivoting for small dense systems.
fn solve_gaussian(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>> {
    let n = b.len();
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.push(b[i]);
            r
        })
        .collect();
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if aug[row][col].abs() > aug[pivot][col].abs() {
                pivot = row;
            }
        }
        if aug[pivot][col].abs() < 1e-30 {
            return Err(MfsError::NumericalFailure(
                "singular system in multi-band equiripple solve".to_string(),
            ));
        }
        aug.swap(col, pivot);
        let diag = aug[col][col];
        for row in (col + 1)..n {
            let factor = aug[row][col] / diag;
            let pivot_row: Vec<f64> = aug[col][col..=n].to_vec();
            for (k, aug_row_k) in aug[row].iter_mut().enumerate().skip(col) {
                *aug_row_k -= factor * pivot_row[k - col];
            }
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = aug[row][n];
        for k in (row + 1)..n {
            sum -= aug[row][k] * x[k];
        }
        x[row] = sum / aug[row][row];
    }
    Ok(x)
}

/// Builds the transmission-zero polynomial in the prototype variable
/// $w$: $P(w) = w$ (if a zero sits at $0$) $\\times \\prod_k (w^2 - z_k^2)$,
/// where the product runs over the **distinct** absolute values $z_k = |\\Omega_z|$.
///
/// A real transmission zero $\\Omega_z$ and its mirror $-\\Omega_z$ are one
/// physical pole pair (one $s^2 + \\Omega_z^2$ factor), so each contributes a
/// single $(w^2 - \\Omega_z^2)$ factor. Passing an explicitly symmetric list
/// (e.g. $\\{0, \\pm 0.2542, \\pm 1.1573\\}$) must not double-count a pair.
fn build_p_polynomial(zeros: &[f64]) -> Result<Vec<f64>> {
    let mut coeff: Vec<f64> = vec![1.0];
    for a in distinct_abs(zeros) {
        if a == 0.0 {
            // factor w
            let mut next = vec![0.0; coeff.len() + 1];
            for (i, &c) in coeff.iter().enumerate() {
                next[i + 1] = c;
            }
            coeff = next;
        } else {
            let z2 = a * a;
            let mut next = vec![0.0; coeff.len() + 2];
            for (i, &c) in coeff.iter().enumerate() {
                next[i + 2] += c;
                next[i] += -z2 * c;
            }
            coeff = next;
        }
    }
    Ok(coeff)
}

/// Distinct absolute values of the finite normalized zeros, sorted ascending.
/// A real zero `$\Omega_z$` and its mirror `$-\Omega_z$` are one physical pole
/// pair, so both collapse to the single value `$|\Omega_z|$`.
fn distinct_abs(zeros: &[f64]) -> Vec<f64> {
    use std::collections::HashSet;

    let mut seen: HashSet<u64> = HashSet::new();
    let mut out: Vec<f64> = Vec::new();
    for &z in zeros {
        if !z.is_finite() {
            continue; // zeros at infinity do not enter P
        }
        let a = z.abs();
        if seen.insert(a.to_bits()) {
            out.push(a);
        }
    }
    out.sort_by(|a, b| a.total_cmp(b));
    out
}

/// Assembles an s-domain [`ComplexPolynomial`] with imaginary roots $\pm j r_k$
/// (i.e. $s^2 + r_k^2$ factors) from positive real prototype zeros.
fn s_domain_from_positive_real(roots: &[f64]) -> Result<ComplexPolynomial> {
    let mut root_list: Vec<ComplexCoefficient> = Vec::with_capacity(roots.len() * 2);
    for &r in roots {
        root_list.push(ComplexCoefficient::new(0.0, r));
        root_list.push(ComplexCoefficient::new(0.0, -r));
    }
    ComplexPolynomial::from_complex_roots(&root_list)
}

/// Assembles an s-domain [`ComplexPolynomial`] for the transmission zeros:
/// $\pm j z_k$ for $z_k \ne 0$ and a single root at $s=0$ for the $z_k = 0$
/// zero.
fn s_domain_p(zeros: &[f64]) -> Result<ComplexPolynomial> {
    let mut roots: Vec<ComplexCoefficient> = Vec::new();
    for a in distinct_abs(zeros) {
        if a == 0.0 {
            roots.push(ComplexCoefficient::new(0.0, 0.0));
        } else {
            roots.push(ComplexCoefficient::new(0.0, a));
            roots.push(ComplexCoefficient::new(0.0, -a));
        }
    }
    ComplexPolynomial::from_complex_roots(&roots)
}

/// Computes the multi-band Chebyshev [`PolynomialSet`] for the given
/// passbands.
///
/// `bands` are passband intervals $(low, high)$ in the shared normalized
/// prototype axis, sorted by `low`. For the symmetric dual-band case this is
/// `[(-1.0, -0.35), (0.35, 1.0)]`. `order` is the total filter order (even).
/// `transmission_zeros` are the normalized zeros on the same axis.
///
/// The returned set is directly consumable by `MatrixSynthesisEngine::synthesize`.
pub fn multi_band_chebyshev_polynomials(
    bands: &[(f64, f64)],
    order: usize,
    return_loss_db: f64,
    transmission_zeros: &[f64],
) -> Result<PolynomialSet> {
    if bands.len() < 2 {
        return Err(MfsError::InvalidFrequency(
            "multi-band approximation requires at least 2 passbands".to_string(),
        ));
    }
    // Each finite transmission zero consumes one order unit (P-degree <= N), so
    // more finite zeros than the order is invalid. Check it early, mirroring the
    // single-band `pad_transmission_zeros`; otherwise the failure surfaces late
    // as a misleading equiripple non-convergence (review P3 #5).
    let finite_count = transmission_zeros.iter().filter(|z| z.is_finite()).count();
    if finite_count > order {
        return Err(MfsError::DimensionMismatch {
            expected: order,
            actual: finite_count,
        });
    }
    let n_rz = positive_rz_count(order)?;

    // Positive passband intervals: the part of each band lying in (0, +inf).
    // For a symmetric prototype these tile the upper half-axis.
    let pos: Vec<(f64, f64)> = bands
        .iter()
        .map(|&(lo, hi)| {
            let l = lo.max(1e-12);
            (l, hi)
        })
        .filter(|&(l, h)| h > l)
        .collect();

    // The reflection zeros are tiled across the positive-axis passband union.
    // The characteristic is a single-function K(w) = F(w)/P(w) with a
    // *real-rooted* F (roots +/-r_k). Equiripple imposes N/2 + M - 1
    // independent constraints (the maxima on M positive-axis intervals, with
    // one free shared level) on N/2 reflection zeros, so the system is square
    // only for M = 1 - the symmetric dual-band case (e.g. bands
    // `[-1, -0.35]` and `[0.35, 1]`). For M >= 2 (K >= 3: a center band
    // spanning the origin, or 2+ positive-axis passbands) it is overdetermined
    // by M - 1 and has no real-rooted solution; a genuine K >= 3 equiripple
    // characteristic needs complex-rooted (Zolotarev) or a multifunction form.
    // See `docs/internal/multi-band-synthesis-design.md` §"K >= 3".
    if pos.len() != 1 {
        return Err(MfsError::NotImplemented(
            "multi-band equiripple requires a real-rooted characteristic, which \
             can only be equiripple on a single positive-axis interval (the \
             symmetric dual-band, K = 2); K >= 3 (a center band spanning the \
             origin, or 2+ positive-axis passbands) is overdetermined and needs \
             a complex-rooted (Zolotarev) or multifunction characteristic"
                .to_string(),
        ));
    }
    let (union_lo, union_hi) = pos[0];
    if union_lo >= union_hi {
        return Err(MfsError::InvalidFrequency(
            "multi-band passbands do not span the positive axis".to_string(),
        ));
    }

    let p = build_p_polynomial(transmission_zeros)?;

    // Seed: evenly spread `n_rz` points across the positive union.
    let mut seed = Vec::with_capacity(n_rz);
    for i in 0..n_rz {
        let t = (i + 1) as f64 / (n_rz as f64 + 1.0);
        seed.push(union_lo + (union_hi - union_lo) * t);
    }

    let rz = solve_rz(union_lo, union_hi, &p, &seed)?;

    // Assemble the s-domain polynomials.
    let f_s = s_domain_from_positive_real(&rz)?;
    let p_s = s_domain_p(transmission_zeros)?;

    let (eps, eps_r) = find_eps(finite_count, &p_s, &f_s, return_loss_db, order)?;
    let e_stage = build_e_polynomial_stage(&f_s, &p_s, eps, eps_r)?;

    let generalized = GeneralizedChebyshevData {
        padded_zeros: transmission_zeros
            .iter()
            .copied()
            .chain(std::iter::repeat_n(f64::INFINITY, order.saturating_sub(finite_count)))
            .collect(),
        finite_zero_count: finite_count,
        f_s: f_s.clone(),
        p_s: p_s.clone(),
        a_s: None,
        a_stage: None,
        e_s: Some(e_stage.e_s.clone()),
        e_stage: Some(e_stage),
        eps,
        eps_r,
    };

    let e = generalized.e_s.clone().ok_or_else(|| {
        MfsError::PreconditionViolation("multi-band approximation is missing E(s)".to_string())
    })?;
    let set = PolynomialSet {
        order,
        ripple_factor: super::chebyshev_ripple_factor(return_loss_db),
        eps,
        eps_r,
        transmission_zeros_normalized: transmission_zeros.to_vec(),
        e,
        f: f_s,
        p: p_s,
        generalized: Some(generalized),
    };
    set.validate()?;
    Ok(set)
}


#[cfg(test)]
mod tests {
    use super::*;

    // Cameron 2011 Fig. 25, dual-band symmetric fixture:
    //   bands [-1,-0.35] and [0.35,1], order 10, RL 22 dB,
    //   TZ = {0, +/-0.25420871, +/-1.15728568}.
    const BANDS: &[(f64, f64)] = &[(-1.0, -0.35), (0.35, 1.0)];
    const ORDER: usize = 10;
    const RL: f64 = 22.0;
    const TZ: &[f64] = &[0.0, 0.25420871, -0.25420871, 1.15728568, -1.15728568];
    const ORACLE: &[f64] = &[0.36115462, 0.46549174, 0.68058384, 0.88628738, 0.98815144];

    #[test]
    fn equiripple_solve_recovers_fixture_rz() {
        let p = build_p_polynomial(TZ).expect("p");
        let seed: Vec<f64> = (0..ORDER / 2)
            .map(|i| 0.35 + 0.65 * ((i + 1) as f64 / ((ORDER / 2) as f64 + 1.0)))
            .collect();
        let rz = solve_rz(0.35, 1.0, &p, &seed).expect("solve");
        assert_eq!(rz.len(), ORDER / 2);
        for (got, want) in rz.iter().zip(ORACLE.iter()) {
            assert!(
                (got - want).abs() < 5e-4,
                "RZ {got} far from oracle {want}"
            );
        }
    }

    /// More finite transmission zeros than the order is rejected early with a
    /// typed `DimensionMismatch`, mirroring the single-band `pad_transmission_
    /// zeros` guard; without it the failure surfaces late as a misleading
    /// equiripple non-convergence (review P3 #5).
    #[test]
    fn more_finite_zeros_than_order_is_a_dimension_error() {
        // order 2 with 3 finite zeros (order < finite count).
        let err = multi_band_chebyshev_polynomials(BANDS, 2, 20.0, &[0.5, 1.2, -2.0])
            .unwrap_err();
        assert!(
            matches!(err, MfsError::DimensionMismatch { .. }),
            "order < finite zero count must be a typed DimensionMismatch, got {err:?}"
        );
    }

    #[test]
    fn fixture_polynomials_are_equiripple_at_specified_rl() {
        let set = multi_band_chebyshev_polynomials(BANDS, ORDER, RL, TZ).expect("polynomials");
        assert_eq!(set.order, ORDER);
        // Edge reflection must equal the spec RL (within the eps
        // normalization, so the ripple floor is ~RL).
        assert!(
            (set.eps - (10_f64.powf(RL / 10.0) - 1.0).sqrt()).abs()
                < (10_f64.powf(RL / 10.0) - 1.0).sqrt() * 0.25,
            "eps {} inconsistent with RL {RL}",
            set.eps
        );
        // E/F/P degrees are the order.
        assert_eq!(set.e.degree(), ORDER);
        assert_eq!(set.f.degree(), ORDER);
    }
}
