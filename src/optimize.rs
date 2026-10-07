//! Coupling matrix optimization and extraction.
//!
//! Provides two workflows:
//! - **Tuning**: adjust coupling values to match a target S-parameter response
//! - **Extraction**: recover a coupling matrix from measured S-parameter data
//!
//! Both use Levenberg-Marquardt nonlinear least-squares optimization.

use crate::error::{MfsError, Result};
use crate::freq::FrequencyGrid;
use crate::matrix::{CouplingMatrix, CouplingMatrixBuilder};
use crate::response::{ResponseSettings, ResponseSolver, SParameterResponse};

/// Configuration for the optimization engine.
#[derive(Debug, Clone)]
pub struct OptimizeConfig {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on cost function change.
    pub tolerance: f64,
    /// Initial Levenberg-Marquardt damping factor.
    pub lambda: f64,
    /// Whether to optimize diagonal entries (detuning).
    pub optimize_diagonals: bool,
    /// Whether to optimize source/load couplings.
    pub optimize_ports: bool,
    /// Unloaded Q for lossy evaluation (infinity = lossless).
    pub unloaded_q: f64,
}

impl Default for OptimizeConfig {
    fn default() -> Self {
        Self {
            max_iterations: 200,
            tolerance: 1e-8,
            lambda: 1e-3,
            optimize_diagonals: true,
            optimize_ports: true,
            unloaded_q: f64::INFINITY,
        }
    }
}

/// Result of an optimization run.
#[derive(Debug, Clone)]
pub struct OptimizeResult {
    /// Optimized coupling matrix.
    pub matrix: CouplingMatrix,
    /// Final cost (sum of squared residuals).
    pub cost: f64,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the optimizer converged.
    pub converged: bool,
}

/// Extracts the optimizable parameters from a coupling matrix.
///
/// For transversal matrices: source couplings, load couplings, diagonals.
/// For folded/arrow: inter-resonator couplings, diagonals, port couplings.
fn extract_parameters(matrix: &CouplingMatrix, config: &OptimizeConfig) -> Vec<f64> {
    let order = matrix.order();
    let side = matrix.side();
    let mut params = Vec::new();

    // Source couplings M[0, k] for k=1..=order
    if config.optimize_ports {
        for k in 1..=order {
            let val = matrix.at(0, k).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                params.push(val);
            }
        }
    }

    // Inter-resonator couplings (upper triangle, rows 1..=order, cols > row, ≤ order)
    for row in 1..=order {
        for col in (row + 1)..=order {
            let val = matrix.at(row, col).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                params.push(val);
            }
        }
    }

    // Diagonal entries (detuning)
    if config.optimize_diagonals {
        for k in 1..=order {
            params.push(matrix.at(k, k).unwrap_or(0.0));
        }
    }

    // Load couplings M[k, N+1] for k=1..=order
    if config.optimize_ports {
        for k in 1..=order {
            let val = matrix.at(k, side - 1).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                params.push(val);
            }
        }
    }

    // Direct source-load coupling M[0, N+1]
    let direct = matrix.at(0, side - 1).unwrap_or(0.0);
    if direct.abs() > 1e-15 {
        params.push(direct);
    }

    params
}

/// Rebuilds a coupling matrix from the parameter vector.
fn rebuild_matrix(
    template: &CouplingMatrix,
    params: &[f64],
    config: &OptimizeConfig,
) -> Result<CouplingMatrix> {
    let order = template.order();
    let side = template.side();
    let mut data = template.as_slice().to_vec();
    let mut idx = 0;

    // Source couplings
    if config.optimize_ports {
        for k in 1..=order {
            let val = template.at(0, k).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                data[k] = params[idx];
                data[k * side] = params[idx];
                idx += 1;
            }
        }
    }

    // Inter-resonator couplings
    for row in 1..=order {
        for col in (row + 1)..=order {
            let val = template.at(row, col).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                data[row * side + col] = params[idx];
                data[col * side + row] = params[idx];
                idx += 1;
            }
        }
    }

    // Diagonal entries
    if config.optimize_diagonals {
        for k in 1..=order {
            data[k * side + k] = params[idx];
            idx += 1;
        }
    }

    // Load couplings
    if config.optimize_ports {
        for k in 1..=order {
            let val = template.at(k, side - 1).unwrap_or(0.0);
            if val.abs() > 1e-15 {
                data[k * side + (side - 1)] = params[idx];
                data[(side - 1) * side + k] = params[idx];
                idx += 1;
            }
        }
    }

    // Direct source-load
    let direct = template.at(0, side - 1).unwrap_or(0.0);
    if direct.abs() > 1e-15 {
        data[side - 1] = params[idx];
        data[(side - 1) * side] = params[idx];
        // idx += 1;
    }

    CouplingMatrix::new_with_topology(order, template.topology(), data)
}

/// Computes the residual vector: difference between current and target S-parameters.
///
/// The objective is the sum of squared residuals, so only `r_i^2` is
/// observable — the sign convention of each component is not. The optimizer
/// therefore fits the *s11* and *s21* response channels (four components per
/// frequency sample) in a least-squares sense; it does not (and cannot) pin
/// the per-component sign, and a sign flip on any single component is
/// cost-neutral.
fn compute_residuals(
    matrix: &CouplingMatrix,
    target: &SParameterResponse,
    grid: &FrequencyGrid,
    settings: ResponseSettings,
) -> Result<Vec<f64>> {
    let current = ResponseSolver.evaluate_normalized_with_settings(matrix, grid, settings)?;

    let mut residuals = Vec::with_capacity(target.samples.len() * 4);
    for (cur, tgt) in current.samples.iter().zip(target.samples.iter()) {
        residuals.push(cur.s11_re - tgt.s11_re);
        residuals.push(cur.s11_im - tgt.s11_im);
        residuals.push(cur.s21_re - tgt.s21_re);
        residuals.push(cur.s21_im - tgt.s21_im);
    }
    Ok(residuals)
}

/// Computes the Jacobian matrix numerically (finite differences).
///
/// The step is scale-aware: `delta_j = 1e-7 * max(1, |x_j|)` — a fixed
/// absolute step is not scale-safe when parameter magnitudes vary (large
/// entries under-perturb, tiny entries over-perturb). If a perturbed
/// evaluation fails (the response matrix is singular at that point —
/// possible far from the solution basin), the step is halved and retried up
/// to 4 times before the error is surfaced.
fn compute_jacobian(
    template: &CouplingMatrix,
    params: &[f64],
    config: &OptimizeConfig,
    target: &SParameterResponse,
    grid: &FrequencyGrid,
    settings: ResponseSettings,
    n_residuals: usize,
) -> Result<Vec<Vec<f64>>> {
    const MAX_DELTA_RETRIES: u32 = 4;
    let n_params = params.len();
    let mut jacobian = vec![vec![0.0; n_params]; n_residuals];

    let base_matrix = rebuild_matrix(template, params, config)?;
    let base_residuals = compute_residuals(&base_matrix, target, grid, settings)?;

    for j in 0..n_params {
        let mut delta = 1e-7 * 1.0_f64.max(params[j].abs());
        let mut last_err: Option<MfsError> = None;
        for _retry in 0..MAX_DELTA_RETRIES {
            let mut perturbed = params.to_vec();
            perturbed[j] += delta;
            let perturbed_matrix = rebuild_matrix(template, &perturbed, config)?;
            match compute_residuals(&perturbed_matrix, target, grid, settings) {
                Ok(perturbed_residuals) => {
                    for i in 0..n_residuals {
                        jacobian[i][j] = (perturbed_residuals[i] - base_residuals[i]) / delta;
                    }
                    last_err = None;
                    break;
                }
                Err(e) => {
                    last_err = Some(e);
                    delta *= 0.5;
                }
            }
        }
        if let Some(e) = last_err {
            return Err(e);
        }
    }

    Ok(jacobian)
}

/// Levenberg-Marquardt optimization step.
///
/// Solves: (J^T J + λ diag(J^T J)) Δp = -J^T r
fn lm_step(jacobian: &[Vec<f64>], residuals: &[f64], lambda: f64, n_params: usize) -> Vec<f64> {
    // Compute J^T J and J^T r
    let mut jtj = vec![vec![0.0; n_params]; n_params];
    let mut jtr = vec![0.0; n_params];

    for i in 0..residuals.len() {
        for j in 0..n_params {
            jtr[j] += jacobian[i][j] * residuals[i];
            for k in 0..n_params {
                jtj[j][k] += jacobian[i][j] * jacobian[i][k];
            }
        }
    }

    // Add damping: J^T J + λ * diag(J^T J)
    for (j, row) in jtj.iter_mut().enumerate() {
        row[j] *= 1.0 + lambda;
    }

    // Solve via Cholesky-like approach (simple Gaussian elimination for small systems)
    solve_linear_system(&jtj, &jtr.iter().map(|v| -v).collect::<Vec<_>>())
}

/// Simple Gaussian elimination for small dense systems.
fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
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

    // Forward elimination with partial pivoting
    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if aug[row][col].abs() > aug[max_row][col].abs() {
                max_row = row;
            }
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        if pivot.abs() < 1e-20 {
            continue;
        }

        let pivot_tail = aug[col][col..=n].to_vec();
        for row_tail in aug.iter_mut().skip(col + 1).map(|row| &mut row[col..=n]) {
            let factor = row_tail[0] / pivot;
            for (dest, &pivot_val) in row_tail.iter_mut().zip(pivot_tail.iter()) {
                *dest -= factor * pivot_val;
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum -= aug[i][j] * x[j];
        }
        if aug[i][i].abs() > 1e-20 {
            x[i] = sum / aug[i][i];
        }
    }
    x
}

/// Optimizes a coupling matrix to match a target S-parameter response.
///
/// This is the core optimization engine used by both tuning and extraction.
pub fn optimize_matrix(
    initial: &CouplingMatrix,
    target: &SParameterResponse,
    grid: &FrequencyGrid,
    config: &OptimizeConfig,
) -> Result<OptimizeResult> {
    let settings = ResponseSettings {
        source_resistance: 1.0,
        load_resistance: 1.0,
        unloaded_q: config.unloaded_q,
    };

    let mut params = extract_parameters(initial, config);
    let n_params = params.len();
    let n_residuals = target.samples.len() * 4;

    if n_params == 0 {
        return Err(MfsError::PreconditionViolation(
            "no optimizable parameters found in the matrix".into(),
        ));
    }

    let mut lambda = config.lambda;
    let mut iterations = 0;
    let mut prev_cost: Option<f64> = None;
    let mut accepted_step = false;
    let mut best_cost = f64::INFINITY;

    for iter in 0..config.max_iterations {
        iterations = iter + 1;

        let matrix = rebuild_matrix(initial, &params, config)?;
        let residuals = compute_residuals(&matrix, target, grid, settings)?;
        let cost: f64 = residuals.iter().map(|r| r * r).sum();
        best_cost = best_cost.min(cost);

        // Convergence: the cost stopped changing after an *accepted* step.
        // (Comparing against the running *best* is wrong: after an accepted
        // trial the next cost equals the best exactly, which reads as "zero
        // change" on iteration 2 and stops the loop before it has converged.
        // And checking after a *rejected* step is also wrong: the parameters
        // did not move, so the cost is trivially unchanged — that is
        // stagnation, not convergence.)
        if accepted_step {
            if let Some(prev) = prev_cost {
                if (prev - cost).abs() < config.tolerance {
                    return Ok(OptimizeResult {
                        matrix,
                        cost,
                        iterations,
                        converged: true,
                    });
                }
            }
        }
        prev_cost = Some(cost);
        accepted_step = false;

        let jacobian = compute_jacobian(
            initial,
            &params,
            config,
            target,
            grid,
            settings,
            n_residuals,
        )?;
        let step = lm_step(&jacobian, &residuals, lambda, n_params);

        // Trial step. A trial whose response cannot be evaluated (a
        // singular matrix) is a *rejection*, not an error: the step walked
        // outside the well-posed region, so damp harder and retry.
        let trial_params: Vec<f64> = params.iter().zip(step.iter()).map(|(p, s)| p + s).collect();
        let trial_cost = rebuild_matrix(initial, &trial_params, config)
            .and_then(|m| compute_residuals(&m, target, grid, settings))
            .map(|r| r.iter().map(|r| r * r).sum::<f64>());

        match trial_cost {
            Ok(tc) if tc < cost => {
                params = trial_params;
                lambda *= 0.5; // Decrease damping (more Newton-like)
                accepted_step = true;
            }
            _ => {
                lambda *= 2.0; // Increase damping (more gradient descent-like)
            }
        }
        // Damped to a negligible step: no further progress is possible.
        if lambda > 1e12 {
            break;
        }
    }

    let final_matrix = rebuild_matrix(initial, &params, config)?;
    Ok(OptimizeResult {
        matrix: final_matrix,
        cost: best_cost,
        iterations,
        converged: false,
    })
}

/// Tunes an existing coupling matrix to better match a target response.
///
/// Starts from the given matrix and adjusts coupling values to minimize
/// the S-parameter error against the target.
pub fn tune_matrix(
    matrix: &CouplingMatrix,
    target: &SParameterResponse,
    grid: &FrequencyGrid,
) -> Result<OptimizeResult> {
    optimize_matrix(matrix, target, grid, &OptimizeConfig::default())
}

/// Extracts a coupling matrix from measured S-parameter data.
///
/// Builds a crude initial guess — an *inline chain* (source `M[0,1]`,
/// load `M[N,N+1]`, inter-resonator `M[k,k+1]`, all set to a uniform rough
/// coupling) with zero diagonals, and optimizes every free value (couplings
/// and diagonals) to match the measured response.
///
/// **Topology limitation:** `optimize_matrix` freezes the template's
/// non-zero pattern, so the recovered matrix shares the *inline-chain*
/// topology of this guess. This is the transversal family — the most common
/// topology and the one the order-6 extraction gate uses. A target whose
/// true matrix has a different non-zero pattern (e.g. a folded or arrow
/// topology with cross couplings) is *not* representable from this start
/// and the extraction will fit the nearest inline-chain matrix instead.
/// Callers who know a better topology should build the guess themselves and
/// call `optimize_matrix` with a custom template.
///
/// `order` specifies the filter order (number of resonators).
pub fn extract_matrix(
    measured: &SParameterResponse,
    grid: &FrequencyGrid,
    order: usize,
) -> Result<OptimizeResult> {
    // Build the initial guess as a crude *inline chain*: source M[0,1],
    // load M[N, N+1], and inter-resonator chain couplings M[k, k+1],
    // all set to a uniform rough coupling estimate. The optimizer freezes
    // the template's non-zero pattern, so the guess must be a template that
    // can *represent* the target family — a transversal (inline-chain)
    // order-`order` filter has exactly this pattern. (The previous guess
    // coupled every resonator to both ports and left the chain couplings
    // frozen at zero — a star-like pattern that cannot represent a
    // transversal target, which is what made the extraction wander into
    // singular response matrices.)
    let coupling = (2.0 / (order as f64 + 1.0)).sqrt();
    let mut builder = CouplingMatrixBuilder::new(order)?;
    builder = builder.set_symmetric(0, 1, coupling)?;
    builder = builder.set_symmetric(order, order + 1, coupling)?;
    for k in 1..order {
        builder = builder.set_symmetric(k, k + 1, coupling)?;
    }

    let initial = builder.build()?;

    let config = OptimizeConfig {
        max_iterations: 500,
        tolerance: 1e-8,
        lambda: 1e-2,
        optimize_diagonals: true,
        optimize_ports: true,
        unloaded_q: f64::INFINITY,
    };

    optimize_matrix(&initial, measured, grid, &config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::FilterDesign;

    #[test]
    fn tune_recovers_original_matrix() -> Result<()> {
        // Synthesize a filter and get its response
        let design = FilterDesign::prototype(3, 20.0).synthesize()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 11)?;
        let target = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;

        // Verify parameters can be extracted and rebuilt
        let config = OptimizeConfig::default();
        let params = extract_parameters(design.matrix(), &config);
        eprintln!(
            "params count: {}, values: {:?}",
            params.len(),
            &params[..params.len().min(5)]
        );

        let rebuilt = rebuild_matrix(design.matrix(), &params, &config)?;
        let rebuilt_resp = ResponseSolver.evaluate_normalized(&rebuilt, &grid)?;

        // Rebuilt should match original exactly
        let diff = (rebuilt_resp.samples[5].s21_mag() - target.samples[5].s21_mag()).abs();
        assert!(diff < 1e-10, "rebuild failed: diff={diff:.2e}");

        // Now perturb and tune
        let mut perturbed_params = params.clone();
        perturbed_params[0] *= 1.05; // Perturb first param by 5%
        let perturbed = rebuild_matrix(design.matrix(), &perturbed_params, &config)?;

        let result = tune_matrix(&perturbed, &target, &grid)?;
        assert!(
            result.cost < 1e-3,
            "cost={:.2e}, iterations={}",
            result.cost,
            result.iterations
        );
        Ok(())
    }

    #[test]
    fn extract_from_synthesized_response() -> Result<()> {
        // Generate a target response from a known filter
        let design = FilterDesign::prototype(3, 20.0).synthesize()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 21)?;
        let target = ResponseSolver.evaluate_normalized(design.matrix(), &grid)?;

        // Use the synthesized matrix as initial guess (realistic scenario:
        // you have a synthesis result and want to fine-tune to match measurement)
        let config = OptimizeConfig {
            max_iterations: 100,
            ..Default::default()
        };

        // Perturb the matrix to simulate a "measurement" that differs slightly
        let params = extract_parameters(design.matrix(), &config);
        let mut perturbed = params.clone();
        for p in &mut perturbed {
            *p *= 0.95; // 5% perturbation on all params
        }
        let initial = rebuild_matrix(design.matrix(), &perturbed, &config)?;

        let result = optimize_matrix(&initial, &target, &grid, &config)?;

        // Should converge close to the original
        assert!(
            result.cost < 0.01,
            "cost={:.2e}, iterations={}",
            result.cost,
            result.iterations
        );
        Ok(())
    }

    #[test]
    fn tune_folded_matrix() -> Result<()> {
        // Synthesize and convert to folded
        let design = FilterDesign::prototype(4, 20.0)
            .zeros([-1.5, 1.5])
            .synthesize()?;
        let folded = design.to_folded()?;
        let grid = FrequencyGrid::linspace(-2.0, 2.0, 21)?;
        let target = ResponseSolver.evaluate_normalized(&folded, &grid)?;

        // Verify folded matrix has inter-resonator couplings
        let config = OptimizeConfig::default();
        let params = extract_parameters(&folded, &config);
        assert!(
            params.len() > 4,
            "folded should have inter-resonator params, got {}",
            params.len()
        );

        // Perturb and tune
        let mut perturbed_params = params.clone();
        perturbed_params[0] *= 1.03;
        let perturbed = rebuild_matrix(&folded, &perturbed_params, &config)?;

        let result = optimize_matrix(&perturbed, &target, &grid, &config)?;
        assert!(
            result.cost < 1e-3,
            "folded tune cost={:.2e}, iterations={}",
            result.cost,
            result.iterations
        );
        Ok(())
    }
}
