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
use crate::response::{ResponseSolver, ResponseSettings, SParameterResponse, ResponseSample};

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
                data[0 * side + k] = params[idx];
                data[k * side + 0] = params[idx];
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
        data[0 * side + (side - 1)] = params[idx];
        data[(side - 1) * side + 0] = params[idx];
        // idx += 1;
    }

    CouplingMatrix::new_with_topology(order, template.topology(), data)
}

/// Computes the residual vector: difference between current and target S-parameters.
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
fn compute_jacobian(
    template: &CouplingMatrix,
    params: &[f64],
    config: &OptimizeConfig,
    target: &SParameterResponse,
    grid: &FrequencyGrid,
    settings: ResponseSettings,
    n_residuals: usize,
) -> Result<Vec<Vec<f64>>> {
    let n_params = params.len();
    let delta = 1e-7;
    let mut jacobian = vec![vec![0.0; n_params]; n_residuals];

    let base_matrix = rebuild_matrix(template, params, config)?;
    let base_residuals = compute_residuals(&base_matrix, target, grid, settings)?;

    for j in 0..n_params {
        let mut perturbed = params.to_vec();
        perturbed[j] += delta;
        let perturbed_matrix = rebuild_matrix(template, &perturbed, config)?;
        let perturbed_residuals = compute_residuals(&perturbed_matrix, target, grid, settings)?;

        for i in 0..n_residuals {
            jacobian[i][j] = (perturbed_residuals[i] - base_residuals[i]) / delta;
        }
    }

    Ok(jacobian)
}

/// Levenberg-Marquardt optimization step.
///
/// Solves: (J^T J + λ diag(J^T J)) Δp = -J^T r
fn lm_step(
    jacobian: &[Vec<f64>],
    residuals: &[f64],
    lambda: f64,
    n_params: usize,
) -> Vec<f64> {
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
    for j in 0..n_params {
        jtj[j][j] *= 1.0 + lambda;
    }

    // Solve via Cholesky-like approach (simple Gaussian elimination for small systems)
    solve_linear_system(&jtj, &jtr.iter().map(|v| -v).collect::<Vec<_>>())
}

/// Simple Gaussian elimination for small dense systems.
fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut aug: Vec<Vec<f64>> = a.iter().enumerate()
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

        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for k in col..=n {
                aug[row][k] -= factor * aug[col][k];
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
    let mut best_cost = f64::INFINITY;
    let mut iterations = 0;

    for iter in 0..config.max_iterations {
        iterations = iter + 1;

        let matrix = rebuild_matrix(initial, &params, config)?;
        let residuals = compute_residuals(&matrix, target, grid, settings)?;
        let cost: f64 = residuals.iter().map(|r| r * r).sum();

        if (best_cost - cost).abs() < config.tolerance && iter > 0 {
            let final_matrix = rebuild_matrix(initial, &params, config)?;
            return Ok(OptimizeResult {
                matrix: final_matrix,
                cost,
                iterations,
                converged: true,
            });
        }

        let jacobian = compute_jacobian(initial, &params, config, target, grid, settings, n_residuals)?;
        let step = lm_step(&jacobian, &residuals, lambda, n_params);

        // Trial step
        let trial_params: Vec<f64> = params.iter().zip(step.iter()).map(|(p, s)| p + s).collect();
        let trial_matrix = rebuild_matrix(initial, &trial_params, config)?;
        let trial_residuals = compute_residuals(&trial_matrix, target, grid, settings)?;
        let trial_cost: f64 = trial_residuals.iter().map(|r| r * r).sum();

        if trial_cost < cost {
            params = trial_params;
            best_cost = trial_cost;
            lambda *= 0.5; // Decrease damping (more Newton-like)
        } else {
            lambda *= 2.0; // Increase damping (more gradient descent-like)
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
/// Starts from an initial guess (e.g., from synthesis) and optimizes
/// all coupling values to match the measured response.
///
/// `order` specifies the filter order (number of resonators).
pub fn extract_matrix(
    measured: &SParameterResponse,
    grid: &FrequencyGrid,
    order: usize,
) -> Result<OptimizeResult> {
    // Build initial guess: estimate couplings from the response bandwidth
    // For a Chebyshev filter, the source/load coupling ≈ 1/sqrt(g0*g1)
    // and chain couplings ≈ 1/sqrt(gk*g(k+1)). Use uniform approximation.
    let coupling = (2.0 / (order as f64 + 1.0)).sqrt();
    let mut builder = CouplingMatrixBuilder::new(order)?;

    // Source couplings (transversal: all resonators couple to source)
    for k in 1..=order {
        builder = builder.set_symmetric(0, k, coupling / (k as f64).sqrt())?;
    }
    // Load couplings
    for k in 1..=order {
        builder = builder.set_symmetric(k, order + 1, coupling / ((order + 1 - k) as f64).sqrt())?;
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
        eprintln!("params count: {}, values: {:?}", params.len(), &params[..params.len().min(5)]);

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
        assert!(result.cost < 1e-3,
            "cost={:.2e}, iterations={}", result.cost, result.iterations);
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
        assert!(result.cost < 0.01,
            "cost={:.2e}, iterations={}", result.cost, result.iterations);
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
        assert!(params.len() > 4, "folded should have inter-resonator params, got {}", params.len());

        // Perturb and tune
        let mut perturbed_params = params.clone();
        perturbed_params[0] *= 1.03;
        let perturbed = rebuild_matrix(&folded, &perturbed_params, &config)?;

        let result = optimize_matrix(&perturbed, &target, &grid, &config)?;
        assert!(result.cost < 1e-3,
            "folded tune cost={:.2e}, iterations={}", result.cost, result.iterations);
        Ok(())
    }
}
