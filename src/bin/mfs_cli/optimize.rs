//! §D.5: the `optimize` and `extract` subcommands.
//!
//! Pure move out of the monolithic `src/bin/mfs_cli.rs` — no behavior
//! change; the shared `load_target` / `parse_matrix_json` helpers live
//! in [`crate::target`].
//!
use crate::Cli;
use crate::io::{read_input, CliRunError, OutputFormat};
use crate::target::{load_target, parse_matrix_json};
/// Re-evaluates a matrix on the grid defined by the target's frequency samples,
/// for display next to the target.
fn recovered_response(
    matrix: &mfs::matrix::CouplingMatrix,
    target: &mfs::response::SParameterResponse,
) -> mfs::Result<mfs::response::SParameterResponse> {
    use mfs::freq::FrequencyGrid;
    use mfs::response::ResponseSolver;
    let n = target.samples.len();
    let first = target.samples.first().map(|s| s.frequency_hz).unwrap_or(-2.0);
    let last = target.samples.last().map(|s| s.frequency_hz).unwrap_or(2.0);
    let grid = FrequencyGrid::linspace(first, last, n.max(2))?;
    ResponseSolver.evaluate_normalized(matrix, &grid)
}

/// Renders the optimized / extracted matrix as JSON.
fn optimize_result_json(
    order: usize,
    result: &mfs::optimize::OptimizeResult,
    target: &mfs::response::SParameterResponse,
) -> Result<serde_json::Value, CliRunError> {
    let side = result.matrix.side();
    let rows: Vec<Vec<f64>> = (0..side)
        .map(|r| (0..side).map(|c| result.matrix.at(r, c).unwrap_or(0.0)).collect())
        .collect();
    let mut value = serde_json::json!({
        "order": order,
        "topology": format!("{:?}", result.matrix.topology()),
        "matrix": rows,
        "data": rows.iter().flatten().copied().collect::<Vec<f64>>(),
        "cost": result.cost,
        "iterations": result.iterations,
        "converged": result.converged,
    });
    if let Ok(rec) = recovered_response(&result.matrix, target) {
        let samples: Vec<serde_json::Value> = rec
            .samples
            .iter()
            .map(|s| {
                serde_json::json!({
                    "frequency_hz": s.frequency_hz,
                    "s11_re": s.s11_re,
                    "s11_im": s.s11_im,
                    "s21_re": s.s21_re,
                    "s21_im": s.s21_im,
                })
            })
            .collect();
        value["recovered"] = serde_json::json!({
            "samples": samples,
        });
    }
    Ok(value)
}

/// Runs the `optimize` subcommand.
pub(crate) fn run_optimize(
    cli: &Cli,
    order: usize,
    target: Option<&std::path::Path>,
    config: &mfs::optimize::OptimizeConfig,
    output: Option<&std::path::Path>,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::optimize::optimize_matrix;
    // The matrix always comes from --input/stdin (JSON); only the *target*
    // may be a .s2p. The documented "matrix from input + Touchstone target"
    // usage requires reading the matrix even when the target is a .s2p
    // (mirrors `run_extract`, review P2-1).
    let json_input = read_input(cli)?;
    let target_resp = load_target(target, &json_input)?;
    let matrix = parse_matrix_json(&json_input)?;
    if matrix.order() != order {
        return Err(CliRunError::InputParse(format!(
            "--order {order} does not match the input matrix order {}",
            matrix.order()
        )));
    }
    let result = optimize_matrix(&matrix, &target_resp, &response_grid(&target_resp)?, config)
        .map_err(CliRunError::Pipeline)?;
    let value = optimize_result_json(order, &result, &target_resp)?;
    if let Some(out) = output {
        std::fs::write(out, serde_json::to_string_pretty(&value)?)
            .map_err(|e| CliRunError::InputRead(format!("failed to write '{}': {e}", out.display())))?;
        eprintln!("Optimized matrix written to: {}", out.display());
    }
    match format {
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Table => {
            println!("── Optimized Matrix ──");
            println!("  Order:     {order}");
            println!("  Topology:  {:?}", result.matrix.topology());
            println!("  Cost:      {:.3e}", result.cost);
            println!("  Iters:     {}", result.iterations);
            println!("  Converged: {}", result.converged);
            let side = result.matrix.side();
            if side <= 10 {
                println!("  Matrix:");
                for i in 0..side {
                    let row: Vec<String> = (0..side)
                        .map(|j| format!("{:>8.4}", result.matrix.at(i, j).unwrap_or(0.0)))
                        .collect();
                    println!("    [{}]", row.join(" "));
                }
            }
        }
    }
    Ok(())
}

/// Builds the grid the optimizer evaluates on, from the target response.
///
/// §C.7: the optimizer evaluates the matrix with `evaluate_normalized`, so the
/// target's `frequency_hz` samples are treated **directly as the normalized
/// grid ω** (the library contract — see `mfs::optimize` tests, which build the
/// grid on `[-2, 2]`). There is no `--normalized` flag: the grid is *always*
/// normalized. A target whose frequency range spans more than 1e6 (i.e. it
/// looks like physical Hz, e.g. a raw Touchstone `.s2p`) is rejected with a
/// typed error pointing at the normalized contract, rather than being
/// silently evaluated as a meaningless ω ∈ [~1e9, ~2e9] grid.
fn response_grid(target: &mfs::response::SParameterResponse) -> Result<mfs::freq::FrequencyGrid, CliRunError> {
    use mfs::freq::FrequencyGrid;
    let n = target.samples.len().max(2);
    let first = target.samples.first().map(|s| s.frequency_hz).unwrap_or(-2.0);
    let last = target.samples.last().map(|s| s.frequency_hz).unwrap_or(2.0);
    // §C.7: reject physical-Hz-scale grids. The normalized bandpass window is
    // |ω| ≲ 3, so a span > 1e6 cannot be a normalized ω.
    if (last - first).abs() > 1e6 {
        return Err(CliRunError::Pipeline(
            mfs::error::MfsError::PreconditionViolation(format!(
                "target frequency range spans {span:.3e} (first {first:.3e}, last {last:.3e}) — this \
                 looks like physical Hz, but the optimizer evaluates on the NORMALIZED grid ω \
                 (|ω| ≲ 3). Provide the target on the normalized grid (e.g. ω ∈ [-2, 2]); there is \
                 no --normalized flag, the grid is always normalized.",
                span = (last - first).abs(),
            )),
        ));
    }
    FrequencyGrid::linspace(first, last, n).map_err(CliRunError::Pipeline)
}

/// Runs the `extract` subcommand.
pub(crate) fn run_extract(
    cli: &Cli,
    order: usize,
    target: Option<&std::path::Path>,
    output: Option<&std::path::Path>,
    format: &OutputFormat,
) -> Result<(), CliRunError> {
    use mfs::optimize::extract_matrix;
    let json_input = read_input(cli)?;
    let target_resp = load_target(target, &json_input)?;
    let result = extract_matrix(&target_resp, &response_grid(&target_resp)?, order)
        .map_err(CliRunError::Pipeline)?;
    let value = optimize_result_json(order, &result, &target_resp)?;
    if let Some(out) = output {
        std::fs::write(out, serde_json::to_string_pretty(&value)?)
            .map_err(|e| CliRunError::InputRead(format!("failed to write '{}': {e}", out.display())))?;
        eprintln!("Extracted matrix written to: {}", out.display());
    }
    match format {
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Table => {
            println!("── Extracted Matrix ──");
            println!("  Order:     {order}");
            println!("  Topology:  {:?}", result.matrix.topology());
            println!("  Cost:      {:.3e}", result.cost);
            println!("  Iters:     {}", result.iterations);
            println!("  Converged: {}", result.converged);
            let side = result.matrix.side();
            if side <= 10 {
                println!("  Matrix:");
                for i in 0..side {
                    let row: Vec<String> = (0..side)
                        .map(|j| format!("{:>8.4}", result.matrix.at(i, j).unwrap_or(0.0)))
                        .collect();
                    println!("    [{}]", row.join(" "));
                }
            }
        }
    }
    Ok(())
}
#[cfg(test)]
mod optimizer_tests {
    use super::*;
    use crate::target::parse_response_json;
    use mfs::design::FilterDesign;
    use mfs::freq::FrequencyGrid;
    use mfs::response::ResponseSolver;

    /// Scales the first `count` non-zero inter-resonator couplings (upper
    /// triangle, off-diagonal) of `matrix` by `factor`, mirroring the
    /// integration-test perturbation helper.
    fn perturb(matrix: &mfs::matrix::CouplingMatrix, count: usize, factor: f64) -> mfs::matrix::CouplingMatrix {
        let side = matrix.side();
        let order = matrix.order();
        let mut data = matrix.as_slice().to_vec();
        let mut done = 0;
        for row in 1..=order {
            for col in (row + 1)..=order {
                let idx = row * side + col;
                if data[idx].abs() > 1e-9 {
                    data[idx] *= factor;
                    data[col * side + row] *= factor;
                    done += 1;
                    if done == count {
                        return mfs::matrix::CouplingMatrix::new_with_topology(
                            order,
                            matrix.topology(),
                            data,
                        )
                        .expect("rebuild");
                    }
                }
            }
        }
        mfs::matrix::CouplingMatrix::new_with_topology(order, matrix.topology(), data).expect("rebuild")
    }

    /// Worst relative deviation between two matrices over entries with |orig| > 1e-9.
    fn worst_rel(a: &mfs::matrix::CouplingMatrix, b: &mfs::matrix::CouplingMatrix) -> f64 {
        let side = a.side();
        let mut w = 0.0_f64;
        for i in 0..side {
            for j in i..side {
                let o = a.at(i, j).unwrap_or(0.0);
                let r = b.at(i, j).unwrap_or(0.0);
                if o.abs() > 1e-9 {
                    w = w.max((o - r).abs() / o.abs());
                }
            }
        }
        w
    }

    #[test]
    fn parse_matrix_json_flat_and_nested() {
        // Flat row-major with order + topology. Order-1 (side-3) ladder:
        // source->res 0.7, res->load 0.7 (symmetric).
        let m = parse_matrix_json(r#"{"order":1,"data":[0,0.7,0,0.7,0,0.7,0,0.7,0]}"#).unwrap();
        assert_eq!(m.order(), 1);
        assert!((m.at(0, 1).unwrap() - 0.7).abs() < 1e-12);
        assert_eq!(m.topology(), mfs::matrix::MatrixTopology::Transversal);

        // Folded topology honored.
        let m = parse_matrix_json(r#"{"order":1,"topology":"folded","data":[0,0.7,0,0.7,0,0.7,0,0.7,0]}"#)
            .unwrap();
        assert_eq!(m.topology(), mfs::matrix::MatrixTopology::Folded);

        // Nested rows, order inferred (side 4 => order 2).
        let m = parse_matrix_json(r#"[ [0,0,0,0],[0,1,1,0],[0,1,2,0],[0,0,0,0] ]"#).unwrap();
        assert_eq!(m.order(), 2);
        assert!((m.at(1, 2).unwrap() - 1.0).abs() < 1e-12);

        // Order / side mismatch is a typed dimension error.
        let err = parse_matrix_json(r#"{"order":2,"data":[0,0.7,0,0.7,0,0.7,0,0.7,0]}"#).unwrap_err();
        assert!(matches!(err, mfs::error::MfsError::DimensionMismatch { .. }));
    }

    #[test]
    fn parse_matrix_json_rejects_one_sided_off_diagonal() {
        // Order-2, side-4: source couplings 0.8 / 0.5, chain 0.6, load 0.7 / 0.4,
        // but the chain pair is deliberately asymmetric (0.6 vs 0.9).
        let data = [
            0.0, 0.8, 0.5, 0.0,
            0.8, 0.0, 0.6, 0.7,
            0.5, 0.9, 0.0, 0.4,
            0.0, 0.7, 0.4, 0.0,
        ];
        let json = serde_json::json!({
            "order": 2,
            "data": data,
        });
        let err = parse_matrix_json(&json.to_string()).unwrap_err();
        assert!(
            matches!(err, mfs::error::MfsError::PreconditionViolation(_)),
            "a one-sided off-diagonal user matrix must be rejected, got {err:?}"
        );
        assert!(
            err.to_string().contains("asymmetric coupling matrix"),
            "the error must say the matrix is asymmetric: {err}"
        );
    }

    #[test]
    fn parse_response_json_samples() {
        let r = parse_response_json(
            r#"{"samples":[
               {"frequency_hz":-1.0,"s11_re":0.1,"s11_im":-0.2,"s21_re":0.9,"s21_im":0.0},
               {"frequency_hz":1.0,"s11_re":0.1,"s11_im":0.2,"s21_re":0.9,"s21_im":0.0}
            ]}"#,
        )
        .unwrap();
        assert_eq!(r.samples.len(), 2);
        assert!((r.samples[0].frequency_hz - (-1.0)).abs() < 1e-12);
        assert!((r.samples[1].s21_im - 0.0).abs() < 1e-12);
    }

    #[test]
    fn optimize_subcommand_recovers_perturbed_matrix() {
        // Order-6 folded, perturb 4 couplings x1.07, recover via the CLI core.
        let original = FilterDesign::prototype(6, 20.0)
            .synthesize()
            .unwrap()
            .to_folded()
            .unwrap();
        let g = FrequencyGrid::linspace(-2.0, 2.0, 41).unwrap();
        let target = ResponseSolver.evaluate_normalized(&original, &g).unwrap();
        let perturbed = perturb(&original, 4, 1.07);

        let config = mfs::optimize::OptimizeConfig {
            max_iterations: 500,
            tolerance: 1e-10,
            lambda: 1e-2,
            ..Default::default()
        };
        let result = mfs::optimize::optimize_matrix(&perturbed, &target, &g, &config).unwrap();
        assert!(result.converged);
        assert!(worst_rel(&original, &result.matrix) < 1e-4);

        // The JSON round-trip through parse_matrix_json must preserve the
        // perturbed matrix so a real CLI invocation sees the same start.
        let side = perturbed.side();
        let nested: Vec<Vec<f64>> = (0..side)
            .map(|r| (0..side).map(|c| perturbed.at(r, c).unwrap_or(0.0)).collect())
            .collect();
        let json = format!(
            r#"{{"order": 6, "topology": "folded", "data": {}}}"#,
            serde_json::to_string(&nested).unwrap()
        );
        let parsed = parse_matrix_json(&json).unwrap();
        assert_eq!(parsed.topology(), mfs::matrix::MatrixTopology::Folded);
        assert!(worst_rel(&perturbed, &parsed) < 1e-12);

    }

    #[test]
    fn extract_subcommand_recovers_transversal() {
        // Order-6 transversal target; extract from a crude guess via the CLI core.
        let original = FilterDesign::prototype(6, 20.0).synthesize().unwrap().matrix().clone();
        let g = FrequencyGrid::linspace(-2.0, 2.0, 41).unwrap();
        let target = ResponseSolver.evaluate_normalized(&original, &g).unwrap();
        let result = mfs::optimize::extract_matrix(&target, &g, 6).unwrap();
        assert!(result.converged);

        let rec = ResponseSolver.evaluate_normalized(&result.matrix, &g).unwrap();
        let mut max_dev = 0.0_f64;
        for (a, b) in rec.samples.iter().zip(target.samples.iter()) {
            max_dev = max_dev
                .max((a.s11_re - b.s11_re).abs())
                .max((a.s11_im - b.s11_im).abs())
                .max((a.s21_re - b.s21_re).abs())
                .max((a.s21_im - b.s21_im).abs());
        }
        assert!(max_dev < 1e-3);
    }

    /// §C.7 regression: the optimizer/extract grid is the target's
    /// `frequency_hz` samples used *directly as the normalized grid ω*. A
    /// target on a physical-Hz scale (range span > 1e6) must be rejected with
    /// a typed error rather than silently evaluated on a meaningless ω ∈
    /// [~1e9, ~2e9] grid; a normalized grid (span ≲ 3) must pass.
    #[test]
    fn response_grid_rejects_physical_hz_and_accepts_normalized() {
        // Physical-Hz-scale target (a raw .s2p grid, e.g. 1 GHz .. 2 GHz):
        // span 1e9 >> 1e6, must be rejected.
        let physical = parse_response_json(
            r#"{"samples":[
               {"frequency_hz":1.0e9,"s21_re":0.0,"s21_im":0.0},
               {"frequency_hz":2.0e9,"s21_re":0.0,"s21_im":0.0}
            ]}"#,
        )
        .unwrap();
        let err = response_grid(&physical).unwrap_err();
        assert!(
            matches!(err, CliRunError::Pipeline(mfs::error::MfsError::PreconditionViolation(_))),
            "a physical-Hz target must be rejected with a typed precondition error, got {err:?}"
        );
        assert!(
            err.to_string().contains("NORMALIZED"),
            "the error must point at the normalized-grid contract: {err}"
        );

        // Normalized-grid target (ω ∈ [-2, 2], span 4 ≪ 1e6): must pass.
        let normalized = parse_response_json(
            r#"{"samples":[
               {"frequency_hz":-2.0,"s21_re":0.0,"s21_im":0.0},
               {"frequency_hz":2.0,"s21_re":0.0,"s21_im":0.0}
            ]}"#,
        )
        .unwrap();
        assert!(
            response_grid(&normalized).is_ok(),
            "a normalized-grid target (span 4) must be accepted"
        );
    }
}
