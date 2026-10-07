//! Integration tests for the `mfs_cli` binary.
//!
//! These tests invoke the compiled CLI binary using `std::process::Command` and verify:
//! - JSON file input → JSON stdout output
//! - stdin input → JSON stdout output
//! - `--stage` flag produces partial output
//! - Invalid input produces structured error on stderr
//!
//! The binary requires the `cli` feature to build, so these tests are gated behind
//! `#[cfg(feature = "cli")]`. Run with:
//!
//! ```sh
//! cargo test --features cli --test cli_integration
//! ```

#![cfg(feature = "cli")]

use std::io::Write;
use std::process::{Command, Stdio};

/// Returns the path to the `mfs_cli` binary built by cargo.
fn cli_bin() -> std::path::PathBuf {
    // `cargo test` places test binaries in target/debug/deps, but the actual
    // binary is in target/debug (or target/release). We use the CARGO_BIN_EXE
    // env var if available (set by cargo for integration tests with [[bin]]),
    // otherwise construct the path manually.
    let bin_name = if cfg!(windows) {
        "mfs_cli.exe"
    } else {
        "mfs_cli"
    };

    // Try the env var first (available since Rust 1.43 for integration tests)
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_mfs_cli") {
        return std::path::PathBuf::from(path);
    }

    // Fallback: construct path relative to the test binary location
    let mut path = std::env::current_exe()
        .expect("cannot determine test binary path")
        .parent()
        .expect("no parent directory")
        .parent()
        .expect("no grandparent directory")
        .to_path_buf();
    path.push(bin_name);
    path
}

/// The `mfs_cli` binary path, asserted to exist. A missing binary is a
/// harness failure (the `cli` feature must build it), not a silent skip.
fn require_bin() -> std::path::PathBuf {
    let bin = cli_bin();
    assert!(
        bin.exists(),
        "mfs_cli binary not found at {} (build with `cargo build --features cli`)",
        bin.display()
    );
    bin
}

/// A valid minimal synthesis request JSON string.
fn minimal_request_json() -> &'static str {
    r#"{"order": 4, "return_loss_db": 20.0}"#
}

/// A valid full synthesis request JSON string with all optional fields.
fn full_request_json() -> &'static str {
    r#"{
        "order": 4,
        "return_loss_db": 20.0,
        "transmission_zeros": [-2.0, 1.5],
        "topology": "Folded",
        "grid": { "start": -3.0, "stop": 3.0, "points": 51 }
    }"#
}

/// An invalid request JSON (order = 0 triggers validation error).
fn invalid_request_json() -> &'static str {
    r#"{"order": 0, "return_loss_db": 20.0}"#
}

// ---------------------------------------------------------------------------
// Test: JSON file input → JSON stdout output
// ---------------------------------------------------------------------------

#[test]
fn test_json_file_input_produces_json_stdout() {
    let bin = require_bin();

    // Write a temporary JSON input file
    let tmp_dir = std::env::temp_dir();
    let input_file = tmp_dir.join("mfs_cli_test_input.json");
    std::fs::write(&input_file, minimal_request_json()).expect("failed to write temp input file");

    let output = Command::new(&bin)
        .arg("--input")
        .arg(&input_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("failed to execute mfs_cli");

    // Clean up temp file
    let _ = std::fs::remove_file(&input_file);

    // Verify successful exit
    assert!(
        output.status.success(),
        "CLI exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify stdout is valid JSON
    let stdout = String::from_utf8(output.stdout).expect("stdout is not valid UTF-8");
    let value: serde_json::Value = serde_json::from_str(&stdout).expect("stdout is not valid JSON");

    // Verify expected top-level fields are present
    assert!(
        value.get("metadata").is_some(),
        "output should contain 'metadata' field"
    );
    assert!(
        value.get("spec").is_some(),
        "output should contain 'spec' field"
    );
    assert!(
        value.get("polynomials").is_some(),
        "output should contain 'polynomials' field"
    );
    assert!(
        value.get("matrix").is_some(),
        "output should contain 'matrix' field"
    );
}

// ---------------------------------------------------------------------------
// Test: stdin input → JSON stdout output
// ---------------------------------------------------------------------------

#[test]
fn test_stdin_input_produces_json_stdout() {
    let bin = require_bin();

    let mut child = Command::new(&bin)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn mfs_cli");

    // Write JSON to stdin
    {
        let stdin = child.stdin.as_mut().expect("failed to open stdin");
        stdin
            .write_all(minimal_request_json().as_bytes())
            .expect("failed to write to stdin");
    }

    let output = child.wait_with_output().expect("failed to wait on mfs_cli");

    // Verify successful exit
    assert!(
        output.status.success(),
        "CLI exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify stdout is valid JSON
    let stdout = String::from_utf8(output.stdout).expect("stdout is not valid UTF-8");
    let value: serde_json::Value = serde_json::from_str(&stdout).expect("stdout is not valid JSON");

    // Verify expected fields
    assert!(
        value.get("metadata").is_some(),
        "output should contain 'metadata' field"
    );
    assert!(
        value.get("polynomials").is_some(),
        "output should contain 'polynomials' field"
    );
    assert!(
        value.get("matrix").is_some(),
        "output should contain 'matrix' field"
    );
}

// ---------------------------------------------------------------------------
// Test: --stage flag produces partial output
// ---------------------------------------------------------------------------

#[test]
fn test_stage_flag_produces_partial_output() {
    let bin = require_bin();

    // Write a temporary JSON input file
    let tmp_dir = std::env::temp_dir();
    let input_file = tmp_dir.join("mfs_cli_test_stage_input.json");
    std::fs::write(&input_file, full_request_json()).expect("failed to write temp input file");

    let output = Command::new(&bin)
        .arg("--input")
        .arg(&input_file)
        .arg("--stage")
        .arg("approximation")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("failed to execute mfs_cli");

    // Clean up temp file
    let _ = std::fs::remove_file(&input_file);

    // Verify successful exit
    assert!(
        output.status.success(),
        "CLI exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify stdout is valid JSON
    let stdout = String::from_utf8(output.stdout).expect("stdout is not valid UTF-8");
    let value: serde_json::Value = serde_json::from_str(&stdout).expect("stdout is not valid JSON");

    // Verify it contains approximation-specific fields
    assert_eq!(
        value.get("stage").and_then(|v| v.as_str()),
        Some("approximation"),
        "output should indicate the 'approximation' stage"
    );
    assert!(
        value.get("order").is_some(),
        "approximation output should contain 'order'"
    );
    assert!(
        value.get("epsilon").is_some(),
        "approximation output should contain 'epsilon'"
    );

    // Verify it does NOT contain full-pipeline fields (matrix, response)
    assert!(
        value.get("matrix").is_none(),
        "partial output should not contain 'matrix' (full pipeline field)"
    );
    assert!(
        value.get("response").is_none(),
        "partial output should not contain 'response' (full pipeline field)"
    );
}

// ---------------------------------------------------------------------------
// Test: Invalid input produces structured error on stderr
// ---------------------------------------------------------------------------

#[test]
fn test_invalid_input_produces_structured_error_on_stderr() {
    let bin = require_bin();

    let mut child = Command::new(&bin)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn mfs_cli");

    // Write invalid JSON to stdin (order = 0 is invalid)
    {
        let stdin = child.stdin.as_mut().expect("failed to open stdin");
        stdin
            .write_all(invalid_request_json().as_bytes())
            .expect("failed to write to stdin");
    }

    let output = child.wait_with_output().expect("failed to wait on mfs_cli");

    // Verify non-zero exit code
    assert!(
        !output.status.success(),
        "CLI should exit with non-zero status for invalid input"
    );

    // Verify stderr contains a structured JSON error object
    let stderr = String::from_utf8(output.stderr).expect("stderr is not valid UTF-8");
    let error_value: serde_json::Value =
        serde_json::from_str(stderr.trim()).expect("stderr should contain valid JSON error object");

    // Verify error object has required fields
    assert!(
        error_value.get("error_type").is_some(),
        "error object should contain 'error_type' field"
    );
    assert!(
        error_value.get("message").is_some(),
        "error object should contain 'message' field"
    );

    // Verify the error_type is meaningful
    let error_type = error_value["error_type"].as_str().unwrap_or("");
    assert!(!error_type.is_empty(), "error_type should not be empty");
}

// ---------------------------------------------------------------------------
// Test: `optimize` subcommand — perturb → recover, end-to-end through the binary
// ---------------------------------------------------------------------------

/// Builds a JSON input string for the `optimize` subcommand: a perturbed
/// order-4 folded matrix (flat `data`) plus its recorded target response
/// (`samples`, re/im channels). Returns the perturbation factor applied.
fn build_optimize_input_json() -> String {
    let original = mfs::design::FilterDesign::prototype(4, 20.0)
        .synthesize()
        .expect("synthesize")
        .to_folded()
        .expect("fold");
    let g = mfs::freq::FrequencyGrid::linspace(-2.0, 2.0, 41).expect("grid");
    let target = mfs::response::ResponseSolver.evaluate_normalized(&original, &g).expect("target");

    // Perturb 3 inter-resonator couplings by 5 %.
    let side = original.side();
    let order = original.order();
    let mut data = original.as_slice().to_vec();
    let mut done = 0;
    'outer: for row in 1..=order {
        for col in (row + 1)..=order {
            let idx = row * side + col;
            if data[idx].abs() > 1e-9 {
                data[idx] *= 1.05;
                data[col * side + row] *= 1.05;
                done += 1;
                if done == 3 {
                    break 'outer;
                }
            }
        }
    }
    let perturbed =
        mfs::matrix::CouplingMatrix::new_with_topology(order, original.topology(), data)
            .expect("perturbed");

    let samples: Vec<serde_json::Value> = target
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
    serde_json::json!({
        "order": order,
        "topology": "folded",
        "data": perturbed.as_slice(),
        "samples": samples,
    })
    .to_string()
}

#[test]
fn test_optimize_subcommand_recovers_perturbed_matrix() {
    let bin = require_bin();

    let original = mfs::design::FilterDesign::prototype(4, 20.0)
        .synthesize()
        .expect("synthesize")
        .to_folded()
        .expect("fold");
    let input_file = std::env::temp_dir().join("mfs_cli_test_optimize_input.json");
    std::fs::write(&input_file, build_optimize_input_json()).expect("write input");

    let output = Command::new(&bin)
        .arg("optimize")
        .arg("--order")
        .arg("4")
        .arg("--max-iterations")
        .arg("500")
        .arg("--tolerance")
        .arg("1e-10")
        .arg("--lambda")
        .arg("1e-2")
        .arg("--input")
        .arg(&input_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("execute optimize");
    let _ = std::fs::remove_file(&input_file);

    assert!(
        output.status.success(),
        "optimize exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("stdout UTF-8");
    let value: serde_json::Value = serde_json::from_str(&stdout).expect("stdout JSON");

    // Physical assertions: the recovered matrix must match the original
    // (the perturbed matrix was a good guess of it) within 1e-4 relative,
    // and the optimizer must report convergence.
    assert!(value["converged"].as_bool().unwrap_or(false), "optimizer must converge");
    let side = original.side();
    let recovered_rows: Vec<Vec<f64>> =
        serde_json::from_value(value["matrix"].clone()).expect("matrix is nested rows");
    assert_eq!(
        recovered_rows.len(),
        side,
        "recovered matrix must have {side} rows"
    );
    let mut worst = 0.0_f64;
    for (i, row) in recovered_rows.iter().enumerate() {
        for (j, r) in row.iter().enumerate().skip(i) {
            let o = original.at(i, j).unwrap_or(0.0);
            if o.abs() > 1e-9 {
                worst = worst.max((o - r).abs() / o.abs());
            }
        }
    }
    assert!(worst < 1e-4, "recovered matrix too far: worst_rel={worst}");
    // The result must carry the re-evaluated response for inspection.
    assert!(value["recovered"]["samples"].is_array(), "recovered samples present");
}

/// E2E: a user-supplied matrix JSON whose off-diagonal pair is asymmetric
/// (`M[1][2]=0.6`, `M[2][1]=0.9`) must be **rejected** by the CLI before the
/// optimizer runs — it must not flow into the response solver, which assumes
/// symmetry. The binary exits non-zero with a `precondition_violation` whose
/// message names the matrix. This pins the `new_with_topology` invariant at
/// the user-input boundary (the only path that can feed asymmetric data).
#[test]
fn test_optimize_rejects_asymmetric_user_matrix() {
    let bin = require_bin();

    let input = serde_json::json!({
        "order": 2,
        "topology": "transversal",
        "data": [
            0.0, 0.8, 0.5, 0.0,
            0.8, 0.0, 0.6, 0.7,
            0.5, 0.9, 0.0, 0.4,
            0.0, 0.7, 0.4, 0.0,
        ],
        "samples": [
            { "frequency_hz": -1.0, "s11_re": 0.1, "s11_im": -0.2, "s21_re": 0.9, "s21_im": 0.0 },
            { "frequency_hz": 1.0,  "s11_re": 0.1, "s11_im": 0.2,  "s21_re": 0.9, "s21_im": 0.0 },
        ],
    });
    let input_file = std::env::temp_dir().join("mfs_cli_test_asym_input.json");
    std::fs::write(&input_file, input.to_string()).expect("write input");

    let output = Command::new(&bin)
        .arg("optimize")
        .arg("--order")
        .arg("2")
        .arg("--input")
        .arg(&input_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("execute optimize");
    let _ = std::fs::remove_file(&input_file);

    // The asymmetric matrix must be rejected before optimization.
    assert!(
        !output.status.success(),
        "an asymmetric user matrix must not reach the optimizer"
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr UTF-8");
    let err: serde_json::Value =
        serde_json::from_str(stderr.trim()).expect("stderr should be a JSON error object");
    assert_eq!(
        err["error_type"].as_str().unwrap_or(""),
        "precondition_violation",
        "an asymmetric user matrix must be a precondition violation, got: {err}"
    );
    assert!(
        err["message"].as_str().unwrap_or("").contains("asymmetric coupling matrix"),
        "the error must name the asymmetry: {err}"
    );
}

// ---------------------------------------------------------------------------
// Test: `extract` subcommand — target → recovered response, end-to-end
// ---------------------------------------------------------------------------

#[test]
fn test_extract_subcommand_recovers_response() {
    let bin = require_bin();

    let original = mfs::design::FilterDesign::prototype(4, 20.0)
        .synthesize()
        .expect("synthesize");
    let g = mfs::freq::FrequencyGrid::linspace(-2.0, 2.0, 41).expect("grid");
    let target = mfs::response::ResponseSolver.evaluate_normalized(original.matrix(), &g).expect("target");

    let samples: Vec<serde_json::Value> = target
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
    let input_file = std::env::temp_dir().join("mfs_cli_test_extract_input.json");
    std::fs::write(&input_file, serde_json::json!({ "samples": samples }).to_string()).expect("write input");

    let output = Command::new(&bin)
        .arg("extract")
        .arg("--order")
        .arg("4")
        .arg("--input")
        .arg(&input_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("execute extract");
    let _ = std::fs::remove_file(&input_file);

    assert!(
        output.status.success(),
        "extract exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("stdout UTF-8");
    let value: serde_json::Value = serde_json::from_str(&stdout).expect("stdout JSON");

    // Physical assertion: the recovered matrix must reproduce the target
    // response (all four S channels) within 1e-3 point-wise.
    assert!(value["converged"].as_bool().unwrap_or(false), "extract must converge");
    let rec_samples = &value["recovered"]["samples"];
    assert!(rec_samples.is_array(), "recovered samples present");
    let n = rec_samples.as_array().expect("samples array").len();
    assert_eq!(n, target.samples.len(), "sample count matches");
    let mut max_dev = 0.0_f64;
    for (i, s) in rec_samples.as_array().unwrap().iter().enumerate() {
        let t = &target.samples[i];
        max_dev = max_dev
            .max((s["s11_re"].as_f64().unwrap() - t.s11_re).abs())
            .max((s["s11_im"].as_f64().unwrap() - t.s11_im).abs())
            .max((s["s21_re"].as_f64().unwrap() - t.s21_re).abs())
            .max((s["s21_im"].as_f64().unwrap() - t.s21_im).abs());
    }
    assert!(max_dev < 1e-3, "recovered response too far: max_dev={max_dev}");
}
