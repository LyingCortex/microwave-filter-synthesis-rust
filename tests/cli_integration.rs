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
    let bin_name = if cfg!(windows) { "mfs_cli.exe" } else { "mfs_cli" };

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
    let bin = cli_bin();
    if !bin.exists() {
        eprintln!("Skipping test: binary not found at {}", bin.display());
        return;
    }

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
    let value: serde_json::Value =
        serde_json::from_str(&stdout).expect("stdout is not valid JSON");

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
    let bin = cli_bin();
    if !bin.exists() {
        eprintln!("Skipping test: binary not found at {}", bin.display());
        return;
    }

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
    let value: serde_json::Value =
        serde_json::from_str(&stdout).expect("stdout is not valid JSON");

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
    let bin = cli_bin();
    if !bin.exists() {
        eprintln!("Skipping test: binary not found at {}", bin.display());
        return;
    }

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
    let value: serde_json::Value =
        serde_json::from_str(&stdout).expect("stdout is not valid JSON");

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
    let bin = cli_bin();
    if !bin.exists() {
        eprintln!("Skipping test: binary not found at {}", bin.display());
        return;
    }

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
    assert!(
        !error_type.is_empty(),
        "error_type should not be empty"
    );
}
