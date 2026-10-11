//! Closed-loop `save_context` ↔ `--resume` end-to-end gate (review §D.6 / P3 #11).
//!
//! The round-1 P3 #11 e2e (full-pipeline vs `mfs_cli --resume` matching
//! S-samples) was a temporary `examples/` harness that got deleted. This is its
//! permanent home: a black-box integration test closing the loop §D.6 asks for
//! (`save_context`/`load_context` ↔ CLI), on top of the P2-9 fix (a snapshot
//! lacks the transform/response artifacts, so resume must **re-run** those
//! stages, not skip them on metadata).
//!
//! Flow:
//! 1. Full pipeline: `mfs_cli --input <request.json>` (order-4 Folded, RL 20 dB,
//!    transmission zeros at −2.0 and 1.5) → JSON stdout with a `response` array.
//! 2. Save: `mfs::pipeline::save_context` on an in-process full run of the same
//!    request → a context-snapshot JSON file. (The CLI has no `--save` flag; the
//!    public `save_context` API is the "save" half of the loop, exactly as in the
//!    deleted P3 #11 harness and the P2-9 regression test.)
//! 3. Resume: `mfs_cli --resume <ctx.json>` (no `--stage`, i.e. the
//!    `run_remaining_stages` driver) → JSON stdout of the completed context.
//! 4. Assert the resumed `response` reproduces the full-pipeline run: all
//!    non-zero samples match within 1e-9 dB, and the transmission-zero samples
//!    are deep S21 nulls (the P3 #11 edge).
//!
//! The binary requires the `cli` feature to build, so this test is gated behind
//! `#[cfg(feature = "cli")]`. Run with:
//!
//! ```sh
//! cargo test --features cli --test cli_save_resume_e2e
//! ```

#![cfg(feature = "cli")]

use std::process::{Command, Stdio};

/// Design request: order-4 Folded prototype, RL 20 dB, transmission zeros at
/// −2.0 and 1.5. The grid (−3..3, 13 points, step 0.5) places both zeros
/// exactly on grid points, so the transmission-zero samples are evaluated at
/// the zeros themselves (the P3 #11 edge).
const REQUEST_JSON: &str = r#"{
    "order": 4,
    "return_loss_db": 20.0,
    "transmission_zeros": [-2.0, 1.5],
    "topology": "Folded",
    "grid": { "start": -3.0, "stop": 3.0, "points": 13 }
}"#;

/// The transmission-zero frequencies, matched against sample `frequency_hz` to
/// identify the null samples. Must match `REQUEST_JSON`.
const ZEROS: [f64; 2] = [-2.0, 1.5];

/// A deep S21 null threshold in dB. A transmission zero drives S21 to 0
/// (numerically a ~1e-17 residual, i.e. s21_db ≈ −310…−340); the nearest
/// non-zero sample is ≥ −30 dB. −60 dB (|S21| < 1e-3) cleanly separates the two.
const DEEP_NULL_DB: f64 = -60.0;

/// Returns the path to the `mfs_cli` binary built by cargo (same helper shape
/// as `tests/cli_integration.rs`).
fn cli_bin() -> std::path::PathBuf {
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

/// Runs `mfs_cli --input <request.json>` (full pipeline) and returns the
/// parsed stdout JSON.
fn run_full_pipeline_cli(bin: &std::path::Path, request_file: &std::path::Path) -> serde_json::Value {
    let output = Command::new(bin)
        .arg("--input")
        .arg(request_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("failed to execute mfs_cli (full pipeline)");
    assert!(
        output.status.success(),
        "full-pipeline CLI exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("stdout is not valid UTF-8");
    serde_json::from_str(&stdout).expect("full-pipeline stdout is not valid JSON")
}

/// Runs `mfs_cli --resume <ctx.json>` (no `--stage`: the `run_remaining_stages`
/// driver) and returns the parsed stdout JSON.
fn run_resume_cli(bin: &std::path::Path, ctx_file: &std::path::Path) -> serde_json::Value {
    let output = Command::new(bin)
        .arg("--resume")
        .arg(ctx_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("failed to execute mfs_cli (--resume)");
    assert!(
        output.status.success(),
        "resume CLI exited with error. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("stdout is not valid UTF-8");
    serde_json::from_str(&stdout).expect("resume stdout is not valid JSON")
}

/// True if the sample frequency is one of the design's transmission zeros.
fn is_zero_sample(frequency_hz: f64) -> bool {
    ZEROS.iter().any(|z| (frequency_hz - z).abs() < 1e-9)
}

/// Reads the `s*_db` field of a response-sample JSON object.
fn db_of(sample: &serde_json::Value, field: &str) -> f64 {
    sample
        .get(field)
        .and_then(|v| v.as_f64())
        .unwrap_or_else(|| panic!("response sample missing finite {field}: {sample}"))
}

/// Closed-loop e2e gate: full pipeline (CLI) vs save_context + `--resume`
/// (CLI). The resumed response must reproduce the full-pipeline S-samples —
/// within 1e-9 dB at every non-zero sample, and the transmission-zero samples
/// must be deep S21 nulls in both runs (the P3 #11 edge: pre-fix the nulls
/// differed ~11 dB in the serialized dB column).
#[test]
fn test_save_resume_e2e_reproduces_full_pipeline_response() {
    let bin = require_bin();

    let tmp_dir = std::env::temp_dir();
    let input_file = tmp_dir.join("mfs_cli_test_save_resume_request.json");
    let ctx_file = tmp_dir.join("mfs_cli_test_save_resume_context.json");
    std::fs::write(&input_file, REQUEST_JSON).expect("failed to write temp request file");
    let _ = std::fs::remove_file(&ctx_file);

    // (1) Full pipeline through the CLI binary.
    let full = run_full_pipeline_cli(&bin, &input_file);
    let full_response = full
        .get("response")
        .and_then(|v| v.as_array())
        .expect("full-pipeline output should contain a 'response' sample array");
    assert_eq!(full_response.len(), 13, "grid has 13 points");

    // (2) Save: the public `save_context` API on an in-process full run of the
    //     same request (the "save" half of the §D.6 closed loop — the CLI has
    //     no `--save` flag; this mirrors the deleted P3 #11 harness and the
    //     P2-9 regression test).
    let request: mfs::pipeline::SynthesisRequest =
        serde_json::from_str(REQUEST_JSON).expect("request JSON must deserialize");
    let ctx = mfs::pipeline::run_full_pipeline(request).expect("in-process full pipeline");
    mfs::pipeline::save_context(&ctx, &ctx_file).expect("save_context must succeed");

    // (3) Resume through the CLI binary (run_remaining_stages driver).
    let resumed = run_resume_cli(&bin, &ctx_file);

    // The spec (including the transmission zeros) must have round-tripped
    // through the snapshot — the zeros drive the transform pattern and the
    // response nulls.
    assert_eq!(
        resumed.get("spec").and_then(|s| s.get("transmission_zeros")),
        Some(&serde_json::json!([-2.0, 1.5])),
        "transmission zeros must survive the save/resume round trip"
    );

    // The resume driver must have re-run the stages whose artifacts the
    // snapshot does not restore (P2-9: transform + response, not skipped on
    // metadata), so all four stages are completed in the resumed context.
    assert_eq!(
        resumed
            .get("metadata")
            .and_then(|m| m.get("stages_executed")),
        Some(&serde_json::json!([
            "approximation",
            "matrix_synthesis",
            "topology_transform",
            "response_evaluation",
        ])),
        "resume must complete all four stages"
    );
    assert_eq!(
        resumed
            .get("transform")
            .and_then(|t| t.get("topology"))
            .and_then(|t| t.as_str()),
        Some("Folded"),
        "resumed context must carry the Folded transform artifact"
    );

    let resumed_response = resumed
        .get("response")
        .and_then(|v| v.as_array())
        .expect("resumed output should contain a 'response' sample array");
    assert_eq!(
        resumed_response.len(),
        full_response.len(),
        "sample count must match the full-pipeline run"
    );

    // (4) S-sample equality. At every non-zero sample the serialized dB must
    //     agree within 1e-9. At a transmission zero, the S21 magnitude is a
    //     machine-epsilon residual (~1e-17) and `20·log10` amplifies the tiny
    //     run-to-run residual into a large dB swing, so the exact dB at the
    //     null is not bit-reproducible; assert the physical property instead —
    //     a deep S21 null with near-total S11 reflection, in BOTH runs.
    let mut zero_samples_seen = 0usize;
    for (i, (f, r)) in full_response.iter().zip(resumed_response.iter()).enumerate() {
        let f_freq = f
            .get("frequency_hz")
            .and_then(|v| v.as_f64())
            .expect("frequency_hz present");
        let r_freq = r
            .get("frequency_hz")
            .and_then(|v| v.as_f64())
            .expect("frequency_hz present");
        assert_eq!(f_freq, r_freq, "sample {i}: frequency must match");

        if is_zero_sample(f_freq) {
            zero_samples_seen += 1;
            for (s, label) in [(f, "full"), (r, "resume")] {
                let s21 = db_of(s, "s21_db");
                let s11 = db_of(s, "s11_db");
                assert!(
                    s21 <= DEEP_NULL_DB,
                    "{label}: expected a deep S21 null at the transmission zero \
                     (f={f_freq}), got s21_db={s21}"
                );
                assert!(
                    s11.abs() <= 0.1,
                    "{label}: expected near-total S11 reflection at the \
                     transmission zero (f={f_freq}), got s11_db={s11}"
                );
            }
        } else {
            for field in ["s11_db", "s21_db"] {
                let f_db = db_of(f, field);
                let r_db = db_of(r, field);
                assert!(
                    (f_db - r_db).abs() <= 1e-9,
                    "sample {i} (f={f_freq}): {field} differs \
                     full={f_db} resume={r_db} (|Δ|={} > 1e-9)",
                    (f_db - r_db).abs()
                );
            }
        }
    }

    // Both transmission-zero grid points must actually have been exercised —
    // this is what makes the gate cover the P3 #11 null edge.
    assert_eq!(
        zero_samples_seen,
        2,
        "both transmission zeros must lie on the grid and be exercised"
    );

    // Clean up temp files
    let _ = std::fs::remove_file(&input_file);
    let _ = std::fs::remove_file(&ctx_file);
}
