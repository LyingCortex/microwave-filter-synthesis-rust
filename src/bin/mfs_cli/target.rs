//! §D.5: shared target/matrix JSON helpers for `optimize` and `extract`.
//!
//! The two optimizer subcommands share `parse_response_json`,
//! `parse_matrix_json`, and `load_target`, so they live in one module.
//!
//! Pure move out of the monolithic `src/bin/mfs_cli.rs` — no behavior
//! change.
//!
use crate::io::CliRunError;
// ---------------------------------------------------------------------------
// Optimizer subcommands (D3-opt): `optimize` and `extract`
// ---------------------------------------------------------------------------

/// Parses a measured / target S-parameter response from JSON.
///
/// Accepts `{"samples":[...]}` or a bare array, each sample carrying
/// `frequency_hz` plus any of `s11_re/s11_im` and `s21_re/s21_im` (missing
/// channels default to 0). Returns a `PreconditionViolation` on malformed
/// input so the CLI reports it as a precondition error, not a crash.
pub(crate) fn parse_response_json(json: &str) -> mfs::Result<mfs::response::SParameterResponse> {
    use mfs::response::{ResponseSample, SParameterResponse};
    let v: serde_json::Value = serde_json::from_str(json)
        .map_err(|e| mfs::error::MfsError::PreconditionViolation(format!("invalid target JSON: {e}")))?;
    let arr = v
        .get("samples")
        .and_then(|s| s.as_array())
        .or_else(|| v.as_array())
        .ok_or_else(|| {
            mfs::error::MfsError::PreconditionViolation(
                "target JSON must contain a `samples` array".into(),
            )
        })?;
    let samples = arr
        .iter()
        .map(|s| {
            let get = |k: &str| s.get(k).and_then(|x| x.as_f64()).unwrap_or(0.0);
            Ok(ResponseSample {
                frequency_hz: get("frequency_hz"),
                normalized_omega: 0.0,
                group_delay: None,
                s11_re: get("s11_re"),
                s11_im: get("s11_im"),
                s21_re: get("s21_re"),
                s21_im: get("s21_im"),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(SParameterResponse { samples })
}

/// Parses a coupling matrix from JSON.
///
/// Accepts `{"order": N, "data": [row-major]}` / `{"data": [[rows]]}` /
/// `{"order": N, "matrix": [flat]}` / `{"matrix": [[rows]]}`. The topology is
/// taken from the input if present (one of `transversal`/`folded`/`arrow`),
/// else defaulted to transversal.
pub(crate) fn parse_matrix_json(json: &str) -> mfs::Result<mfs::matrix::CouplingMatrix> {
    use mfs::error::MfsError;
    use mfs::matrix::{CouplingMatrix, MatrixTopology};
    let v: serde_json::Value = serde_json::from_str(json)
        .map_err(|e| MfsError::PreconditionViolation(format!("invalid matrix JSON: {e}")))?;

    // Locate the raw matrix value: an object with "data"/"matrix", a bare
    // object, or a bare array.
    let (raw, order_hint, topo_hint) = if v.is_object() {
        (
            v.get("data")
                .cloned()
                .or_else(|| v.get("matrix").cloned()),
            v.get("order").and_then(|x| x.as_u64()).map(|n| n as usize),
            v.get("topology").and_then(|t| t.as_str()),
        )
    } else {
        (Some(v.clone()), None, None)
    };
    let raw = raw.ok_or_else(|| {
        MfsError::PreconditionViolation(
            "matrix JSON must be {{order, data}} or an array (flat row-major or nested rows)"
                .into(),
        )
    })?;

    // Accept either a flat row-major vector or nested rows.
    let rows: Vec<Vec<f64>> = if let Some(a) = raw.as_array() {
        match a.first() {
            Some(inner) if inner.is_array() => {
                // Nested rows: array of arrays.
                a.iter()
                    .map(|r| {
                        r.as_array().ok_or_else(|| {
                            MfsError::PreconditionViolation("matrix rows must be arrays".into())
                        })?
                        .iter()
                        .map(|x| x.as_f64())
                        .collect::<Option<Vec<f64>>>()
                        .ok_or_else(|| {
                            MfsError::PreconditionViolation("matrix entries must be numbers".into())
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            _ => {
                // Flat row-major vector: reshape into a square matrix.
                let flat: Vec<f64> = a
                    .iter()
                    .map(|x| x.as_f64())
                    .collect::<Option<Vec<f64>>>()
                    .ok_or_else(|| {
                        MfsError::PreconditionViolation(
                            "matrix `data` must be a flat array of numbers".into(),
                        )
                    })?;
                let side = (flat.len() as f64).sqrt().round() as usize;
                if side * side != flat.len() || side < 3 {
                    return Err(MfsError::DimensionMismatch {
                        expected: side * side,
                        actual: flat.len(),
                    });
                }
                flat.chunks_exact(side).map(|c| c.to_vec()).collect()
            }
        }
    } else {
        return Err(MfsError::PreconditionViolation(
            "matrix JSON must be an array or an object with `data`/`matrix`".into(),
        ));
    };

    let side = rows.len();
    for r in &rows {
        if r.len() != side {
            return Err(MfsError::DimensionMismatch {
                expected: side,
                actual: r.len(),
            });
        }
    }
    if side < 3 {
        return Err(MfsError::PreconditionViolation(
            "matrix side must be >= 3 (order >= 1)".into(),
        ));
    }

    let order = match order_hint {
        Some(o) => {
            if o + 2 != side {
                return Err(MfsError::DimensionMismatch {
                    expected: (o + 2) * (o + 2),
                    actual: side * side,
                });
            }
            o
        }
        None => side - 2,
    };
    let topology = match topo_hint {
        Some("folded") => MatrixTopology::Folded,
        Some("arrow") => MatrixTopology::Arrow,
        Some("transversal") | None => MatrixTopology::Transversal,
        Some(other) => {
            return Err(MfsError::PreconditionViolation(format!(
                "unknown matrix topology '{other}' (expected transversal|folded|arrow)"
            )));
        }
    };
    CouplingMatrix::new_with_topology(order, topology, rows.into_iter().flatten().collect())
}

/// Loads a target response from a Touchstone file (`--target` pointing at a
/// `.s2p`) or from the JSON string otherwise.
pub(crate) fn load_target(
    target: Option<&std::path::Path>,
    json: &str,
) -> Result<mfs::response::SParameterResponse, CliRunError> {
    match target {
        Some(path) if path.extension().and_then(|e| e.to_str()) == Some("s2p") => {
            mfs::touchstone::read_touchstone(path).map_err(|e| CliRunError::Touchstone(e.to_string()))
        }
        Some(path) => Err(CliRunError::InputParse(format!(
            "--target must be a Touchstone (.s2p) file, got '{}'",
            path.display()
        ))),
        None => parse_response_json(json).map_err(CliRunError::Pipeline),
    }
}
