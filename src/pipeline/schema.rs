//! JSON Schema generation for the synthesis pipeline's structured I/O.
//!
//! Provides [`describe_schema()`] which returns a JSON string containing the
//! JSON Schema that describes both the input format (`SynthesisRequest`) and
//! the output format (`SynthesisResponse`).

use serde_json::{json, Value};

/// Returns a JSON string containing the JSON Schema for the pipeline's
/// input (`SynthesisRequest`) and output (`SynthesisResponse`) formats.
///
/// The schema is constructed programmatically and describes:
/// - Required input fields (`order`, `return_loss_db`)
/// - Optional input fields with their types and defaults
/// - The full output structure including all artifacts and metadata
///
/// # Example
///
/// ```
/// let schema_json = mfs::pipeline::describe_schema();
/// let schema: serde_json::Value = serde_json::from_str(&schema_json).unwrap();
/// assert!(schema["input"].is_object());
/// assert!(schema["output"].is_object());
/// ```
pub fn describe_schema() -> String {
    let schema = json!({
        "$schema": "http://json-schema.org/draft-07/schema#",
        "title": "MFS Pipeline Schema",
        "description": "JSON Schema for the MFS synthesis pipeline input and output formats.",
        "version": env!("CARGO_PKG_VERSION"),
        "input": build_input_schema(),
        "output": build_output_schema()
    });

    serde_json::to_string_pretty(&schema).expect("schema serialization should never fail")
}

/// Builds the JSON Schema for `SynthesisRequest` (pipeline input).
fn build_input_schema() -> Value {
    json!({
        "type": "object",
        "title": "SynthesisRequest",
        "description": "Input specification for a filter synthesis pipeline run.",
        "required": ["order", "return_loss_db"],
        "properties": {
            "order": {
                "type": "integer",
                "minimum": 1,
                "description": "Number of resonators in the synthesized network."
            },
            "return_loss_db": {
                "type": "number",
                "exclusiveMinimum": 0.0,
                "description": "Minimum passband return loss in dB."
            },
            "transmission_zeros": {
                "type": "array",
                "items": { "type": "number" },
                "description": "Transmission zeros in normalized low-pass prototype coordinates."
            },
            "unloaded_q": {
                "type": "number",
                "exclusiveMinimum": 0.0,
                "description": "Unloaded Q factor for lossy synthesis (optional)."
            },
            "topology": {
                "type": "string",
                "enum": ["transversal", "folded", "arrow"],
                "description": "Requested output topology for the coupling matrix."
            },
            "mapping": {
                "type": "object",
                "description": "Frequency mapping configuration for physical-frequency evaluation.",
                "properties": {
                    "kind": {
                        "type": "string",
                        "enum": ["bandpass", "lowpass"],
                        "description": "Mapping type."
                    },
                    "center_hz": {
                        "type": "number",
                        "description": "Center frequency in Hz (for bandpass mappings)."
                    },
                    "bandwidth_hz": {
                        "type": "number",
                        "description": "Bandwidth in Hz (for bandpass mappings)."
                    },
                    "cutoff_hz": {
                        "type": "number",
                        "description": "Cutoff frequency in Hz (for lowpass mappings)."
                    }
                }
            },
            "grid": {
                "type": "object",
                "description": "Frequency grid configuration for response evaluation.",
                "properties": {
                    "start": {
                        "type": "number",
                        "description": "Start frequency in Hz."
                    },
                    "stop": {
                        "type": "number",
                        "description": "Stop frequency in Hz."
                    },
                    "points": {
                        "type": "integer",
                        "minimum": 2,
                        "description": "Number of evaluation points."
                    }
                },
                "required": ["start", "stop", "points"]
            }
        },
        "additionalProperties": false
    })
}

/// Builds the JSON Schema for `SynthesisResponse` (pipeline output).
fn build_output_schema() -> Value {
    json!({
        "type": "object",
        "title": "SynthesisResponse",
        "description": "Output of a completed synthesis pipeline run.",
        "required": ["metadata"],
        "properties": {
            "metadata": {
                "type": "object",
                "description": "Execution metadata and diagnostics.",
                "required": ["version", "stages_executed", "stage_timings_ms", "warnings"],
                "properties": {
                    "version": {
                        "type": "string",
                        "description": "Library version that produced this output."
                    },
                    "stages_executed": {
                        "type": "array",
                        "items": { "type": "string" },
                        "description": "Names of stages that were executed."
                    },
                    "stage_timings_ms": {
                        "type": "array",
                        "items": { "type": "number" },
                        "description": "Execution time in milliseconds for each completed stage."
                    },
                    "warnings": {
                        "type": "array",
                        "items": { "type": "string" },
                        "description": "Non-fatal warnings accumulated during execution."
                    }
                }
            },
            "spec": {
                "type": "object",
                "description": "Echo of the input specification.",
                "properties": {
                    "order": { "type": "integer" },
                    "return_loss_db": { "type": "number" },
                    "transmission_zeros": {
                        "type": "array",
                        "items": { "type": "number" }
                    }
                }
            },
            "polynomials": {
                "type": "object",
                "description": "Approximation stage output: prototype polynomials.",
                "properties": {
                    "order": { "type": "integer" },
                    "epsilon": { "type": "number" },
                    "transmission_zeros_normalized": {
                        "type": "array",
                        "items": { "type": "number" }
                    }
                }
            },
            "matrix": {
                "type": "object",
                "description": "Matrix synthesis stage output: coupling matrix.",
                "properties": {
                    "order": { "type": "integer" },
                    "topology": {
                        "type": "string",
                        "enum": ["transversal", "folded", "arrow"]
                    },
                    "data": {
                        "type": "array",
                        "items": { "type": "number" },
                        "description": "Row-major dense matrix entries including source/load."
                    }
                }
            },
            "transform": {
                "type": "object",
                "description": "Topology transform stage output.",
                "properties": {
                    "topology": {
                        "type": "string",
                        "enum": ["transversal", "folded", "arrow"]
                    },
                    "pattern_verified": { "type": "boolean" },
                    "response_invariant": { "type": "boolean" }
                }
            },
            "response": {
                "type": "object",
                "description": "Response evaluation stage output: S-parameter samples.",
                "properties": {
                    "samples": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "frequency_hz": { "type": "number" },
                                "s11_db": { "type": "number" },
                                "s21_db": { "type": "number" }
                            }
                        }
                    }
                }
            }
        }
    })
}

/// Validates a JSON value against the input schema by checking required fields
/// and basic type constraints. Returns a list of validation error messages.
///
/// This is a lightweight validator (not a full JSON Schema validator) that checks:
/// - Presence of required fields
/// - Basic type correctness for known fields
pub fn validate_request(request: &Value) -> Vec<String> {
    let mut errors = Vec::new();

    // Check top-level is an object
    let obj = match request.as_object() {
        Some(o) => o,
        None => {
            errors.push("request must be a JSON object".to_string());
            return errors;
        }
    };

    // Required fields
    if !obj.contains_key("order") {
        errors.push("missing required field: order".to_string());
    } else if let Some(order) = obj.get("order") {
        if !order.is_u64() && !order.is_i64() {
            errors.push("field 'order' must be an integer".to_string());
        } else if let Some(v) = order.as_u64() {
            if v == 0 {
                errors.push("field 'order' must be >= 1".to_string());
            }
        } else if let Some(v) = order.as_i64() {
            if v <= 0 {
                errors.push("field 'order' must be >= 1".to_string());
            }
        }
    }

    if !obj.contains_key("return_loss_db") {
        errors.push("missing required field: return_loss_db".to_string());
    } else if let Some(rl) = obj.get("return_loss_db") {
        if !rl.is_number() {
            errors.push("field 'return_loss_db' must be a number".to_string());
        } else if let Some(v) = rl.as_f64() {
            if v <= 0.0 {
                errors.push("field 'return_loss_db' must be > 0".to_string());
            }
        }
    }

    // Optional field type checks
    if let Some(tz) = obj.get("transmission_zeros") {
        if !tz.is_array() {
            errors.push("field 'transmission_zeros' must be an array".to_string());
        } else if let Some(arr) = tz.as_array() {
            for (i, item) in arr.iter().enumerate() {
                if !item.is_number() {
                    errors.push(format!(
                        "field 'transmission_zeros[{i}]' must be a number"
                    ));
                }
            }
        }
    }

    if let Some(uq) = obj.get("unloaded_q") {
        if !uq.is_number() {
            errors.push("field 'unloaded_q' must be a number".to_string());
        } else if let Some(v) = uq.as_f64() {
            if v <= 0.0 {
                errors.push("field 'unloaded_q' must be > 0".to_string());
            }
        }
    }

    if let Some(topo) = obj.get("topology") {
        if let Some(s) = topo.as_str() {
            if !["transversal", "folded", "arrow"].contains(&s) {
                errors.push(format!(
                    "field 'topology' must be one of: transversal, folded, arrow (got '{s}')"
                ));
            }
        } else {
            errors.push("field 'topology' must be a string".to_string());
        }
    }

    if let Some(grid) = obj.get("grid") {
        if !grid.is_object() {
            errors.push("field 'grid' must be an object".to_string());
        } else if let Some(g) = grid.as_object() {
            if !g.contains_key("start") {
                errors.push("field 'grid.start' is required".to_string());
            }
            if !g.contains_key("stop") {
                errors.push("field 'grid.stop' is required".to_string());
            }
            if !g.contains_key("points") {
                errors.push("field 'grid.points' is required".to_string());
            } else if let Some(p) = g.get("points") {
                if let Some(v) = p.as_u64() {
                    if v < 2 {
                        errors.push("field 'grid.points' must be >= 2".to_string());
                    }
                } else {
                    errors.push("field 'grid.points' must be an integer".to_string());
                }
            }
        }
    }

    errors
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn describe_schema_returns_valid_json() {
        let schema_str = describe_schema();
        let schema: Value =
            serde_json::from_str(&schema_str).expect("describe_schema must return valid JSON");

        // Top-level structure
        assert!(schema["input"].is_object(), "schema must have 'input' key");
        assert!(schema["output"].is_object(), "schema must have 'output' key");
        assert_eq!(schema["version"].as_str().unwrap(), env!("CARGO_PKG_VERSION"));
    }

    #[test]
    fn input_schema_declares_required_fields() {
        let schema_str = describe_schema();
        let schema: Value = serde_json::from_str(&schema_str).unwrap();

        let input = &schema["input"];
        let required = input["required"].as_array().expect("input must have 'required'");
        let required_strs: Vec<&str> = required.iter().map(|v| v.as_str().unwrap()).collect();

        assert!(required_strs.contains(&"order"));
        assert!(required_strs.contains(&"return_loss_db"));
    }

    #[test]
    fn output_schema_includes_metadata_fields() {
        let schema_str = describe_schema();
        let schema: Value = serde_json::from_str(&schema_str).unwrap();

        let output = &schema["output"];
        let metadata_props = &output["properties"]["metadata"]["properties"];

        assert!(metadata_props["version"].is_object());
        assert!(metadata_props["stages_executed"].is_object());
        assert!(metadata_props["stage_timings_ms"].is_object());
        assert!(metadata_props["warnings"].is_object());
    }

    #[test]
    fn output_schema_includes_all_artifacts() {
        let schema_str = describe_schema();
        let schema: Value = serde_json::from_str(&schema_str).unwrap();

        let output_props = &schema["output"]["properties"];

        assert!(output_props["spec"].is_object(), "output must include 'spec'");
        assert!(output_props["polynomials"].is_object(), "output must include 'polynomials'");
        assert!(output_props["matrix"].is_object(), "output must include 'matrix'");
        assert!(output_props["transform"].is_object(), "output must include 'transform'");
        assert!(output_props["response"].is_object(), "output must include 'response'");
    }

    #[test]
    fn validate_request_accepts_known_good_request() {
        let request = json!({
            "order": 4,
            "return_loss_db": 20.0,
            "transmission_zeros": [-2.0, 1.5],
            "topology": "folded",
            "mapping": {
                "kind": "bandpass",
                "center_hz": 6.75e9,
                "bandwidth_hz": 300e6
            },
            "grid": {
                "start": 6.0e9,
                "stop": 7.5e9,
                "points": 201
            }
        });

        let errors = validate_request(&request);
        assert!(
            errors.is_empty(),
            "known-good request should validate without errors, got: {errors:?}"
        );
    }

    #[test]
    fn validate_request_rejects_missing_required_fields() {
        let request = json!({});
        let errors = validate_request(&request);

        assert!(errors.iter().any(|e| e.contains("order")));
        assert!(errors.iter().any(|e| e.contains("return_loss_db")));
    }

    #[test]
    fn validate_request_rejects_invalid_order() {
        let request = json!({ "order": 0, "return_loss_db": 20.0 });
        let errors = validate_request(&request);
        assert!(errors.iter().any(|e| e.contains("order") && e.contains(">= 1")));
    }

    #[test]
    fn validate_request_rejects_invalid_return_loss() {
        let request = json!({ "order": 4, "return_loss_db": -5.0 });
        let errors = validate_request(&request);
        assert!(errors.iter().any(|e| e.contains("return_loss_db") && e.contains("> 0")));
    }

    #[test]
    fn validate_request_rejects_invalid_topology() {
        let request = json!({ "order": 4, "return_loss_db": 20.0, "topology": "invalid" });
        let errors = validate_request(&request);
        assert!(errors.iter().any(|e| e.contains("topology")));
    }

    #[test]
    fn validate_request_accepts_minimal_request() {
        let request = json!({ "order": 4, "return_loss_db": 20.0 });
        let errors = validate_request(&request);
        assert!(
            errors.is_empty(),
            "minimal request (order + return_loss_db) should be valid, got: {errors:?}"
        );
    }

    #[test]
    fn validate_request_rejects_non_object() {
        let request = json!([1, 2, 3]);
        let errors = validate_request(&request);
        assert!(errors.iter().any(|e| e.contains("JSON object")));
    }

    #[test]
    fn validate_request_rejects_invalid_grid_points() {
        let request = json!({
            "order": 4,
            "return_loss_db": 20.0,
            "grid": { "start": 1.0, "stop": 2.0, "points": 1 }
        });
        let errors = validate_request(&request);
        assert!(errors.iter().any(|e| e.contains("grid.points") && e.contains(">= 2")));
    }
}
