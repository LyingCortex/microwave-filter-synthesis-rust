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
                "enum": ["Transversal", "Folded", "Arrow", "Wheel"],
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
                        "items": { "type": ["number", "null"] },
                        "description": "Normalized zeros; a zero at infinity is null."
                    },
                    "unloaded_q": {
                        "type": ["number", "null"],
                        "description": "Unloaded Q (finite only; lossless is omitted)."
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
                        "items": { "type": ["number", "null"] },
                        "description": "Normalized zeros; a zero at infinity is null."
                    }
                }
            },
            "matrix": {
                "type": "object",
                "description": "Matrix synthesis stage output: coupling matrix.",
                "properties": {
                    "order": { "type": "integer" },
                    "side": { "type": "integer" },
                    "topology": {
                        "type": "string",
                        "enum": ["Transversal", "Folded", "Arrow", "Wheel"]
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
                        "enum": ["Transversal", "Folded", "Arrow", "Wheel"]
                    },
                    "pattern_verified": { "type": "boolean" },
                    "notes": {
                        "type": "array",
                        "items": { "type": "string" },
                        "description": "Notes describing the transform backend used."
                    }
                }
            },
            "response": {
                "type": "array",
                "description": "Response evaluation stage output: S-parameter samples \
                                 (a direct array of samples, matching the \
                                 SynthesisContext serializer).",
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
    })
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
        assert!(
            schema["output"].is_object(),
            "schema must have 'output' key"
        );
        assert_eq!(
            schema["version"].as_str().unwrap(),
            env!("CARGO_PKG_VERSION")
        );
    }

    #[test]
    fn input_schema_declares_required_fields() {
        let schema_str = describe_schema();
        let schema: Value = serde_json::from_str(&schema_str).unwrap();

        let input = &schema["input"];
        let required = input["required"]
            .as_array()
            .expect("input must have 'required'");
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

        assert!(
            output_props["spec"].is_object(),
            "output must include 'spec'"
        );
        assert!(
            output_props["polynomials"].is_object(),
            "output must include 'polynomials'"
        );
        assert!(
            output_props["matrix"].is_object(),
            "output must include 'matrix'"
        );
        assert!(
            output_props["transform"].is_object(),
            "output must include 'transform'"
        );
        assert!(
            output_props["response"].is_object(),
            "output must include 'response'"
        );
    }

}
