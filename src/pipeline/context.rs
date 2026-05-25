//! Pipeline context: accumulates stage artifacts across the synthesis workflow.

use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::approx::PolynomialSet;
use crate::matrix::CouplingMatrix;
use crate::response::SParameterResponse;
use crate::spec::FilterSpec;
use crate::transform::{TopologyKind, TransformOutcome};

/// Options controlling pipeline behavior and stage configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineOptions {
    /// Requested output topology for the matrix transform stage.
    pub topology: Option<TopologyKind>,
    /// Frequency grid configuration for response evaluation.
    pub grid: Option<GridConfig>,
    /// Frequency mapping configuration (bandpass, lowpass, etc.).
    pub mapping: Option<MappingConfig>,
    /// Tolerance used for response invariance checks.
    pub response_tolerance: Option<ResponseToleranceConfig>,
}

impl Default for PipelineOptions {
    fn default() -> Self {
        Self {
            topology: None,
            grid: None,
            mapping: None,
            response_tolerance: None,
        }
    }
}

/// Frequency mapping configuration for physical-frequency evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MappingConfig {
    /// Mapping type: "bandpass", "lowpass", etc.
    pub kind: String,
    /// Center frequency in Hz (for bandpass mappings).
    pub center_hz: Option<f64>,
    /// Bandwidth in Hz (for bandpass mappings).
    pub bandwidth_hz: Option<f64>,
    /// Cutoff frequency in Hz (for lowpass mappings).
    pub cutoff_hz: Option<f64>,
}

/// Frequency grid configuration for response evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridConfig {
    /// Start frequency in Hz.
    pub start: f64,
    /// Stop frequency in Hz.
    pub stop: f64,
    /// Number of evaluation points.
    pub points: usize,
}

/// Serializable tolerance configuration for response invariance checks.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResponseToleranceConfig {
    /// Maximum allowed S11 magnitude deviation.
    pub s11_magnitude: Option<f64>,
    /// Maximum allowed S21 magnitude deviation.
    pub s21_magnitude: Option<f64>,
    /// Maximum allowed group delay deviation.
    pub group_delay: Option<f64>,
}

/// Metadata tracking pipeline execution progress and diagnostics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineMetadata {
    /// Library version that produced this context.
    pub version: String,
    /// Names of stages that have been executed so far.
    pub stages_executed: Vec<String>,
    /// Execution time in milliseconds for each completed stage.
    pub stage_timings_ms: Vec<f64>,
    /// Non-fatal warnings accumulated during pipeline execution.
    pub warnings: Vec<String>,
}

impl Default for PipelineMetadata {
    fn default() -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            stages_executed: Vec::new(),
            stage_timings_ms: Vec::new(),
            warnings: Vec::new(),
        }
    }
}

/// Central context object that accumulates all stage artifacts across the
/// synthesis pipeline. Each stage writes its output into the context without
/// consuming or invalidating earlier artifacts.
///
/// Implements `Serialize` with a custom implementation that outputs metadata,
/// spec summary, and stage completion status. Domain types that do not yet
/// implement `Serialize`/`Deserialize` are represented as summary objects.
#[derive(Debug, Clone)]
pub struct SynthesisContext {
    /// Filter specification driving this synthesis run.
    pub spec: FilterSpec,
    /// Pipeline-level options and stage configuration.
    pub options: PipelineOptions,
    /// Approximation stage output: prototype polynomials.
    pub polynomials: Option<Arc<PolynomialSet>>,
    /// Matrix synthesis stage output: canonical coupling matrix.
    pub matrix: Option<Arc<CouplingMatrix>>,
    /// Topology transform stage output.
    pub transform: Option<TransformOutcome>,
    /// Response evaluation stage output.
    pub response: Option<SParameterResponse>,
    /// Execution metadata and diagnostics.
    pub metadata: PipelineMetadata,
}

impl SynthesisContext {
    /// Creates a new pipeline context from a filter specification.
    pub fn new(spec: FilterSpec) -> Self {
        Self {
            spec,
            options: PipelineOptions::default(),
            polynomials: None,
            matrix: None,
            transform: None,
            response: None,
            metadata: PipelineMetadata::default(),
        }
    }

    /// Creates a new pipeline context with explicit options.
    pub fn with_options(spec: FilterSpec, options: PipelineOptions) -> Self {
        Self {
            spec,
            options,
            polynomials: None,
            matrix: None,
            transform: None,
            response: None,
            metadata: PipelineMetadata::default(),
        }
    }

    /// Returns a reference to the polynomial set if the approximation stage has completed.
    pub fn polynomials(&self) -> Option<&PolynomialSet> {
        self.polynomials.as_deref()
    }

    /// Returns a reference to the coupling matrix if the matrix synthesis stage has completed.
    pub fn matrix(&self) -> Option<&CouplingMatrix> {
        self.matrix.as_deref()
    }

    /// Returns a reference to the transform outcome if the topology transform stage has completed.
    pub fn transform(&self) -> Option<&TransformOutcome> {
        self.transform.as_ref()
    }

    /// Returns a reference to the S-parameter response if the response evaluation stage has completed.
    pub fn response(&self) -> Option<&SParameterResponse> {
        self.response.as_ref()
    }
}

impl Serialize for SynthesisContext {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;

        // Count fields: metadata + spec + options + optional stages
        let mut field_count = 3; // metadata, spec, options always present
        if self.polynomials.is_some() {
            field_count += 1;
        }
        if self.matrix.is_some() {
            field_count += 1;
        }
        if self.transform.is_some() {
            field_count += 1;
        }
        if self.response.is_some() {
            field_count += 1;
        }

        let mut state = serializer.serialize_struct("SynthesisContext", field_count)?;
        state.serialize_field("metadata", &self.metadata)?;

        // Serialize spec as a summary
        let spec_summary = SpecSummary {
            order: self.spec.order,
            return_loss_db: self.spec.return_loss_db,
            transmission_zeros: self.spec.transmission_zeros.iter().map(|tz| tz.value).collect(),
            unloaded_q: self.spec.unloaded_q,
        };
        state.serialize_field("spec", &spec_summary)?;
        state.serialize_field("options", &self.options)?;

        // Serialize polynomials summary if present
        if let Some(ref polys) = self.polynomials {
            let poly_summary = PolynomialsSummary {
                order: polys.order,
                epsilon: polys.eps,
                epsilon_r: polys.eps_r,
                transmission_zeros_normalized: polys.transmission_zeros_normalized.clone(),
            };
            state.serialize_field("polynomials", &poly_summary)?;
        }

        // Serialize matrix summary if present
        if let Some(ref matrix) = self.matrix {
            let matrix_summary = MatrixSummary {
                order: matrix.order(),
                side: matrix.side(),
                topology: matrix.topology(),
                data: matrix.as_slice().to_vec(),
            };
            state.serialize_field("matrix", &matrix_summary)?;
        }

        // Serialize transform summary if present
        if let Some(ref transform) = self.transform {
            let transform_summary = TransformSummary {
                topology: transform.topology,
                pattern_verified: transform.report.pattern_verified,
                notes: transform.report.notes.clone(),
            };
            state.serialize_field("transform", &transform_summary)?;
        }

        // Serialize response if present
        if let Some(ref response) = self.response {
            let samples: Vec<ResponseSampleSummary> = response
                .samples
                .iter()
                .map(|s| ResponseSampleSummary {
                    frequency_hz: s.frequency_hz,
                    s11_db: 10.0 * (s.s11_re.powi(2) + s.s11_im.powi(2)).log10(),
                    s21_db: 10.0 * (s.s21_re.powi(2) + s.s21_im.powi(2)).log10(),
                })
                .collect();
            state.serialize_field("response", &samples)?;
        }

        state.end()
    }
}

/// Serializable summary of the filter specification.
#[derive(Serialize)]
struct SpecSummary {
    order: usize,
    return_loss_db: f64,
    transmission_zeros: Vec<f64>,
    unloaded_q: Option<f64>,
}

/// Serializable summary of the polynomial set.
#[derive(Serialize)]
struct PolynomialsSummary {
    order: usize,
    epsilon: f64,
    epsilon_r: f64,
    transmission_zeros_normalized: Vec<f64>,
}

/// Serializable summary of the coupling matrix.
#[derive(Serialize)]
struct MatrixSummary {
    order: usize,
    side: usize,
    topology: TopologyKind,
    data: Vec<f64>,
}

/// Serializable summary of the transform outcome.
#[derive(Serialize)]
struct TransformSummary {
    topology: TopologyKind,
    pattern_verified: bool,
    notes: Vec<String>,
}

/// Serializable response sample with dB magnitudes.
#[derive(Serialize)]
struct ResponseSampleSummary {
    frequency_hz: f64,
    s11_db: f64,
    s21_db: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::FilterSpec;

    #[test]
    fn context_starts_with_no_artifacts() {
        let spec = FilterSpec::new(4, 20.0).expect("valid spec");
        let ctx = SynthesisContext::new(spec);

        assert!(ctx.polynomials().is_none());
        assert!(ctx.matrix().is_none());
        assert!(ctx.transform().is_none());
        assert!(ctx.response().is_none());
        assert!(ctx.metadata.stages_executed.is_empty());
    }

    #[test]
    fn pipeline_metadata_defaults_to_current_version() {
        let metadata = PipelineMetadata::default();
        assert_eq!(metadata.version, env!("CARGO_PKG_VERSION"));
        assert!(metadata.stages_executed.is_empty());
        assert!(metadata.stage_timings_ms.is_empty());
        assert!(metadata.warnings.is_empty());
    }

    #[test]
    fn pipeline_options_serializes_to_json() {
        let options = PipelineOptions {
            topology: None,
            grid: Some(GridConfig {
                start: 6.0e9,
                stop: 7.5e9,
                points: 201,
            }),
            mapping: Some(MappingConfig {
                kind: "bandpass".to_string(),
                center_hz: Some(6.75e9),
                bandwidth_hz: Some(300.0e6),
                cutoff_hz: None,
            }),
            response_tolerance: None,
        };

        let json = serde_json::to_string(&options).expect("serialization should succeed");
        let deserialized: PipelineOptions =
            serde_json::from_str(&json).expect("deserialization should succeed");

        assert_eq!(deserialized.grid.as_ref().unwrap().points, 201);
        assert_eq!(
            deserialized.mapping.as_ref().unwrap().kind,
            "bandpass"
        );
    }

    #[test]
    fn pipeline_metadata_serializes_round_trip() {
        let metadata = PipelineMetadata {
            version: "0.1.0".to_string(),
            stages_executed: vec!["approximation".to_string(), "matrix_synthesis".to_string()],
            stage_timings_ms: vec![1.2, 0.8],
            warnings: vec!["some warning".to_string()],
        };

        let json = serde_json::to_string(&metadata).expect("serialization should succeed");
        let deserialized: PipelineMetadata =
            serde_json::from_str(&json).expect("deserialization should succeed");

        assert_eq!(deserialized.version, "0.1.0");
        assert_eq!(deserialized.stages_executed.len(), 2);
        assert_eq!(deserialized.stage_timings_ms.len(), 2);
        assert_eq!(deserialized.warnings.len(), 1);
    }
}
