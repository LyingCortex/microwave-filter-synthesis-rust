use std::fs;
use std::path::PathBuf;

use serde::Deserialize;

use crate::error::{MfsError, Result};
use crate::freq::{BandPassMapping, FrequencyGrid};
use crate::spec::{FilterSpec, TransmissionZero};

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct FilterDatabaseDocument {
    pub case: Vec<FilterDatabaseCase>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct FilterDatabaseCase {
    pub case_id: String,
    pub specification: CaseSpecification,
    pub mathematical_model: Option<MathematicalModel>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct CaseSpecification {
    pub filter_order: usize,
    pub center_freq: ScalarWithUnit,
    pub start_freq: ScalarWithUnit,
    pub stop_freq: ScalarWithUnit,
    pub return_loss: ScalarWithUnit,
    pub normalized_transmission_zeros: Vec<NormalizedTransmissionZero>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct ScalarWithUnit {
    pub value: f64,
    pub unit: String,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct NormalizedTransmissionZero {
    pub value: ComplexValue,
    pub unit: String,
    pub r#type: String,
    pub domain: String,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct ComplexValue {
    pub re: f64,
    pub im: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct MathematicalModel {
    pub polynomial_coefficients: PolynomialCoefficients,
    pub singularities: Option<Singularities>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct PolynomialCoefficients {
    #[serde(rename = "E")]
    pub e: Vec<ComplexValue>,
    #[serde(rename = "F")]
    pub f: Vec<ComplexValue>,
    #[serde(rename = "P")]
    pub p: Vec<ComplexValue>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct Singularities {
    #[serde(rename = "epsilon_R")]
    pub epsilon_r: Option<f64>,
    pub epsilon: Option<f64>,
    pub reflection_zeros: Option<Vec<ComplexValue>>,
    pub reflection_poles: Option<Vec<ComplexValue>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EndToEndFixture {
    pub case_id: String,
    pub spec: FilterSpec,
    pub mapping: BandPassMapping,
    pub grid: FrequencyGrid,
}

impl ScalarWithUnit {
    pub fn as_hz(&self) -> Result<f64> {
        let scale = match self.unit.as_str() {
            "Hz" => 1.0,
            "kHz" => 1.0e3,
            "MHz" => 1.0e6,
            "GHz" => 1.0e9,
            other => {
                return Err(MfsError::Unsupported(format!(
                    "unsupported frequency unit in fixture: {other}"
                )))
            }
        };
        Ok(self.value * scale)
    }

    pub fn as_db(&self) -> Result<f64> {
        if self.unit == "dB" {
            Ok(self.value)
        } else {
            Err(MfsError::Unsupported(format!(
                "unsupported return-loss unit in fixture: {}",
                self.unit
            )))
        }
    }
}

impl FilterDatabaseCase {
    pub fn to_end_to_end_fixture(&self) -> Result<EndToEndFixture> {
        let center_hz = self.specification.center_freq.as_hz()?;
        let start_hz = self.specification.start_freq.as_hz()?;
        let stop_hz = self.specification.stop_freq.as_hz()?;
        let return_loss_db = self.specification.return_loss.as_db()?;
        let prototype_order = self.prototype_order();
        let transmission_zeros = self
            .specification
            .normalized_transmission_zeros
            .iter()
            .map(|zero| {
                if zero.value.re.abs() > 1.0e-12 {
                    return Err(MfsError::Unsupported(
                        "fixture loader currently expects purely imaginary normalized transmission zeros"
                            .to_string(),
                    ));
                }
                Ok(TransmissionZero::normalized(zero.value.im))
            })
            .collect::<Result<Vec<_>>>()?;

        let spec =
            FilterSpec::generalized_chebyshev(prototype_order, return_loss_db)?
                .with_transmission_zeros(transmission_zeros);
        let mapping = BandPassMapping::new(center_hz, stop_hz - start_hz)?;
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, 21)?;

        Ok(EndToEndFixture {
            case_id: self.case_id.clone(),
            spec,
            mapping,
            grid,
        })
    }

    pub fn prototype_order(&self) -> usize {
        self.mathematical_model
            .as_ref()
            .map(|model| model.polynomial_coefficients.f.len().saturating_sub(1))
            .filter(|order| *order > 0)
            .unwrap_or(self.specification.filter_order)
    }
}

pub fn load_filter_database_document(path: impl Into<PathBuf>) -> Result<FilterDatabaseDocument> {
    let path = path.into();
    let raw = fs::read_to_string(&path).map_err(|error| {
        MfsError::Unsupported(format!(
            "failed to read fixture file {}: {error}",
            path.display()
        ))
    })?;

    serde_json::from_str(&raw).map_err(|error| {
        MfsError::Unsupported(format!(
            "failed to deserialize fixture file {}: {error}",
            path.display()
        ))
    })
}

pub fn load_filter_database_case(
    path: impl Into<PathBuf>,
    case_id: &str,
) -> Result<FilterDatabaseCase> {
    let document = load_filter_database_document(path)?;
    document
        .case
        .into_iter()
        .find(|case| case.case_id == case_id)
        .ok_or_else(|| {
            MfsError::Unsupported(format!("fixture case not found in filter_database.json: {case_id}"))
        })
}

pub fn load_filter_database_case_from_repo(case_id: &str) -> Result<FilterDatabaseCase> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("filter_database.json");
    load_filter_database_case(path, case_id)
}

pub fn load_filter_database_end_to_end_fixture(case_id: &str) -> Result<EndToEndFixture> {
    load_filter_database_case_from_repo(case_id)?.to_end_to_end_fixture()
}
