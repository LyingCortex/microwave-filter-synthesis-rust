//! Touchstone (.s2p) file export and import.
//!
//! Provides export of `SParameterResponse` to industry-standard Touchstone v1.0
//! format, and import via the `touchstone` crate.

use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::Path;

use crate::error::{MfsError, Result};
use crate::response::{ResponseSample, SParameterResponse};

/// Frequency unit for Touchstone export.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FreqUnit { Hz, KHz, MHz, GHz }

impl FreqUnit {
    fn label(&self) -> &'static str {
        match self { Self::Hz => "HZ", Self::KHz => "KHZ", Self::MHz => "MHZ", Self::GHz => "GHZ" }
    }
    fn divisor(&self) -> f64 {
        match self { Self::Hz => 1.0, Self::KHz => 1e3, Self::MHz => 1e6, Self::GHz => 1e9 }
    }
}

/// Data format for S-parameter values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataFormat { RI, MA, DB }

impl DataFormat {
    fn label(&self) -> &'static str {
        match self { Self::RI => "RI", Self::MA => "MA", Self::DB => "DB" }
    }
}

/// Touchstone file version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TouchstoneVersion { V1, V2 }

/// Configuration for Touchstone export.
#[derive(Debug, Clone)]
pub struct TouchstoneConfig {
    pub freq_unit: FreqUnit,
    pub format: DataFormat,
    pub impedance: f64,
    pub version: TouchstoneVersion,
    pub comments: Vec<String>,
}

impl Default for TouchstoneConfig {
    fn default() -> Self {
        Self {
            freq_unit: FreqUnit::GHz,
            format: DataFormat::RI,
            impedance: 50.0,
            version: TouchstoneVersion::V1,
            comments: Vec::new(),
        }
    }
}

/// Formats an `SParameterResponse` as a Touchstone string (v1.0 or v2.0).
///
/// Assumes symmetric reciprocal 2-port (S12=S21, S22=S11).
pub fn to_touchstone_string(response: &SParameterResponse, config: &TouchstoneConfig) -> Result<String> {
    if response.samples.is_empty() {
        return Err(MfsError::PreconditionViolation("cannot export empty response".into()));
    }
    if config.impedance <= 0.0 {
        return Err(MfsError::PreconditionViolation(
            format!("impedance must be > 0, got {}", config.impedance),
        ));
    }

    let mut out = String::new();

    // V2.0 header
    if config.version == TouchstoneVersion::V2 {
        writeln!(out, "[Version] 2.0").unwrap();
    }

    // Comments
    for c in &config.comments {
        writeln!(out, "! {c}").unwrap();
    }

    // V2.0 port/frequency declarations
    if config.version == TouchstoneVersion::V2 {
        writeln!(out, "[Number of Ports] 2").unwrap();
        writeln!(out, "[Number of Frequencies] {}", response.samples.len()).unwrap();
    }

    // Option line
    writeln!(out, "# {} S {} R {}", config.freq_unit.label(), config.format.label(), config.impedance).unwrap();

    // V2.0 network data marker
    if config.version == TouchstoneVersion::V2 {
        writeln!(out, "[Network Data]").unwrap();
    }

    // Data lines
    let div = config.freq_unit.divisor();
    for s in &response.samples {
        let freq = s.frequency_hz / div;
        let line = format_data_line(freq, s, config.format);
        writeln!(out, "{line}").unwrap();
    }

    // V2.0 end marker
    if config.version == TouchstoneVersion::V2 {
        writeln!(out, "[End]").unwrap();
    }

    Ok(out)
}

/// Writes an `SParameterResponse` to a Touchstone .s2p file.
pub fn write_touchstone(
    response: &SParameterResponse,
    config: &TouchstoneConfig,
    path: impl AsRef<Path>,
) -> Result<()> {
    let content = to_touchstone_string(response, config)?;
    fs::write(path, content).map_err(|e| MfsError::PreconditionViolation(
        format!("failed to write Touchstone file: {e}"),
    ))
}

/// Reads a Touchstone file using the `touchstone` crate and returns S-parameter data.
pub fn read_touchstone(path: impl AsRef<Path>) -> Result<SParameterResponse> {
    let path_str = path.as_ref().to_str()
        .ok_or_else(|| MfsError::PreconditionViolation("invalid file path".into()))?;
    let network = touchstone::Network::new(path_str.to_string());

    if network.rank != 2 {
        return Err(MfsError::PreconditionViolation(
            format!("expected 2-port network, got {}-port", network.rank),
        ));
    }

    let s11_data = network.s_ri(1, 1);
    let s21_data = network.s_ri(2, 1);

    let samples = s11_data.iter().zip(s21_data.iter())
        .map(|(s11, s21)| ResponseSample {
            frequency_hz: s11.frequency,
            normalized_omega: 0.0,
            group_delay: 0.0,
            s11_re: s11.s_ri.0,
            s11_im: s11.s_ri.1,
            s21_re: s21.s_ri.0,
            s21_im: s21.s_ri.1,
        })
        .collect();

    Ok(SParameterResponse { samples })
}

fn format_data_line(freq: f64, s: &ResponseSample, fmt: DataFormat) -> String {
    match fmt {
        DataFormat::RI => {
            // S11_re S11_im S21_re S21_im S12_re S12_im S22_re S22_im
            format!("{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}",
                freq,
                s.s11_re, s.s11_im,
                s.s21_re, s.s21_im,
                s.s21_re, s.s21_im,  // S12 = S21 (reciprocal)
                s.s11_re, s.s11_im,  // S22 = S11 (symmetric)
            )
        }
        DataFormat::MA => {
            let s11_mag = s.s11_mag();
            let s11_ang = s.s11_im.atan2(s.s11_re).to_degrees();
            let s21_mag = s.s21_mag();
            let s21_ang = s.s21_im.atan2(s.s21_re).to_degrees();
            format!("{:.9e} {:.9e} {:.4} {:.9e} {:.4} {:.9e} {:.4} {:.9e} {:.4}",
                freq,
                s11_mag, s11_ang,
                s21_mag, s21_ang,
                s21_mag, s21_ang,
                s11_mag, s11_ang,
            )
        }
        DataFormat::DB => {
            let s11_db = s.s11_db();
            let s11_ang = s.s11_im.atan2(s.s11_re).to_degrees();
            let s21_db = s.s21_db();
            let s21_ang = s.s21_im.atan2(s.s21_re).to_degrees();
            format!("{:.9e} {:.6} {:.4} {:.6} {:.4} {:.6} {:.4} {:.6} {:.4}",
                freq,
                s11_db, s11_ang,
                s21_db, s21_ang,
                s21_db, s21_ang,
                s11_db, s11_ang,
            )
        }
    }
}
