# Design Document: Touchstone Export

## Overview

Adds a `touchstone` module to MFS that serializes/parses S-parameter data in
IEEE 370 Touchstone format (.s2p). Integrates with `FilterDesign` and Python bindings.

## Architecture

```
src/
  touchstone/
    mod.rs          Public API: TouchstoneExporter, TouchstoneParser, enums
    writer.rs       Formatting logic (v1.0 and v2.0)
    parser.rs       Parsing logic (v1.0 and v2.0)
    types.rs        FreqUnit, DataFormat, TouchstoneConfig, TouchstoneFile
```

## Data Model

```rust
/// Frequency unit for Touchstone files.
pub enum FreqUnit { Hz, KHz, MHz, GHz }

/// Complex number representation format.
pub enum DataFormat { RI, MA, DB }

/// Touchstone file version.
pub enum TouchstoneVersion { V1, V2 }

/// Configuration for Touchstone export.
pub struct TouchstoneConfig {
    pub version: TouchstoneVersion,    // default: V1
    pub freq_unit: FreqUnit,           // default: GHz
    pub format: DataFormat,            // default: RI
    pub impedance: f64,                // default: 50.0
    pub comments: Vec<String>,         // optional user comments
}

/// Parsed Touchstone file data.
pub struct TouchstoneFile {
    pub config: TouchstoneConfig,
    pub frequencies: Vec<f64>,         // in Hz (internal)
    pub s11: Vec<(f64, f64)>,          // (re, im)
    pub s21: Vec<(f64, f64)>,
    pub s12: Vec<(f64, f64)>,
    pub s22: Vec<(f64, f64)>,
}
```

## Key Functions

```rust
// Core export
pub fn to_touchstone_string(
    response: &SParameterResponse,
    config: &TouchstoneConfig,
) -> Result<String>;

pub fn write_touchstone(
    response: &SParameterResponse,
    config: &TouchstoneConfig,
    path: impl AsRef<Path>,
) -> Result<()>;

// Core parse
pub fn parse_touchstone(content: &str) -> Result<TouchstoneFile>;

// FilterDesign integration
impl FilterDesign {
    pub fn to_touchstone(&self) -> TouchstoneBuilder;
    pub fn save_touchstone(&self, path: impl AsRef<Path>) -> Result<()>;
}

// Builder for configuration
pub struct TouchstoneBuilder { ... }
impl TouchstoneBuilder {
    pub fn freq_unit(self, unit: FreqUnit) -> Self;
    pub fn format(self, fmt: DataFormat) -> Self;
    pub fn impedance(self, z0: f64) -> Self;
    pub fn version(self, v: TouchstoneVersion) -> Self;
    pub fn comment(self, text: impl Into<String>) -> Self;
    pub fn build(self) -> Result<String>;
    pub fn save(self, path: impl AsRef<Path>) -> Result<()>;
}
```

## Touchstone v1.0 Format

```
! MFS v0.2.0 - 6th order bandpass filter
! Center: 6.75 GHz, BW: 300 MHz, RL: 23 dB
# GHz S RI R 50
6.000000  -0.9500  0.0200  0.0100  0.3000  0.0100  0.3000  -0.9500  0.0200
6.100000  -0.8800  0.0500  0.0300  0.4500  0.0300  0.4500  -0.8800  0.0500
...
```

## Touchstone v2.0 Format

```
[Version] 2.0
! MFS v0.2.0 - 6th order bandpass filter
[Number of Ports] 2
[Number of Frequencies] 201
# GHz S RI R 50
[Network Data]
6.000000  -0.9500  0.0200  0.0100  0.3000  0.0100  0.3000  -0.9500  0.0200
...
[End]
```

## Reciprocal Network Handling

For symmetric filters (default): S12 = S21, S22 = S11.
For asymmetric filters: S22 computed from response if available, otherwise = S11.

## Python Integration

```python
# PyFilterDesign methods
def to_touchstone(self, freq_unit="GHz", format="RI", impedance=50.0, version=1):
    """Returns Touchstone-formatted string."""

def save_touchstone(self, path, freq_unit="GHz", format="RI", impedance=50.0, version=1):
    """Writes Touchstone file to disk."""
```

## Error Handling

- Invalid impedance (≤ 0): `MfsError::PreconditionViolation`
- File write failure: `MfsError::PreconditionViolation` with OS error message
- Parse failure: `MfsError::PreconditionViolation` with line number and description
- Missing frequency data: `MfsError::PreconditionViolation`
