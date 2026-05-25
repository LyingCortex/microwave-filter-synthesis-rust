# Requirements Document

## Introduction

This feature adds Touchstone (.s2p) file export to the MFS microwave filter synthesis library. Touchstone is the industry-standard format (IEEE 370 / IBIS Touchstone specification) for exchanging S-parameter data between EDA tools such as Keysight ADS, Ansys HFSS, and CST Studio Suite. The export module will serialize `SParameterResponse` data into compliant Touchstone v1.0 and v2.0 files, enabling seamless integration of synthesized filter responses into downstream simulation and verification workflows.

## Glossary

- **Touchstone_Exporter**: The module responsible for serializing S-parameter response data into Touchstone-formatted text
- **Touchstone_Parser**: The module responsible for reading Touchstone-formatted text back into S-parameter response data
- **SParameterResponse**: The existing MFS struct containing a vector of `ResponseSample` values representing the synthesized filter response
- **Option_Line**: The Touchstone header line beginning with `#` that specifies frequency unit, parameter type, data format, and reference impedance
- **Frequency_Unit**: The unit for frequency values in the Touchstone file; one of Hz, kHz, MHz, or GHz
- **Data_Format**: The representation format for S-parameter complex values; one of RI (real/imaginary), MA (magnitude/angle), or DB (dB/angle)
- **Reference_Impedance**: The characteristic impedance used for normalization, defaulting to 50 ohms
- **FilterDesign**: The existing MFS high-level API struct that synthesizes coupling matrices and evaluates S-parameter responses
- **Port_Count**: The number of ports in the network; for filters this is 2 (two-port), producing .s2p files

## Requirements

### Requirement 1: Touchstone v1.0 Export

**User Story:** As a microwave engineer, I want to export synthesized S-parameter data to Touchstone v1.0 format, so that I can import the data into any EDA tool for further analysis and verification.

#### Acceptance Criteria

1. WHEN an SParameterResponse is provided, THE Touchstone_Exporter SHALL produce a string conforming to Touchstone v1.0 specification for 2-port networks
2. THE Touchstone_Exporter SHALL write an Option_Line containing the frequency unit, parameter type (S), data format, and reference impedance
3. WHEN the Data_Format is RI, THE Touchstone_Exporter SHALL write each data row as: frequency S11_real S11_imag S21_real S21_imag S12_real S12_imag S22_real S22_imag
4. WHEN the Data_Format is MA, THE Touchstone_Exporter SHALL write each data row with magnitude (linear) and angle (degrees) pairs for each S-parameter
5. WHEN the Data_Format is DB, THE Touchstone_Exporter SHALL write each data row with magnitude (dB) and angle (degrees) pairs for each S-parameter
6. THE Touchstone_Exporter SHALL write frequency values in ascending order
7. THE Touchstone_Exporter SHALL prefix comment lines with an exclamation mark character

### Requirement 2: Touchstone v2.0 Export

**User Story:** As a microwave engineer, I want to export to Touchstone v2.0 format, so that I can use the extended metadata capabilities of the newer specification when my tools support it.

#### Acceptance Criteria

1. WHEN Touchstone v2.0 format is selected, THE Touchstone_Exporter SHALL write a `[Version] 2.0` header line
2. WHEN Touchstone v2.0 format is selected, THE Touchstone_Exporter SHALL write a `[Number of Ports] 2` declaration
3. WHEN Touchstone v2.0 format is selected, THE Touchstone_Exporter SHALL write a `[Number of Frequencies]` declaration matching the sample count
4. WHEN Touchstone v2.0 format is selected, THE Touchstone_Exporter SHALL enclose network data within `[Network Data]` and `[End]` markers

### Requirement 3: Frequency Unit Configuration

**User Story:** As a microwave engineer, I want to choose the frequency unit for the exported file, so that the data matches the conventions used by my specific EDA tool or project.

#### Acceptance Criteria

1. THE Touchstone_Exporter SHALL support Hz, kHz, MHz, and GHz as Frequency_Unit options
2. WHEN a Frequency_Unit is specified, THE Touchstone_Exporter SHALL divide the internal frequency values (stored in Hz) by the appropriate scaling factor (1, 1e3, 1e6, or 1e9)
3. WHEN no Frequency_Unit is specified, THE Touchstone_Exporter SHALL default to GHz

### Requirement 4: Data Format Configuration

**User Story:** As a microwave engineer, I want to choose the complex number representation format, so that the exported data is in the format expected by my analysis workflow.

#### Acceptance Criteria

1. THE Touchstone_Exporter SHALL support RI, MA, and DB as Data_Format options
2. WHEN no Data_Format is specified, THE Touchstone_Exporter SHALL default to RI
3. WHEN converting to MA format, THE Touchstone_Exporter SHALL compute magnitude as sqrt(re² + im²) and angle as atan2(im, re) in degrees
4. WHEN converting to DB format, THE Touchstone_Exporter SHALL compute magnitude as 20*log10(sqrt(re² + im²)) and angle as atan2(im, re) in degrees

### Requirement 5: Reference Impedance Configuration

**User Story:** As a microwave engineer, I want to specify the reference impedance, so that the exported data is normalized to the correct system impedance for my design.

#### Acceptance Criteria

1. THE Touchstone_Exporter SHALL include the Reference_Impedance value in the Option_Line
2. WHEN no Reference_Impedance is specified, THE Touchstone_Exporter SHALL default to 50 ohms
3. IF a Reference_Impedance of zero or negative value is provided, THEN THE Touchstone_Exporter SHALL return an error

### Requirement 6: Reciprocal Network Handling

**User Story:** As a microwave engineer, I want the exporter to correctly populate all four 2-port S-parameters, so that the exported file is physically valid for a reciprocal passive filter.

#### Acceptance Criteria

1. THE Touchstone_Exporter SHALL set S12 equal to S21 for all frequency points (reciprocal network assumption)
2. THE Touchstone_Exporter SHALL set S22 equal to S11 for all frequency points (symmetric network assumption)
3. WHERE a non-symmetric filter is specified, THE Touchstone_Exporter SHALL allow S22 to be provided independently from S11

### Requirement 7: File Write Support

**User Story:** As a microwave engineer, I want to write the Touchstone data directly to a file, so that I can save the exported data without manual string handling.

#### Acceptance Criteria

1. WHEN a file path is provided, THE Touchstone_Exporter SHALL write the Touchstone-formatted data to the specified file
2. THE Touchstone_Exporter SHALL use the `.s2p` file extension for 2-port data
3. IF the file cannot be written, THEN THE Touchstone_Exporter SHALL return an error describing the failure reason

### Requirement 8: Touchstone Parsing (Round-Trip)

**User Story:** As a microwave engineer, I want to parse Touchstone files back into S-parameter data, so that I can verify exported data and import reference measurements.

#### Acceptance Criteria

1. WHEN a valid Touchstone v1.0 string is provided, THE Touchstone_Parser SHALL parse the Option_Line and extract frequency unit, data format, and reference impedance
2. WHEN a valid Touchstone v1.0 string is provided, THE Touchstone_Parser SHALL parse all data rows into ResponseSample values
3. WHEN a valid Touchstone v2.0 string is provided, THE Touchstone_Parser SHALL parse the version header, port count, and network data section
4. IF an invalid or malformed Touchstone string is provided, THEN THE Touchstone_Parser SHALL return a descriptive error indicating the line number and nature of the problem
5. FOR ALL valid SParameterResponse values, exporting to Touchstone then parsing the result SHALL produce an equivalent SParameterResponse within floating-point tolerance (round-trip property)

### Requirement 9: FilterDesign API Integration

**User Story:** As a developer using MFS, I want to export Touchstone files directly from a FilterDesign object, so that I can go from synthesis to file export in a single fluent call chain.

#### Acceptance Criteria

1. WHEN a FilterDesign with band-pass parameters is available, THE FilterDesign SHALL provide a method to export the S-parameter response as a Touchstone string
2. WHEN exporting from FilterDesign, THE Touchstone_Exporter SHALL use the design's frequency range and point count as defaults
3. THE FilterDesign export method SHALL accept optional configuration for frequency unit, data format, reference impedance, and Touchstone version

### Requirement 10: Python Binding Integration

**User Story:** As a Python user of MFS, I want to export Touchstone files from Python, so that I can integrate filter synthesis into Python-based EDA workflows and Jupyter notebooks.

#### Acceptance Criteria

1. THE PyFilterDesign class SHALL expose a `to_touchstone` method that returns the Touchstone-formatted string
2. THE PyFilterDesign class SHALL expose a `save_touchstone` method that writes the Touchstone data to a file path
3. THE `to_touchstone` method SHALL accept optional keyword arguments for frequency unit, data format, reference impedance, and version
4. WHEN called from Python, THE `to_touchstone` method SHALL raise a ValueError if the design lacks band-pass parameters and no frequency range is provided

### Requirement 11: Comment and Metadata Support

**User Story:** As a microwave engineer, I want the exported Touchstone file to include descriptive comments, so that I can identify the filter design parameters when reviewing the file later.

#### Acceptance Criteria

1. THE Touchstone_Exporter SHALL write a comment header including the filter order, return loss, center frequency, and bandwidth
2. THE Touchstone_Exporter SHALL write a comment indicating the MFS library version that generated the file
3. WHERE custom comments are provided, THE Touchstone_Exporter SHALL include the custom comments in the file header
