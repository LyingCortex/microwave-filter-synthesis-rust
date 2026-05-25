# Requirements Document

## Introduction

This feature addresses the numerical failure in the MFS synthesis engine where the residue expansion algorithm (`src/synthesis/residues.rs`) fails with `NumericalFailure("y12 residue is unexpectedly complex in the current real-valued synthesis path")` for 11 test cases. The root cause is that the `real_part_if_almost_real` guard rejects residues with non-negligible imaginary parts, which legitimately arise for certain filter configurations (order 3 with 2 transmission zeros, all-pole filters of order 2 and 3). The fix extends the synthesis path to properly handle complex residue pairs, complex-valued polynomial division, and the transversal matrix construction from complex residue data.

## Glossary

- **Synthesis_Engine**: The `MatrixSynthesisEngine` struct and its associated functions that orchestrate coupling matrix construction from approximation polynomials.
- **Residue_Expansion**: The partial-fraction decomposition algorithm in `residues.rs` that decomposes Y-parameter rational functions into simple pole-residue pairs.
- **Transversal_Builder**: The `build_transversal_from_residues` function that constructs a transversal coupling matrix from pole-residue data.
- **Admittance_Polynomials**: The Y-parameter polynomial numerators (y11, y12, y22) and their shared denominator derived from the generalized Chebyshev E, F, and P polynomials.
- **Complex_Residue_Pair**: A pair of residues at conjugate poles where the residue values themselves are complex conjugates; these must be combined to produce real coupling matrix entries.
- **Coupling_Matrix**: The dense N+2 × N+2 real-valued matrix (including source and load nodes) representing the filter network.
- **DurandKerner_Solver**: The iterative root-finding algorithm (`DurandKernerRootSolver`) used to locate polynomial roots in the complex plane.
- **PolynomialSet**: The bundle of approximation polynomials (E, F, P) with metadata (eps, eps_r, order, transmission zeros) that serves as input to synthesis.

## Requirements

### Requirement 1: Complex Residue Recognition and Pairing

**User Story:** As a filter designer, I want the residue expansion to correctly identify and pair complex-conjugate residues, so that the synthesis engine can handle all valid filter configurations without numerical failure.

#### Acceptance Criteria

1. WHEN the Residue_Expansion computes a residue with imaginary magnitude exceeding 1e-6, THE Residue_Expansion SHALL identify it as part of a Complex_Residue_Pair by locating its conjugate pole and residue within the expansion.
2. WHEN two poles in the Residue_Expansion are complex conjugates (imaginary parts equal in magnitude and opposite in sign within tolerance 1e-8), THE Residue_Expansion SHALL pair their residues as a Complex_Residue_Pair.
3. THE Residue_Expansion SHALL preserve all residue data (both real and complex) without discarding imaginary components before passing results to the Transversal_Builder.
4. IF the Residue_Expansion encounters a complex residue without a matching conjugate pair, THEN THE Residue_Expansion SHALL return an `MfsError::NumericalFailure` with a descriptive message indicating the unpaired residue.

### Requirement 2: Transversal Matrix Construction from Complex Residues

**User Story:** As a filter designer, I want the transversal matrix builder to correctly derive real-valued coupling entries from complex residue pairs, so that the resulting Coupling_Matrix is physically valid.

#### Acceptance Criteria

1. WHEN the Transversal_Builder receives a Complex_Residue_Pair at conjugate poles ±jω, THE Transversal_Builder SHALL combine the pair to produce real-valued source coupling, load coupling, and diagonal entries for the corresponding resonator rows.
2. THE Transversal_Builder SHALL produce a Coupling_Matrix where all entries are real-valued (imaginary residual below 1e-10) for any valid PolynomialSet input.
3. WHEN the Transversal_Builder receives residues that are already real-valued (imaginary magnitude below 1e-6), THE Transversal_Builder SHALL use the existing real-valued extraction path without modification.
4. THE Transversal_Builder SHALL produce a Coupling_Matrix of dimension (order+2) × (order+2) with source-to-resonator and resonator-to-load couplings that are non-zero for all resonator indices.

### Requirement 3: All-Pole Filter Synthesis Support

**User Story:** As a filter designer, I want to synthesize coupling matrices for all-pole filters (no finite transmission zeros), so that basic Chebyshev filter designs work through the full synthesis pipeline.

#### Acceptance Criteria

1. WHEN a PolynomialSet has zero finite transmission zeros and order N ≥ 2, THE Synthesis_Engine SHALL produce a valid Coupling_Matrix of order N without returning a NumericalFailure error.
2. WHEN a PolynomialSet represents an all-pole filter of order 2 with 20 dB return loss, THE Synthesis_Engine SHALL produce a Coupling_Matrix where source and load couplings are symmetric (equal within tolerance 1e-6).
3. WHEN a PolynomialSet represents an all-pole filter of order 3 with 20 dB return loss, THE Synthesis_Engine SHALL produce a Coupling_Matrix where the diagonal entries of all resonators are zero (within tolerance 1e-6).

### Requirement 4: Order-3 Filter with Two Transmission Zeros

**User Story:** As a filter designer, I want to synthesize coupling matrices for order-3 filters with 2 finite transmission zeros, so that asymmetric response designs are supported.

#### Acceptance Criteria

1. WHEN a PolynomialSet has order 3 and 2 finite transmission zeros (e.g., normalized at ±1.5), THE Synthesis_Engine SHALL produce a valid Coupling_Matrix of order 3 without returning a NumericalFailure error.
2. WHEN a PolynomialSet has order 3 and 2 finite transmission zeros, THE Synthesis_Engine SHALL produce a Coupling_Matrix with non-zero source coupling at position (0,1) and non-zero load coupling at position (3,4).
3. WHEN a PolynomialSet has order 3 and 2 finite transmission zeros, THE Synthesis_Engine SHALL produce a Coupling_Matrix that can be successfully transformed to Folded topology without error.

### Requirement 5: Downstream Operation Compatibility

**User Story:** As a filter designer, I want all downstream operations (bandpass scaling, external Q extraction, triplet/quadruplet/trisection extraction) to work with matrices produced by the fixed synthesis path, so that the full design workflow is unblocked.

#### Acceptance Criteria

1. WHEN the Synthesis_Engine produces a Coupling_Matrix from a previously-failing configuration, THE Coupling_Matrix SHALL be accepted by the bandpass scaling operation (`denormalize_bandpass`) without error.
2. WHEN the Synthesis_Engine produces a Coupling_Matrix from a previously-failing configuration, THE Coupling_Matrix SHALL yield positive source and load external Q values through `denormalize_bandpass_with_external_q`.
3. WHEN the Synthesis_Engine produces a Coupling_Matrix from a filter with finite transmission zeros, THE Coupling_Matrix SHALL be accepted by triplet extraction (`synthesize_triplet`) and trisection extraction (`synthesize_trisection`) without error.

### Requirement 6: Numerical Stability of Polynomial Root Finding

**User Story:** As a filter designer, I want the root-finding algorithm to converge reliably for all denominator polynomials encountered during residue expansion, so that synthesis does not fail due to root-solver divergence.

#### Acceptance Criteria

1. WHEN the DurandKerner_Solver is applied to the denominator polynomial of the Admittance_Polynomials, THE DurandKerner_Solver SHALL converge (residual below 1e-12) within 128 iterations for all filter orders from 2 to 8.
2. IF the DurandKerner_Solver fails to converge for a denominator polynomial, THEN THE Synthesis_Engine SHALL return an `MfsError::NumericalFailure` with a message identifying the polynomial degree and maximum residual achieved.
3. WHEN the DurandKerner_Solver finds roots of the denominator polynomial, THE DurandKerner_Solver SHALL produce roots where each root satisfies |P(root)| / |leading_coefficient| < 1e-8.

### Requirement 7: Backward Compatibility with Working Configurations

**User Story:** As a filter designer, I want the existing working filter configurations (order 4 with 2 transmission zeros) to continue producing identical results after the fix, so that no regressions are introduced.

#### Acceptance Criteria

1. WHEN a PolynomialSet has order 4 and 2 finite transmission zeros at ±1.5 with 20 dB return loss, THE Synthesis_Engine SHALL produce a Coupling_Matrix numerically identical (within tolerance 1e-10) to the matrix produced by the current implementation.
2. THE Residue_Expansion SHALL produce residue values identical (within tolerance 1e-10) to the current implementation for all configurations where residues are already real-valued.
3. WHEN the Transversal_Builder processes purely real residues, THE Transversal_Builder SHALL follow the same source/load coupling sign conventions as the current implementation.

### Requirement 8: Admittance Polynomial Construction Correctness

**User Story:** As a filter designer, I want the admittance polynomial construction to correctly form y11, y12, y22 numerators and the shared denominator from E, F, P polynomials, so that the residue expansion receives mathematically correct input.

#### Acceptance Criteria

1. THE Admittance_Polynomials SHALL satisfy the degree constraint: degree(denominator) = order for all valid PolynomialSet inputs.
2. THE Admittance_Polynomials SHALL satisfy degree(y12) ≤ number of finite transmission zeros for all valid PolynomialSet inputs.
3. FOR ALL valid PolynomialSet inputs, parsing the Admittance_Polynomials into residues then reconstructing the rational function SHALL produce a polynomial that matches the original y12 numerator within tolerance 1e-8 (round-trip property).
4. WHEN the PolynomialSet has purely imaginary E(s) coefficients on odd powers and purely real coefficients on even powers, THE Admittance_Polynomials denominator SHALL have purely imaginary coefficients on odd powers and purely real coefficients on even powers.
