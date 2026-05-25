# Implementation Plan: Synthesis Numerical Fix

## Overview

This plan implements the fix for the numerical failure in the residue expansion algorithm. The approach adds residue classification, conjugate pairing, and a new transversal builder path for complex residue pairs while preserving the existing real-residue path unchanged. Implementation is localized to `src/synthesis/residues.rs` with property-based tests in `tests/synthesis_properties.rs` and regression tests in `tests/synthesis_regression.rs`.

## Tasks

- [x] 1. Add residue classification and conjugate pairing infrastructure
  - [x] 1.1 Implement `ResidueClassification` enum and `ClassifiedResidues` struct
    - Add `ResidueClassification` enum with `Real { index }` and `ComplexPair { index_a, index_b }` variants
    - Add `ClassifiedResidues` struct with `real_indices: Vec<usize>` and `conjugate_pairs: Vec<(usize, usize)>` fields
    - Implement `classify_residues` function that iterates residues, detects complex ones (imaginary magnitude > 1e-6), matches conjugate poles (within 1e-8), and returns error for unpaired residues
    - _Requirements: 1.1, 1.2, 1.4_

  - [x] 1.2 Write property test for residue classification (Property 1)
    - **Property 1: Complex conjugate residues are correctly paired**
    - **Validates: Requirements 1.1, 1.2**
    - Generate random `FilterSpec` (order 2–8, return loss 15–25 dB, 0 to order-1 TZs with magnitudes in [1.1, 3.0])
    - Assert that every complex residue is paired with exactly one conjugate partner satisfying `|Im(p_a) + Im(p_b)| < 1e-8` and `|Re(p_a) - Re(p_b)| < 1e-8`

  - [x] 1.3 Implement `combine_conjugate_pair` function
    - Implement the function that takes a complex pole and three complex residues (y11, y12, y22) and produces four real values (diagonal, source_coupling, load_coupling, cross_coupling)
    - Diagonal entry: `-Im(pole)` (pole location on imaginary axis)
    - Source coupling: `sqrt(2 * Re(r11))`
    - Load coupling: derived from `r12 / sqrt(r11)` or `sqrt(r22)` depending on magnitudes
    - Return error if combined result has non-real residual above 1e-10
    - _Requirements: 2.1, 2.2_

- [x] 2. Modify `build_transversal_from_residues` to support complex residue pairs
  - [x] 2.1 Integrate classification into the transversal builder
    - Call `classify_residues` at the start of `build_transversal_from_residues`
    - Route real residues through the existing extraction logic (preserving `real_part_if_almost_real` for those)
    - Route complex-conjugate pairs through `combine_conjugate_pair`
    - Ensure the builder handles both paths producing correct matrix entries at the right indices
    - _Requirements: 1.3, 2.1, 2.3, 2.4_

  - [x] 2.2 Write property test for transversal matrix realness (Property 3)
    - **Property 3: Transversal matrix entries are real-valued**
    - **Validates: Requirements 2.1, 2.2**
    - Generate random `FilterSpec` (order 2–8, 0 to order-1 TZs)
    - Assert that `build_transversal_from_residues` produces a `CouplingMatrix` where every entry is real (stored as `f64`)

  - [x] 2.3 Write property test for matrix structural invariants (Property 4)
    - **Property 4: Matrix structural invariants**
    - **Validates: Requirements 2.4, 4.2**
    - Generate random `FilterSpec` (order 2–8)
    - Assert matrix dimension is (N+2)×(N+2), source-to-resonator and resonator-to-load couplings each have at least one non-zero value > 1e-12

- [x] 3. Checkpoint - Verify core logic compiles and existing tests pass
  - Ensure all tests pass, ask the user if questions arise.

- [x] 4. Add unit tests for specific filter configurations
  - [x] 4.1 Write unit tests for all-pole and order-3 with 2 TZ configurations
    - Order-2 all-pole 20 dB: verify symmetric source/load couplings (within 1e-6)
    - Order-3 all-pole 20 dB: verify zero diagonal entries (within 1e-6)
    - Order-3 with 2 TZs at ±1.5: verify synthesis succeeds and produces valid matrix
    - Add tests in `src/synthesis/residues.rs` `#[cfg(test)]` module
    - _Requirements: 3.1, 3.2, 3.3, 4.1, 4.2_

  - [x] 4.2 Write property test for synthesis success across configurations (Property 5)
    - **Property 5: Synthesis succeeds for all valid filter configurations**
    - **Validates: Requirements 3.1, 4.1**
    - Generate random `FilterSpec` (order 2–8, return loss 15–25 dB, 0 to N-1 TZs with magnitudes in [1.1, 3.0])
    - Assert `MatrixSynthesisEngine::synthesize` returns `Ok(CouplingMatrix)` without `NumericalFailure`

  - [x] 4.3 Write property test for topology transformation (Property 6)
    - **Property 6: Topology transformation succeeds on synthesized matrices**
    - **Validates: Requirements 4.3**
    - Generate random `FilterSpec` (order 3–8, at least 1 TZ)
    - Assert synthesized matrix can be transformed to `MatrixTopology::Folded` without error

- [x] 5. Add backward compatibility and admittance polynomial tests
  - [x] 5.1 Write backward compatibility unit tests
    - Order-4 with 2 TZs at ±1.5, 20 dB: capture reference matrix values from current implementation and assert new implementation matches within 1e-10
    - Verify real-residue path selection: assert existing path is used when residues are real
    - Verify coupling sign conventions match current implementation
    - Add tests in `src/synthesis/residues.rs` `#[cfg(test)]` module
    - _Requirements: 7.1, 7.2, 7.3_

  - [x] 5.2 Write property test for backward compatibility (Property 8)
    - **Property 8: Backward compatibility for real-residue configurations**
    - **Validates: Requirements 7.1, 7.2, 7.3**
    - Generate random `FilterSpec` (order ≥ 4, 1 to order-2 TZs) — configurations that produce real residues
    - Assert the fixed path produces a matrix numerically identical (within 1e-10) to the current implementation

  - [x] 5.3 Write property test for admittance polynomial degree constraints (Property 9)
    - **Property 9: Admittance polynomial degree constraints**
    - **Validates: Requirements 8.1, 8.2**
    - Generate random `FilterSpec` (order 2–8, 0 to order-1 TZs)
    - Assert `degree(denominator) == N` and `degree(y12) <= K` (number of finite TZs)

  - [x] 5.4 Write property test for admittance polynomial coefficient parity (Property 10)
    - **Property 10: Admittance polynomial coefficient parity**
    - **Validates: Requirements 8.4**
    - Generate random `FilterSpec` where E(s) has the expected parity structure
    - Assert denominator polynomial preserves parity: purely imaginary coefficients on odd powers, purely real on even powers

- [x] 6. Checkpoint - Ensure all unit and property tests pass
  - Ensure all tests pass, ask the user if questions arise.

- [x] 7. Add integration and regression tests
  - [x] 7.1 Create `tests/synthesis_properties.rs` with remaining property tests
    - **Property 2: Residue data preservation through the pipeline**
    - **Validates: Requirements 1.3, 8.3**
    - Generate random `FilterSpec`, compute residue expansion, reconstruct rational function at 10 random imaginary-axis points, assert match within 1e-8
    - **Property 7: Root solver accuracy**
    - **Validates: Requirements 6.1, 6.3**
    - Generate random `FilterSpec` (order 2–8), compute admittance denominator, find roots, assert `|P(root)| / |leading_coefficient| < 1e-8` for every root
    - Add `proptest` as a dev-dependency in `Cargo.toml` if not already present
    - _Requirements: 1.3, 6.1, 6.3, 8.3_

  - [x] 7.2 Create `tests/synthesis_regression.rs` with regression tests for previously-failing configurations
    - Test all 11 previously-failing configurations through the full synthesis pipeline
    - Verify bandpass scaling succeeds on each (`denormalize_bandpass`)
    - Verify external Q extraction produces positive values (`denormalize_bandpass_with_external_q`)
    - Verify triplet/trisection extraction works for configs with TZs
    - _Requirements: 5.1, 5.2, 5.3_

  - [x] 7.3 Write unit test for unpaired complex residue error
    - Construct a synthetic scenario where a complex residue has no conjugate match
    - Assert `MfsError::NumericalFailure` is returned with descriptive message
    - _Requirements: 1.4_

  - [x] 7.4 Write unit test for root solver convergence error format
    - Verify error message format when root solver fails to converge (if testable via mock or extreme input)
    - _Requirements: 6.2_

- [x] 8. Final checkpoint - Ensure all tests pass
  - Ensure all tests pass, ask the user if questions arise.

## Notes

- Tasks marked with `*` are optional and can be skipped for faster MVP
- Each task references specific requirements for traceability
- Checkpoints ensure incremental validation
- Property tests validate universal correctness properties from the design document
- Unit tests validate specific examples and edge cases
- The `proptest` crate must be added as a dev-dependency in `Cargo.toml`
- All property tests should use minimum 100 iterations and be tagged with `// Feature: synthesis-numerical-fix, Property {N}: {title}`
- The implementation is localized to `src/synthesis/residues.rs` — no other source files need modification

## Task Dependency Graph

```json
{
  "waves": [
    { "id": 0, "tasks": ["1.1"] },
    { "id": 1, "tasks": ["1.2", "1.3"] },
    { "id": 2, "tasks": ["2.1"] },
    { "id": 3, "tasks": ["2.2", "2.3", "4.1"] },
    { "id": 4, "tasks": ["4.2", "4.3", "5.1"] },
    { "id": 5, "tasks": ["5.2", "5.3", "5.4"] },
    { "id": 6, "tasks": ["7.1", "7.2"] },
    { "id": 7, "tasks": ["7.3", "7.4"] }
  ]
}
```
