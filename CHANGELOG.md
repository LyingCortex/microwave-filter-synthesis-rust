# Changelog

All notable changes to this project should be documented in this file.

The format is loosely based on Keep a Changelog.

## [Unreleased]

### Fixed

- **Odd-order generalized-Chebyshev synthesis produced spurious passband
  transmission zeros and destroyed return loss.** Two compounding defects in
  the E-polynomial and admittance stages:
  - `build_e_polynomial_stage`
    (`src/approx/generalized_chebyshev_helpers.rs`) now normalizes F(w) and
    P(w) to monic after the `s = jw` transform, applies the required `j`
    factor in Cameron's `E(w) = F(w)/eps_r + j·P(w)/eps`, and re-normalizes
    E(s) to monic after the `w → s` transform.
  - `synthesize_admittance_polynomials`
    (`src/synthesis/residues.rs`) now applies the global `(-1)^N` factor to
    the `E(-s)*` / `F(-s)*` mirror terms; without it the admittance
    denominator's leading coefficient cancels for odd N and the poles drift
    off the imaginary axis.
- The order-3 E(s) fixture anchor
  (`src/fixtures/mod.rs::cameron_order3_generalized_pipeline_exact_case`)
  carried a stale capture of the pre-fix un-normalized E(s) gauge (leading
  coefficient `j`); updated to the monic gauge used by the reference
  implementation.

### Added

- Open-source project metadata and repository scaffolding
- Real lossless response solving with group-delay extraction
- Generalized Chebyshev helper routines derived from the Python core

### Changed

- Approximation, matrix, and response modules were split into directory-based
  layouts
- Filter specification semantics were expanded with filter class,
  approximation family, and performance specification types
