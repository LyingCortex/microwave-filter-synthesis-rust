# First GitHub Release Notes Draft

## Summary

This release publishes the initial open-source version of `mfs`, a Rust core
library for microwave filter synthesis.

The current release focuses on architecture, typed domain modeling, early
generalized Chebyshev helper support, coupling-matrix artifacts, and lossless
response evaluation.

## Included In This Release

- typed filter specifications and frequency-plan APIs
- transmission-zero normalization helpers
- validated polynomial and coupling-matrix artifacts
- generalized Chebyshev helper routines derived from the Python prototype
- lossless response solving with group-delay extraction
- examples, tests, and contributor-facing project scaffolding

## Current Limitations

- the main approximation path still uses placeholder `e/f` coefficient
  generation in some cases
- generalized Chebyshev helper outputs are not yet fully integrated into the
  primary approximation artifact model
- topology transforms are not implemented yet
- coupling-matrix recovery and optimization stages are still in progress

## Intended Audience

This release is best suited for:

- contributors
- researchers exploring the architecture
- users evaluating the Rust rewrite direction

It should not yet be presented as a fully complete production synthesis stack.
