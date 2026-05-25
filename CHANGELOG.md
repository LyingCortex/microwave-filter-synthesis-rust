# Changelog

## [0.1.0] — 2025-05-25

### Added

- `FilterDesign` high-level API with `bandpass()` and `prototype()` entry points
- Generalized Chebyshev (Cameron) polynomial synthesis for arbitrary order
- Coupling matrix synthesis via residue expansion
- Topology transforms: Transversal → Folded, Arrow
- S-parameter response evaluation with two backends:
  - LU-based solver (exact, O(N³) per frequency point)
  - Pole-expansion solver (fast, O(N) per frequency point, auto-verified)
- Adaptive root-finding with three-level fallback:
  - Durand-Kerner (fastest, degree ≤ ~28)
  - Aberth-Ehrlich (cubic convergence, degree 25-40)
  - Companion Matrix eigenvalues (most robust, any degree)
- In-place Givens similarity rotations (O(N) per rotation vs O(N³))
- Band-pass scaling and external Q computation
- Section extraction: triplet, quadruplet, trisection
- Python bindings via PyO3 (optional `python` feature)
- 164 unit and integration tests
- High-order stability verified up to order 30

### Performance

- Matrix topology transforms: ~500× faster than naive implementation
- Frequency response (pole expansion): ~100-500× faster than LU for large sweeps
- Numerical precision: machine epsilon (2.2e-16) symmetry error at order 30
