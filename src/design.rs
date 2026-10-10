//! High-level filter design API.
//!
//! `FilterDesign` is the primary entry point for this library.
//!
//! # Quick Start
//!
//! ```rust
//! use mfs::prelude::*;
//!
//! // Band-pass filter: 6 resonators, 23 dB RL, center 6.75 GHz, BW 300 MHz
//! let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
//!     .zeros_hz([6.4e9, 6.5e9, 7.0e9])
//!     .synthesize()?;
//!
//! // Normalized prototype (no physical frequencies)
//! let design = FilterDesign::prototype(4, 20.0)
//!     .zeros([-1.5, 1.5])
//!     .synthesize()?;
//! # Ok::<(), MfsError>(())
//! ```

use crate::approx::{generalized_chebyshev_polynomials, PolynomialSet};
use crate::error::Result;
use crate::freq::{BandPassMapping, FrequencyGrid, FrequencyMapping};
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::response::{ResponseSolver, SParameterResponse};
use crate::spec::FilterSpec;
use crate::synthesis::{build_filter_spec, MatrixSynthesisEngine};

// ═══════════════════════════════════════════════════════════════════════════════
// Builder
// ═══════════════════════════════════════════════════════════════════════════════

/// Builder for specifying and synthesizing a filter design.
#[derive(Debug, Clone)]
pub struct FilterDesignBuilder {
    order: usize,
    return_loss_db: f64,
    center_hz: Option<f64>,
    bandwidth_hz: Option<f64>,
    zeros_normalized: Vec<f64>,
    zeros_hz: Vec<f64>,
    unloaded_q: Option<f64>,
}

impl FilterDesignBuilder {
    /// Adds transmission zeros in normalized prototype coordinates.
    ///
    /// Use this when working with a normalized prototype (no physical frequencies).
    /// Positive = upper stopband, negative = lower stopband.
    pub fn zeros(mut self, zeros: impl IntoIterator<Item = f64>) -> Self {
        self.zeros_normalized.extend(zeros);
        self
    }

    /// Adds transmission zeros in physical Hz.
    ///
    /// Requires that the builder was created with `bandpass()` (center/bandwidth known).
    /// The zeros are automatically normalized during `synthesize()`.
    pub fn zeros_hz(mut self, zeros: impl IntoIterator<Item = f64>) -> Self {
        self.zeros_hz.extend(zeros);
        self
    }

    /// Adds a single transmission zero in physical Hz.
    pub fn zero_hz(mut self, freq_hz: f64) -> Self {
        self.zeros_hz.push(freq_hz);
        self
    }

    /// Adds a single normalized transmission zero.
    pub fn zero(mut self, value: f64) -> Self {
        self.zeros_normalized.push(value);
        self
    }

    /// Sets the unloaded Q factor.
    pub fn unloaded_q(mut self, q: f64) -> Self {
        self.unloaded_q = Some(q);
        self
    }

    /// Runs the synthesis pipeline. Returns a completed `FilterDesign`.
    pub fn synthesize(self) -> Result<FilterDesign> {
        // Validate band-pass parameters at this common choke point. Both
        // `bandpass()` and `bandpass_edges()` (which *derives* center/bandwidth
        // from the passband edges — inverted edges yield bw < 0, non-positive
        // edges yield a NaN center) land here, so a single `BandPassMapping::new`
        // check rejects an invalid mapping with a typed error up front instead
        // of storing it and failing later at `.response()`/`.scale()`.
        if let (Some(center), Some(bw)) = (self.center_hz, self.bandwidth_hz) {
            BandPassMapping::new(center, bw)?;
        }
        // Normalize Hz zeros if band-pass parameters are available
        let mut all_zeros = self.zeros_normalized;
        if !self.zeros_hz.is_empty() {
            let center = self.center_hz.ok_or_else(|| {
                crate::error::MfsError::PreconditionViolation(
                    "zeros_hz requires band-pass parameters (use FilterDesign::bandpass)".into(),
                )
            })?;
            let bw = self.bandwidth_hz.ok_or_else(|| {
                crate::error::MfsError::PreconditionViolation(
                    "zeros_hz requires a band-pass bandwidth (use \
                     FilterDesign::bandpass)".into(),
                )
            })?;
            let mapping = BandPassMapping::new(center, bw)?;
            for hz in &self.zeros_hz {
                let sample = mapping.map_hz_to_normalized(*hz)?;
                all_zeros.push(sample.omega);
            }
        }

        let spec = build_filter_spec(self.order, self.return_loss_db, &all_zeros, self.unloaded_q)?;

        let polynomials = generalized_chebyshev_polynomials(&spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;

        Ok(FilterDesign {
            spec,
            polynomials,
            matrix,
            center_hz: self.center_hz,
            bandwidth_hz: self.bandwidth_hz,
        })
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// FilterDesign
// ═══════════════════════════════════════════════════════════════════════════════

/// A synthesized filter design.
///
/// Holds the coupling matrix and provides methods for topology transformation,
/// S-parameter evaluation, and band-pass scaling.
#[derive(Debug, Clone)]
pub struct FilterDesign {
    spec: FilterSpec,
    polynomials: PolynomialSet,
    matrix: CouplingMatrix,
    center_hz: Option<f64>,
    bandwidth_hz: Option<f64>,
}

impl FilterDesign {
    // ─── Constructors ────────────────────────────────────────────────────

    /// Band-pass filter design with physical frequency parameters.
    ///
    /// `center_hz` is the geometric center frequency used in the lowpass-to-bandpass
    /// frequency transformation: `Ω = (1/FBW) × (f/f0 - f0/f)` where `f0 = center_hz`
    /// and `FBW = bandwidth_hz / center_hz`.
    ///
    /// For the passband edges to map exactly to Ω = ±1, `center_hz` must be
    /// the geometric mean of the band edges: `center_hz = √(f_low × f_high)`.
    /// The bandwidth is then `bandwidth_hz = f_high - f_low`.
    ///
    /// If you have band edge frequencies, prefer [`Self::bandpass_edges()`] which
    /// computes the correct geometric center automatically.
    ///
    /// Transmission zeros can be specified in Hz via `.zeros_hz(...)`.
    ///
    /// ```rust
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
    ///     .zeros_hz([6.4e9, 7.0e9])
    ///     .synthesize()?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn bandpass(
        order: usize,
        return_loss_db: f64,
        center_hz: f64,
        bandwidth_hz: f64,
    ) -> FilterDesignBuilder {
        FilterDesignBuilder {
            order,
            return_loss_db,
            center_hz: Some(center_hz),
            bandwidth_hz: Some(bandwidth_hz),
            zeros_normalized: Vec::new(),
            zeros_hz: Vec::new(),
            unloaded_q: None,
        }
    }

    /// Band-pass filter design specified by passband edge frequencies.
    ///
    /// This is an alternative to `bandpass()` for engineers who think in terms
    /// of passband start/stop rather than center/bandwidth.
    ///
    /// The passband edges must be positive and ordered (`f_start_hz <
    /// f_stop_hz`): the center is derived as their geometric mean
    /// `√(f_start·f_stop)` and the bandwidth as `f_stop − f_start`. Inverted
    /// edges (bandwidth < 0) or a non-positive edge (non-finite center) are
    /// rejected with a typed [`crate::error::MfsError::InvalidFrequency`] at
    /// [`FilterDesignBuilder::synthesize`].
    /// ```rust
    /// # use mfs::prelude::*;
    /// // 6.6 GHz to 6.9 GHz passband → center=√(6.6e9×6.9e9), BW=300 MHz
    /// let d = FilterDesign::bandpass_edges(6, 23.0, 6.6e9, 6.9e9)
    ///     .zeros_hz([6.4e9, 7.0e9])
    ///     .synthesize()?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn bandpass_edges(
        order: usize,
        return_loss_db: f64,
        f_start_hz: f64,
        f_stop_hz: f64,
    ) -> FilterDesignBuilder {
        // Geometric mean ensures Ω(f_start)=-1 and Ω(f_stop)=+1 exactly.
        let (center_hz, bandwidth_hz) =
            crate::freq::edges_to_center_bandwidth(f_start_hz, f_stop_hz);
        FilterDesignBuilder {
            order,
            return_loss_db,
            center_hz: Some(center_hz),
            bandwidth_hz: Some(bandwidth_hz),
            zeros_normalized: Vec::new(),
            zeros_hz: Vec::new(),
            unloaded_q: None,
        }
    }

    /// Normalized prototype design (no physical frequencies).
    ///
    /// Transmission zeros are in normalized low-pass prototype coordinates.
    ///
    /// ```rust
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::prototype(4, 20.0)
    ///     .zeros([-1.5, 1.5])
    ///     .synthesize()?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn prototype(order: usize, return_loss_db: f64) -> FilterDesignBuilder {
        FilterDesignBuilder {
            order,
            return_loss_db,
            center_hz: None,
            bandwidth_hz: None,
            zeros_normalized: Vec::new(),
            zeros_hz: Vec::new(),
            unloaded_q: None,
        }
    }

    /// Backward-compatible alias for `prototype()`.
    pub fn chebyshev(order: usize, return_loss_db: f64) -> FilterDesignBuilder {
        Self::prototype(order, return_loss_db)
    }

    /// Creates a design from a pre-built `FilterSpec`.
    pub fn from_spec(spec: &FilterSpec) -> Result<Self> {
        let polynomials = generalized_chebyshev_polynomials(spec)?;
        let matrix = MatrixSynthesisEngine.synthesize(&polynomials)?;
        Ok(Self {
            spec: spec.clone(),
            polynomials,
            matrix,
            center_hz: None,
            bandwidth_hz: None,
        })
    }

    /// Creates a design from pre-computed polynomials and coupling matrix.
    ///
    /// Used internally by the duplexer synthesis to construct filter designs
    /// from modified polynomial sets.
    pub(crate) fn from_polynomials_and_matrix(
        order: usize,
        return_loss_db: f64,
        center_hz: f64,
        bandwidth_hz: f64,
        polynomials: PolynomialSet,
        matrix: CouplingMatrix,
    ) -> Result<Self> {
        // Build a minimal spec matching the design parameters. Propagate a
        // validation failure (order 0 / non-positive return loss) rather than
        // silently substituting a 20 dB spec.
        let spec = FilterSpec::new(order, return_loss_db)?;

        Ok(Self {
            spec,
            polynomials,
            matrix,
            center_hz: Some(center_hz),
            bandwidth_hz: Some(bandwidth_hz),
        })
    }

    // ─── Basic info ──────────────────────────────────────────────────────

    /// Filter order (number of resonators).
    pub fn order(&self) -> usize {
        self.spec.order
    }

    /// The transversal coupling matrix.
    pub fn matrix(&self) -> &CouplingMatrix {
        &self.matrix
    }

    /// Returns normalized transmission zeros (Omega values).
    pub fn transmission_zeros_normalized(&self) -> Vec<f64> {
        self.spec
            .transmission_zeros
            .iter()
            .filter(|tz| tz.value.is_finite())
            .map(|tz| tz.value)
            .collect()
    }
    /// The prototype polynomials (E, F, P).
    pub fn polynomials(&self) -> &PolynomialSet {
        &self.polynomials
    }

    /// The filter specification.
    pub fn spec(&self) -> &FilterSpec {
        &self.spec
    }

    /// Center frequency in Hz (if band-pass design).
    pub fn center_hz(&self) -> Option<f64> {
        self.center_hz
    }

    /// Bandwidth in Hz (if band-pass design).
    pub fn bandwidth_hz(&self) -> Option<f64> {
        self.bandwidth_hz
    }

    /// Lower passband edge in Hz (if band-pass design).
    ///
    /// The band-pass mapping uses a geometric-mean center, so the Ω=±1 edges
    /// are `√(c²+bw²)−c` / `√(c²+bw²)+c` (the `band_edge_pair` convention),
    /// **not** the arithmetic-mean `c ± bw/2` (review P2-2).
    pub fn f_start_hz(&self) -> Option<f64> {
        match (self.center_hz, self.bandwidth_hz) {
            (Some(c), Some(bw)) => Some(crate::freq::band_edge_pair(c, bw).0),
            _ => None,
        }
    }

    /// Upper passband edge in Hz (if band-pass design).
    ///
    /// See [`Self::f_start_hz`] for the geometric-mean edge convention.
    pub fn f_stop_hz(&self) -> Option<f64> {
        match (self.center_hz, self.bandwidth_hz) {
            (Some(c), Some(bw)) => Some(crate::freq::band_edge_pair(c, bw).1),
            _ => None,
        }
    }

    // ─── Topology ────────────────────────────────────────────────────────

    /// Converts to folded coupling matrix.
    pub fn to_folded(&self) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, MatrixTopology::Folded)
    }

    /// Converts to arrow coupling matrix.
    pub fn to_arrow(&self) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, MatrixTopology::Arrow)
    }

    /// Converts to Cascaded Trisection (CT) topology.
    ///
    /// A single condition-at-load-tail + peel-to-source (`extract_triplet`)
    /// reveals the full trisection cascade: `zeros[0]` conditions the peel and
    /// becomes the source-side trisection's zero; the remaining zeros are
    /// already encoded in the arrow form and are used only for the count
    /// validation (`zeros.len() <= order - 2`).
    pub fn to_ct(&self, zeros: &[f64]) -> Result<CouplingMatrix> {
        self.matrix.to_cascaded_trisections(zeros)
    }

    /// Converts to Cascaded Quadruplet (CQ) topology.
    ///
    /// Each zero pair is realized as one 4-resonator quadruplet, tiled
    /// non-overlappingly from the source side. Requires `order >= 4 *
    /// zero_pairs.len()` (use [`to_ct`](Self::to_ct) for denser configurations).
    pub fn to_cq(&self, zero_pairs: &[(f64, f64)]) -> Result<CouplingMatrix> {
        self.matrix.to_cascaded_quadruplets(zero_pairs)
    }

    /// Converts to a specific topology.
    pub fn to_topology(&self, topology: MatrixTopology) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, topology)
    }

    // ─── S-parameter response ────────────────────────────────────────────

    /// Evaluates S-parameters on a normalized prototype frequency grid.
    ///
    /// Uses the fast pole-expansion method (O(N) per point) by default.
    /// Automatically verifies against the LU method on first call and falls
    /// back to LU if results don't match.
    ///
    /// `start` and `stop` are normalized frequencies (0 = center of passband).
    pub fn response_normalized(
        &self,
        start: f64,
        stop: f64,
        points: usize,
    ) -> Result<SParameterResponse> {
        use crate::response::pole_expansion::PoleExpansionData;

        let grid = FrequencyGrid::linspace(start, stop, points)?;

        // Try pole expansion first (fast path). Verify it against LU at a few
        // sample points before trusting it: the fast path is numerically
        // unreliable near a network pole, so any error (or a check grid that
        // compared no valid points) falls back to LU rather than returning a
        // wrong/unverified response (P3 #1, #9).
        if let Ok(pole_data) = PoleExpansionData::from_matrix(&self.matrix) {
            let check_grid = FrequencyGrid::linspace(start, stop, 5.min(points))?;
            let lu_check = ResponseSolver.evaluate_normalized(&self.matrix, &check_grid)?;
            let pole_check = match pole_data.evaluate_normalized(&check_grid) {
                Ok(r) => r,
                Err(_) => return ResponseSolver.evaluate_normalized(&self.matrix, &grid),
            };

            // Deep-stopband check points (|S21| at the noise floor) are skipped,
            // so count how many were actually compared; an all-skipped grid
            // means the fast path was never verified.
            let mut max_diff = 0.0_f64;
            let mut compared = 0usize;
            for (p, l) in pole_check.samples.iter().zip(lu_check.samples.iter()) {
                let p_mag = p.s21_mag();
                let l_mag = l.s21_mag();
                if l_mag > 1e-10 {
                    let rel_diff = ((p_mag - l_mag) / l_mag).abs();
                    max_diff = max_diff.max(rel_diff);
                    compared += 1;
                }
            }

            if compared > 0 && max_diff < 1e-3 {
                // Verified fast path; a near-singular full-grid point still
                // falls back to LU.
                return pole_data
                    .evaluate_normalized(&grid)
                    .or_else(|_| ResponseSolver.evaluate_normalized(&self.matrix, &grid));
            }
        }

        // Fallback to LU-based evaluation
        ResponseSolver.evaluate_normalized(&self.matrix, &grid)
    }

    /// Evaluates S-parameters on the physical frequency grid.
    ///
    /// Uses the center/bandwidth from `bandpass()`. If this is a prototype design,
    /// you must provide explicit parameters via `response_bandpass()`.
    ///
    /// ```rust
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let r = d.response(6.5e9, 7.0e9, 201)?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn response(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<SParameterResponse> {
        let center = self.center_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "response() requires band-pass parameters; use response_bandpass() or \
                 create design with FilterDesign::bandpass()"
                    .into(),
            )
        })?;
        let bw = self.bandwidth_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "response() requires a band-pass bandwidth; use \
                 response_bandpass() or FilterDesign::bandpass()".into(),
            )
        })?;
        self.response_bandpass(center, bw, start_hz, stop_hz, points)
    }

    /// Evaluates S-parameters over a default frequency range.
    ///
    /// For bandpass designs: covers center ± 1.5×bandwidth (201 points).
    /// For prototype designs: covers normalized ω = -3 to +3 (201 points).
    ///
    /// This is a convenience method when you just want a quick look at the response.
    ///
    /// ```rust
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let r = d.response_default()?;
    /// assert_eq!(r.samples.len(), 201);
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn response_default(&self) -> Result<SParameterResponse> {
        self.default_response()
    }

    /// Evaluates S-parameters with explicit band-pass parameters.
    ///
    /// Uses pole expansion (fast) with automatic LU fallback verification.
    pub fn response_bandpass(
        &self,
        center_hz: f64,
        bandwidth_hz: f64,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<SParameterResponse> {
        use crate::response::pole_expansion::PoleExpansionData;

        let mapping = BandPassMapping::new(center_hz, bandwidth_hz)?;
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;

        // Try pole expansion first (fast path). Verify it against LU at a few
        // sample points before trusting it: the fast path is numerically
        // unreliable near a network pole, so any error (or a check grid that
        // compared no valid points) falls back to LU rather than returning a
        // wrong/unverified response (P3 #1, #9).
        if let Ok(pole_data) = PoleExpansionData::from_matrix(&self.matrix) {
            let n_check = 5.min(points);
            let check_grid = FrequencyGrid::linspace(start_hz, stop_hz, n_check)?;
            let lu_check = ResponseSolver.evaluate(&self.matrix, &check_grid, &mapping)?;
            let pole_check = match pole_data.evaluate_bandpass(&check_grid, &mapping) {
                Ok(r) => r,
                Err(_) => return ResponseSolver.evaluate(&self.matrix, &grid, &mapping),
            };

            // Deep-stopband check points (|S21| at the noise floor) are skipped,
            // so count how many were actually compared; an all-skipped grid
            // means the fast path was never verified.
            let mut max_diff = 0.0_f64;
            let mut compared = 0usize;
            for (p, l) in pole_check.samples.iter().zip(lu_check.samples.iter()) {
                let p_mag = p.s21_mag();
                let l_mag = l.s21_mag();
                if l_mag > 1e-10 {
                    let rel_diff = ((p_mag - l_mag) / l_mag).abs();
                    max_diff = max_diff.max(rel_diff);
                    compared += 1;
                }
            }

            if compared > 0 && max_diff < 1e-3 {
                // Verified fast path; a near-singular full-grid point still
                // falls back to LU.
                return pole_data
                    .evaluate_bandpass(&grid, &mapping)
                    .or_else(|_| ResponseSolver.evaluate(&self.matrix, &grid, &mapping));
            }
        }

        // Fallback to LU
        ResponseSolver.evaluate(&self.matrix, &grid, &mapping)
    }

    /// Evaluates S-parameters of any matrix on a normalized grid.
    pub fn eval(
        &self,
        matrix: &CouplingMatrix,
        start: f64,
        stop: f64,
        points: usize,
    ) -> Result<SParameterResponse> {
        let grid = FrequencyGrid::linspace(start, stop, points)?;
        ResponseSolver.evaluate_normalized(matrix, &grid)
    }

    /// Evaluates S-parameters with a custom mapping and grid.
    pub fn eval_mapped(
        &self,
        matrix: &CouplingMatrix,
        grid: &FrequencyGrid,
        mapping: &impl FrequencyMapping,
    ) -> Result<SParameterResponse> {
        ResponseSolver.evaluate(matrix, grid, mapping)
    }

    // ─── Lossy response ──────────────────────────────────────────────────

    /// Evaluates S-parameters with finite unloaded Q (lossy response).
    ///
    /// This predicts the actual filter response including insertion loss
    /// caused by resonator dissipation.
    ///
    /// ```rust
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let lossy = d.response_lossy(6.5e9, 7.0e9, 201, 3000.0)?;
    /// // lossy.samples[100].s21_db() will show insertion loss
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn response_lossy(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
        unloaded_q: f64,
    ) -> Result<SParameterResponse> {
        use crate::response::ResponseSettings;

        // P2-4: the backend drops a non-finite or non-positive Q silently
        // (lossless response for a "lossy" request); reject it here. ∞ = lossless.
        validate_unloaded_q(unloaded_q)?;
        let center = self.center_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "response_lossy() requires band-pass parameters".into(),
            )
        })?;
        let bw = self.bandwidth_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "response_lossy() requires a band-pass bandwidth; create \
                 the design with FilterDesign::bandpass()".into(),
            )
        })?;
        let mapping = BandPassMapping::new(center, bw)?;
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;

        // In normalized prototype domain, the dissipation per resonator is:
        //   δ = f₀ / (Qu × BW)
        // This converts physical Qu to the normalized loss term.
        let normalized_dissipation_q = unloaded_q * bw / center;

        let settings = ResponseSettings {
            source_resistance: 1.0,
            load_resistance: 1.0,
            unloaded_q: normalized_dissipation_q,
        };
        ResponseSolver.evaluate_with_settings(&self.matrix, &grid, &mapping, settings)
    }

    /// Evaluates lossy response on a normalized frequency grid.
    pub fn response_lossy_normalized(
        &self,
        start: f64,
        stop: f64,
        points: usize,
        unloaded_q: f64,
    ) -> Result<SParameterResponse> {
        use crate::response::ResponseSettings;

        // P2-4: reject non-finite / non-positive Q (∞ = lossless is allowed).
        validate_unloaded_q(unloaded_q)?;
        let grid = FrequencyGrid::linspace(start, stop, points)?;
        let settings = ResponseSettings {
            source_resistance: 1.0,
            load_resistance: 1.0,
            unloaded_q,
        };
        ResponseSolver.evaluate_normalized_with_settings(&self.matrix, &grid, settings)
    }

    // ─── Band-pass scaling ───────────────────────────────────────────────

    /// Scales the transversal matrix to physical band-pass units (Hz).
    ///
    /// Uses stored center/bandwidth if available, otherwise requires explicit params.
    pub fn scale(&self) -> Result<CouplingMatrix> {
        let center = self.center_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "scale() requires band-pass parameters; use scale_to() or \
                 create design with FilterDesign::bandpass()"
                    .into(),
            )
        })?;
        let bw = self.bandwidth_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "scale() requires a band-pass bandwidth; use scale_to() or \
                 create the design with FilterDesign::bandpass()".into(),
            )
        })?;
        self.scale_to(center, bw)
    }

    /// Scales the transversal matrix with explicit band-pass parameters.
    pub fn scale_to(&self, center_hz: f64, bandwidth_hz: f64) -> Result<CouplingMatrix> {
        let mapping = BandPassMapping::new(center_hz, bandwidth_hz)?;
        self.matrix.denormalize_bandpass(&mapping)
    }

    /// Scales any matrix to physical band-pass units.
    pub fn scale_matrix(
        &self,
        matrix: &CouplingMatrix,
        center_hz: f64,
        bandwidth_hz: f64,
    ) -> Result<CouplingMatrix> {
        let mapping = BandPassMapping::new(center_hz, bandwidth_hz)?;
        matrix.denormalize_bandpass(&mapping)
    }

    // ─── Touchstone export ───────────────────────────────────────────────

    /// Saves the S-parameter response to a Touchstone .s2p file.
    ///
    /// Uses default config (GHz, RI format, 50Ω). For band-pass designs,
    /// automatically computes the response over the design bandwidth.
    ///
    /// ```no_run
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// d.save_touchstone("filter.s2p")?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn save_touchstone(&self, path: impl AsRef<std::path::Path>) -> Result<()> {
        use crate::touchstone::{write_touchstone, TouchstoneConfig};

        let response = self.default_response()?;
        let config = TouchstoneConfig {
            comments: self.auto_comments(),
            ..TouchstoneConfig::default()
        };
        write_touchstone(&response, &config, path)
    }

    /// Returns a `TouchstoneBuilder` for configuring the export.
    ///
    /// ```no_run
    /// # use mfs::prelude::*;
    /// # use mfs::touchstone::{FreqUnit, DataFormat};
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// d.to_touchstone()
    ///     .freq_unit(FreqUnit::MHz)
    ///     .format(DataFormat::DB)
    ///     .save("filter.s2p")?;
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn to_touchstone(&self) -> TouchstoneBuilder<'_> {
        TouchstoneBuilder {
            design: self,
            config: crate::touchstone::TouchstoneConfig::default(),
            start_hz: None,
            stop_hz: None,
            points: None,
            range_normalized: false,
        }
    }

    /// Builds a live `touchstone::Network` (via `NetworkBuilder`) from the
    /// default response of this design. Unlike the text export, the returned
    /// network supports resampling, cascading, and Y/Z/ABCD conversion.
    /// Uses a 50 Ω reference impedance.
    ///
    /// ```no_run
    /// # use mfs::prelude::*;
    /// let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let net = d.build_network()?;
    /// assert_eq!(net.rank, 2);
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn build_network(&self) -> Result<touchstone::Network> {
        let response = self.default_response()?;
        crate::touchstone::build_network(&response, 50.0)
    }

    /// Cascades this design with `other` (port 2 → port 1) into a single
    /// 2-port `touchstone::Network`, both referenced to 50 Ω.
    ///
    /// Both designs are evaluated on their default response grid (passband
    /// edges padded by 0.5·bw, 201 points), so the call returns a typed error
    /// when the two grids differ — use `build_network` + `touchstone::cascade`
    /// with explicit `resample` for non-default grids.
    ///
    /// ```no_run
    /// # use mfs::prelude::*;
    /// let a = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let b = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    /// let net = a.cascade(&b)?;
    /// assert_eq!(net.rank, 2);
    /// # Ok::<(), MfsError>(())
    /// ```
    pub fn cascade(&self, other: &FilterDesign) -> Result<touchstone::Network> {
        let net_a = self.build_network()?;
        let net_b = other.build_network()?;
        crate::touchstone::cascade(&net_a, &net_b)
    }

    pub(crate) fn default_response(&self) -> Result<SParameterResponse> {
        if let (Some(center), Some(bw)) = (self.center_hz, self.bandwidth_hz) {
            // Default range: the passband (geometric-mean `band_edge_pair`
            // edges) padded by half a bandwidth on each side. This mirrors the
            // former `center ± 1.5·bw` window for narrowband designs, but is
            // derived from the *edges* so a wideband design (fbw > 2/3, where
            // `center − 1.5·bw < 0`) still yields a positive window instead of
            // a grid that fails in `map_hz_to_normalized` (review P2-3).
            let (lower_edge, upper_edge) = crate::freq::band_edge_pair(center, bw);
            let pad = 0.5 * bw;
            let start = (lower_edge - pad).max(1e-9 * center);
            let stop = upper_edge + pad;
            self.response_bandpass(center, bw, start, stop, 201)
        } else {
            self.response_normalized(-3.0, 3.0, 201)
        }
    }

    /// Generates auto-comment lines for Touchstone export.
    pub fn auto_comments(&self) -> Vec<String> {
        let mut comments = vec![
            format!(
                "MFS v{} - Generalized Chebyshev filter",
                env!("CARGO_PKG_VERSION")
            ),
            format!(
                "Order: {}, Return Loss: {:.1} dB",
                self.order(),
                self.spec.return_loss_db
            ),
        ];
        if let (Some(c), Some(bw)) = (self.center_hz, self.bandwidth_hz) {
            comments.push(format!(
                "Center: {:.6} GHz, BW: {:.3} MHz",
                c / 1e9,
                bw / 1e6
            ));
        }
        comments
    }
}

/// Validates a user-supplied unloaded Q for the lossy-response APIs.
///
/// The response backend silently drops a non-finite or non-positive Q (it
/// only applies the dissipation term when `is_finite() && > 0.0`), which
/// would return a *lossless* response for a "lossy" request with no warning
/// (review P2-4). `f64::INFINITY` is the documented lossless sentinel and is
/// allowed through; everything else must be a finite positive value.
fn validate_unloaded_q(unloaded_q: f64) -> Result<()> {
    if unloaded_q == f64::INFINITY {
        return Ok(());
    }
    if !unloaded_q.is_finite() || unloaded_q <= 0.0 {
        return Err(crate::error::MfsError::InvalidFrequency(format!(
            "unloaded_q must be a positive finite value (or f64::INFINITY for \
             lossless), got {unloaded_q}"
        )));
    }
    Ok(())
}

/// Builder for configuring Touchstone export.
pub struct TouchstoneBuilder<'a> {
    design: &'a FilterDesign,
    config: crate::touchstone::TouchstoneConfig,
    start_hz: Option<f64>,
    stop_hz: Option<f64>,
    points: Option<usize>,
    /// When set, `start_hz`/`stop_hz`/`points` are a *normalized ω* range
    /// (review P3 #4: a physical-Hz `freq_range` on a prototype design must
    /// not be silently reinterpreted as ω).
    range_normalized: bool,
}

impl<'a> TouchstoneBuilder<'a> {
    /// Sets the frequency unit.
    pub fn freq_unit(mut self, unit: crate::touchstone::FreqUnit) -> Self {
        self.config.freq_unit = unit;
        self
    }

    /// Sets the data format.
    pub fn format(mut self, fmt: crate::touchstone::DataFormat) -> Self {
        self.config.format = fmt;
        self
    }

    /// Sets the reference impedance.
    pub fn impedance(mut self, z0: f64) -> Self {
        self.config.impedance = z0;
        self
    }

    /// Sets the Touchstone version (V1 or V2).
    pub fn version(mut self, v: crate::touchstone::TouchstoneVersion) -> Self {
        self.config.version = v;
        self
    }

    /// Sets a **physical Hz** frequency range for the export.
    ///
    /// Only meaningful for band-pass designs (which carry a center/bandwidth).
    /// A prototype design has no physical frequency axis, so `build`/`save`
    /// return a typed error rather than silently treating these values as
    /// normalized ω (review P3 #4). Use [`Self::freq_range_normalized`] for a
    /// prototype.
    pub fn freq_range(mut self, start_hz: f64, stop_hz: f64, points: usize) -> Self {
        self.start_hz = Some(start_hz);
        self.stop_hz = Some(stop_hz);
        self.points = Some(points);
        self.range_normalized = false;
        self
    }

    /// Sets a **normalized ω** frequency range for the export.
    ///
    /// This is the explicit way to sweep a prototype design (no center/bw);
    /// it also works on a band-pass design if a prototype-domain sweep is
    /// wanted. The exported frequency column carries the ω values.
    pub fn freq_range_normalized(
        mut self,
        start_omega: f64,
        stop_omega: f64,
        points: usize,
    ) -> Self {
        self.start_hz = Some(start_omega);
        self.stop_hz = Some(stop_omega);
        self.points = Some(points);
        self.range_normalized = true;
        self
    }

    /// Adds a comment line.
    pub fn comment(mut self, text: impl Into<String>) -> Self {
        self.config.comments.push(text.into());
        self
    }

    /// Builds the Touchstone string.
    pub fn build(self) -> Result<String> {
        let response = self.get_response()?;
        let mut config = self.config;
        config.comments.extend(self.design.auto_comments());
        crate::touchstone::to_touchstone_string(&response, &config)
    }

    /// Saves to a file.
    pub fn save(self, path: impl AsRef<std::path::Path>) -> Result<()> {
        let response = self.get_response()?;
        let mut config = self.config;
        config.comments.extend(self.design.auto_comments());
        crate::touchstone::write_touchstone(&response, &config, path)
    }

    fn get_response(&self) -> Result<SParameterResponse> {
        if let (Some(start), Some(stop), Some(pts)) = (self.start_hz, self.stop_hz, self.points) {
            if self.range_normalized {
                // Explicit normalized-ω range: valid for any design.
                return self.design.response_normalized(start, stop, pts);
            }
            // Physical-Hz range: requires band-pass parameters. A prototype
            // design has no physical frequency axis, so error out instead of
            // silently reinterpreting Hz as ω (review P3 #4).
            match (self.design.center_hz, self.design.bandwidth_hz) {
                (Some(center), Some(bw)) => {
                    self.design.response_bandpass(center, bw, start, stop, pts)
                }
                _ => Err(crate::error::MfsError::PreconditionViolation(
                    "freq_range(start_hz, stop_hz) takes physical Hz and requires a \
                     band-pass design (center/bandwidth); this is a prototype design \
                     with no physical frequency axis. Use freq_range_normalized(\
                     start_omega, stop_omega, points) to sweep in the normalized \
                     prototype domain."
                        .into(),
                )),
            }
        } else {
            self.design.default_response()
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prototype_all_pole() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;
        assert_eq!(d.order(), 4);
        Ok(())
    }

    /// A physical-Hz `freq_range` on a prototype design (no center/bw) must
    /// be rejected, not silently reinterpreted as normalized ω (review P3
    /// #4); `freq_range_normalized` is the explicit prototype-domain path.
    #[test]
    fn touchstone_freq_range_hz_on_prototype_is_rejected() -> Result<()> {
        use crate::error::MfsError;
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;
        let err = d
            .to_touchstone()
            .freq_range(6.5e9, 7.0e9, 51)
            .build()
            .unwrap_err();
        assert!(
            matches!(err, MfsError::PreconditionViolation(_)),
            "a physical-Hz range on a prototype design must be a typed error, got {err:?}"
        );
        // The explicit normalized-ω range must work and sweep that domain.
        let r = d
            .to_touchstone()
            .freq_range_normalized(-2.0, 2.0, 51)
            .build()?;
        assert!(r.contains("# "), "touchstone output empty");
        Ok(())
    }

    #[test]
    fn prototype_with_zeros() -> Result<()> {
        let d = FilterDesign::prototype(6, 23.0)
            .zeros([-2.0, -1.2, 1.5])
            .synthesize()?;
        assert_eq!(d.order(), 6);
        Ok(())
    }

    #[test]
    fn bandpass_with_hz_zeros() -> Result<()> {
        let d = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
            .zeros_hz([6.4e9, 6.5e9, 7.0e9])
            .synthesize()?;
        assert_eq!(d.order(), 6);
        assert_eq!(d.center_hz(), Some(6.75e9));
        assert_eq!(d.bandwidth_hz(), Some(300e6));
        Ok(())
    }

    #[test]
    fn bandpass_edges_equivalent() -> Result<()> {
        // bandpass_edges(f1, f2) uses geometric center sqrt(f1*f2) and bw = f2 - f1
        // (so that Omega(f1) = -1 and Omega(f2) = +1 exactly), and must equal
        // bandpass(geometric_center, bw).
        let f_start: f64 = 6.6e9;
        let f_stop: f64 = 6.9e9;
        let center = (f_start * f_stop).sqrt();
        let bw = f_stop - f_start;

        let d1 = FilterDesign::bandpass(4, 20.0, center, bw).synthesize()?;
        let d2 = FilterDesign::bandpass_edges(4, 20.0, f_start, f_stop).synthesize()?;
        assert_eq!(d2.center_hz(), Some(center));
        assert_eq!(d2.bandwidth_hz(), Some(bw));
        // P2-2: the edges must round-trip to the *input* edges (geometric-
        // mean convention), not the arithmetic-mean `c ± bw/2`.
        assert!((d2.f_start_hz().unwrap() - f_start).abs() < 1e-6);
        assert!((d2.f_stop_hz().unwrap() - f_stop).abs() < 1e-6);
        // Coupling matrices should be identical
        let m1 = d1.matrix();
        let m2 = d2.matrix();
        for i in 0..m1.side() {
            for j in 0..m2.side() {
                let v1 = m1.at(i, j).unwrap_or(0.0);
                let v2 = m2.at(i, j).unwrap_or(0.0);
                assert!(
                    (v1 - v2).abs() < 1e-12,
                    "mismatch at ({i},{j}): {v1} vs {v2}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn response_default_works() -> Result<()> {
        let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
        let r = d.response_default()?;
        assert_eq!(r.samples.len(), 201);
        // P2-3: the window is the geometric passband edges padded by 0.5·bw,
        // which stays positive even for wideband designs (fbw > 2/3).
        let (lower_edge, upper_edge) = crate::freq::band_edge_pair(6.75e9, 300e6);
        let expected_start = (lower_edge - 0.5 * 300e6).max(1e-9 * 6.75e9);
        let expected_stop = upper_edge + 0.5 * 300e6;
        assert!((r.samples[0].frequency_hz - expected_start).abs() < 1.0);
        assert!((r.samples[200].frequency_hz - expected_stop).abs() < 1.0);
        Ok(())
    }

    /// Regression for review P2-3: a wideband design (fbw > 2/3, where the
    /// old `center − 1.5·bw` window went below 0 Hz) must still produce a
    /// valid default response with an all-positive window.
    #[test]
    fn response_default_wideband_positive_window() -> Result<()> {
        // fbw = bw/center = 1.5 > 2/3, so the old window went negative.
        let d = FilterDesign::bandpass(4, 20.0, 1.0e9, 1.5e9).synthesize()?;
        let r = d.response_default()?;
        assert!(r.samples[0].frequency_hz > 0.0, "window must stay positive");
        assert!(r.samples[0].frequency_hz <= 1.0e9, "window must start in/near the passband region");
        Ok(())
    }

    #[test]
    fn bandpass_response_shortcut() -> Result<()> {
        let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
        let r = d.response(6.5e9, 7.0e9, 101)?;
        assert_eq!(r.samples.len(), 101);
        Ok(())
    }

    #[test]
    fn bandpass_scale_shortcut() -> Result<()> {
        let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
        let scaled = d.scale()?;
        assert_eq!(scaled.order(), 4);
        Ok(())
    }

    #[test]
    fn topologies() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0)
            .zeros([1.5, -1.5])
            .synthesize()?;
        let folded = d.to_folded()?;
        let arrow = d.to_arrow()?;
        assert_eq!(folded.order(), 4);
        assert_eq!(arrow.order(), 4);
        Ok(())
    }

    #[test]
    fn normalized_response_power_conservation() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;
        let r = d.response_normalized(-2.0, 2.0, 41)?;
        assert_eq!(r.samples.len(), 41);
        let c = &r.samples[20];
        let power = c.s11_mag().powi(2) + c.s21_mag().powi(2);
        assert!((power - 1.0).abs() < 1e-9);
        Ok(())
    }

    #[test]
    fn unloaded_q() -> Result<()> {
        let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
            .unloaded_q(3000.0)
            .synthesize()?;
        assert_eq!(d.spec().unloaded_q, Some(3000.0));
        Ok(())
    }

    #[test]
    fn chebyshev_alias() -> Result<()> {
        let d = FilterDesign::chebyshev(4, 20.0)
            .zeros([1.5, -1.5])
            .synthesize()?;
        assert_eq!(d.order(), 4);
        Ok(())
    }

    #[test]
    fn from_spec() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?.with_normalized_transmission_zeros(vec![2.0]);
        let d = FilterDesign::from_spec(&spec)?;
        assert_eq!(d.order(), 3);
        Ok(())
    }

    #[test]
    fn eval_transformed() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;
        let folded = d.to_folded()?;
        let r = d.eval(&folded, -2.0, 2.0, 21)?;
        assert_eq!(r.samples.len(), 21);
        Ok(())
    }

    #[test]
    fn pole_expansion_matches_lu_all_pole() -> Result<()> {
        use crate::response::pole_expansion::PoleExpansionData;

        let d = FilterDesign::prototype(4, 20.0).synthesize()?;
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 101)?;

        // LU-based reference
        let lu_resp = ResponseSolver.evaluate_normalized(d.matrix(), &grid)?;
        // Pole expansion from matrix
        let pole_data = PoleExpansionData::from_matrix(d.matrix())?;
        let pole_resp = pole_data.evaluate_normalized(&grid)?;

        for (lu, pole) in lu_resp.samples.iter().zip(pole_resp.samples.iter()) {
            let lu_s21 = lu.s21_mag();
            let pole_s21 = pole.s21_mag();
            if lu_s21 > 1e-10 {
                let rel_err = ((pole_s21 - lu_s21) / lu_s21).abs();
                assert!(rel_err < 1e-6, "S21 mismatch at ω={}: LU={lu_s21:.6e}, pole={pole_s21:.6e}, rel_err={rel_err:.2e}",
                    lu.normalized_omega);
            }
        }
        Ok(())
    }

    #[test]
    fn pole_expansion_matches_lu_with_zeros() -> Result<()> {
        use crate::response::pole_expansion::PoleExpansionData;

        let d = FilterDesign::prototype(6, 23.0)
            .zeros([-2.0, -1.2, 1.5])
            .synthesize()?;
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 101)?;

        let lu_resp = ResponseSolver.evaluate_normalized(d.matrix(), &grid)?;
        let pole_data = PoleExpansionData::from_matrix(d.matrix())?;
        let pole_resp = pole_data.evaluate_normalized(&grid)?;

        for (lu, pole) in lu_resp.samples.iter().zip(pole_resp.samples.iter()) {
            let lu_s21 = lu.s21_mag();
            let pole_s21 = pole.s21_mag();
            if lu_s21 > 1e-10 {
                let rel_err = ((pole_s21 - lu_s21) / lu_s21).abs();
                assert!(rel_err < 1e-4, "S21 mismatch at ω={}: LU={lu_s21:.6e}, pole={pole_s21:.6e}, rel_err={rel_err:.2e}",
                    lu.normalized_omega);
            }
        }
        Ok(())
    }

    #[test]
    fn pole_expansion_matches_lu_high_order() -> Result<()> {
        use crate::response::pole_expansion::PoleExpansionData;

        let d = FilterDesign::prototype(12, 20.0)
            .zeros([-1.3, 1.3, -1.8, 1.8])
            .synthesize()?;
        let grid = FrequencyGrid::linspace(-3.0, 3.0, 201)?;

        let lu_resp = ResponseSolver.evaluate_normalized(d.matrix(), &grid)?;
        let pole_data = PoleExpansionData::from_matrix(d.matrix())?;
        let pole_resp = pole_data.evaluate_normalized(&grid)?;

        for (lu, pole) in lu_resp.samples.iter().zip(pole_resp.samples.iter()) {
            let lu_s21 = lu.s21_mag();
            let pole_s21 = pole.s21_mag();
            if lu_s21 > 1e-8 {
                let rel_err = ((pole_s21 - lu_s21) / lu_s21).abs();
                assert!(rel_err < 1e-3, "S21 mismatch at ω={:.3}: LU={lu_s21:.6e}, pole={pole_s21:.6e}, rel_err={rel_err:.2e}",
                    lu.normalized_omega);
            }
        }
        Ok(())
    }

    #[test]
    fn lossy_response_has_insertion_loss() -> Result<()> {
        let d = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;

        let lossless = d.response(6.5e9, 7.0e9, 21)?;
        let lossy = d.response_lossy(6.5e9, 7.0e9, 21, 200.0)?; // Qu=200 (moderate loss)

        // At center frequency, lossy S21 should be lower than lossless
        let center_idx = 10;
        let lossless_s21 = lossless.samples[center_idx].s21_db();
        let lossy_s21 = lossy.samples[center_idx].s21_db();

        // Lossy should have more insertion loss (more negative dB)
        // For order 4, Qu=200, BW/f0=0.044: IL ≈ 4.343*4/(200*0.044) ≈ 2.0 dB
        assert!(lossy_s21 < lossless_s21 - 0.5,
            "lossy S21 ({lossy_s21:.2} dB) should be noticeably less than lossless ({lossless_s21:.2} dB)");

        // Power should NOT be conserved for lossy (|S11|² + |S21|² < 1)
        let s = &lossy.samples[center_idx];
        let power = s.s11_mag().powi(2) + s.s21_mag().powi(2);
        assert!(
            power < 0.99,
            "lossy power sum should be < 1, got {power:.4}"
        );

        Ok(())
    }

    #[test]
    fn lossy_normalized_response() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;

        let lossless = d.response_normalized(-2.0, 2.0, 41)?;
        let lossy = d.response_lossy_normalized(-2.0, 2.0, 41, 5.0)?; // Very low normalized Q

        // Center point: lossy should have more insertion loss
        let center = 20;
        assert!(
            lossy.samples[center].s21_db() < lossless.samples[center].s21_db() - 0.1,
            "lossy={:.2} dB, lossless={:.2} dB",
            lossy.samples[center].s21_db(),
            lossless.samples[center].s21_db()
        );
        Ok(())
    }

    #[test]
    fn infinite_q_equals_lossless() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;

        let lossless = d.response_normalized(-2.0, 2.0, 11)?;
        let inf_q = d.response_lossy_normalized(-2.0, 2.0, 11, f64::INFINITY)?;

        for (a, b) in lossless.samples.iter().zip(inf_q.samples.iter()) {
            assert!((a.s21_mag() - b.s21_mag()).abs() < 1e-10);
        }
        Ok(())
    }
}
