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
use crate::synthesis::MatrixSynthesisEngine;

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
        // Normalize Hz zeros if band-pass parameters are available
        let mut all_zeros = self.zeros_normalized;
        if !self.zeros_hz.is_empty() {
            let center = self.center_hz.ok_or_else(|| {
                crate::error::MfsError::PreconditionViolation(
                    "zeros_hz requires band-pass parameters (use FilterDesign::bandpass)".into(),
                )
            })?;
            let bw = self.bandwidth_hz.unwrap_or(0.0);
            let mapping = BandPassMapping::new(center, bw)?;
            for hz in &self.zeros_hz {
                let sample = mapping.map_hz_to_normalized(*hz)?;
                all_zeros.push(sample.omega);
            }
        }

        let mut spec = FilterSpec::new(self.order, self.return_loss_db)?;
        if !all_zeros.is_empty() {
            spec = spec.with_normalized_transmission_zeros(all_zeros);
        }
        if let Some(q) = self.unloaded_q {
            spec = spec.with_unloaded_q(q);
        }

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
    /// This is the most common entry point. Transmission zeros can be specified
    /// in Hz via `.zeros_hz(...)`.
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

    // ─── Basic info ──────────────────────────────────────────────────────

    /// Filter order (number of resonators).
    pub fn order(&self) -> usize { self.spec.order }

    /// The transversal coupling matrix.
    pub fn matrix(&self) -> &CouplingMatrix { &self.matrix }

    /// The prototype polynomials (E, F, P).
    pub fn polynomials(&self) -> &PolynomialSet { &self.polynomials }

    /// The filter specification.
    pub fn spec(&self) -> &FilterSpec { &self.spec }

    /// Center frequency in Hz (if band-pass design).
    pub fn center_hz(&self) -> Option<f64> { self.center_hz }

    /// Bandwidth in Hz (if band-pass design).
    pub fn bandwidth_hz(&self) -> Option<f64> { self.bandwidth_hz }

    // ─── Topology ────────────────────────────────────────────────────────

    /// Converts to folded coupling matrix.
    pub fn to_folded(&self) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, MatrixTopology::Folded)
    }

    /// Converts to arrow coupling matrix.
    pub fn to_arrow(&self) -> Result<CouplingMatrix> {
        MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, MatrixTopology::Arrow)
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

        // Try pole expansion first (fast path)
        if let Ok(pole_data) = PoleExpansionData::from_matrix(&self.matrix) {
            let pole_response = pole_data.evaluate_normalized(&grid)?;

            // Verify against LU at a few sample points
            let check_grid = FrequencyGrid::linspace(start, stop, 5.min(points))?;
            let lu_check = ResponseSolver.evaluate_normalized(&self.matrix, &check_grid)?;
            let pole_check = pole_data.evaluate_normalized(&check_grid)?;

            let mut max_diff = 0.0_f64;
            for (p, l) in pole_check.samples.iter().zip(lu_check.samples.iter()) {
                let p_mag = p.s21_mag();
                let l_mag = l.s21_mag();
                if l_mag > 1e-10 {
                    let rel_diff = ((p_mag - l_mag) / l_mag).abs();
                    max_diff = max_diff.max(rel_diff);
                }
            }

            // If pole expansion matches LU within 0.1% relative error, use it
            if max_diff < 1e-3 {
                return Ok(pole_response);
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
        let bw = self.bandwidth_hz.unwrap_or(0.0);
        self.response_bandpass(center, bw, start_hz, stop_hz, points)
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

        // Try pole expansion first (fast path)
        if let Ok(pole_data) = PoleExpansionData::from_matrix(&self.matrix) {
            let pole_response = pole_data.evaluate_bandpass(&grid, &mapping)?;

            // Verify at a few sample points
            let n_check = 5.min(points);
            let check_grid = FrequencyGrid::linspace(start_hz, stop_hz, n_check)?;
            let lu_check = ResponseSolver.evaluate(&self.matrix, &check_grid, &mapping)?;
            let pole_check = pole_data.evaluate_bandpass(&check_grid, &mapping)?;

            let mut max_diff = 0.0_f64;
            for (p, l) in pole_check.samples.iter().zip(lu_check.samples.iter()) {
                let p_mag = p.s21_mag();
                let l_mag = l.s21_mag();
                if l_mag > 1e-10 {
                    let rel_diff = ((p_mag - l_mag) / l_mag).abs();
                    max_diff = max_diff.max(rel_diff);
                }
            }

            if max_diff < 1e-3 {
                return Ok(pole_response);
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

        let center = self.center_hz.ok_or_else(|| {
            crate::error::MfsError::PreconditionViolation(
                "response_lossy() requires band-pass parameters".into(),
            )
        })?;
        let bw = self.bandwidth_hz.unwrap_or(0.0);
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
        let bw = self.bandwidth_hz.unwrap_or(0.0);
        self.scale_to(center, bw)
    }

    /// Scales the transversal matrix with explicit band-pass parameters.
    pub fn scale_to(
        &self,
        center_hz: f64,
        bandwidth_hz: f64,
    ) -> Result<CouplingMatrix> {
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
        let mut config = TouchstoneConfig::default();
        config.comments = self.auto_comments();
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
        }
    }

    pub(crate) fn default_response(&self) -> Result<SParameterResponse> {
        if let (Some(center), Some(bw)) = (self.center_hz, self.bandwidth_hz) {
            let start = center - bw;
            let stop = center + bw;
            self.response_bandpass(center, bw, start, stop, 201)
        } else {
            self.response_normalized(-3.0, 3.0, 201)
        }
    }

    /// Generates auto-comment lines for Touchstone export.
    pub fn auto_comments(&self) -> Vec<String> {
        let mut comments = vec![
            format!("MFS v{} - Generalized Chebyshev filter", env!("CARGO_PKG_VERSION")),
            format!("Order: {}, Return Loss: {:.1} dB", self.order(), self.spec.return_loss_db),
        ];
        if let (Some(c), Some(bw)) = (self.center_hz, self.bandwidth_hz) {
            comments.push(format!("Center: {:.6} GHz, BW: {:.3} MHz", c / 1e9, bw / 1e6));
        }
        comments
    }
}

/// Builder for configuring Touchstone export.
pub struct TouchstoneBuilder<'a> {
    design: &'a FilterDesign,
    config: crate::touchstone::TouchstoneConfig,
    start_hz: Option<f64>,
    stop_hz: Option<f64>,
    points: Option<usize>,
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

    /// Sets the frequency range for the export.
    pub fn freq_range(mut self, start_hz: f64, stop_hz: f64, points: usize) -> Self {
        self.start_hz = Some(start_hz);
        self.stop_hz = Some(stop_hz);
        self.points = Some(points);
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
            if let (Some(center), Some(bw)) = (self.design.center_hz, self.design.bandwidth_hz) {
                self.design.response_bandpass(center, bw, start, stop, pts)
            } else {
                self.design.response_normalized(start, stop, pts)
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
        let spec = FilterSpec::new(3, 20.0)?
            .with_normalized_transmission_zeros(vec![2.0]);
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
        assert!(power < 0.99, "lossy power sum should be < 1, got {power:.4}");

        Ok(())
    }

    #[test]
    fn lossy_normalized_response() -> Result<()> {
        let d = FilterDesign::prototype(4, 20.0).synthesize()?;

        let lossless = d.response_normalized(-2.0, 2.0, 41)?;
        let lossy = d.response_lossy_normalized(-2.0, 2.0, 41, 5.0)?; // Very low normalized Q

        // Center point: lossy should have more insertion loss
        let center = 20;
        assert!(lossy.samples[center].s21_db() < lossless.samples[center].s21_db() - 0.1,
            "lossy={:.2} dB, lossless={:.2} dB",
            lossy.samples[center].s21_db(), lossless.samples[center].s21_db());
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
