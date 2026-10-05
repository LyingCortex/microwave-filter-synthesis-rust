//! Python bindings via PyO3.
//!
//! Activated with `cargo build --features python`.
//! Install with `maturin develop` or `pip install .`
//!
//! Python usage:
//! ```python
//! import mfs
//!
//! # Band-pass filter design
//! design = mfs.bandpass(order=6, rl=23.0, center=6.75e9, bw=300e6,
//!                       zeros=[6.4e9, 6.5e9, 7.0e9], q=3000.0)
//!
//! # Get coupling matrices
//! m = design.matrix()          # transversal (numpy 2D array)
//! m = design.folded()          # folded topology
//! m = design.arrow()           # arrow topology
//!
//! # S-parameter response
//! freq, s21, s11 = design.response(6.0e9, 7.5e9, 201)
//!
//! # Scaled matrix
//! m = design.scaled_matrix()
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::design::FilterDesign;
use crate::matrix::CouplingMatrix;

/// Convert a CouplingMatrix to a nested list for Python.
fn matrix_to_pylist(py: Python<'_>, matrix: &CouplingMatrix) -> PyObject {
    use pyo3::types::PyList;
    let side = matrix.side();
    let rows: Vec<_> = (0..side)
        .map(|row| {
            let row_data: Vec<f64> = (0..side)
                .map(|col| matrix.at(row, col).unwrap_or(0.0))
                .collect();
            PyList::new(py, &row_data).unwrap().into()
        })
        .collect::<Vec<PyObject>>();
    PyList::new(py, &rows).unwrap().into()
}

/// A synthesized filter design exposed to Python.
#[pyclass(name = "FilterDesign")]
#[derive(Clone)]
struct PyFilterDesign {
    inner: FilterDesign,
}

#[pymethods]
impl PyFilterDesign {
    /// Filter order.
    #[getter]
    fn order(&self) -> usize {
        self.inner.order()
    }

    /// Center frequency in Hz (None for prototype designs).
    #[getter]
    fn center(&self) -> Option<f64> {
        self.inner.center_hz()
    }

    /// Bandwidth in Hz (None for prototype designs).
    #[getter]
    fn bw(&self) -> Option<f64> {
        self.inner.bandwidth_hz()
    }

    /// Lower passband edge in Hz (None for prototype designs).
    #[getter]
    fn f_start(&self) -> Option<f64> {
        self.inner.f_start_hz()
    }

    /// Upper passband edge in Hz (None for prototype designs).
    #[getter]
    fn f_stop(&self) -> Option<f64> {
        self.inner.f_stop_hz()
    }

    /// Returns the transversal coupling matrix as a 2D list.
    fn matrix(&self, py: Python<'_>) -> PyObject {
        matrix_to_pylist(py, self.inner.matrix())
    }

    /// Returns the folded coupling matrix as a 2D list.
    fn folded(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self
            .inner
            .to_folded()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns the arrow coupling matrix as a 2D list.
    fn arrow(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self
            .inner
            .to_arrow()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns CT (Cascaded Trisection) topology matrix.
    /// positions: list of [from_node, to_node] pairs (1-indexed resonator numbers).
    #[pyo3(signature = (positions))]
    fn ct(&self, py: Python<'_>, positions: Vec<[usize; 2]>) -> PyResult<PyObject> {
        let arrow = self
            .inner
            .to_arrow()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let zeros = self.inner.transmission_zeros_normalized();
        let pos: Vec<(usize, usize)> = positions.iter().map(|p| (p[0], p[1])).collect();
        let m = arrow
            .extract_cascade_trisections(&zeros, &pos)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns CQ (Cascaded Quadruplet) topology matrix.
    /// z_locs: list of [from, to] pairs. span=2 for trisection, span=3 for quadruplet diagonal.
    #[pyo3(signature = (z_locs))]
    fn cq(&self, py: Python<'_>, z_locs: Vec<[usize; 2]>) -> PyResult<PyObject> {
        let arrow = self
            .inner
            .to_arrow()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let zeros = self.inner.transmission_zeros_normalized();
        let m = arrow
            .extract_cascade_quadruplets(&zeros, &z_locs)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns BOX topology matrix.
    /// spans: list of [start, end] pairs with end-start=3, one per zero.
    #[pyo3(signature = (spans))]
    fn box_topology(&self, py: Python<'_>, spans: Vec<[usize; 2]>) -> PyResult<PyObject> {
        let arrow = self
            .inner
            .to_arrow()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let zeros = self.inner.transmission_zeros_normalized();
        let m = arrow
            .extract_box_sections(&zeros, &spans)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns the band-pass scaled matrix as a 2D list.
    fn scaled_matrix(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self
            .inner
            .scale()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Computes S-parameter response.
    ///
    /// Returns (frequencies, s21_db, s11_db) as three lists.
    ///
    /// If start/stop are None, uses default range (center ± 1.5*bw for bandpass,
    /// or -3 to +3 for prototype).
    ///
    /// For bandpass designs: `response(start_hz, stop_hz, points)`
    /// For prototype designs: `response(start, stop, points)`
    #[pyo3(signature = (start=None, stop=None, points=201))]
    fn response(
        &self,
        start: Option<f64>,
        stop: Option<f64>,
        points: usize,
    ) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let resp = match (start, stop) {
            (Some(s), Some(e)) => {
                if self.inner.center_hz().is_some() {
                    self.inner.response(s, e, points)
                } else {
                    self.inner.response_normalized(s, e, points)
                }
            }
            _ => {
                // Use default range
                self.inner.response_default()
            }
        };
        let resp = resp.map_err(|e| PyValueError::new_err(e.to_string()))?;

        let freq: Vec<f64> = resp.samples.iter().map(|s| s.frequency_hz).collect();
        let s21: Vec<f64> = resp.samples.iter().map(|s| s.s21_db()).collect();
        let s11: Vec<f64> = resp.samples.iter().map(|s| s.s11_db()).collect();
        Ok((freq, s21, s11))
    }

    /// Computes full S-parameter response with all fields.
    ///
    /// Returns a list of dicts with keys: freq, s21_db, s11_db, s21_mag, s11_mag,
    /// s21_phase, s11_phase, group_delay (float, or None where the solver path
    /// does not compute group delay).
    #[pyo3(signature = (start=None, stop=None, points=201))]
    fn response_full(
        &self,
        py: Python<'_>,
        start: Option<f64>,
        stop: Option<f64>,
        points: usize,
    ) -> PyResult<PyObject> {
        use pyo3::types::PyDict;

        let resp = match (start, stop) {
            (Some(s), Some(e)) => {
                if self.inner.center_hz().is_some() {
                    self.inner.response(s, e, points)
                } else {
                    self.inner.response_normalized(s, e, points)
                }
            }
            _ => self.inner.response_default(),
        };
        let resp = resp.map_err(|e| PyValueError::new_err(e.to_string()))?;

        let result: Vec<_> = resp
            .samples
            .iter()
            .map(|s| {
                let dict = PyDict::new(py);
                dict.set_item("freq", s.frequency_hz).unwrap();
                dict.set_item("s21_db", s.s21_db()).unwrap();
                dict.set_item("s11_db", s.s11_db()).unwrap();
                dict.set_item("s21_mag", s.s21_mag()).unwrap();
                dict.set_item("s11_mag", s.s11_mag()).unwrap();
                dict.set_item("s21_phase", s.s21_phase_deg()).unwrap();
                dict.set_item("s11_phase", s.s11_phase_deg()).unwrap();
                dict.set_item("group_delay", s.group_delay).unwrap();
                dict.into()
            })
            .collect::<Vec<PyObject>>();

        Ok({
            use pyo3::types::PyList;
            PyList::new(py, &result)?.into()
        })
    }

    fn __repr__(&self) -> String {
        let zeros_info = if self.inner.spec().transmission_zeros.is_empty() {
            "all-pole".to_string()
        } else {
            format!("{} zeros", self.inner.spec().transmission_zeros.len())
        };
        match (self.inner.center_hz(), self.inner.bandwidth_hz()) {
            (Some(c), Some(bw)) => format!(
                "FilterDesign(order={}, rl={:.1}dB, center={:.4}GHz, bw={:.1}MHz, {})",
                self.inner.order(),
                self.inner.spec().return_loss_db,
                c / 1e9,
                bw / 1e6,
                zeros_info,
            ),
            _ => format!(
                "FilterDesign(order={}, rl={:.1}dB, prototype, {})",
                self.inner.order(),
                self.inner.spec().return_loss_db,
                zeros_info,
            ),
        }
    }

    /// Returns the Touchstone-formatted string.
    #[pyo3(signature = (freq_unit="GHz", format="RI", impedance=50.0, version=1))]
    fn to_touchstone(
        &self,
        freq_unit: &str,
        format: &str,
        impedance: f64,
        version: u8,
    ) -> PyResult<String> {
        use crate::touchstone::{DataFormat, FreqUnit, TouchstoneConfig, TouchstoneVersion};

        let fu = match freq_unit.to_uppercase().as_str() {
            "HZ" => FreqUnit::Hz,
            "KHZ" => FreqUnit::KHz,
            "MHZ" => FreqUnit::MHz,
            "GHZ" | _ => FreqUnit::GHz,
        };
        let df = match format.to_uppercase().as_str() {
            "MA" => DataFormat::MA,
            "DB" => DataFormat::DB,
            "RI" | _ => DataFormat::RI,
        };
        let ver = if version >= 2 {
            TouchstoneVersion::V2
        } else {
            TouchstoneVersion::V1
        };

        let response = self
            .inner
            .default_response()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut config = TouchstoneConfig {
            freq_unit: fu,
            format: df,
            impedance,
            version: ver,
            comments: Vec::new(),
        };
        config.comments = self.inner.auto_comments();
        crate::touchstone::to_touchstone_string(&response, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Saves the Touchstone file to disk.
    #[pyo3(signature = (path, freq_unit="GHz", format="RI", impedance=50.0, version=1))]
    fn save_touchstone(
        &self,
        path: &str,
        freq_unit: &str,
        format: &str,
        impedance: f64,
        version: u8,
    ) -> PyResult<()> {
        use crate::touchstone::{DataFormat, FreqUnit, TouchstoneConfig, TouchstoneVersion};

        let fu = match freq_unit.to_uppercase().as_str() {
            "HZ" => FreqUnit::Hz,
            "KHZ" => FreqUnit::KHz,
            "MHZ" => FreqUnit::MHz,
            _ => FreqUnit::GHz,
        };
        let df = match format.to_uppercase().as_str() {
            "MA" => DataFormat::MA,
            "DB" => DataFormat::DB,
            _ => DataFormat::RI,
        };
        let ver = if version >= 2 {
            TouchstoneVersion::V2
        } else {
            TouchstoneVersion::V1
        };

        let response = self
            .inner
            .default_response()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut config = TouchstoneConfig {
            freq_unit: fu,
            format: df,
            impedance,
            version: ver,
            comments: Vec::new(),
        };
        config.comments = self.inner.auto_comments();
        crate::touchstone::write_touchstone(&response, &config, path)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Computes lossy S-parameter response with finite unloaded Q.
    ///
    /// Returns (frequencies, s21_db, s11_db).
    #[pyo3(signature = (start, stop, points=201, q=1000.0))]
    fn response_lossy(
        &self,
        start: f64,
        stop: f64,
        points: usize,
        q: f64,
    ) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let resp = if self.inner.center_hz().is_some() {
            self.inner.response_lossy(start, stop, points, q)
        } else {
            self.inner.response_lossy_normalized(start, stop, points, q)
        };
        let resp = resp.map_err(|e| PyValueError::new_err(e.to_string()))?;

        let freq: Vec<f64> = resp.samples.iter().map(|s| s.frequency_hz).collect();
        let s21: Vec<f64> = resp.samples.iter().map(|s| s.s21_db()).collect();
        let s11: Vec<f64> = resp.samples.iter().map(|s| s.s11_db()).collect();
        Ok((freq, s21, s11))
    }

    /// Tunes the coupling matrix to match a target response.
    ///
    /// `target_s21` and `target_s11` are lists of complex values (re, im)
    /// pairs at the same frequency points as `freqs`.
    ///
    /// The optimizer rebuilds its evaluation grid as
    /// `linspace(freqs[0], freqs[-1], len(freqs))`, so `freqs` must be
    /// uniformly ascending; non-uniform `freqs` are mis-interpreted (the
    /// targets would be matched at the wrong normalized frequencies).
    ///
    /// Returns a dict with: matrix (2D list), cost (float), iterations (int), converged (bool).
    #[pyo3(signature = (freqs, target_s21_re, target_s21_im, target_s11_re, target_s11_im, max_iter=200))]
    fn tune(
        &self,
        py: Python<'_>,
        freqs: Vec<f64>,
        target_s21_re: Vec<f64>,
        target_s21_im: Vec<f64>,
        target_s11_re: Vec<f64>,
        target_s11_im: Vec<f64>,
        max_iter: usize,
    ) -> PyResult<PyObject> {
        use crate::freq::FrequencyGrid;
        use crate::optimize::{optimize_matrix, OptimizeConfig};
        use crate::response::{ResponseSample, SParameterResponse};
        use pyo3::types::PyDict;

        let n = freqs.len();
        for (name, v) in [
            ("target_s21_re", &target_s21_re),
            ("target_s21_im", &target_s21_im),
            ("target_s11_re", &target_s11_re),
            ("target_s11_im", &target_s11_im),
        ] {
            if v.len() != n {
                return Err(PyValueError::new_err(format!(
                    "{name} must have the same length as freqs ({n})"
                )));
            }
        }

        let target = SParameterResponse {
            samples: freqs
                .iter()
                .enumerate()
                .map(|(i, &f)| ResponseSample {
                    frequency_hz: f,
                    normalized_omega: f,
                    group_delay: None,
                    s11_re: target_s11_re[i],
                    s11_im: target_s11_im[i],
                    s21_re: target_s21_re[i],
                    s21_im: target_s21_im[i],
                })
                .collect(),
        };

        let grid = FrequencyGrid::linspace(
            *freqs.first().unwrap_or(&-2.0),
            *freqs.last().unwrap_or(&2.0),
            freqs.len(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let config = OptimizeConfig {
            max_iterations: max_iter,
            ..Default::default()
        };
        let result = optimize_matrix(self.inner.matrix(), &target, &grid, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let dict = PyDict::new(py);
        dict.set_item("matrix", matrix_to_pylist(py, &result.matrix))
            .unwrap();
        dict.set_item("cost", result.cost).unwrap();
        dict.set_item("iterations", result.iterations).unwrap();
        dict.set_item("converged", result.converged).unwrap();
        Ok(dict.into())
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Module-level functions
// ═══════════════════════════════════════════════════════════════════════════════

/// Creates a band-pass filter design.
///
/// Args:
///     order: Number of resonators
///     rl: Return loss in dB
///     center: Center frequency in Hz
///     bw: Bandwidth in Hz
///     zeros: Transmission zero frequencies in Hz (optional)
///     q: Unloaded Q factor (optional)
///
/// Returns:
///     FilterDesign object
#[pyfunction]
#[pyo3(signature = (order, rl, center, bw, zeros=None, q=None))]
fn bandpass(
    order: usize,
    rl: f64,
    center: f64,
    bw: f64,
    zeros: Option<Vec<f64>>,
    q: Option<f64>,
) -> PyResult<PyFilterDesign> {
    let mut builder = FilterDesign::bandpass(order, rl, center, bw);
    if let Some(z) = zeros {
        builder = builder.zeros_hz(z);
    }
    if let Some(q_val) = q {
        builder = builder.unloaded_q(q_val);
    }
    let inner = builder
        .synthesize()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyFilterDesign { inner })
}

/// Creates a band-pass filter design from passband edge frequencies.
///
/// This is an alternative to `bandpass()` for specifying the filter by
/// its start and stop frequencies instead of center and bandwidth.
///
/// Args:
///     order: Number of resonators
///     rl: Return loss in dB
///     f_start: Lower passband edge in Hz
///     f_stop: Upper passband edge in Hz
///     zeros: Transmission zero frequencies in Hz (optional)
///     q: Unloaded Q factor (optional)
///
/// Returns:
///     FilterDesign object
#[pyfunction]
#[pyo3(signature = (order, rl, f_start, f_stop, zeros=None, q=None))]
fn bandpass_edges(
    order: usize,
    rl: f64,
    f_start: f64,
    f_stop: f64,
    zeros: Option<Vec<f64>>,
    q: Option<f64>,
) -> PyResult<PyFilterDesign> {
    let mut builder = FilterDesign::bandpass_edges(order, rl, f_start, f_stop);
    if let Some(z) = zeros {
        builder = builder.zeros_hz(z);
    }
    if let Some(q_val) = q {
        builder = builder.unloaded_q(q_val);
    }
    let inner = builder
        .synthesize()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyFilterDesign { inner })
}

/// Creates a normalized prototype filter design.
///
/// Args:
///     order: Number of resonators
///     rl: Return loss in dB
///     zeros: Normalized transmission zeros (optional)
///     q: Unloaded Q factor (optional)
///
/// Returns:
///     FilterDesign object
#[pyfunction]
#[pyo3(signature = (order, rl, zeros=None, q=None))]
fn prototype(
    order: usize,
    rl: f64,
    zeros: Option<Vec<f64>>,
    q: Option<f64>,
) -> PyResult<PyFilterDesign> {
    let mut builder = FilterDesign::prototype(order, rl);
    if let Some(z) = zeros {
        builder = builder.zeros(z);
    }
    if let Some(q_val) = q {
        builder = builder.unloaded_q(q_val);
    }
    let inner = builder
        .synthesize()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyFilterDesign { inner })
}

/// Extracts a coupling matrix from a measured S-parameter response.
///
/// Builds a crude inline-chain initial guess for an `order`-resonator
/// transversal filter and optimizes it against the measured response. The
/// recovered matrix shares the inline-chain topology of the guess (the
/// optimizer freezes the template's non-zero pattern), so this is the right
/// entry point for transversal designs; a folded or arrow target is not
/// representable from this start.
///
/// Args:
///     freqs: List of sample frequencies (the values used as the normalized
///         prototype grid the optimizer evaluates on).
///     s21_re: Real parts of the measured S21 at each frequency.
///     s21_im: Imaginary parts of the measured S21 at each frequency.
///     s11_re: Real parts of the measured S11 at each frequency.
///     s11_im: Imaginary parts of the measured S11 at each frequency.
///     order: Filter order (number of resonators).
///
/// Returns:
///     dict with keys: matrix (2D list), cost (float), iterations (int),
///     converged (bool).
#[pyfunction]
#[pyo3(signature = (freqs, s21_re, s21_im, s11_re, s11_im, order))]
fn extract(
    py: Python<'_>,
    freqs: Vec<f64>,
    s21_re: Vec<f64>,
    s21_im: Vec<f64>,
    s11_re: Vec<f64>,
    s11_im: Vec<f64>,
    order: usize,
) -> PyResult<PyObject> {
    use crate::freq::FrequencyGrid;
    use crate::optimize::extract_matrix;
    use crate::response::{ResponseSample, SParameterResponse};
    use pyo3::types::PyDict;

    let n = freqs.len();
    if n < 2 {
        return Err(PyValueError::new_err("at least two frequency samples required"));
    }
    for (name, v) in [
        ("s21_re", &s21_re),
        ("s21_im", &s21_im),
        ("s11_re", &s11_re),
        ("s11_im", &s11_im),
    ] {
        if v.len() != n {
            return Err(PyValueError::new_err(format!(
                "{name} must have the same length as freqs ({n})"
            )));
        }
    }

    let target = SParameterResponse {
        samples: freqs
            .iter()
            .enumerate()
            .map(|(i, &f)| ResponseSample {
                frequency_hz: f,
                normalized_omega: f,
                group_delay: None,
                s11_re: s11_re[i],
                s11_im: s11_im[i],
                s21_re: s21_re[i],
                s21_im: s21_im[i],
            })
            .collect(),
    };

    let grid = FrequencyGrid::linspace(freqs[0], *freqs.last().unwrap(), n)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let result = extract_matrix(&target, &grid, order)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let dict = PyDict::new(py);
    dict.set_item("matrix", matrix_to_pylist(py, &result.matrix))
        .unwrap();
    dict.set_item("cost", result.cost).unwrap();
    dict.set_item("iterations", result.iterations).unwrap();
    dict.set_item("converged", result.converged).unwrap();
    Ok(dict.into())
}

/// Python module definition.
// ── Duplexer bindings ───────────────────────────────────────────────
use crate::duplexer::Duplexer;

/// Python wrapper for Duplexer.
#[pyclass(name = "Duplexer")]
struct PyDuplexer {
    inner: Duplexer,
}

#[pymethods]
impl PyDuplexer {
    /// Compute duplexer 3-port response.
    ///
    /// Returns dict with keys: frequencies, tx_s21_db, rx_s21_db,
    /// tx_s11_db, rx_s11_db, antenna_s11_db, isolation_db.
    #[pyo3(signature = (start, stop, points=401))]
    fn response(&self, start: f64, stop: f64, points: usize) -> PyResult<PyObject> {
        let resp = self
            .inner
            .response(start, stop, points)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Python::with_gil(|py| {
            use pyo3::types::PyDict;
            let d = PyDict::new(py);
            d.set_item("frequencies", resp.frequencies.clone())?;
            d.set_item("tx_s21_db", resp.tx_s21_db.clone())?;
            d.set_item("rx_s21_db", resp.rx_s21_db.clone())?;
            d.set_item("tx_s11_db", resp.tx_s11_db.clone())?;
            d.set_item("rx_s11_db", resp.rx_s11_db.clone())?;
            d.set_item("antenna_s11_db", resp.antenna_s11_db.clone())?;
            d.set_item("isolation_db", resp.isolation_db.clone())?;
            Ok(d.into())
        })
    }

    /// Get Tx filter coupling matrix (folded).
    fn tx_matrix(&self, py: Python<'_>) -> PyObject {
        matrix_to_pylist(py, self.inner.tx.matrix())
    }

    /// Get Rx filter coupling matrix (folded).
    fn rx_matrix(&self, py: Python<'_>) -> PyObject {
        matrix_to_pylist(py, self.inner.rx.matrix())
    }

    fn __repr__(&self) -> String {
        format!(
            "Duplexer(tx_order={}, rx_order={})",
            self.inner.tx.order(),
            self.inner.rx.order()
        )
    }
}

/// Create a duplexer from Tx and Rx bandpass specifications.
///
/// Args:
///     tx_order, tx_rl, tx_center, tx_bw: Tx filter parameters
///     rx_order, rx_rl, rx_center, rx_bw: Rx filter parameters
///     tx_zeros, rx_zeros: Optional transmission zeros (Hz)
///
/// Returns:
///     Duplexer object
#[pyfunction]
#[pyo3(signature = (tx_order, tx_rl, tx_center, tx_bw, rx_order, rx_rl, rx_center, rx_bw, tx_zeros=None, rx_zeros=None))]
fn duplexer(
    tx_order: usize,
    tx_rl: f64,
    tx_center: f64,
    tx_bw: f64,
    rx_order: usize,
    rx_rl: f64,
    rx_center: f64,
    rx_bw: f64,
    tx_zeros: Option<Vec<f64>>,
    rx_zeros: Option<Vec<f64>>,
) -> PyResult<PyDuplexer> {
    let mut tx_builder = FilterDesign::bandpass(tx_order, tx_rl, tx_center, tx_bw);
    if let Some(z) = tx_zeros {
        tx_builder = tx_builder.zeros_hz(z);
    }
    let tx = tx_builder
        .synthesize()
        .map_err(|e| PyValueError::new_err(format!("Tx: {e}")))?;

    let mut rx_builder = FilterDesign::bandpass(rx_order, rx_rl, rx_center, rx_bw);
    if let Some(z) = rx_zeros {
        rx_builder = rx_builder.zeros_hz(z);
    }
    let rx = rx_builder
        .synthesize()
        .map_err(|e| PyValueError::new_err(format!("Rx: {e}")))?;

    let dpx = Duplexer::new(tx, rx).map_err(|e| PyValueError::new_err(format!("{e}")))?;
    Ok(PyDuplexer { inner: dpx })
}

// ── Multiplexer bindings ─────────────────────────────────────────────

use crate::multiplexer::{Multiplexer, MultiplexerSynthesis};

/// Python wrapper for an N-channel Multiplexer.
#[pyclass(name = "Multiplexer")]
struct PyMultiplexer {
    inner: Multiplexer,
}

#[pymethods]
impl PyMultiplexer {
    /// Number of channels.
    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }

    /// Compute the (N+1)-port response.
    ///
    /// Returns dict with keys: frequencies, antenna_s11_db, channel_s21_db,
    /// channel_s11_db, isolation_db.
    #[pyo3(signature = (start, stop, points=401))]
    fn response(&self, start: f64, stop: f64, points: usize) -> PyResult<PyObject> {
        let resp = self
            .inner
            .response(start, stop, points)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Python::with_gil(|py| {
            use pyo3::types::PyDict;
            let d = PyDict::new(py);
            d.set_item("frequencies", resp.frequencies.clone())?;
            d.set_item("antenna_s11_db", resp.antenna_s11_db.clone())?;
            d.set_item("channel_s21_db", resp.channel_s21_db.clone())?;
            d.set_item("channel_s11_db", resp.channel_s11_db.clone())?;
            d.set_item("isolation_db", resp.isolation_db.clone())?;
            Ok(d.into())
        })
    }

    /// Get the i-th channel's folded coupling matrix as a 2D list.
    fn channel_folded(&self, py: Python<'_>, i: usize) -> PyResult<PyObject> {
        let m = self
            .inner
            .channel_folded(i)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    fn __repr__(&self) -> String {
        format!("Multiplexer(channels={})", self.inner.len())
    }
}

/// Create an N-channel multiplexer from repeated channel specs.
///
/// Args:
///     channels: List of channel specs, each formatted as
///         `"order:return_loss_db:center_hz:bandwidth_hz"`.
///     max_iterations: Maximum number of convergence steps (default 40).
///     tol: Convergence tolerance on root displacement (default 1e-8).
///
/// Returns:
///     Multiplexer object
#[pyfunction]
#[pyo3(signature = (channels, max_iterations=40, tol=1e-8))]
fn multiplexer(channels: Vec<String>, max_iterations: usize, tol: f64) -> PyResult<PyMultiplexer> {
    let mut specs = Vec::with_capacity(channels.len());
    for raw in &channels {
        let parts: Vec<&str> = raw.split(':').collect();
        if parts.len() != 4 {
            return Err(PyValueError::new_err(format!(
                "channel spec '{raw}' must be formatted as order:return_loss_db:center_hz:bandwidth_hz"
            )));
        }
        let order = parts[0].trim().parse::<usize>().map_err(|_| {
            PyValueError::new_err(format!("invalid order '{}' in '{raw}'", parts[0]))
        })?;
        let return_loss_db = parts[1].trim().parse::<f64>().map_err(|_| {
            PyValueError::new_err(format!("invalid return loss '{}' in '{raw}'", parts[1]))
        })?;
        let center_hz = parts[2].trim().parse::<f64>().map_err(|_| {
            PyValueError::new_err(format!(
                "invalid center frequency '{}' in '{raw}'",
                parts[2]
            ))
        })?;
        let bandwidth_hz = parts[3].trim().parse::<f64>().map_err(|_| {
            PyValueError::new_err(format!("invalid bandwidth '{}' in '{raw}'", parts[3]))
        })?;
        specs.push(crate::duplexer::ChannelSpec {
            order,
            return_loss_db,
            zeros: Vec::new(),
            center_hz,
            bandwidth_hz,
        });
    }

    let result = MultiplexerSynthesis::new()
        .channels(specs)
        .max_iterations(max_iterations)
        .convergence_tol(tol)
        .synthesize()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let inner =
        Multiplexer::new(result.channels).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyMultiplexer { inner })
}

#[pymodule]
pub fn mfs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFilterDesign>()?;
    m.add_class::<PyDuplexer>()?;
    m.add_class::<PyMultiplexer>()?;
    m.add_function(wrap_pyfunction!(bandpass, m)?)?;
    m.add_function(wrap_pyfunction!(bandpass_edges, m)?)?;
    m.add_function(wrap_pyfunction!(prototype, m)?)?;
    m.add_function(wrap_pyfunction!(duplexer, m)?)?;
    m.add_function(wrap_pyfunction!(multiplexer, m)?)?;
    m.add_function(wrap_pyfunction!(extract, m)?)?;
    Ok(())
}
