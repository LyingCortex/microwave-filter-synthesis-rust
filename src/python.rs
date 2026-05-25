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

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

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
            PyList::new_bound(py, &row_data).into()
        })
        .collect::<Vec<PyObject>>();
    PyList::new_bound(py, &rows).into()
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

    /// Returns the transversal coupling matrix as a 2D list.
    fn matrix(&self, py: Python<'_>) -> PyObject {
        matrix_to_pylist(py, self.inner.matrix())
    }

    /// Returns the folded coupling matrix as a 2D list.
    fn folded(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self.inner.to_folded().map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns the arrow coupling matrix as a 2D list.
    fn arrow(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self.inner.to_arrow().map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Returns the band-pass scaled matrix as a 2D list.
    fn scaled_matrix(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = self.inner.scale().map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(matrix_to_pylist(py, &m))
    }

    /// Computes S-parameter response.
    ///
    /// Returns (frequencies, s21_db, s11_db) as three lists.
    ///
    /// For bandpass designs: `response(start_hz, stop_hz, points)`
    /// For prototype designs: `response_normalized(start, stop, points)`
    #[pyo3(signature = (start, stop, points=201))]
    fn response(&self, start: f64, stop: f64, points: usize) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let resp = if self.inner.center_hz().is_some() {
            self.inner.response(start, stop, points)
        } else {
            self.inner.response_normalized(start, stop, points)
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
    /// s21_phase, s11_phase, group_delay.
    #[pyo3(signature = (start, stop, points=201))]
    fn response_full(&self, py: Python<'_>, start: f64, stop: f64, points: usize) -> PyResult<PyObject> {
        use pyo3::types::{PyDict, PyList};

        let resp = if self.inner.center_hz().is_some() {
            self.inner.response(start, stop, points)
        } else {
            self.inner.response_normalized(start, stop, points)
        };
        let resp = resp.map_err(|e| PyValueError::new_err(e.to_string()))?;

        let result: Vec<_> = resp.samples.iter().map(|s| {
            let dict = PyDict::new_bound(py);
            dict.set_item("freq", s.frequency_hz).unwrap();
            dict.set_item("s21_db", s.s21_db()).unwrap();
            dict.set_item("s11_db", s.s11_db()).unwrap();
            dict.set_item("s21_mag", s.s21_mag()).unwrap();
            dict.set_item("s11_mag", s.s11_mag()).unwrap();
            dict.set_item("s21_phase", s.s21_phase_deg()).unwrap();
            dict.set_item("s11_phase", s.s11_phase_deg()).unwrap();
            dict.set_item("group_delay", s.group_delay).unwrap();
            dict.into()
        }).collect::<Vec<PyObject>>();

        Ok(PyList::new_bound(py, &result).into())
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
    fn to_touchstone(&self, freq_unit: &str, format: &str, impedance: f64, version: u8) -> PyResult<String> {
        use crate::touchstone::{FreqUnit, DataFormat, TouchstoneVersion, TouchstoneConfig};

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
        let ver = if version >= 2 { TouchstoneVersion::V2 } else { TouchstoneVersion::V1 };

        let response = self.inner.default_response()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut config = TouchstoneConfig { freq_unit: fu, format: df, impedance, version: ver, comments: Vec::new() };
        config.comments = self.inner.auto_comments();
        crate::touchstone::to_touchstone_string(&response, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Saves the Touchstone file to disk.
    #[pyo3(signature = (path, freq_unit="GHz", format="RI", impedance=50.0, version=1))]
    fn save_touchstone(&self, path: &str, freq_unit: &str, format: &str, impedance: f64, version: u8) -> PyResult<()> {
        use crate::touchstone::{FreqUnit, DataFormat, TouchstoneVersion, TouchstoneConfig};

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
        let ver = if version >= 2 { TouchstoneVersion::V2 } else { TouchstoneVersion::V1 };

        let response = self.inner.default_response()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut config = TouchstoneConfig { freq_unit: fu, format: df, impedance, version: ver, comments: Vec::new() };
        config.comments = self.inner.auto_comments();
        crate::touchstone::write_touchstone(&response, &config, path)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Computes lossy S-parameter response with finite unloaded Q.
    ///
    /// Returns (frequencies, s21_db, s11_db).
    #[pyo3(signature = (start, stop, points=201, q=1000.0))]
    fn response_lossy(&self, start: f64, stop: f64, points: usize, q: f64) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {
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
    /// target_s21 and target_s11 are lists of complex values (re, im) pairs
    /// at the same frequency points as the grid.
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
        use pyo3::types::{PyDict, PyList};
        use crate::response::{ResponseSample, SParameterResponse};
        use crate::freq::FrequencyGrid;
        use crate::optimize::{tune_matrix, OptimizeConfig, optimize_matrix};

        if freqs.len() != target_s21_re.len() {
            return Err(PyValueError::new_err("all arrays must have the same length"));
        }

        let target = SParameterResponse {
            samples: freqs.iter().enumerate().map(|(i, &f)| ResponseSample {
                frequency_hz: f,
                normalized_omega: f,
                group_delay: 0.0,
                s11_re: target_s11_re[i],
                s11_im: target_s11_im[i],
                s21_re: target_s21_re[i],
                s21_im: target_s21_im[i],
            }).collect(),
        };

        let grid = FrequencyGrid::linspace(
            *freqs.first().unwrap_or(&-2.0),
            *freqs.last().unwrap_or(&2.0),
            freqs.len(),
        ).map_err(|e| PyValueError::new_err(e.to_string()))?;

        let config = OptimizeConfig { max_iterations: max_iter, ..Default::default() };
        let result = optimize_matrix(self.inner.matrix(), &target, &grid, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let dict = PyDict::new_bound(py);
        dict.set_item("matrix", matrix_to_pylist(py, &result.matrix)).unwrap();
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
    let inner = builder.synthesize().map_err(|e| PyValueError::new_err(e.to_string()))?;
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
    let inner = builder.synthesize().map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyFilterDesign { inner })
}

/// Python module definition.
#[pymodule]
pub fn mfs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFilterDesign>()?;
    m.add_function(wrap_pyfunction!(bandpass, m)?)?;
    m.add_function(wrap_pyfunction!(prototype, m)?)?;
    Ok(())
}
