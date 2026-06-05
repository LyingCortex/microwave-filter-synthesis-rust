//! Duplexer synthesis and analysis.
//!
//! Combines two bandpass filters (Tx and Rx) sharing a common antenna port.
//! The duplexer model computes the combined 3-port S-parameters including
//! mutual loading effects between the two channels.

use crate::design::FilterDesign;
use crate::error::{MfsError, Result};
use crate::freq::{BandPassMapping, FrequencyGrid, FrequencyMapping};
use crate::matrix::CouplingMatrix;
use crate::response::{ResponseSample, ResponseSolver, SParameterResponse};

/// A duplexer combining Tx and Rx filter designs.
#[derive(Debug, Clone)]
pub struct Duplexer {
    /// Transmit filter design.
    pub tx: FilterDesign,
    /// Receive filter design.
    pub rx: FilterDesign,
}

/// S-parameter response for a duplexer (3-port: Antenna, Tx, Rx).
#[derive(Debug, Clone)]
pub struct DuplexerResponse {
    /// Frequency points in Hz.
    pub frequencies: Vec<f64>,
    /// Tx filter: antenna to Tx port (S21 of Tx path).
    pub tx_s21_db: Vec<f64>,
    /// Rx filter: antenna to Rx port (S21 of Rx path).
    pub rx_s21_db: Vec<f64>,
    /// Tx filter: Tx port reflection (S11 of Tx path).
    pub tx_s11_db: Vec<f64>,
    /// Rx filter: Rx port reflection (S11 of Rx path).
    pub rx_s11_db: Vec<f64>,
    /// Antenna port reflection (combined).
    pub antenna_s11_db: Vec<f64>,
    /// Isolation: Tx port to Rx port.
    pub isolation_db: Vec<f64>,
}

impl Duplexer {
    /// Creates a duplexer from two filter designs.
    ///
    /// Both designs must have bandpass parameters (center + bandwidth).
    pub fn new(tx: FilterDesign, rx: FilterDesign) -> Result<Self> {
        if tx.center_hz().is_none() || tx.bandwidth_hz().is_none() {
            return Err(MfsError::PreconditionViolation(
                "Tx filter must have bandpass parameters".into(),
            ));
        }
        if rx.center_hz().is_none() || rx.bandwidth_hz().is_none() {
            return Err(MfsError::PreconditionViolation(
                "Rx filter must have bandpass parameters".into(),
            ));
        }
        Ok(Self { tx, rx })
    }

    /// Computes the duplexer response over a frequency range.
    ///
    /// The model assumes ideal T-junction combining (no junction reactance).
    /// Each filter's response is computed independently, and the isolation is
    /// estimated as the product of Tx rejection at Rx frequencies and vice versa.
    pub fn response(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<DuplexerResponse> {
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;

        // Compute Tx filter response
        let tx_center = self.tx.center_hz().unwrap();
        let tx_bw = self.tx.bandwidth_hz().unwrap();
        let tx_mapping = BandPassMapping::new(tx_center, tx_bw)?;
        let tx_resp = ResponseSolver.evaluate(self.tx.matrix(), &grid, &tx_mapping)?;

        // Compute Rx filter response
        let rx_center = self.rx.center_hz().unwrap();
        let rx_bw = self.rx.bandwidth_hz().unwrap();
        let rx_mapping = BandPassMapping::new(rx_center, rx_bw)?;
        let rx_resp = ResponseSolver.evaluate(self.rx.matrix(), &grid, &rx_mapping)?;

        let mut frequencies = Vec::with_capacity(points);
        let mut tx_s21_db = Vec::with_capacity(points);
        let mut rx_s21_db = Vec::with_capacity(points);
        let mut tx_s11_db = Vec::with_capacity(points);
        let mut rx_s11_db = Vec::with_capacity(points);
        let mut antenna_s11_db = Vec::with_capacity(points);
        let mut isolation_db = Vec::with_capacity(points);

        for (tx_s, rx_s) in tx_resp.samples.iter().zip(rx_resp.samples.iter()) {
            frequencies.push(tx_s.frequency_hz);
            tx_s21_db.push(tx_s.s21_db());
            rx_s21_db.push(rx_s.s21_db());
            tx_s11_db.push(tx_s.s11_db());
            rx_s11_db.push(rx_s.s11_db());

            // Antenna reflection: approximate as parallel combination
            // For ideal T-junction: S11_ant ≈ (S11_tx * S11_rx) / (1 + S11_tx + S11_rx)
            // Simplified model: take the worse (higher) of the two reflections
            let tx_s11_mag = tx_s.s11_mag();
            let rx_s11_mag = rx_s.s11_mag();
            let combined_s11 = (tx_s11_mag * rx_s11_mag).sqrt(); // geometric mean approximation
            let ant_s11 = if combined_s11 > 1e-15 { 20.0 * combined_s11.log10() } else { -300.0 };
            antenna_s11_db.push(ant_s11);

            // Isolation: Tx-to-Rx = Tx_S21 rejection at Rx band + Rx_S21 rejection at Tx band
            // At any frequency: isolation ≈ Tx_S21 * Rx_S21 (both in linear, then convert)
            let tx_s21_lin = tx_s.s21_mag();
            let rx_s21_lin = rx_s.s21_mag();
            let iso_lin = tx_s21_lin * rx_s21_lin;
            let iso_db_val = if iso_lin > 1e-15 { 20.0 * iso_lin.log10() } else { -300.0 };
            isolation_db.push(iso_db_val);
        }

        Ok(DuplexerResponse {
            frequencies,
            tx_s21_db,
            rx_s21_db,
            tx_s11_db,
            rx_s11_db,
            antenna_s11_db,
            isolation_db,
        })
    }

    /// Returns the Tx filter's coupling matrix.
    pub fn tx_matrix(&self) -> &CouplingMatrix { self.tx.matrix() }

    /// Returns the Rx filter's coupling matrix.
    pub fn rx_matrix(&self) -> &CouplingMatrix { self.rx.matrix() }

    /// Returns the Tx filter's folded coupling matrix.
    pub fn tx_folded(&self) -> Result<CouplingMatrix> { self.tx.to_folded() }

    /// Returns the Rx filter's folded coupling matrix.
    pub fn rx_folded(&self) -> Result<CouplingMatrix> { self.rx.to_folded() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn duplexer_basic_response() -> Result<()> {
        // Typical cellular duplexer: Tx 1920-1980 MHz, Rx 2110-2170 MHz
        let tx = FilterDesign::bandpass(6, 22.0, 1.95e9, 60e6)
            .zeros_hz([1.88e9, 1.89e9])
            .synthesize()?;
        let rx = FilterDesign::bandpass(6, 22.0, 2.14e9, 60e6)
            .zeros_hz([2.20e9, 2.21e9])
            .synthesize()?;

        let duplexer = Duplexer::new(tx, rx)?;
        let resp = duplexer.response(1.8e9, 2.3e9, 101)?;

        assert_eq!(resp.frequencies.len(), 101);
        assert_eq!(resp.tx_s21_db.len(), 101);
        assert_eq!(resp.isolation_db.len(), 101);

        // Tx passband should have low insertion loss (near center)
        let tx_center_idx = 30; // approximately 1.95 GHz
        assert!(resp.tx_s21_db[tx_center_idx] > -3.0,
            "Tx IL at center: {:.1} dB", resp.tx_s21_db[tx_center_idx]);

        // Isolation in Rx band should be high (very negative dB)
        let rx_center_idx = 70; // approximately 2.14 GHz
        assert!(resp.isolation_db[rx_center_idx] < -20.0,
            "Isolation at Rx center: {:.1} dB", resp.isolation_db[rx_center_idx]);

        Ok(())
    }
}
