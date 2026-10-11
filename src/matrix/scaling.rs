use crate::error::{MfsError, Result};
use crate::freq::{BandPassMapping, FrequencyMapping};

use super::core::{BandPassScaledCouplingMatrix, CouplingMatrix};

impl CouplingMatrix {
    /// Scales a normalized band-pass matrix into physical-frequency units.
    ///
    /// Internal resonator couplings and diagonal terms are converted into Hz.
    /// Source/load couplings remain as couplings, matching the Python helper's
    /// bandwidth-scaled mode rather than the external-Q representation.
    pub fn denormalize_bandpass(&self, mapping: &BandPassMapping) -> Result<Self> {
        let center_hz = mapping.center_hz();
        let bandwidth_hz = mapping.bandwidth_hz();
        let side = self.side();
        let mut data = vec![0.0; side * side];

        for row in 0..side {
            for col in 0..side {
                let value = self.get(row, col);
                let scaled = if row == col {
                    if row == 0 || row == side - 1 {
                        0.0
                    } else {
                        denormalize_resonator_frequency(value, center_hz, bandwidth_hz)
                    }
                } else {
                    value * bandwidth_hz
                };
                data[row * side + col] = scaled;
            }
        }

        Self::new_with_topology(self.order(), self.topology(), data)
    }

    /// Converts a normalized matrix into a physical band-pass matrix plus port Q values.
    pub fn denormalize_bandpass_with_external_q(
        &self,
        mapping: &BandPassMapping,
    ) -> Result<BandPassScaledCouplingMatrix> {
        let mut matrix_hz = self.denormalize_bandpass(mapping)?;
        let fractional_bw = mapping.bandwidth_hz() / mapping.center_hz();
        let source_coupling = self.get(0, 1);
        let load_coupling = self.get(self.order(), self.side() - 1);
        let source_external_q =
            external_q_from_normalized_coupling(source_coupling, fractional_bw)?;
        let load_external_q = external_q_from_normalized_coupling(load_coupling, fractional_bw)?;

        matrix_hz.set_entry(0, 1, source_external_q);
        matrix_hz.set_entry(1, 0, source_external_q);
        matrix_hz.set_entry(self.order(), self.side() - 1, load_external_q);
        matrix_hz.set_entry(self.side() - 1, self.order(), load_external_q);

        Ok(BandPassScaledCouplingMatrix {
            matrix_hz,
            source_external_q,
            load_external_q,
        })
    }

    /// Converts a physical band-pass matrix back into normalized units.
    ///
    /// This expects source/load entries to still be couplings, not external Q values.
    pub fn normalize_bandpass(&self, mapping: &BandPassMapping) -> Result<Self> {
        let center_hz = mapping.center_hz();
        let bandwidth_hz = mapping.bandwidth_hz();
        let side = self.side();
        let mut data = vec![0.0; side * side];

        for row in 0..side {
            for col in 0..side {
                let value = self.get(row, col);
                let normalized = if row == col {
                    if row == 0 || row == side - 1 {
                        0.0
                    } else {
                        normalize_resonator_frequency(value, center_hz, bandwidth_hz)?
                    }
                } else {
                    value / bandwidth_hz
                };
                data[row * side + col] = normalized;
            }
        }

        Self::new_with_topology(self.order(), self.topology(), data)
    }

    /// Converts a physical band-pass matrix that stores external Q values back into normalized form.
    pub fn normalize_bandpass_with_external_q(&self, mapping: &BandPassMapping) -> Result<Self> {
        let mut normalized = self.normalize_bandpass(mapping)?;
        let fractional_bw = mapping.bandwidth_hz() / mapping.center_hz();
        let source_q = self.get(0, 1);
        let load_q = self.get(self.order(), self.side() - 1);

        let source_coupling = normalized_coupling_from_external_q(source_q, fractional_bw)?;
        let load_coupling = normalized_coupling_from_external_q(load_q, fractional_bw)?;
        normalized.set_entry(0, 1, source_coupling);
        normalized.set_entry(1, 0, source_coupling);
        normalized.set_entry(self.order(), self.side() - 1, load_coupling);
        normalized.set_entry(self.side() - 1, self.order(), load_coupling);

        Ok(normalized)
    }
}

/// Inverse band-pass map for a resonator rest frequency: the positive-
/// frequency branch, exactly matching
/// [`BandPassMapping::map_normalized_to_hz`] (the same quadratic, same branch).
/// The historical form subtracted the linear term, which is the *negative-
/// frequency* branch — `denorm(n) == freq.denorm(-n)` — mirroring every
/// physical detuning sign (review P1-2).
fn denormalize_resonator_frequency(normalized: f64, center_hz: f64, bandwidth_hz: f64) -> f64 {
    let t = normalized * bandwidth_hz;
    (t + (t * t + 4.0 * center_hz * center_hz).sqrt()) / 2.0
}

/// Normalized detuning of a resonator rest frequency: exactly
/// [`BandPassMapping::map_hz_to_normalized`], so the two cannot drift.
fn normalize_resonator_frequency(
    physical_hz: f64,
    center_hz: f64,
    bandwidth_hz: f64,
) -> Result<f64> {
    if !physical_hz.is_finite() || physical_hz <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "physical resonator frequency must be > 0, got {physical_hz}"
        )));
    }

    let mapping = BandPassMapping::new(center_hz, bandwidth_hz)?;
    Ok(mapping.map_hz_to_normalized(physical_hz)?.omega)
}

fn external_q_from_normalized_coupling(coupling: f64, fractional_bw: f64) -> Result<f64> {
    if !coupling.is_finite() || coupling.abs() <= 1e-12 {
        return Err(MfsError::InvalidFrequency(
            "normalized source/load coupling must be non-zero when converting to external Q"
                .to_string(),
        ));
    }

    Ok(1.0 / (coupling * coupling * fractional_bw))
}

fn normalized_coupling_from_external_q(external_q: f64, fractional_bw: f64) -> Result<f64> {
    if !external_q.is_finite() || external_q <= 0.0 {
        return Err(MfsError::InvalidFrequency(format!(
            "external Q must be > 0, got {external_q}"
        )));
    }

    Ok((1.0 / (external_q * fractional_bw)).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The resonator rest-frequency denorm/normalize must agree with the
    /// canonical `BandPassMapping` inverse/forward maps (review P1-2): the
    /// denorm is the *positive-frequency* branch, not its sign-mirrored twin.
    #[test]
    fn resonator_frequency_maps_agree_with_bandpass_mapping() {
        let center_hz = 6.75e9;
        let bandwidth_hz = 300.0e6;
        let mapping = BandPassMapping::new(center_hz, bandwidth_hz).unwrap();

        for &normalized in &[-3.0, -1.5, -0.75, -0.1, 0.0, 0.1, 0.75, 1.5, 3.0] {
            let denorm =
                denormalize_resonator_frequency(normalized, center_hz, bandwidth_hz);
            let canonical = mapping
                .map_normalized_to_hz(crate::freq::NormalizedSample { omega: normalized })
                .unwrap();
            assert!(
                (denorm - canonical).abs() <= 1e-9 * canonical.max(1.0),
                "denorm({normalized}) = {denorm} vs {canonical}"
            );

            let re_norm = normalize_resonator_frequency(denorm, center_hz, bandwidth_hz)
                .unwrap();
            assert!(
                (re_norm - normalized).abs() <= 1e-9,
                "normalize({denorm}) = {re_norm} vs {normalized}"
            );

            // Cross-direction: forward map of the denormed frequency must
            // return the original normalized detuning.
            let forward = mapping.map_hz_to_normalized(denorm).unwrap().omega;
            assert!(
                (forward - normalized).abs() <= 1e-9,
                "forward({denorm}) = {forward} vs {normalized}"
            );
        }
    }

    /// §C.4: the external-Q convention `Qe = 1/(M01²·Δ)` (vs the "standard"
    /// `1/(M01·Δ)` from some textbooks) is anchored to a physical order-1
    /// bandpass. In this code's coupling convention the source/load coupling
    /// `M01` is `sqrt(rate)`, so the *fractional* external coupling *rate*
    /// is `M01²·Δ` and the external Q is its reciprocal. This is verified
    /// against the measured physical external Q of a single-resonator,
    /// narrowband (Δ = 0.01, so the half-power points sit inside the
    /// passband) bandpass: for a symmetric, lossless, critically-coupled
    /// resonator, the loaded Q is `f0/BW_3dB` and `Qse = Qle = 2·Q_L`.
    ///
    /// The measurement (RL = 20 dB, f0 = 1 GHz, Δ = 0.01) gives
    /// `Qse_meas = 2f0/BW_3dB = 20.10` while `1/(M01²·Δ) = 20.10`
    /// (ratio 1.0000) and `1/(M01·Δ) = 44.83` (ratio 0.448, off by M01 = 2.23).
    /// The `g²` form is the physically correct one in this code's convention;
    /// the "standard" `1/(M01·Δ)` would be off by the coupling magnitude.
    /// (A maximally-flat order-1 bandpass has a return-loss-free, flat
    /// passband, so its -3 dB points coincide with the half-power points.)
    #[test]
    fn external_q_convention_matches_physical_order1_bandpass() -> Result<()> {
        use crate::design::FilterDesign;
        let f0 = 1.0e9;
        let bw = 1.0e7; // Δ = 0.01, narrowband
        let delta = bw / f0;
        let design = FilterDesign::bandpass(1, 20.0, f0, bw).synthesize()?;
        let m = design.matrix();
        let m01 = m.at(0, 1).ok_or_else(|| {
            MfsError::PreconditionViolation(
                "order-1 bandpass matrix has no source coupling".into(),
            )
        })?;
        let mapping = BandPassMapping::new(f0, bw)?;

        // The physical external Q measured from the response.
        let grid = crate::freq::FrequencyGrid::linspace(f0 * 0.90, f0 * 1.10, 20001)?;
        let solver = crate::response::ResponseSolver;
        let resp = solver.evaluate(m, &grid, &mapping)?;
        let mag2 = |i: usize| {
            resp.samples[i].s21_re * resp.samples[i].s21_re + resp.samples[i].s21_im.powi(2)
        };
        let peak = (0..resp.samples.len()).map(mag2).fold(0.0_f64, f64::max);
        let level = 0.5 * peak;
        let center = resp.samples.len() / 2;
        let mut lo: Option<f64> = None;
        let mut hi: Option<f64> = None;
        for i in 0..center {
            if lo.is_none() && mag2(i) < level && mag2(i + 1) >= level {
                lo = Some(resp.samples[i + 1].frequency_hz);
            }
        }
        for i in center..resp.samples.len() - 1 {
            if hi.is_none() && mag2(i) >= level && mag2(i + 1) < level {
                hi = Some(resp.samples[i + 1].frequency_hz);
            }
        }
        let (lo, hi) = match (lo, hi) {
            (Some(l), Some(h)) => (l, h),
            _ => {
                return Err(MfsError::NumericalFailure(format!(
                    "-3dB points not found in [0.9,1.1]f0 (peak={peak:.3e})"
                )))
            }
        };
        let qse_meas = 2.0 * f0 / (hi - lo); // Qse = Qle = 2·Q_L (symmetric)

        // The code's convention `1/(M01²·Δ)` must agree with the measurement.
        let qse_code = 1.0 / (m01 * m01 * delta);
        assert!(
            ((qse_meas - qse_code).abs() / qse_meas) < 1.0e-2,
            "code external Q 1/(M01²·Δ)={qse_code:.4} disagrees with measured \
             Qse=2f0/BW3dB={qse_meas:.4} (rel diff {})",
            (qse_meas - qse_code).abs() / qse_meas
        );
        // The "standard" 1/(M01·Δ) must NOT agree (it would be off by M01).
        let qse_std = 1.0 / (m01 * delta);
        assert!(
            ((qse_meas - qse_std).abs() / qse_meas) > 1.0,
            "the 'standard' 1/(M01·Δ)={qse_std:.4} should be off by M01 from \
             the measured Qse={qse_meas:.4}; if this ever passes the convention \
             has changed and the g² anchor must be re-derived"
        );
        Ok(())
    }
}
