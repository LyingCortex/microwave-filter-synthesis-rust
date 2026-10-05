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
}
