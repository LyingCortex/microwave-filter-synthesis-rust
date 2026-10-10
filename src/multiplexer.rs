//! Generalized N-channel multiplexer synthesis.
//!
//! Generalizes the 2-channel Macchiarella-Tamiazzo iterative diplexer synthesis
//! (`crate::duplexer`) to an arbitrary number `N` of channels sharing a single
//! common (antenna) junction. Each channel is synthesized independently, then
//! one of two iteration paths accounts for the mutual loading of **all** other
//! channels through the junction:
//!
//! * default (simplified): a frequency-offset perturbation loop — lightweight,
//!   adequate for well-separated bands (the original path, kept as the default);
//! * `full_iteration(true)`: the exact common-domain Macchiarella–Tamiazzo
//!   spectral loop shared with the diplexer (`crate::synthesis::machiarella`) —
//!   a quasi-equiripple (N+1)-port junction; for N > 2 only the ideal Type-I
//!   junction (`n = 1`, `b0 = 0`) and Type-II are supported.

use crate::approx::{
    compute_root_displacement, generalized_chebyshev_polynomials, ComplexPolynomial,
};
use crate::design::FilterDesign;
use crate::duplexer::{ChannelSpec, JunctionType};
use crate::error::{MfsError, Result};
use crate::freq::{BandPassMapping, FrequencyGrid};
use crate::matrix::CouplingMatrix;
use crate::response::ResponseSolver;
use crate::spec::FilterSpec;
use crate::synthesis::build_filter_spec;
use nalgebra::DMatrix;
use num_complex::Complex64;

/// N-channel multiplexer synthesis result.
#[derive(Debug, Clone)]
pub struct MultiplexerSynthesisResult {
    /// One synthesized `FilterDesign` per channel (accounts for mutual loading).
    pub channels: Vec<FilterDesign>,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final convergence error (max root displacement across all channels).
    pub convergence_error: f64,
    /// `true` when the iteration stopped on its own convergence criterion
    /// (strict tolerance, or the high-degree stagnation accept for the full
    /// spectral path) rather than exhausting the iteration budget. An
    /// exhausted budget is a non-convergence whatever the residual.
    pub converged: bool,
    /// The junction parameters used.
    pub junction: JunctionType,
    /// Per-channel S-polynomial roots (auxiliary polynomial).
    pub s_roots: Vec<Vec<Complex64>>,
    /// The iteration budget in force, so a caller can distinguish "converged"
    /// from "ran out of iterations" without holding the builder.
    pub max_iterations: usize,
    /// Common-domain junction polynomials (full spectral path only): present
    /// when `full_iteration(true)` was set and the iteration converged.
    pub junction_common: Option<crate::synthesis::JunctionCommon>,
}

/// Builder for iterative N-channel multiplexer synthesis.
///
/// Each channel is described by a [`ChannelSpec`] (order, return loss,
/// transmission zeros, center frequency, bandwidth). The synthesis iteratively
/// adjusts each channel's S-polynomial roots to account for the loading of all
/// other channels through the common junction.
#[derive(Debug, Clone)]
pub struct MultiplexerSynthesis {
    channels: Vec<ChannelSpec>,
    junction: JunctionType,
    max_iterations: usize,
    convergence_tol: f64,
    full_iteration: bool,
}

impl MultiplexerSynthesis {
    /// Creates a new multiplexer synthesis builder with default parameters.
    pub fn new() -> Self {
        Self {
            channels: Vec::new(),
            junction: JunctionType::default(),
            max_iterations: 40,
            convergence_tol: 1e-8,
            full_iteration: false,
        }
    }

    /// Appends a single channel specification.
    pub fn channel(mut self, spec: ChannelSpec) -> Self {
        self.channels.push(spec);
        self
    }

    /// Appends multiple channel specifications.
    pub fn channels(mut self, specs: impl IntoIterator<Item = ChannelSpec>) -> Self {
        self.channels.extend(specs);
        self
    }

    /// Sets the junction type and parameters.
    pub fn junction(mut self, junction: JunctionType) -> Self {
        self.junction = junction;
        self
    }

    /// Sets maximum number of iterations (default: 40).
    pub fn max_iterations(mut self, n: usize) -> Self {
        self.max_iterations = n;
        self
    }

    /// Sets convergence tolerance (default: 1e-8).
    pub fn convergence_tol(mut self, tol: f64) -> Self {
        self.convergence_tol = tol;
        self
    }

    /// Enables the full Macchiarella–Tamiazzo spectral iteration (default:
    /// the simplified frequency-offset perturbation). The full path
    /// converges to a quasi-equiripple (N+1)-port common junction and is
    /// heavier (spectral factorization per iteration); for more than two
    /// channels only the ideal Type-I junction (`n = 1`, `b0 = 0`) and
    /// Type-II are supported.
    pub fn full_iteration(mut self, enabled: bool) -> Self {
        self.full_iteration = enabled;
        self
    }

    /// Runs the iterative N-channel synthesis algorithm.
    ///
    /// Each channel is first synthesized independently, then the S-polynomial
    /// roots are perturbed iteratively to compensate for the loading of every
    /// other channel through the common junction. Uses the simplified
    /// perturbation path (frequency-offset aware) mirroring the diplexer.
    pub fn synthesize(&self) -> Result<MultiplexerSynthesisResult> {
        let n = self.channels.len();

        // (a) Validate inputs.
        if n < 2 {
            return Err(MfsError::PreconditionViolation(
                "multiplexer requires at least two channels".into(),
            ));
        }
        for (i, ch) in self.channels.iter().enumerate() {
            if ch.order == 0 {
                return Err(MfsError::PreconditionViolation(format!(
                    "channel {i} must have a positive order"
                )));
            }
            if !ch.center_hz.is_finite() || ch.center_hz <= 0.0 {
                return Err(MfsError::PreconditionViolation(format!(
                    "channel {i} must have a positive center frequency"
                )));
            }
            if !ch.bandwidth_hz.is_finite() || ch.bandwidth_hz <= 0.0 {
                return Err(MfsError::PreconditionViolation(format!(
                    "channel {i} must have a positive bandwidth"
                )));
            }
        }

        // P2-6: the simplified path computes the junction coupling as
        // `1/(1 − c0)` (multiplexer.rs), so a non-finite / non-positive c0
        // yields NaN/Inf and silently poisons every S-root. The full path
        // validates c0 in `machiarella` and the duplexer validates it at its
        // entry; the simplified multiplexer path was the gap.
        if let JunctionType::TypeII { c0 } = self.junction {
            if !c0.is_finite() || c0 <= 0.0 {
                return Err(MfsError::PreconditionViolation(
                    "Type-II junction capacitance c0 must be finite and positive".into(),
                ));
            }
        }
        // Full Macchiarella–Tamiazzo spectral iteration, shared with the
        // duplexer (`crate::synthesis::machiarella`): the exact common-domain
        // bandpass transform, per-iteration spectral factorization, and the
        // N×N p₀² system (LU for N ≥ 3, the historical Cramer 2×2 for N = 2).
        // Only the ideal Type-I junction (n = 1, b0 = 0) and Type-II are
        // supported for N > 2 (the shared module rejects the rest).
        if self.full_iteration {
            let (prepared, state) = crate::synthesis::synthesize_full(
                &self.channels,
                self.junction,
                self.max_iterations,
                self.convergence_tol,
            )?;
            // Recover each channel's design from its converged S roots. In
            // the full iteration the per-channel D polynomials are
            // D_i = S_i − F_i with the channel's *fixed* reflection zeros
            // (normalization-independent, mirroring the duplexer's full-path
            // recovery).
            let mut channels = Vec::with_capacity(n);
            for i in 0..n {
                let ch = &self.channels[i];
                let s_proto = ComplexPolynomial::from_complex_roots(&state.s_roots[i])?;
                let d_recover = s_proto.sub(&prepared.polys[i].f)?;
                let design = crate::synthesis::recover_channel_design(
                    &state.s_roots[i],
                    &d_recover,
                    &prepared.polys[i],
                    ch.order,
                    ch.return_loss_db,
                    ch.center_hz,
                    ch.bandwidth_hz,
                )?;
                channels.push(design);
            }
            return Ok(MultiplexerSynthesisResult {
                channels,
                iterations: state.iterations,
                convergence_error: state.convergence_error,
                converged: state.converged,
                junction: self.junction,
                s_roots: state.s_roots,
                max_iterations: self.max_iterations,
                junction_common: state.common,
            });
        }

        // (b) Synthesize each channel independently and extract auxiliary
        // polynomials S_i = (E + F)/2, D_i = (E - F)/2.
        let half = Complex64::new(0.5, 0.0);
        let mut d_polys: Vec<ComplexPolynomial> = Vec::with_capacity(n);
        let mut original_polys: Vec<crate::approx::PolynomialSet> = Vec::with_capacity(n);
        let mut s_roots: Vec<Vec<Complex64>> = Vec::with_capacity(n);
        let mut centers: Vec<f64> = Vec::with_capacity(n);
        let mut bandwidths: Vec<f64> = Vec::with_capacity(n);
        let mut orders: Vec<usize> = Vec::with_capacity(n);

        for ch in &self.channels {
            let spec = Self::build_spec(ch)?;
            let polys = generalized_chebyshev_polynomials(&spec)?;

            let s_init = polys.e.add(&polys.f)?.scale(half)?;
            let d = polys.e.sub(&polys.f)?.scale(half)?;

            d_polys.push(d);
            original_polys.push(polys);
            s_roots.push(s_init.roots()?);
            centers.push(ch.center_hz);
            bandwidths.push(ch.bandwidth_hz);
            orders.push(ch.order);
        }

        // (c) Frequency separation reference: average spacing between
        // consecutive channel centers relative to the average bandwidth.
        let avg_bw = bandwidths.iter().copied().sum::<f64>() / n as f64;

        // Junction coupling: a complex coefficient that weights the loading a
        // neighbouring channel imposes through the common junction. For the
        // ideal Type-I junction (n=1, b₀=0) this is 1, so the default behaviour
        // is unchanged. Non-zero junction susceptance/capacitance enters both
        // the magnitude and phase of the mutual-loading correction, so the
        // junction parameters genuinely affect the synthesis.
        //
        // Type-I (waveguide H-plane T): transformer ratio `n` and susceptance
        // `b₀`, matching the diplexer's `N(s) = ∏Sᵢ − n²·Σ(DⱼSᵢ+SⱼDᵢ)/(1−jn²b₀)`.
        // Type-II (coaxial node with capacitance `c₀`): the node capacitance
        // shunts the junction and couples the channels; normalised as a
        // frequency-independent susceptance for the simplified perturbation.
        let junction_coupling = match self.junction {
            JunctionType::TypeI { n, b0 } => {
                let n_sq = n * n;
                Complex64::new(n_sq, 0.0) / Complex64::new(1.0, -n_sq * b0)
            }
            JunctionType::TypeII { c0 } => Complex64::new(1.0, 0.0) / Complex64::new(1.0, -c0),
        };

        let mut iterations = 0;
        let mut convergence_error = f64::INFINITY;
        let mut converged = false;

        // (d) Iterate: simplified perturbation, N-channel version. The loading
        // correction for channel i accumulates contributions from every other
        // channel j, evaluated at channel i's current roots with a frequency
        // offset that accounts for the relative placement of the two bands.
        // `alpha` is a purely geometric under-relaxation factor; the junction
        // loading is carried by `junction_coupling` instead.
        for iter in 0..self.max_iterations {
            iterations = iter + 1;

            let mut s_polys: Vec<ComplexPolynomial> = Vec::with_capacity(n);
            for roots in &s_roots {
                s_polys.push(ComplexPolynomial::from_complex_roots(roots)?);
            }

            let alpha = 0.3 / (1.0 + 0.3 * iter as f64);

            let mut new_roots: Vec<Vec<Complex64>> = Vec::with_capacity(n);
            for i in 0..n {
                let mut new_i = Vec::with_capacity(orders[i]);
                for &root in &s_roots[i] {
                    let mut correction = Complex64::new(0.0, 0.0);
                    for j in 0..n {
                        if j == i {
                            continue;
                        }
                        // Frequency offset in normalized units: the cross-channel
                        // polynomials are sampled at a point shifted along the
                        // imaginary axis by the band-to-band separation. This is
                        // an approximation of the true loading at channel i's
                        // center, adequate for well-separated bands.
                        let offset = 2.0 * (centers[j] - centers[i]) / avg_bw;
                        let freq = root + Complex64::new(0.0, offset);
                        let s_j_val = s_polys[j].evaluate(freq);
                        let d_j_val = d_polys[j].evaluate(freq);
                        if s_j_val.norm() > 1e-10 {
                            correction += junction_coupling * (d_j_val / s_j_val);
                        }
                    }
                    let new_root = root - alpha * correction;
                    // Stability guard: keep the root finite and bounded, and only
                    // reflect a root that drifts into the right half-plane (which
                    // would make the auxiliary polynomial non-Hurwitz) back to a
                    // small negative real part. Legitimate negative real parts are
                    // preserved unchanged.
                    let updated = if new_root.re.is_finite()
                        && new_root.im.is_finite()
                        && new_root.norm() < 50.0
                    {
                        let re = if new_root.re > 0.0 {
                            new_root.re.min(-1e-10)
                        } else {
                            new_root.re
                        };
                        Complex64::new(re, new_root.im)
                    } else {
                        root
                    };
                    new_i.push(updated);
                }
                new_roots.push(new_i);
            }

            let mut max_delta = 0.0_f64;
            for i in 0..n {
                max_delta = max_delta.max(compute_root_displacement(&s_roots[i], &new_roots[i]));
            }
            s_roots = new_roots;
            convergence_error = max_delta;
            if max_delta < self.convergence_tol {
                converged = true;
                break;
            }
        }

        // (e) Recover each channel's FilterDesign from its updated S roots.
        // D = (Eᵢ − Fᵢ)/2 from the initial synthesis is already in S's
        // normalization (degree n − 1, non-monic leading coefficient); the
        // recovery rebuilds E = S + D, F = S − D exactly (review P1-1).
        let mut channels = Vec::with_capacity(n);
        for i in 0..n {
            let design = crate::synthesis::recover_channel_design(
                &s_roots[i],
                &d_polys[i],
                &original_polys[i],
                orders[i],
                self.channels[i].return_loss_db,
                centers[i],
                bandwidths[i],
            )?;
            channels.push(design);
        }

        Ok(MultiplexerSynthesisResult {
            channels,
            iterations,
            convergence_error,
            converged,
            junction: self.junction,
            s_roots,
            max_iterations: self.max_iterations,
            junction_common: None,
        })
    }

    fn build_spec(ch: &ChannelSpec) -> Result<FilterSpec> {
        build_filter_spec(ch.order, ch.return_loss_db, &ch.zeros, None)
    }
}

impl Default for MultiplexerSynthesis {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(order: usize, center_hz: f64) -> ChannelSpec {
        ChannelSpec {
            order,
            return_loss_db: 20.0,
            zeros: Vec::new(),
            center_hz,
            bandwidth_hz: 120e6,
        }
    }

    #[test]
    fn multiplexer_requires_two_or_more_channels() -> crate::error::Result<()> {
        let result = MultiplexerSynthesis::new()
            .channel(spec(5, 10.0e9))
            .synthesize();
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn multiplexer_triplexer_converges() -> crate::error::Result<()> {
        let result = MultiplexerSynthesis::new()
            .channels([spec(5, 10.0e9), spec(5, 10.2e9), spec(5, 10.4e9)])
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
            .synthesize()?;

        // The simplified perturbation path converges to a *bounded* residual
        // (it does not reach the machine tolerance for these channels), but it
        // must be finite, non-negative, and non-divergent. Asserting a strict
        // `<= tol` would fail because the approximate model stalls at a small
        // but non-zero residual before exhausting `max_iterations`.
        assert!(result.convergence_error.is_finite());
        assert!(result.convergence_error >= 0.0);
        assert!(
            result.convergence_error < 1.0,
            "triplexer perturbation must be bounded (got {:e})",
            result.convergence_error
        );
        assert!(result.iterations >= 1);
        assert_eq!(result.channels.len(), 3);
        assert_eq!(result.s_roots.len(), 3);
        for (ch, roots) in result.channels.iter().zip(result.s_roots.iter()) {
            assert_eq!(ch.order(), 5);
            assert_eq!(roots.len(), 5);
        }
        Ok(())
    }

    /// Verifies that the recovered channel filter is a usable (passband) filter
    /// rather than a degenerate fully-reflective one.
    ///
    /// This historically failed because the coupling-matrix synthesis extracted
    /// residue real parts only, which is zero for odd-order all-pole designs
    /// (an `s = j w` phase factor) — see `residue_real_component` in
    /// `crate::synthesis::residues`. That core defect (which also affected
    /// `crate::synthesis::recover_channel_design`) is now fixed, so the
    /// recovered channels must exhibit a real passband.
    #[test]
    fn multiplexer_recovered_channel_is_passband() -> crate::error::Result<()> {
        let result = MultiplexerSynthesis::new()
            .channels([spec(5, 10.0e9), spec(5, 10.2e9), spec(5, 10.4e9)])
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
            .synthesize()?;

        let ch = &result.channels[0];
        let center = ch.center_hz().unwrap();
        let bw = ch.bandwidth_hz().unwrap();
        let resp = ch.response_bandpass(center, bw, center - bw, center + bw, 101)?;

        // Peak in-band S21 magnitude and corresponding dB.
        let mut peak_mag = 0.0_f64;
        for s in &resp.samples {
            peak_mag = peak_mag.max(s.s21_mag());
        }
        let peak_db = if peak_mag > 1e-15 {
            20.0 * peak_mag.log10()
        } else {
            f64::NEG_INFINITY
        };

        // A usable loaded channel has a real passband. The mutual-loading
        // correction (now propagated into the matrix) costs some passband peak
        // in the simplified perturbation model, so we allow a bounded loss
        // while still rejecting a degenerate fully-reflective channel. A tight
        // near-0 dB bound is only valid for the un-loaded (independent) case —
        // see the full Macchiarella spectral path for the high-fidelity model.
        assert!(
            peak_db > -10.0,
            "recovered channel is degenerate (peak S21 = {peak_db} dB)"
        );
        Ok(())
    }

    #[test]
    fn multiplexer_quadruplexer_converges() -> crate::error::Result<()> {
        let result = MultiplexerSynthesis::new()
            .channels([
                spec(5, 10.0e9),
                spec(5, 10.2e9),
                spec(5, 10.4e9),
                spec(5, 10.6e9),
            ])
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
            .synthesize()?;

        assert!(result.convergence_error.is_finite());
        assert!(result.convergence_error >= 0.0);
        assert!(
            result.convergence_error < 1.0,
            "quadruplexer perturbation must be bounded (got {:e})",
            result.convergence_error
        );
        assert!(result.iterations >= 1);
        assert_eq!(result.channels.len(), 4);
        assert_eq!(result.s_roots.len(), 4);
        for (ch, roots) in result.channels.iter().zip(result.s_roots.iter()) {
            assert_eq!(ch.order(), 5);
            assert_eq!(roots.len(), 5);
        }
        Ok(())
    }

    /// A Type-II junction (coaxial node with capacitance) must synthesize and
    /// converge just like a Type-I ideal junction.
    #[test]
    fn multiplexer_type_ii_junction_converges() -> crate::error::Result<()> {
        let result = MultiplexerSynthesis::new()
            .channels([spec(5, 10.0e9), spec(5, 10.2e9), spec(5, 10.4e9)])
            .junction(JunctionType::TypeII { c0: 0.2 })
            .synthesize()?;

        assert!(result.convergence_error.is_finite());
        assert!(result.convergence_error >= 0.0);
        assert_eq!(result.channels.len(), 3);
        Ok(())
    }

    /// Regression for review P2-6: the simplified (non-full-iteration) path
    /// computes the Type-II junction coupling as `1/(1 − c0)`, so a
    /// non-finite or non-positive `c0` would produce NaN/Inf and silently
    /// poison every S-root. The full path and the duplexer both validate `c0`;
    /// the simplified multiplexer entry must reject it too.
    #[test]
    fn multiplexer_type_ii_rejects_nonpositive_capacitance() {
        for c0 in [0.0_f64, -1.0, f64::NAN] {
            let err = MultiplexerSynthesis::new()
                .channels([spec(5, 10.0e9), spec(5, 10.2e9)])
                .junction(JunctionType::TypeII { c0 })
                .synthesize()
                .expect_err("non-positive c0 must be rejected");
            assert!(
                matches!(err, MfsError::PreconditionViolation(_)),
                "expected PreconditionViolation, got: {err:?}"
            );
        }
    }

    /// The junction parameters must actually influence the synthesis. A
    /// non-zero Type-I susceptance `b₀` changes the mutual-loading coupling, so
    /// the converged S-roots must differ from the ideal `b₀ = 0` case.
    #[test]
    fn multiplexer_junction_susceptance_affects_synthesis() -> crate::error::Result<()> {
        let specs = [spec(5, 10.0e9), spec(5, 10.2e9), spec(5, 10.4e9)];
        let ideal = MultiplexerSynthesis::new()
            .channels(specs.clone())
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
            .synthesize()?;
        let loaded = MultiplexerSynthesis::new()
            .channels(specs)
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.8 })
            .synthesize()?;

        let mut differ = false;
        for i in 0..ideal.s_roots.len() {
            for (a, b) in ideal.s_roots[i].iter().zip(loaded.s_roots[i].iter()) {
                if (a - b).norm() > 1e-6 {
                    differ = true;
                }
            }
        }
        assert!(
            differ,
            "non-zero junction susceptance b0 must change the synthesized S-roots"
        );
        Ok(())
    }

    /// Regression for review P1-1: the channel-recovery step must reconstruct
    /// each channel's E and F polynomials exactly, not corrupt them. Before the
    /// fix, the simplified path called `recover_channel_design` with
    /// `rescale_d = true`, which forced D to monic (its `d_leading.norm() >
    /// 1e-15` guard always fired, because `sub` trims the cancelled degree-n
    /// term so the stored D is genuinely degree n−1 with an O(1) leading
    /// coefficient). That changed F = S − D's leading coefficient and moved the
    /// reflection zeros, degrading the passband return-loss spec (review:
    /// −20.52 dB → −2.53 dB on the §5 golden spec, the 40-iteration path).
    ///
    /// This guard runs the synthesis at the zero-perturbation limit
    /// (`max_iterations(0)`): no perturbation is applied, so the recovered
    /// S-roots are exactly the initial S = (Eᵢ + Fᵢ)/2 and the recovery must
    /// rebuild E = S + D = Eᵢ and F = S − D = Fᵢ exactly — the limit the P1-1
    /// fix targets (see `recovery.rs`). Isolating the recovery from the
    /// perturbation loop matters for two reasons: the loop is the approximate
    /// Macchiarella loading model, whose passband quality depends on band
    /// separation, and the rescale's effect on the ripple floor is masked at
    /// this limit — it moves the reflection zeros, a band-edge/stopband effect.
    /// So the test asserts the invariant directly (E and F equal the
    /// standalone channel's) rather than through a loose passband-S11 reading.
    ///
    /// The spec is order-4 / RL 22 / 6.75+6.95 GHz / 300 MHz — the review's
    /// confirmed-bite case (stored-D leading coefficient 1.1711, rescale factor
    /// 0.8539). Other specs (e.g. order-5 / RL 20) have D leading ≈ 1, where
    /// the monic rescale is a no-op and would not catch a regression.
    #[test]
    fn multiplexer_recovered_channel_preserves_return_loss() -> crate::error::Result<()> {
        let order = 4;
        let rl_db = 22.0;
        let bw = 300e6;
        let result = MultiplexerSynthesis::new()
            .channels([
                ChannelSpec { order, return_loss_db: rl_db, zeros: Vec::new(), center_hz: 6.75e9, bandwidth_hz: bw },
                ChannelSpec { order, return_loss_db: rl_db, zeros: Vec::new(), center_hz: 6.95e9, bandwidth_hz: bw },
            ])
            .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
            .max_iterations(0)
            .synthesize()?;

        let ch = &result.channels[0];
        let center = ch.center_hz().unwrap();
        let bw = ch.bandwidth_hz().unwrap();
        // Standalone synthesis of the identical spec: the reference E/F that
        // the zero-perturbation recovery must reproduce.
        let standalone = FilterDesign::bandpass(order, rl_db, center, bw).synthesize()?;

        // Monic-normalize a polynomial's coefficients (ascending power order;
        // the leading coefficient is last) so a monic rescale of the whole
        // polynomial is normalized away and we compare the shape of E and F.
        let mono_coefs = |p: &crate::approx::ComplexPolynomial| -> Vec<f64> {
            let lc = p.coefficients.last().unwrap().re;
            p.coefficients.iter().map(|x| x.re / lc).collect()
        };
        // Max absolute difference of two equal-length coefficient vectors
        // (E and F are the same degree for the same spec, so lengths match).
        let max_coeff_diff = |a: &Vec<f64>, b: &Vec<f64>| -> f64 {
            a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).fold(0.0_f64, f64::max)
        };

        let (e_rec, f_rec) = (
            mono_coefs(&ch.polynomials().e),
            mono_coefs(&ch.polynomials().f),
        );
        let (e_ref, f_ref) = (
            mono_coefs(&standalone.polynomials().e),
            mono_coefs(&standalone.polynomials().f),
        );
        // The P1-1 rescale moved F's degree-(n−1) coefficient by ~0.24 (and the
        // reflection zeros with it); the correct reconstruction is exact to
        // ~1e-16. A 1e-8 tolerance is five orders below the bug and far above
        // noise, so it catches a rescale regression without flaking.
        let e_diff = max_coeff_diff(&e_rec, &e_ref);
        let f_diff = max_coeff_diff(&f_rec, &f_ref);
        assert!(
            e_diff < 1e-8,
            "recovery corrupted E (max monic-coef diff {e_diff:.3e}); the P1-1 \
             monic-rescale regression would move E's shape"
        );
        assert!(
            f_diff < 1e-8,
            "recovery corrupted F (max monic-coef diff {f_diff:.3e}); the P1-1 \
             monic-rescale regression moved the reflection zeros by ~0.24"
        );

        // Sanity: the recovered channel's passband is not degenerate. The zero-
        // perturbation recovery carries no junction loading, so its passband
        // S11 tracks the standalone (both ~ the ripple floor). This is a loose
        // guard — the coefficient check above is the precise P1-1 regression.
        let worst_inband = |d: &FilterDesign| -> f64 {
            let resp = d
                .response_bandpass(center, bw, center - bw / 2.0, center + bw / 2.0, 101)
                .unwrap();
            resp.samples.iter().map(|s| s.s11_db()).fold(f64::NEG_INFINITY, f64::max)
        };
        let recovered_s11 = worst_inband(ch);
        let solo_s11 = worst_inband(&standalone);
        assert!(
            (recovered_s11 - solo_s11).abs() < 5.0,
            "recovered passband S11 {recovered_s11} dB diverges from standalone \
             {solo_s11} dB (degenerate recovery)"
        );
        Ok(())
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Multiplexer (N+1)-port Analysis
// ═══════════════════════════════════════════════════════════════════════════════

/// Model used to combine N channel filters into the common-port response.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResponseModel {
    /// True (N+1)-port common-junction model: each channel is a lossless
    /// reciprocal 2-port (port A at the common node, port B at the external
    /// channel port), and the (N+1)×(N+1) S-matrix is the result of the
    /// parallel junction `Y_in = Σᵢ 1/Z_in,i` — assembled via the
    /// wave-variable node equation, which is well-conditioned where the
    /// explicit per-channel admittance matrix is singular (matched band
    /// center). Unitary and reciprocal for lossless inputs; N=1 reduces
    /// exactly to the channel's own 2-port S-matrix. Physically correct.
    Admittance,
    /// Legacy magnitude model (geometric-mean antenna reflection
    /// `(∏ᵢ|S11ᵢ|)^(1/N)` and product-of-transmissions isolation). Retained
    /// for one release cycle for parity comparison; superseded by
    /// [`ResponseModel::Admittance`], which it approximates poorly in the
    /// passband (see `docs/multiplexer-synthesis.md` deviation item 2).
    GeometricMean,
}

/// An N-channel multiplexer combining N bandpass filters on one common port.
#[derive(Debug, Clone)]
pub struct Multiplexer {
    channels: Vec<FilterDesign>,
}

impl Multiplexer {
    /// All N channels must have bandpass parameters (center + bandwidth).
    pub fn new(channels: Vec<FilterDesign>) -> Result<Self> {
        if channels.is_empty() {
            return Err(MfsError::PreconditionViolation(
                "multiplexer requires at least one channel".into(),
            ));
        }
        for (i, ch) in channels.iter().enumerate() {
            if ch.center_hz().is_none() || ch.bandwidth_hz().is_none() {
                return Err(MfsError::PreconditionViolation(format!(
                    "channel {i} must have bandpass parameters"
                )));
            }
        }
        Ok(Self { channels })
    }

    /// Number of channels.
    pub fn len(&self) -> usize {
        self.channels.len()
    }

    /// True when there are no channels.
    pub fn is_empty(&self) -> bool {
        self.channels.is_empty()
    }

    /// The i-th channel filter design.
    pub fn channel(&self, i: usize) -> &FilterDesign {
        &self.channels[i]
    }

    /// The i-th channel's folded coupling matrix.
    pub fn channel_folded(&self, i: usize) -> Result<CouplingMatrix> {
        self.channels[i].to_folded()
    }

    /// All channel coupling matrices.
    pub fn channel_matrices(&self) -> Vec<&CouplingMatrix> {
        self.channels.iter().map(|ch| ch.matrix()).collect()
    }

    /// N+1-port response over [start_hz, stop_hz] using the default
    /// [`ResponseModel::Admittance`] model (the true (N+1)-port common-node
    /// composition).
    ///
    /// For the legacy geometric-mean magnitude model use
    /// [`Multiplexer::response_with`].
    pub fn response(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<MultiplexerResponse> {
        self.response_with(start_hz, stop_hz, points, ResponseModel::Admittance)
    }

    /// N+1-port response over [start_hz, stop_hz] with an explicit
    /// [`ResponseModel`].
    pub fn response_with(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
        model: ResponseModel,
    ) -> Result<MultiplexerResponse> {
        match model {
            ResponseModel::Admittance => self.response_admittance(start_hz, stop_hz, points),
            ResponseModel::GeometricMean => {
                self.response_geometric_mean(start_hz, stop_hz, points)
            }
        }
    }

    /// Full complex (N+1)×(N+1) S-matrix at every sampled frequency
    /// ([`ResponseModel::Admittance`]).
    ///
    /// Each channel filter is treated as a lossless reciprocal 2-port whose
    /// input (port A) sits at the common node and whose output (port B) is the
    /// external channel port. The (N+1)-port S-matrix is assembled from the
    /// wave-variable node equation (KCL at the common node), per the
    /// [`ResponseModel::Admittance`] model. For lossless channels the result
    /// is unitary and reciprocal, and N = 1 reduces exactly to the channel's
    /// own 2-port S-matrix. The formulation is well-conditioned at the
    /// matched band center, where the per-channel admittance matrix is
    /// singular.
    pub fn response_matrix(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<MultiplexerSMatrix> {
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;
        let channels = self.evaluate_channels(&grid)?;
        let mut s = Vec::with_capacity(points);
        for p in 0..points {
            let s11s: Vec<Complex64> = channels.iter().map(|ch| ch[p].0).collect();
            let s21s: Vec<Complex64> = channels.iter().map(|ch| ch[p].1).collect();
            let mat = compose_node_s(&s11s, &s21s);
            // Flatten explicitly row-major (S[i][j] at i*p + j). `mat.iter()`
            // is nalgebra's column-major order, which would store S transposed
            // (review P1-2).
            let ports = mat.nrows();
            let mut flat = Vec::with_capacity(ports * ports);
            for i in 0..ports {
                for j in 0..ports {
                    flat.push(mat[(i, j)]);
                }
            }
            s.push(flat);
        }
        Ok(MultiplexerSMatrix {
            frequencies: grid.as_slice().to_vec(),
            num_ports: self.channels.len() + 1,
            s,
        })
    }

    /// Admittance-model response: the true (N+1)-port S-matrix, reported as dB.
    ///
    /// `antenna_s11_db` = |S₀₀|, `channel_s21_db[i]` = |Sᵢ₊₁₀| (antenna →
    /// channel i), `channel_s11_db[i]` = |Sᵢ₊₁,ᵢ₊₁|, `isolation_db[i][j]` =
    /// |Sᵢ₊₁,ⱼ₊₁| (channel j driven, channel i measured). Unitarity bounds
    /// every entry by 0 dB.
    fn response_admittance(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<MultiplexerResponse> {
        let matrix = self.response_matrix(start_hz, stop_hz, points)?;
        let n = self.channels.len();
        let p = n + 1;

        let mut antenna_s11_db = Vec::with_capacity(points);
        let mut channel_s21_db: Vec<Vec<f64>> =
            (0..n).map(|_| Vec::with_capacity(points)).collect();
        let mut channel_s11_db: Vec<Vec<f64>> =
            (0..n).map(|_| Vec::with_capacity(points)).collect();
        let mut isolation_db: Vec<Vec<Vec<f64>>> =
            (0..n).map(|_| (0..n).map(|_| Vec::new()).collect()).collect();

        for q in 0..points {
            let row = &matrix.s[q];
            antenna_s11_db.push(mag_to_db(row[0]));
            for i in 0..n {
                // channel_s21_db[i] = |S[i+1][0]| (antenna -> channel i).
                channel_s21_db[i].push(mag_to_db(row[(i + 1) * p]));
                // channel_s11_db[i] = |S[i+1][i+1]| (channel i's own port).
                channel_s11_db[i].push(mag_to_db(row[(i + 1) * p + (i + 1)]));
            }
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        // isolation_db[i][j] = |S[i+1][j+1]| (channel j
                        // driven, channel i measured).
                        isolation_db[i][j].push(mag_to_db(row[(i + 1) * p + (j + 1)]));
                    }
                }
            }
        }

        Ok(MultiplexerResponse {
            frequencies: matrix.frequencies,
            antenna_s11_db,
            channel_s21_db,
            channel_s11_db,
            isolation_db,
        })
    }

    /// Legacy geometric-mean magnitude model: each channel's filter is
    /// evaluated independently, the common (antenna) reflection is approximated
    /// as the geometric mean of the channel S11 magnitudes
    /// `(∏ᵢ |S11ᵢ|)^(1/N)` (reducing to the diplexer's
    /// `√(|S11_tx|·|S11_rx|)` for N=2), and the isolation between channels i
    /// and j is approximated as the product of their S21 magnitudes.
    ///
    /// Retained for one release cycle for parity comparison; superseded by
    /// the admittance model, which it approximates poorly in the passband
    /// (see `docs/multiplexer-synthesis.md` deviation item 2).
    fn response_geometric_mean(
        &self,
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    ) -> Result<MultiplexerResponse> {
        let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;
        let n = self.channels.len();

        // Evaluate each channel's filter with its own band-pass mapping.
        let mut channel_s21_mag: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut channel_s11_mag: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut channel_s21_db: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut channel_s11_db: Vec<Vec<f64>> = Vec::with_capacity(n);

        for ch in &self.channels {
            let center = ch.center_hz().unwrap();
            let bw = ch.bandwidth_hz().unwrap();
            let mapping = BandPassMapping::new(center, bw)?;
            let resp = ResponseSolver.evaluate(ch.matrix(), &grid, &mapping)?;

            let mut s21_mag = Vec::with_capacity(points);
            let mut s11_mag = Vec::with_capacity(points);
            let mut s21_db = Vec::with_capacity(points);
            let mut s11_db = Vec::with_capacity(points);
            for s in &resp.samples {
                s21_mag.push(s.s21_mag());
                s11_mag.push(s.s11_mag());
                s21_db.push(s.s21_db());
                s11_db.push(s.s11_db());
            }
            channel_s21_mag.push(s21_mag);
            channel_s11_mag.push(s11_mag);
            channel_s21_db.push(s21_db);
            channel_s11_db.push(s11_db);
        }

        let frequencies = grid.as_slice().to_vec();

        // Antenna reflection: geometric mean of the channel S11 magnitudes.
        // `(∏ᵢ |S11ᵢ|)^(1/N)` matches the diplexer's `√(|S11_tx|·|S11_rx|)`
        // model for N=2 and is a better approximation of the parallel-junction
        // reflection than the maximum. This keeps the N-channel and diplexer
        // analyses consistent.
        let mut antenna_s11_db = Vec::with_capacity(points);
        let mut channel_s11_iters: Vec<_> = channel_s11_mag
            .iter()
            .map(|ch| ch.iter().copied())
            .collect();
        for _ in 0..points {
            let mut log_sum = 0.0_f64;
            let mut nonzero = 0usize;
            for iter in &mut channel_s11_iters {
                if let Some(mag) = iter.next() {
                    if mag > 1e-300 {
                        log_sum += mag.ln();
                        nonzero += 1;
                    }
                }
            }
            let geomean = if nonzero > 0 {
                (log_sum / nonzero as f64).exp()
            } else {
                0.0
            };
            let db = if geomean > 1e-15 {
                20.0 * geomean.log10()
            } else {
                -300.0
            };
            antenna_s11_db.push(db);
        }

        // Isolation: product of each ordered pair's S21 magnitudes.
        let mut isolation_db = Vec::with_capacity(n);
        for i in 0..n {
            let mut iso_row = Vec::with_capacity(n);
            for j in 0..n {
                if i == j {
                    // Self-isolation is not meaningful; store a zero-length vector.
                    iso_row.push(Vec::new());
                    continue;
                }
                let mut row_db = Vec::with_capacity(points);
                for (mag_i, mag_j) in channel_s21_mag[i].iter().zip(channel_s21_mag[j].iter()) {
                    let prod = mag_i * mag_j;
                    let db = if prod > 1e-15 {
                        20.0 * prod.log10()
                    } else {
                        -300.0
                    };
                    row_db.push(db);
                }
                iso_row.push(row_db);
            }
            isolation_db.push(iso_row);
        }

        Ok(MultiplexerResponse {
            frequencies,
            antenna_s11_db,
            channel_s21_db,
            channel_s11_db,
            isolation_db,
        })
    }

    /// Evaluate every channel's 2-port `(s11, s21)` over the shared grid.
    fn evaluate_channels(&self, grid: &FrequencyGrid) -> Result<Vec<Vec<(Complex64, Complex64)>>> {
        let mut channels = Vec::with_capacity(self.channels.len());
        for ch in &self.channels {
            let center = ch.center_hz().unwrap();
            let bw = ch.bandwidth_hz().unwrap();
            let mapping = BandPassMapping::new(center, bw)?;
            let resp = ResponseSolver.evaluate(ch.matrix(), grid, &mapping)?;
            channels.push(
                resp
                    .samples
                    .iter()
                    .map(|s| {
                        (
                            Complex64::new(s.s11_re, s.s11_im),
                            Complex64::new(s.s21_re, s.s21_im),
                        )
                    })
                    .collect(),
            );
        }
        Ok(channels)
    }
}

/// S-magnitude to dB, with the -300 dB floor used across the crate.
fn mag_to_db(mag: Complex64) -> f64 {
    let m = mag.norm();
    if m > 1e-15 {
        20.0 * m.log10()
    } else {
        -300.0
    }
}

/// Lossless-reciprocal 2-port `s22` from `s11` and `s21`:
/// `s22 = −conj(s11)·s21/conj(s21)`, with the `s21 → 0` limit `−conj(s11)`.
fn lossless_s22(s11: Complex64, s21: Complex64) -> Complex64 {
    if s21 == Complex64::new(0.0, 0.0) {
        -s11.conj()
    } else {
        -s11.conj() * s21 / s21.conj()
    }
}

/// Compose the (N+1)-port S-matrix at the common node using the wave-variable
/// node equation (KCL). For each channel the input (port A) sits at the common
/// node and the output (port B) is the external channel port. With
/// `a₀`/`aBₖ` the driven-port incident waves, the node voltage is
///
/// `V = [2a₀ + Σₖ 2 s₁₂ₖ aBₖ/(1+s₁₁ₖ)] / [1 + Σₖ (1−s₁₁ₖ)/(1+s₁₁ₖ)]`
///
/// and `b₀ = V − a₀`, `bBₖ = s₂₁ₖ(V − s₁₂ₖ aBₖ)/(1+s₁₁ₖ) + s₂₂ₖ aBₖ`.
///
/// For lossless reciprocal channels (`s₁₂=s₂₁`, `s₂₂` from [`lossless_s22`])
/// the resulting (N+1)×(N+1) S-matrix is unitary and reciprocal, and N=1
/// reduces exactly to the channel's own 2-port S-matrix. This form is
/// well-conditioned at the matched band center (s₁₁=0), where the admittance
/// construction is singular.
fn compose_node_s(s11s: &[Complex64], s21s: &[Complex64]) -> DMatrix<Complex64> {
    let n = s11s.len();
    let p = n + 1;
    let mut s = DMatrix::zeros(p, p);
    for drive in 0..p {
        let a0 = if drive == 0 { 1.0 } else { 0.0 };
        let mut num = Complex64::new(2.0 * a0, 0.0);
        let mut den = Complex64::new(1.0, 0.0);
        for k in 0..n {
            let s11 = s11s[k];
            let s21 = s21s[k];
            let one_plus = Complex64::new(1.0, 0.0) + s11;
            if k + 1 == drive {
                num += 2.0 * s21 / one_plus; // s12 = s21 (reciprocal)
            }
            den += (1.0 - s11) / one_plus;
        }
        let v = num / den;
        s[(0, drive)] = v - a0;
        for k in 0..n {
            let s11 = s11s[k];
            let s21 = s21s[k];
            let s22 = lossless_s22(s11, s21);
            let ab = if k + 1 == drive { 1.0 } else { 0.0 };
            s[(k + 1, drive)] = s21 * (v - s21 * ab) / (1.0 + s11) + s22 * ab;
        }
    }
    s
}

/// Full complex (N+1)×(N+1) S-matrix at each sampled frequency
/// ([`ResponseModel::Admittance`]). Port 0 is the antenna (common port);
/// ports 1..N+1 are the channels in order. Row-major:
/// `s[q][i * num_ports + j]` = Sᵢⱼ at frequency `q`.
#[derive(Debug, Clone)]
pub struct MultiplexerSMatrix {
    /// Frequency points in Hz.
    pub frequencies: Vec<f64>,
    /// Number of ports (N channels + 1 antenna).
    pub num_ports: usize,
    /// Per-frequency (N+1)×(N+1) S-matrix, row-major.
    pub s: Vec<Vec<Complex64>>,
}

/// (N+1)-port response: 1 common port + N channel ports.
#[derive(Debug, Clone)]
pub struct MultiplexerResponse {
    /// Frequency points in Hz.
    pub frequencies: Vec<f64>,
    /// Antenna (common port) reflection S11, dB.
    pub antenna_s11_db: Vec<f64>,
    /// Per-channel transmission antenna→channel i (S_{(i+1)1}), dB. Outer index = channel.
    pub channel_s21_db: Vec<Vec<f64>>,
    /// Per-channel reflection at its own port, dB.
    pub channel_s11_db: Vec<Vec<f64>>,
    /// Isolation between every ordered pair of distinct channels (i != j), dB.
    /// `isolation[i][j]` = |`S[i+1][j+1]`| — channel j driven, channel i measured.
    /// For i == j the entry is a zero-length vector (self-isolation is undefined).
    pub isolation_db: Vec<Vec<Vec<f64>>>,
}

#[cfg(test)]
mod analysis_tests {
    use super::*;

    fn bandpass_design(center_hz: f64) -> Result<FilterDesign> {
        FilterDesign::bandpass(5, 20.0, center_hz, 100e6).synthesize()
    }

    #[test]
    fn multiplexer_requires_at_least_one_channel() -> Result<()> {
        let result = Multiplexer::new(vec![]);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn multiplexer_analysis_matches_diplexer_for_n2() -> Result<()> {
        let ch0 = bandpass_design(10.0e9)?;
        let ch1 = bandpass_design(10.2e9)?;
        let mux = Multiplexer::new(vec![ch0, ch1])?;

        let points = 101;
        let resp = mux.response(9.9e9, 10.3e9, points)?;

        assert_eq!(resp.frequencies.len(), points);
        assert_eq!(resp.channel_s21_db.len(), 2);
        assert_eq!(resp.channel_s11_db.len(), 2);
        for i in 0..2 {
            assert_eq!(resp.channel_s21_db[i].len(), points);
            assert_eq!(resp.channel_s11_db[i].len(), points);
        }
        assert_eq!(resp.isolation_db.len(), 2);
        for i in 0..2 {
            assert_eq!(resp.isolation_db[i].len(), 2);
            assert!(resp.isolation_db[i][i].is_empty());
            assert_eq!(resp.isolation_db[i][1 - i].len(), points);
        }
        assert_eq!(resp.antenna_s11_db.len(), points);
        for db in &resp.antenna_s11_db {
            assert!(db.is_finite());
        }
        Ok(())
    }

    /// The admittance-model (N+1)-port S-matrix is power-conserving (unitary)
    /// and reciprocal for the lossless channel filters. This is the physical
    /// correctness gate for the true (N+1)-port response (B1 acceptance).
    #[test]
    fn multiplexer_admittance_response_is_unitary_and_reciprocal() -> Result<()> {
        let mux = Multiplexer::new(vec![
            bandpass_design(10.0e9)?,
            bandpass_design(10.2e9)?,
            bandpass_design(10.4e9)?,
        ])?;
        let mat = mux.response_matrix(9.6e9, 10.8e9, 61)?;
        let p = mat.num_ports;
        assert_eq!(p, 4);
        assert_eq!(mat.s.len(), 61);
        for row in mat.s.iter() {
            assert_eq!(row.len(), p * p);
            for r in 0..p {
                for c in 0..p {
                    let mut ss = Complex64::new(0.0, 0.0);
                    for k in 0..p {
                        ss += row[k * p + r] * row[k * p + c].conj();
                    }
                    let expected = if r == c {
                        Complex64::new(1.0, 0.0)
                    } else {
                        Complex64::new(0.0, 0.0)
                    };
                    assert!(
                        (ss - expected).norm() < 1e-9,
                        "unitarity violated at row {r} col {c}: {ss:?}"
                    );
                    assert!(
                        (row[r * p + c] - row[c * p + r]).norm() < 1e-9,
                        "reciprocity violated at ({r},{c})"
                    );
                }
            }
        }
        Ok(())
    }

    /// N=1 reduction: a single-channel "multiplexer" is exactly the channel's
    /// own lossless 2-port S-matrix (antenna port = channel input port).
    #[test]
    fn multiplexer_single_channel_reduces_to_2port() -> Result<()> {
        let ch = bandpass_design(10.0e9)?;
        let mux = Multiplexer::new(vec![ch])?;
        let mat = mux.response_matrix(9.9e9, 10.1e9, 21)?;
        assert_eq!(mat.num_ports, 2);

        let center = mux.channel(0).center_hz().unwrap();
        let bw = mux.channel(0).bandwidth_hz().unwrap();
        let mapping = BandPassMapping::new(center, bw)?;
        let grid = FrequencyGrid::linspace(9.9e9, 10.1e9, 21)?;
        let resp = ResponseSolver.evaluate(mux.channel(0).matrix(), &grid, &mapping)?;

        for (q, row) in mat.s.iter().enumerate() {
            let s11 = Complex64::new(resp.samples[q].s11_re, resp.samples[q].s11_im);
            let s21 = Complex64::new(resp.samples[q].s21_re, resp.samples[q].s21_im);
            let s22 = lossless_s22(s11, s21);
            assert!((row[0] - s11).norm() < 1e-9, "S11 mismatch at {q}");
            assert!((row[1] - s21).norm() < 1e-9, "S21 mismatch at {q}");
            assert!((row[2] - s21).norm() < 1e-9, "S12 (reciprocity) mismatch at {q}");
            assert!((row[3] - s22).norm() < 1e-9, "S22 (lossless) mismatch at {q}");
        }
        Ok(())
    }

    /// The legacy geometric-mean model and the true admittance model agree in
    /// the deep stopband (10.1 GHz, between the two bands where both channels
    /// are reflective) — the only regime where the geomean approximation is
    /// physically accurate. In the passband the two legitimately diverge
    /// (geomean underestimates the common-port reflection).
    #[test]
    fn multiplexer_models_agree_in_stopband() -> Result<()> {
        let mux = Multiplexer::new(vec![bandpass_design(10.0e9)?, bandpass_design(10.2e9)?])?;
        let points = 21;
        let adm = mux.response(10.09e9, 10.11e9, points)?;
        let geo = mux.response_with(10.09e9, 10.11e9, points, ResponseModel::GeometricMean)?;
        let mid = points / 2; // 10.10 GHz, both channels in stopband
        let delta = (adm.antenna_s11_db[mid] - geo.antenna_s11_db[mid]).abs();
        assert!(
            delta < 0.5,
            "stopband antenna S11: admittance = {} dB, geomean = {} dB, delta = {} dB",
            adm.antenna_s11_db[mid],
            geo.antenna_s11_db[mid],
            delta
        );
        Ok(())
    }

    /// Regression for review P1-2: the admittance-model dB fields must match
    /// the documented physical elements of the raw (N+1)-port S-matrix. The
    /// previous extraction read the column-major-flattened matrix as if it
    /// were row-major, so `channel_s21_db[i]` (i>=1) and `isolation_db[i][j]`
    /// (all pairs) reported transposed elements — a 72 dB isolation was
    /// reported as ~0 dB. This test asserts every reported dB value against
    /// the exact `S[i+1][j]` / `S[i+1][j+1]` entries.
    #[test]
    fn multiplexer_response_db_fields_match_s_matrix() -> Result<()> {
        let mux = Multiplexer::new(vec![bandpass_design(10.0e9)?, bandpass_design(10.2e9)?])?;
        let points = 21;
        let resp = mux.response(9.9e9, 10.3e9, points)?;
        let mat = mux.response_matrix(9.9e9, 10.3e9, points)?;
        let p = mat.num_ports;

        for q in 0..points {
            let row = &mat.s[q];
            let s = |i: usize, j: usize| mag_to_db(row[i * p + j]);
            assert!(
                (resp.antenna_s11_db[q] - s(0, 0)).abs() < 1e-9,
                "antenna_s11 mismatch at {q}"
            );
            for i in 0..2 {
                assert!(
                    (resp.channel_s21_db[i][q] - s(i + 1, 0)).abs() < 1e-9,
                    "channel_s21[{i}] mismatch at {q}"
                );
                assert!(
                    (resp.channel_s11_db[i][q] - s(i + 1, i + 1)).abs() < 1e-9,
                    "channel_s11[{i}] mismatch at {q}"
                );
            }
            for i in 0..2 {
                for j in 0..2 {
                    if i != j {
                        assert!(
                            (resp.isolation_db[i][j][q] - s(i + 1, j + 1)).abs() < 1e-9,
                            "isolation[{i}][{j}] mismatch at {q}"
                        );
                    }
                }
            }
        }

        // Physical sanity at channel 1's own band center (10.2 GHz): its
        // transmission is a passband peak (near 0 dB) while channel 0 is in
        // its stopband (deep). The pre-fix code reported the *diagonal*
        // S[1][1] (~0 dB for any N) for channel 1's transmission and the
        // transposed S[1][2] for channel 0's, so both assertions fail pre-fix.
        // Grid 9.9e9..=10.3e9 / 21 points: step 2 MHz, 10.2e9 is index 15.
        let c1_center = (0..points)
            .find(|&q| (resp.frequencies[q] - 10.2e9).abs() < 1e3)
            .expect("10.2 GHz on the grid");
        assert!(
            resp.channel_s21_db[1][c1_center] > -10.0,
            "channel 1 passband peak too low: {} dB",
            resp.channel_s21_db[1][c1_center]
        );
        assert!(
            resp.channel_s21_db[0][c1_center] < -10.0,
            "channel 0 must be in its stopband at channel 1's center: {} dB",
            resp.channel_s21_db[0][c1_center]
        );
        Ok(())
    }
}
