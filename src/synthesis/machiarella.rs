//! Shared Macchiarella–Tamiazzo full-spectral diplexer/multiplexer iteration.
//!
//! Generalizes the converged 2-channel full-spectral loop (originally in
//! `crate::duplexer`) to an arbitrary number of channels sharing one common
//! (antenna) junction. All polynomials live in ONE shared bandpass variable
//! `s = j·Ω(f)` with
//!
//! ```text
//! Ω(f) = (f/f₀ − f₀/f)/b,   f₀ = √(f_L·f_H),  b = (f_H − f_L)/f₀,
//! ```
//!
//! where `f_L`/`f_H` are the outer band edges of the combined channel set, so
//! the two extreme edges sit at `s = ±j`. Each channel's prototype polynomial
//! `P_p(ω)` is composed into the common domain by `P_c(Ω) = P_p(AᵢΩ + Bᵢ)`
//! with `Aᵢ = (Ωᵢ,hi − Ωᵢ,lo)/2`, `Bᵢ = (Ωᵢ,hi + Ωᵢ,lo)/2`, mapping the
//! prototype band `[−1,+1]` onto the channel's common sub-band.
//!
//! The fixed-point iteration (paper §IV, eq. 13–17):
//!
//! 1. `N(Ω)` — fixed numerator built once from the pure-imaginary channel
//!    reflection zeros (roots of each channel's `F`), plus the resonating
//!    node's extra real zero `1/c0` for a Type-II junction.
//! 2. Per iteration, with the running channel characteristic polynomials
//!    `Sᵢ(Ω)` (root state, prototype domain): the transmission polynomial
//!    of channel `i` is `Ptᵢ = Pᵢ · ∏_{j≠i} Sⱼ` where `Pᵢ` is the
//!    finite-zero polynomial (degree 0 for zero-free channels), and the
//!    junction power scalings `p0ᵢ²` are the solution of the N×N linear
//!    system obtained by imposing the return-loss condition
//!    `Σⱼ |Ptⱼ(eᵢ)|²·p0ⱼ² = |N(eᵢ)|²·(1−εᵢ)/εᵢ` at each channel's *outer*
//!    band edge `eᵢ` (the edge farthest from the other channels' centers;
//!    for two channels these are exactly `s = −j` and `s = +j`).
//! 3. `R = N·Ñ + Σᵢ p0ᵢ²·Ptᵢ·Pt̃ᵢ`; `D = spectral_factor(R, deg N)`.
//! 4. `∏ᵢ Sᵢ = (a·N + b·D)/2` (Type-I, `a = 1−j·n²·b0`, `b = 1+j·n²·b0`)
//!    or `(c0/2)·(D_monic − N)` (Type-II); the product roots are assigned
//!    to the channels in contiguous blocks of sizes `nᵢ`, sorted by
//!    increasing imaginary part (paper Step 4) and taken in ascending
//!    channel-center order.
//! 5. Accept when the max root displacement over all channels drops below
//!    the tolerance, or (high-degree stagnation) when the displacement
//!    stops improving and is already small.
//!
//! The N = 2 branch reproduces the historical 2-channel loop **with
//! bit-identical operation order** (root state held as prototype-domain
//! vectors, Cramer's-rule 2×2 p₀ solve, `±j` edge selection by band
//! position, block sizes in lower/upper order) so the existing
//! gate-anchored diplexer results — including the paper §V printed-design
//! golden — are unchanged.

use crate::approx::{
    compose_affine, compute_root_displacement, generalized_chebyshev_polynomials,
    pad_polynomial, ComplexPolynomial, PolynomialSet,
};
use crate::duplexer::{ChannelSpec, JunctionType};
use crate::error::{MfsError, Result};
use crate::freq::FrequencyGrid;
use crate::response::magnitude_to_db;
use crate::spectral::spectral_factor;
use crate::synthesis::build_filter_spec;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// Per-channel state of the fixed-point iteration.
#[derive(Debug, Clone)]
pub struct ChannelState {
    /// Filter order `nᵢ`.
    pub order: usize,
    /// Return loss in dB (εᵢ = 10^(−rlᵢ/10)).
    pub return_loss_db: f64,
    /// Affine scale Aᵢ (prototype→common).
    pub a: f64,
    /// Affine offset Bᵢ (prototype→common, imaginary).
    pub b: f64,
    /// Finite-zero polynomial Pᵢ (prototype domain, degree 0 for zero-free).
    pub p_proto: ComplexPolynomial,
    /// Characteristic polynomial Sᵢ (prototype domain, latest iteration).
    pub s_proto: ComplexPolynomial,
    /// Outer-band return-loss edge in the common domain (±j).
    pub edge: Complex64,
    /// Center frequency in Hz (channel ordering reference).
    pub center_hz: f64,
}

impl ChannelState {
    /// Common-domain finite-zero polynomial.
    fn p_common(&self) -> Result<ComplexPolynomial> {
        compose_affine(&self.p_proto, self.a, Complex64::new(0.0, self.b))
    }

    /// Maps common-domain roots back to the prototype domain, dropping
    /// non-finite and out-of-range roots (the historical stability filter).
    fn to_proto(&self, roots: &[Complex64]) -> Vec<Complex64> {
        roots
            .iter()
            .map(|&r| (r - Complex64::new(0.0, self.b)) / self.a)
            .filter(|r| r.re.is_finite() && r.im.is_finite() && r.norm() < 50.0)
            .collect()
    }
}

/// Common-domain junction polynomials captured from the converged iteration.
#[derive(Debug, Clone)]
pub struct JunctionCommon {
    /// Fixed antenna-reflection numerator N(Ω).
    pub n: ComplexPolynomial,
    /// Junction spectral factor D(Ω).
    pub d: ComplexPolynomial,
    /// Per-channel transmission polynomials Ptᵢ(Ω) = Pᵢ·∏_{j≠i} Sⱼ, in
    /// channel order.
    pub pt: Vec<ComplexPolynomial>,
    /// Per-channel junction power scalings |p0ᵢ|², in channel order.
    pub p0_sq: Vec<f64>,
    /// Combined-band geometric-mean center f₀ = √(f_L·f_H).
    pub f0: f64,
    /// Combined-band fractional width b = (f_H − f_L)/f₀.
    pub b: f64,
}

/// Convergence outcome of the full fixed-point iteration.
#[derive(Debug)]
pub struct FullIterationState {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Max root displacement at capture (strict tolerance or stagnation
    /// accept); `INFINITY` when no capture happened.
    pub convergence_error: f64,
    /// `true` when the loop stopped on its own convergence criterion (strict
    /// tolerance or the high-degree stagnation accept) rather than exhausting
    /// the iteration budget.
    pub converged: bool,
    /// Per-channel converged S roots (prototype domain), in channel order.
    pub s_roots: Vec<Vec<Complex64>>,
    /// Junction polynomials (present when the iteration was captured).
    pub common: Option<JunctionCommon>,
}

const STAGNATION_ACCEPT: f64 = 1e-3;
const STAGNATION_WINDOW: usize = 3;
const STAGNATION_MIN_ITER: usize = 5;

/// Floor for the p₀² power scaling: absorbs the ~0 round-off case (a
/// physically-absent channel) down to a tiny positive value so `sqrt` is
/// well-defined. Legitimate p₀² are far above this (measured floor
/// ~4e-12 on the gated well-separated case).
const P0_SQ_FLOOR: f64 = 1e-15;

/// Rejection threshold for a *genuinely* negative p₀² (review P1-4). A
/// negative beyond this means the p₀² system has no all-positive solution —
/// the channels are too strongly coupled/overlapping for a realizable
/// quasi-equiripple junction. Calibrated 2026-10-06: the smallest *legitimate*
/// p₀² measured on the gated cases is ~4.3e-12 (well-separated N=2), so
/// -1e-12 rejects every observed non-physical negative (O(1e-5) or larger on
/// the overlapping specs) while never touching a valid one. Values in
/// (-1e-12, 0] are treated as round-off and floored, not rejected.
const P0_SQ_REJECTION_TOL: f64 = 1e-12;

/// Common-domain setup for the fixed-point iteration: per-channel affine
/// maps and edges, the fixed numerator N, the combined-band constants, and
/// the independently synthesized per-channel polynomial sets (needed by
/// callers to recover each channel's design).
#[derive(Debug)]
pub struct Prepared {
    /// Per-channel iteration state (prototype-domain S polynomials).
    pub states: Vec<ChannelState>,
    /// Fixed antenna-reflection numerator N(Ω).
    pub n_poly: ComplexPolynomial,
    /// Combined-band geometric-mean center f₀ = √(f_L·f_H).
    pub f0: f64,
    /// Combined-band fractional width b = (f_H − f_L)/f₀.
    pub b: f64,
    /// Per-channel independently synthesized polynomial sets (E, F, P, …).
    pub polys: Vec<PolynomialSet>,
}

/// Sets up the common-domain state for `channels` and builds the fixed
/// numerator N together with the combined-band constants (f₀, b).
///
/// Junction support: the full spectral path implements the *ideal* Type-I
/// junction (`n = 1`, `b0 = 0`) and Type-II (`c0 > 0`) for any N ≥ 2. The
/// 2-channel path additionally supports a non-ideal Type-I verbatim from the
/// historical implementation; that generalization has no paper basis for
/// N > 2 and is rejected with an explicit error rather than invented math.
pub fn prepare(
    channels: &[ChannelSpec],
    junction: JunctionType,
) -> Result<Prepared> {
    let n = channels.len();
    if n < 2 {
        return Err(MfsError::PreconditionViolation(
            "full iteration requires at least two channels".into(),
        ));
    }
    match junction {
        JunctionType::TypeI { n: ratio, b0 }
            if n > 2 && (ratio != 1.0 || b0 != 0.0) =>
        {
            return Err(MfsError::PreconditionViolation(
                "the full spectral iteration supports only the ideal Type-I \
                 junction (n = 1, b0 = 0) and Type-II for more than two \
                 channels".into(),
            ));
        }
        _ => {}
    }
    if let JunctionType::TypeII { c0 } = junction {
        if !c0.is_finite() || c0 <= 0.0 {
            return Err(MfsError::PreconditionViolation(
                "Type-II junction capacitance c0 must be finite and positive".into(),
            ));
        }
    }

    // Independent per-channel synthesis (the same path the historical
    // 2-channel implementation uses).
    let polys: Vec<PolynomialSet> = channels
        .iter()
        .map(|ch| {
            let spec = build_filter_spec(ch.order, ch.return_loss_db, &ch.zeros, None)?;
            generalized_chebyshev_polynomials(&spec)
        })
        .collect::<Result<_>>()?;

    // Combined band limits: min/max over all channel edges.
    let f_lo = channels
        .iter()
        .map(|ch| ch.center_hz - 0.5 * ch.bandwidth_hz)
        .fold(f64::INFINITY, f64::min);
    let f_hi = channels
        .iter()
        .map(|ch| ch.center_hz + 0.5 * ch.bandwidth_hz)
        .fold(f64::NEG_INFINITY, f64::max);
    if !(f_lo.is_finite() && f_hi.is_finite()) || f_hi <= f_lo {
        return Err(MfsError::PreconditionViolation(
            "channels must span a non-empty combined band".into(),
        ));
    }
    let f0 = (f_lo * f_hi).sqrt();
    let b_common = (f_hi - f_lo) / f0;
    let omega = |f: f64| (f / f0 - f0 / f) / b_common;

    // Per-channel affine maps (prototype [−1,+1] → common sub-band).
    let mut centers: Vec<f64> = Vec::with_capacity(n);
    let mut affine: Vec<(f64, f64)> = Vec::with_capacity(n);
    for ch in channels {
        centers.push(ch.center_hz);
        let w_lo = omega(ch.center_hz - 0.5 * ch.bandwidth_hz);
        let w_hi = omega(ch.center_hz + 0.5 * ch.bandwidth_hz);
        affine.push(((w_hi - w_lo) / 2.0, (w_hi + w_lo) / 2.0));
    }

    // Outer-edge selection: for channel i, the band edge farthest (in Hz)
    // from the other channels' centers. For two channels this is exactly
    // s = −j (lower channel) / s = +j (upper channel), matching the
    // historical ±j assignment (its `tx_is_lower` tie-break falls to the
    // lower edge, the same rule as here).
    let mut edges: Vec<Complex64> = Vec::with_capacity(n);
    for i in 0..n {
        let lo_f = centers[i] - 0.5 * channels[i].bandwidth_hz;
        let hi_f = centers[i] + 0.5 * channels[i].bandwidth_hz;
        let mut max_lo = 0.0_f64;
        let mut max_hi = 0.0_f64;
        for (j, c_j) in centers.iter().enumerate() {
            if j == i {
                continue;
            }
            max_lo = max_lo.max((c_j - lo_f).abs());
            max_hi = max_hi.max((c_j - hi_f).abs());
        }
        // On an exact tie the lower edge wins (the historical rule).
        let outer_hi = max_hi > max_lo;
        if n == 2 {
            // Two channels: the outer edges of the combined band sit exactly
            // at s = ±j by construction of the bandpass transform, so the
            // historical hardcoded ±j edges are used bit-for-bit.
            edges.push(if outer_hi {
                Complex64::new(0.0, 1.0)
            } else {
                Complex64::new(0.0, -1.0)
            });
        } else {
            // N ≥ 3: the channel's own outer band edge maps to the real
            // frequency `f_edge` on the common axis, s = j·ω(f_edge)
            // (in general not ±j, which are the combined band's outermost
            // edges).
            let f_edge = if outer_hi { hi_f } else { lo_f };
            edges.push(Complex64::new(0.0, omega(f_edge)));
        }
    }

    // Fixed numerator N: union of the channel reflection zeros mapped to the
    // common domain, plus the Type-II real zero 1/c0.
    let mut n_roots: Vec<Complex64> = Vec::new();
    for i in 0..n {
        let (a, b) = affine[i];
        for r in polys[i].f.roots()? {
            n_roots.push(Complex64::new(0.0, b + a * r.im));
        }
    }
    if let JunctionType::TypeII { c0 } = junction {
        n_roots.push(Complex64::new(1.0 / c0, 0.0));
    }
    let n_poly = ComplexPolynomial::from_complex_roots(&n_roots)?;

    // Initial S state: the channel characteristic polynomials S = (E+F)/2
    // (left-half-plane roots), as in the historical implementation.
    let half = Complex64::new(0.5, 0.0);
    let mut states: Vec<ChannelState> = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = affine[i];
        let s_init = polys[i]
            .e
            .clone()
            .add(&pad_polynomial(&polys[i].f, polys[i].e.coefficients.len())?)?
            .scale(half)?;
        states.push(ChannelState {
            order: channels[i].order,
            return_loss_db: channels[i].return_loss_db,
            a,
            b,
            p_proto: polys[i].p.clone(),
            s_proto: s_init,
            edge: edges[i],
            center_hz: centers[i],
        });
    }

    Ok(Prepared {
        states,
        n_poly,
        f0,
        b: b_common,
        polys,
    })
}

/// Solves the N×N p₀² system from the return-loss conditions.
///
/// Rows (channel i, at its outer edge eᵢ):
/// `Σⱼ |Ptⱼ(eᵢ)|²·p0ⱼ² = |N(eᵢ)|²·(1−εᵢ)/εᵢ`,  εᵢ = 10^(−rlᵢ/10).
///
/// The diagonal entries are the own-band terms (O(1)) and the off-diagonal
/// ones are the other channels' transmission polynomials at the foreign edge
/// (small for separated bands), so the system is diagonally dominant and
/// well-conditioned. For N = 2 the exact historical Cramer's-rule 2×2
/// arithmetic is preserved bit-for-bit (row 0 = channel 0's edge, row 1 =
/// channel 1's edge, columns in channel order); for N ≥ 3 a
/// partial-pivoting LU solve is used.
fn solve_p0(
    states: &[ChannelState],
    n_poly: &ComplexPolynomial,
    pts: &[ComplexPolynomial],
) -> Result<Vec<f64>> {
    let n = states.len();
    if n == 2 {
        // Bit-identical historical 2×2.
        let (e0, e1) = (states[0].edge, states[1].edge);
        let (pt0, pt1) = (&pts[0], &pts[1]);
        let n0 = n_poly.evaluate(e0).norm_sqr();
        let pt00 = pt0.evaluate(e0).norm_sqr();
        let pt10 = pt1.evaluate(e0).norm_sqr();
        let n1 = n_poly.evaluate(e1).norm_sqr();
        let pt01 = pt0.evaluate(e1).norm_sqr();
        let pt11 = pt1.evaluate(e1).norm_sqr();

        let eps0 = 10.0_f64.powf(-states[0].return_loss_db / 10.0);
        let eps1 = 10.0_f64.powf(-states[1].return_loss_db / 10.0);
        let rhs0 = n0 * (1.0 - eps0) / eps0;
        let rhs1 = n1 * (1.0 - eps1) / eps1;

        let det = pt00 * pt11 - pt01 * pt10;
        if det.abs() < 1e-30 {
            // Degenerate (nearly) singular 2×2: distribute the total equally
            // between the two channels. Reject a genuinely negative total.
            let p0_sq = (rhs0 + rhs1) / (pt00 + pt10 + pt01 + pt11);
            if p0_sq < -P0_SQ_REJECTION_TOL {
                return Err(MfsError::PreconditionViolation(format!(
                    "p₀² is non-physical (negative {p0_sq:.3e}); the channels are \
                     too strongly coupled/overlapping for a realizable \
                     quasi-equiripple junction"
                )));
            }
            let clamped = p0_sq.max(P0_SQ_FLOOR);
            return Ok(vec![clamped, clamped]);
        }
        let x0 = (rhs0 * pt11 - rhs1 * pt10) / det;
        let x1 = (pt00 * rhs1 - pt01 * rhs0) / det;
        for (i, v) in [x0, x1].iter().enumerate() {
            if *v < -P0_SQ_REJECTION_TOL {
                return Err(MfsError::PreconditionViolation(format!(
                    "channel {i}'s p₀² is non-physical (negative {v:.3e}); the \
                     channels are too strongly coupled/overlapping for a \
                     realizable quasi-equiripple junction"
                )));
            }
        }
        return Ok(vec![x0.max(P0_SQ_FLOOR), x1.max(P0_SQ_FLOOR)]);
    }

    // N ≥ 3: LU solve with partial pivoting (nalgebra).
    let mut matrix = DMatrix::<f64>::zeros(n, n);
    let mut rhs = DVector::<f64>::zeros(n);
    for i in 0..n {
        let e = states[i].edge;
        let eps = 10.0_f64.powf(-states[i].return_loss_db / 10.0);
        rhs[i] = n_poly.evaluate(e).norm_sqr() * (1.0 - eps) / eps;
        for j in 0..n {
            matrix[(i, j)] = pts[j].evaluate(e).norm_sqr();
        }
    }
    let solution = matrix
        .clone()
        .lu()
        .solve(&rhs)
        .ok_or_else(|| {
            MfsError::NumericalFailure("p₀² system is singular (LU solve failed)".into())
        })?;
    for (i, &v) in solution.iter().enumerate() {
        if v < -P0_SQ_REJECTION_TOL {
            return Err(MfsError::PreconditionViolation(format!(
                "channel {i}'s p₀² is non-physical (negative {v:.3e}); the \
                 channels are too strongly coupled/overlapping for a \
                 realizable quasi-equiripple junction"
            )));
        }
    }
    Ok((0..n).map(|i| solution[i].max(P0_SQ_FLOOR)).collect())
}

/// Computes the junction product polynomial ∏ᵢ Sᵢ from N and D.
///
/// Type-I: `(a·N + b·D)/2` with `a = 1−j·n²·b0`, `b = 1+j·n²·b0` (for the
/// ideal junction a = b = 1, i.e. `∏S = (N + D)/2`, the junction-relation
/// identity `N + D = 2·∏S`). Type-II (resonating node): `(c0/2)·(D_monic −
/// N)` — with `N = P·(s−1/c0) + (1/c0)·cross` and `D = P·(s+1/c0) +
/// (1/c0)·cross` the leading terms cancel in `D − N`, leaving
/// `(2/c0)·∏Sᵢ`. D is re-monic'd first because `spectral_factor`
/// phase-rotates its leading coefficient (≈ 1 + O(ε)); the cancellation must
/// hold to machine precision or a phantom root runs to ~1/ε.
fn compute_s_product(
    junction: JunctionType,
    n_poly: &ComplexPolynomial,
    d_poly: &ComplexPolynomial,
) -> Result<ComplexPolynomial> {
    match junction {
        JunctionType::TypeI { n, b0 } => {
            let n_sq_b0 = n * n * b0;
            let a = Complex64::new(1.0, -n_sq_b0);
            let b = Complex64::new(1.0, n_sq_b0);
            let half = Complex64::new(0.5, 0.0);
            let a_n = n_poly.scale(a)?;
            let b_d = d_poly.scale(b)?;
            a_n.add(&b_d)?.scale(half)
        }
        JunctionType::TypeII { c0 } => {
            let d_lead = d_poly.leading_coefficient();
            // Guard the monic normalization: `spectral.rs::monic` leaves a
            // degenerate factor (|lead| < 1e-15) unnormalized, and dividing by
            // such a tiny coefficient would scale the polynomial across ~200
            // decades and produce catastrophic roots. Surface a typed error
            // instead (review round-2, new defect 1).
            if d_lead.norm() < 1e-15 {
                return Err(MfsError::NumericalFailure(
                    "Type-II junction: spectral factor D(s) has a near-zero \
                     leading coefficient; cannot monic-normalize"
                        .into(),
                ));
            }
            let d_monic = d_poly.scale(Complex64::new(1.0, 0.0) / d_lead)?;
            d_monic.sub(n_poly)?.scale(Complex64::new(c0 / 2.0, 0.0))
        }
    }
}

/// Assigns the ∏Sᵢ product roots to contiguous blocks (paper Step 4):
/// sort by (imaginary, real) and take consecutive runs of the requested
/// sizes. Block order is the caller's band-position ordering (ascending
/// channel center); for two channels that is exactly the historical
/// lower/upper split.
fn assign_blocks(
    product_roots: &[Complex64],
    block_sizes: &[usize],
) -> Result<Vec<Vec<Complex64>>> {
    let total: usize = block_sizes.iter().sum();
    if product_roots.len() < total {
        return Err(MfsError::NumericalFailure(format!(
            "∏Sᵢ product has {} roots, expected {total}",
            product_roots.len()
        )));
    }
    let mut sorted = product_roots.to_vec();
    sorted.sort_by(|x, y| {
        x.im.partial_cmp(&y.im)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| x.re.partial_cmp(&y.re).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut blocks = Vec::with_capacity(block_sizes.len());
    let mut offset = 0usize;
    for &count in block_sizes {
        blocks.push(sorted[offset..offset + count].to_vec());
        offset += count;
    }
    Ok(blocks)
}

/// Runs the full fixed-point iteration to convergence (strict tolerance or
/// high-degree stagnation accept).
///
/// `states` must come from [`prepare`], which also returns the
/// combined-band constants (f₀, b) stamped into the captured
/// [`JunctionCommon`].
pub fn run_full_iteration(
    states_in: Vec<ChannelState>,
    n_poly: &ComplexPolynomial,
    junction: JunctionType,
    f0: f64,
    b: f64,
    max_iterations: usize,
    convergence_tol: f64,
) -> Result<FullIterationState> {
    let n = states_in.len();
    let orders: Vec<usize> = states_in.iter().map(|st| st.order).collect();

    // Root state (prototype domain), exactly as the historical 2-channel
    // implementation tracks it: the S = (E+F)/2 left-half-plane roots.
    let mut s_roots: Vec<Vec<Complex64>> = states_in
        .iter()
        .map(|st| st.s_proto.roots())
        .collect::<Result<Vec<_>>>()?;

    // Block order for the root assignment: ascending channel-center (==
    // ascending band position for the well-separated case; identical to the
    // historical lower/upper split for two channels).
    let mut block_idx: Vec<usize> = (0..n).collect();
    block_idx.sort_by(|&x, &y| {
        states_in[x]
            .center_hz
            .partial_cmp(&states_in[y].center_hz)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let block_sizes: Vec<usize> = block_idx.iter().map(|&i| orders[i]).collect();

    let mut iterations = 0usize;
    let mut convergence_error = f64::INFINITY;
    let mut best_delta = f64::INFINITY;
    let mut stagnant = 0usize;
    let mut conv: Option<(
        ComplexPolynomial,
        ComplexPolynomial,
        Vec<ComplexPolynomial>,
        Vec<f64>,
    )> = None;

    for iter in 0..max_iterations {
        iterations = iter + 1;

        // ── State (prototype) → common domain ───────────────────────────
        let mut s_common: Vec<ComplexPolynomial> = Vec::with_capacity(n);
        let mut p_common: Vec<ComplexPolynomial> = Vec::with_capacity(n);
        for i in 0..n {
            let roots: Vec<Complex64> = s_roots[i]
                .iter()
                .map(|&r| states_in[i].a * r + Complex64::new(0.0, states_in[i].b))
                .collect();
            s_common.push(ComplexPolynomial::from_complex_roots(&roots)?);
            p_common.push(states_in[i].p_common()?);
        }

        // Transmission polynomials: Ptᵢ = Pᵢ · ∏_{j≠i} Sⱼ.
        let mut pts: Vec<ComplexPolynomial> = Vec::with_capacity(n);
        for (i, p) in p_common.iter().enumerate() {
            let mut acc = p.clone();
            for (j, s_j) in s_common.iter().enumerate() {
                if j == i {
                    continue;
                }
                acc = acc.multiply(s_j)?;
            }
            pts.push(acc);
        }

        let p0_sq = solve_p0(&states_in, n_poly, &pts)?;

        // R = N·Ñ + Σᵢ p0ᵢ²·Ptᵢ·Pt̃ᵢ.
        let n_conj = n_poly.alternating_conjugate()?;
        let mut r_poly = n_poly.multiply(&n_conj)?;
        for (pt, scale) in pts.iter().zip(p0_sq.iter()) {
            let pt_conj = pt.alternating_conjugate()?;
            let ptpt = pt.multiply(&pt_conj)?;
            r_poly = r_poly.add(&ptpt.scale(Complex64::new(*scale, 0.0))?)?;
        }

        // D carries the same degree as N (Type-II: the extra real zero is in
        // N, so deg N = Σ nᵢ + 1 covers both junctions).
        let d_poly = spectral_factor(&r_poly, n_poly.degree())?;

        // New ∏Sᵢ product and its roots, assigned in band-position blocks.
        let product_poly = compute_s_product(junction, n_poly, &d_poly)?;
        let product_roots = product_poly.roots()?;
        let blocks = assign_blocks(&product_roots, &block_sizes)?;
        // Map band-position blocks back to channels.
        let mut new_proto: Vec<Vec<Complex64>> = vec![Vec::new(); n];
        for (slot, &ch) in block_idx.iter().enumerate() {
            new_proto[ch] = states_in[ch].to_proto(&blocks[slot]);
        }

        // ── Common → prototype (pure fixed-point update) ────────────────
        let mut max_delta = 0.0_f64;
        for i in 0..n {
            if new_proto[i].len() != orders[i] {
                // §C.6: `to_proto` drops *both* non-finite roots and finite
                // roots with |r| ≥ 50 (the historical stability filter). The
                // old message ("non-finite root set") was misleading for the
                // finite-but-out-of-range case — classify the offending mapped
                // root and report what actually happened, naming it.
                let slot = block_idx
                    .iter()
                    .position(|&c| c == i)
                    .expect("block_idx is a permutation of channels");
                let (a, b) = (states_in[i].a, states_in[i].b);
                let dropped = orders[i].saturating_sub(new_proto[i].len());
                let offending = blocks[slot].iter().find_map(|&r| {
                    let mapped = (r - Complex64::new(0.0, b)) / a;
                    if !mapped.re.is_finite() || !mapped.im.is_finite() {
                        Some(format!("non-finite mapped root {r}"))
                    } else if mapped.norm() >= 50.0 {
                        Some(format!(
                            "finite root {r} (prototype |r| = {:.3} >= 50, dropped)",
                            mapped.norm()
                        ))
                    } else {
                        None
                    }
                });
                let detail = offending
                    .map(|o| format!(", e.g. {o}"))
                    .unwrap_or_default();
                return Err(MfsError::NumericalFailure(format!(
                    "Full iteration: channel {i} (center {} Hz) lost {} of {} \
                     prototype roots during the common->prototype map because \
                     they are non-finite or |r| >= 50{detail}; expected {} roots",
                    states_in[i].center_hz,
                    dropped,
                    orders[i],
                    orders[i],
                )));
            }
            max_delta = max_delta.max(compute_root_displacement(&s_roots[i], &new_proto[i]));
        }
        s_roots = new_proto;
        convergence_error = max_delta;

        // Fast path: strict convergence. Stagnation path: high-degree
        // root-extraction noise floor exceeds the tolerance while the
        // physical fixed point (p0², D) has converged — accept once the
        // displacement stops improving and is already small.
        if max_delta < convergence_tol
            || (iter + 1) >= STAGNATION_MIN_ITER
                && stagnant >= STAGNATION_WINDOW
                && best_delta < STAGNATION_ACCEPT
        {
            // §C.5: captured at convergence, but note the one-iteration lag.
            // `d_poly`/`pts`/`p0_sq` were computed from `s_roots` *before* the
            // update `s_roots = new_proto` (above), so `common.{d,pt,p0_sq}`
            // describe the junction as of iteration k while the returned
            // `s_roots` field is iteration k+1. The two are therefore off by
            // one iteration: `common` (which feeds the junction-level output
            // `junction_response`) lags `s_roots` (which feeds per-channel
            // recovery) by exactly this iteration's root displacement, which is
            // reported as the returned `convergence_error` — < convergence_tol
            // on the strict-convergence path, and O(best_delta <
            // STAGNATION_ACCEPT, 1e-3) when accepted by stagnation. It is not
            // recomputed from the final `s_roots` because doing so would change
            // the captured output (convergence behavior).
            conv = Some((n_poly.clone(), d_poly, pts, p0_sq));
            break;
        }
        if max_delta < best_delta * 0.99 {
            best_delta = max_delta;
            stagnant = 0;
        } else {
            stagnant += 1;
        }
    }

    let converged = conv.is_some();
    let common = conv.map(
        |(n_poly_, d_poly_, pts_, p0_sq_)| JunctionCommon {
            n: n_poly_,
            d: d_poly_,
            pt: pts_,
            p0_sq: p0_sq_,
            f0,
            b,
        },
    );
    Ok(FullIterationState {
        iterations,
        convergence_error,
        converged,
        s_roots,
        common,
    })
}

/// Full fixed-point synthesis for N channels sharing one junction.
///
/// Returns the common-domain setup (per-channel states, fixed N, band
/// constants, and the independently synthesized polynomial sets needed to
/// recover each channel's design) together with the iteration outcome:
/// per-channel converged S roots and, when the iteration was captured, the
/// common-domain junction polynomials.
pub fn synthesize_full(
    channels: &[ChannelSpec],
    junction: JunctionType,
    max_iterations: usize,
    convergence_tol: f64,
) -> Result<(Prepared, FullIterationState)> {
    let prepared = prepare(channels, junction)?;
    let state = run_full_iteration(
        prepared.states.clone(),
        &prepared.n_poly,
        junction,
        prepared.f0,
        prepared.b,
        max_iterations,
        convergence_tol,
    )?;
    Ok((prepared, state))
}

/// Common-junction response over a frequency sweep (the physical observable
/// of the full path): `S_ant = N/D`, `Sᵢ = √p0ᵢ·Ptᵢ/D` in `s = j·Ω(f)`.
#[derive(Debug, Clone)]
pub struct JunctionResponse {
    /// Frequencies (Hz).
    pub frequencies: Vec<f64>,
    /// Common (antenna) port reflection in dB.
    pub antenna_s11_db: Vec<f64>,
    /// Per-channel transmission in dB, channel-major.
    pub channel_s21_db: Vec<Vec<f64>>,
}

/// Evaluates the junction-level response from the common-domain polynomials.
pub fn junction_response(
    common: &JunctionCommon,
    channel_count: usize,
    start_hz: f64,
    stop_hz: f64,
    points: usize,
) -> Result<JunctionResponse> {
    let grid = FrequencyGrid::linspace(start_hz, stop_hz, points)?;
    let mut channel_s21_db = vec![Vec::with_capacity(points); channel_count];
    let mut antenna_s11_db = Vec::with_capacity(points);
    let mut frequencies = Vec::with_capacity(points);

    for &f in grid.as_slice() {
        let omega = (f / common.f0 - common.f0 / f) / common.b;
        let s = Complex64::new(0.0, omega);
        let dd = common.d.evaluate(s);
        if dd.norm() < 1e-300 {
            return Err(MfsError::NumericalFailure(
                "junction response: D vanishes on the sweep".into(),
            ));
        }
        let dn = common.n.evaluate(s);
        frequencies.push(f);
        antenna_s11_db.push(magnitude_to_db(dn.norm() / dd.norm()));
        for (i, pt) in common.pt.iter().enumerate() {
            let p0 = common.p0_sq[i].sqrt();
            channel_s21_db[i].push(magnitude_to_db(p0 * pt.evaluate(s).norm() / dd.norm()));
        }
    }
    Ok(JunctionResponse {
        frequencies,
        antenna_s11_db,
        channel_s21_db,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::ComplexPolynomial;

    #[test]
    fn type_ii_s_product_rejects_degenerate_leading_coefficient() {
        // A near-zero leading coefficient (the degenerate path
        // `spectral.rs::monic` leaves unnormalized) must be rejected rather
        // than divided into — that division would scale the polynomial across
        // ~200 decades and produce catastrophic roots (review round-2, new
        // defect 1).
        let d = ComplexPolynomial::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1e-20, 0.0),
        ])
        .unwrap();
        let n =
            ComplexPolynomial::new(vec![Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)]).unwrap();
        let err = compute_s_product(JunctionType::TypeII { c0: 1.0 }, &n, &d).unwrap_err();
        assert!(
            matches!(err, MfsError::NumericalFailure(_)),
            "expected NumericalFailure, got {err:?}"
        );
        assert!(err.to_string().contains("leading coefficient"));
    }
     #[test]
     fn full_path_tz_at_neighbor_band_edge_converges() {
         // Regression guard: a finite transmission zero at the *neighboring*
         // channel's band edge (B3's zero at B1's lower edge 1925 MHz, B1's
         // zero at B3's upper edge 1885 MHz) places several roots of the
         // junction R polynomial within ΔΩ ≈ 0.01 of the jω axis. The
         // companion-matrix eigenvalues of such a clustered near-axis R carry
         // ~1e-3 forward error, which broke the T-symmetric mirror pairing in
         // `select_lhp_roots` and failed the para-product verification at
         // iteration 1 (typed NumericalFailure from `synthesize_full`). The
         // Newton polish in `companion_spectral_factor` restores the pairing;
         // this test pins the production path on that worst case and asserts
         // convergence plus a real junction-level passband with working TZs.
         use crate::freq::FrequencyMapping;
         fn spec(order: usize, lo: f64, hi: f64, tz: Option<f64>) -> ChannelSpec {
             let (c, bw) = crate::freq::edges_to_center_bandwidth(lo, hi);
             let m = crate::freq::BandPassMapping::new(c, bw).unwrap();
             let zeros = match tz {
                 Some(f) => vec![m.map_hz_to_normalized(f).unwrap().omega],
                 None => Vec::new(),
             };
             ChannelSpec { order, return_loss_db: 20.0, zeros, center_hz: c, bandwidth_hz: bw }
         }
         let junction = JunctionType::TypeI { n: 1.0, b0: 0.0 };
         let (_, state) = synthesize_full(
             &[spec(4, 1805e6, 1890e6, Some(1925e6)),
              spec(7, 1920e6, 2170e6, Some(1885e6)),
              spec(5, 2500e6, 2690e6, None)],
             junction,
             80,
             1e-8,
         )
         .expect("full path must converge with TZs at the neighbor band edges");
         assert!(state.converged, "full path exhausted its iteration budget");

         // The junction response must show a real passband at each band
         // center (the recovered D is a valid spectral factor, not just a
         // non-erroring one) and a deep dip at each transmission zero.
         let common = state.common.as_ref().expect("captured junction state");
         let resp = junction_response(common, 3, 1750e6, 2750e6, 2001).unwrap();
         let nearest = |f: f64| {
             resp.frequencies
                 .iter()
                 .enumerate()
                 .min_by(|(_, a), (_, b)| (*a - f).abs().partial_cmp(&(*b - f).abs()).unwrap())
                 .map(|(i, _)| i)
                 .unwrap()
         };
         let centers = [1.847e9, 2.041e9, 2.593e9];
         for (i, &c) in centers.iter().enumerate() {
             let db = resp.channel_s21_db[i][nearest(c)];
             assert!(db > -1.0, "channel {i} junction S21 at band center: {db:.3} dB");
         }
         let db_b3_tz = resp.channel_s21_db[0][nearest(1925e6)];
         let db_b1_tz = resp.channel_s21_db[1][nearest(1885e6)];
         assert!(db_b3_tz < -10.0, "B3 TZ dip at 1925 MHz: {db_b3_tz:.3} dB");
         assert!(db_b1_tz < -10.0, "B1 TZ dip at 1885 MHz: {db_b1_tz:.3} dB");
         assert!(
             resp.antenna_s11_db.iter().all(|x| x.is_finite()),
             "antenna S11 non-finite"
         );
     }
}
