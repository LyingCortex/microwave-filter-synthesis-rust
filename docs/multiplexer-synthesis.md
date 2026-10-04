# Multiplexer (Multi-Channel) Synthesis: Design & Plan

## Purpose

This document defines the design and implementation plan for extending the MFS
crate from a **2-channel diplexer** to a general **N-channel multiplexer**
(三工器/四工器/多工器). It generalizes the existing Macchiarella-Tamiazzo
iterative diplexer synthesis to an arbitrary number of channels sharing a
single common (antenna) junction.

The plan is written first so the work can be split across multiple agents with
clear, non-overlapping responsibilities and a well-defined integration point.

## Current State

> **Status as of 0.5.0:** the N-channel multiplexer described below has been
> implemented and shipped. The pre-implementation state recorded here is
> kept for context; "Current State" now reflects the code.

- `src/duplexer.rs` (~950 lines) implements the 2-channel (TX/RX) duplexer:
  - `Duplexer` — analysis: combines two `FilterDesign`s, evaluates the 3-port
    response.
  - `DuplexerSynthesis` — Macchiarella-Tamiazzo iterative synthesis.
  - `JunctionType::{TypeI{ n, b0 }, TypeII{ c0 }}` — junction model.
  - `ChannelSpec` — single-channel spec (order, RL, zeros, center, bandwidth).
  - `DuplexerSynthesisResult` — `{ tx, rx, iterations, convergence_error,
    junction, s_tx_roots, s_rx_roots }`.
  - The full spectral-factorization path is opt-in via
    `DuplexerSynthesis::full_iteration(true)` (0.5.0; previously a hardcoded
    `false` flag). It is **not** numerically converged for order ≥ 4 — see
    `docs/architecture.md` Known Limitations and
    `docs/multiplexer-code-review.md`. The active path is the simplified
    frequency-offset perturbation.
- `src/multiplexer.rs` (~700 lines) — the N-channel generalization shipped
  with this design; the post-iteration channel recovery is the shared
  `mfs::synthesis::recover_channel_design` (0.5.0+; previously a
  per-module `recover_filter_design` copy) which writes the perturbed E/F
  back to `generalized.e_s/f_s` (0.5.0 P1 fix), so the mutual loading
  actually propagates into the output coupling matrices.
- `src/spectral.rs` — spectral factorization (`spectral_factor`), used by the
  full-iteration path.
- Existing doc `docs/duplexer-synthesis.md` derives the 2-channel math.

## Design Goal

A general `Multiplexer` with N channels all connected to one common port:

```
                    ┌── Channel 1 filter ── Port 2
Port 1 ── Common ──┼── Channel 2 filter ── Port 3
(antenna)  junction│  ...
                    └── Channel N filter ── Port (N+1)
```

Synthesis must produce N channel `FilterDesign`s whose **combined** common-port
reflection is quasi-equiripple, accounting for the mutual loading of **all**
channels through the junction (not just 2).

## Mathematical Generalization (N Channels)

The existing 2-channel derivation (see `docs/duplexer-synthesis.md`) extends
directly. For channels $i = 1..N$ with auxiliary polynomials
$S_i = (E_i + F_i)/2$, $D_i = (E_i - F_i)/2$:

### Common-port numerator (Type-I ideal, $n=1, b_0=0$)

The diplexer N(s) is $S_1 S_2 - n^2 (D_1 S_2 + S_1 D_2)/(1 - jn^2 b_0)$. The
N-channel form sums each channel's admittance contribution against the product
of the others:

$$N(s) = \prod_{i=1}^{N} S_i(s) \;-\; \sum_{i=1}^{N} \frac{D_i(s)}{S_i(s)}\prod_{j=1}^{N} S_j(s)\cdot w_i(s)$$

where $w_i$ are the junction weighting terms (Type-I transformer ratio $n_i$,
Type-II node capacitance $c_0$). In practice this is implemented as a
polynomial sum of the form:

$$N(s) = \prod_i S_i(s) \;-\; \sum_i \alpha_i\, D_i(s) \prod_{j \neq i} S_j(s)$$

which keeps every term a well-defined product of the N channel polynomials
(no polynomial division needed).

### Unitary condition and p₀ scaling

$$|D(s)|^2 = |n_0|^2 |N(s)|^2 + \sum_{i=1}^{N} |p_{0i}|^2 |P_i(s)|^2,
\qquad P_i(s) = P_{in}(s)\prod_{j \neq i} S_j(s)$$

Imposing return loss at each channel's passband edge ($s = \pm j$) gives an
**N×N linear system** for the $|p_{0i}|^2$ (the diplexer's 2×2 case generalizes
to N×N).

### Iterative root update

Factor $\prod_{i=1}^{N} S_i(s)$ from $(a N(s) + b D(s))/2$, assign roots to the
N channels, under-relax, and iterate until the max root displacement is below
tolerance (diplexer converges in 5-10 iterations; multiplexers typically
converge in 10-30).

### Filter recovery

$F_i(s) = S_i(s) - D_i(s)$, $E_i(s) = S_i(s) + D_i(s)$, with $D_i$ recovered by
interpolation at the roots of $S_i$ — identical to the diplexer, applied to
each channel.

## New API Surface

Add a new module `src/multiplexer.rs` (do **not** disturb the existing
`Duplexer` API; keep `Duplexer` as the N=2 specialization for backwards
compatibility).

```rust
// Channel descriptor shared by synthesis + analysis.
pub struct ChannelSpec {            // reuse / extend existing, add a `label`
    pub label: String,
    pub order: usize,
    pub return_loss_db: f64,
    pub zeros: Vec<f64>,
    pub center_hz: f64,
    pub bandwidth_hz: f64,
}

// N-port analysis.
pub struct Multiplexer {
    channels: Vec<FilterDesign>,
}
impl Multiplexer {
    pub fn new(channels: Vec<FilterDesign>) -> Result<Self>;
    pub fn response(&self, freq: &FrequencyGrid) -> Result<MultiplexerResponse>;
    pub fn channel(&self, i: usize) -> &FilterDesign;
    pub fn channel_folded(&self, i: usize) -> Result<CouplingMatrix>;
}

// (N+1)-port response.
pub struct MultiplexerResponse {
    pub common_reflection: Vec<Complex64>,      // S11 at antenna
    pub channel_transmissions: Vec<Vec<Complex64>>, // S(i+1)1
    pub channel_isolation: Vec<Vec<Vec<Complex64>>>, // S(i+1)(j+1)
}

// Builder + iterative synthesis.
pub struct MultiplexerSynthesis {
    channels: Vec<ChannelSpec>,
    junction: JunctionType,
    max_iterations: usize,
    convergence_tol: f64,
    full_iteration: bool,   // default false (simplified perturbation)
}
impl MultiplexerSynthesis {
    pub fn new() -> Self;
    pub fn channel(mut self, spec: ChannelSpec) -> Self;   // push one
    pub fn channels(mut self, specs: impl IntoIterator<Item = ChannelSpec>) -> Self;
    pub fn junction(mut self, j: JunctionType) -> Self;
    pub fn max_iterations(mut self, n: usize) -> Self;
    pub fn convergence_tol(mut self, tol: f64) -> Self;
    pub fn full_iteration(mut self, enabled: bool) -> Self; // B2: shared spectral path
    pub fn synthesize(&self) -> Result<MultiplexerSynthesisResult>;
}

pub struct MultiplexerSynthesisResult {
    pub channels: Vec<FilterDesign>,
    pub iterations: usize,
    pub convergence_error: f64,
    pub converged: bool,             // B3: stopped on own criterion, not the budget
    pub junction: JunctionType,
    pub s_roots: Vec<Vec<Complex64>>,   // per-channel S-polynomial roots
    pub max_iterations: usize,
    pub junction_common: Option<JunctionCommon>, // B2: full-path common-domain N/D/Pt/p0²
}
```

`DuplexerSynthesis` is retained and can be re-expressed internally as
`MultiplexerSynthesis` with `channels = [tx, rx]` (parity test).

## N-Port Response Evaluation

The 3-port diplexer response uses the combination of two filter
admittances. For N channels, build the common-port input admittance as the
sum of each channel's input admittance seen through the junction:

$$Y_{in}(\omega) = \sum_{i=1}^{N} \frac{1}{Z_{in,i}^{(\text{shifted})}(\omega)}$$

and derive $S_{11} = (1-Y_{in})/(1+Y_{in})$ (scaled), with each
$S_{(i+1)1}$ obtained from the individual channel transmission times the
power split at the common port. Implement by summing normalized
admittance contributions per channel; reuse the existing
`ResponseSolver`/`evaluate_normalized` per channel where possible.

## Phase Breakdown (for parallel agents)

The work is cohesive (all lands in the same crate), so phases are sequenced
around a shared core. Each agent owns a distinct file/surface to avoid
edit conflicts.

| Phase | Agent | Owner / Scope | Depends on |
|-------|-------|---------------|-----------|
| 0 | (plan) | `docs/multiplexer-synthesis.md`, task split | — |
| 1 | **A** | `src/multiplexer.rs`: generalized N-channel core — `ChannelSpec`, generalized N(s) for N channels, N×N p₀ solve, product factorization + root assignment | — |
| 2 | **B** | `src/multiplexer.rs`: `Multiplexer` N+1-port analysis + `MultiplexerResponse` | A |
| 3 | **C** | `src/multiplexer.rs`: `MultiplexerSynthesis` builder + result; `src/lib.rs` re-export; `src/bin/mfs_cli.rs` subcommand; `src/python.rs` bindings | A |
| 4 | **D** | `tests/multiplexer.rs`: diplexer parity (N=2 vs `DuplexerSynthesis`), triplexer/quadruplexer RL & isolation checks | A,B,C |
| 5 | (integrate) | `cargo test` all green; docs cross-reference; CHANGELOG | all |

## Implementation Deviations (as built)

The landed code (`src/multiplexer.rs`, `tests/multiplexer.rs`) deviates from the
API and math described above in the following ways. These are recorded here so
the doc reflects reality, not aspiration.

1. **`ChannelSpec` has no `label` field.** The `ChannelSpec` struct is reused
   from `crate::duplexer` (`duplexer.rs:169`) and only carries
   `{ order, return_loss_db, zeros, center_hz, bandwidth_hz }`. The doc's
   promised `label: String` was not added, to avoid disturbing the shared
   duplexer type.

2. **Response is now a true (N+1)-port model (B1, 2026-10-04); the
   geometric-mean model is retained behind a flag.** `Multiplexer::response`
   defaults to [`ResponseModel::Admittance`], which treats each channel as a
   lossless reciprocal 2-port (input port A at the common node, output port B
   at the external channel port) and composes the full complex (N+1)×(N+1)
   S-matrix by the wave-variable node equation. For lossless channels the
   result is unitary and reciprocal, and N=1 reduces exactly to the channel's
   own 2-port S-matrix (gate tests in `src/multiplexer.rs`). The
   `MultiplexerResponse` struct still stores dB (`f64`) arrays, but they are
   the *true* (N+1)-port S-magnitude values: `antenna_s11_db = |S₀₀|`,
   `channel_s21_db[i] = |Sᵢ₊₁₀|`, `channel_s11_db[i] = |Sᵢ₊₁,ᵢ₊₁|`,
   `isolation_db[i][j] = |Sᵢ₊₁,ⱼ₊₁|`. The full complex matrix is exposed
   directly by `Multiplexer::response_matrix`. The legacy geometric-mean
   magnitude model (antenna = `(∏ᵢ|S11ᵢ|)^(1/N)`, isolation = product of S21
   magnitudes) is kept as `ResponseModel::GeometricMean` via
   `Multiplexer::response_with` for one release cycle for parity comparison;
   it approximates the admittance model well in the deep stopband (both
   channels reflective) but is off by up to ~150 dB in a channel's passband
   (the single in-band channel drives the geometric mean toward zero while
   the other channels' stopband reactances load the common node).

   The `Duplexer::response` (N=2 diplexer) still uses its own geometric-mean
   model and is unchanged.

3. **p₀ N×N linear system — implemented in the full path (B2, 2026-10-04).**
   The doc's §Unitary condition (solving for `|p_{0i}|²` from a
   return-loss constraint per passband edge) is implemented in the shared
   `src/synthesis/machiarella.rs` full-spectral loop (opt-in via
   `MultiplexerSynthesis::full_iteration(true)`, also used by
   `DuplexerSynthesis::full_iteration(true)`): N×N `nalgebra` LU solve with
   partial pivoting for N ≥ 3 (the historical Cramer's-rule 2×2 is preserved
   verbatim for N = 2). The default simplified perturbation path is unchanged
   and performs no N×N solve.

4. **`recover_filter_design` yields fully-reflective channels — RESOLVED.**
   This was previously reported as a KNOWN BUG. Investigation revealed the root
   cause was **not** in `recover_filter_design` at all: `build_transversal_from_residues`
   (`src/synthesis/residues.rs`) extracted coupling values from residue **real parts**
   only, which is zero for odd-order all-pole designs (an `s = j w` phase factor on the
   transfer residue `r12`). This affected every odd-order filter in the crate, not just
   the multiplexer. The fix (`residue_real_component`) recovers the physical coupling
   magnitude/sign for both parities. The diagnostic test
   `multiplexer_recovered_channel_is_passband` is now un-`#[ignore]`d and passes.

5. **Mutual loading now reaches the output (P1, 0.5.0) — the earlier
   no-op limitation is resolved.** Originally, the probe below showed the
   recovered channel matrix **bit-identical** (`max |Δ| = 0`) to independent
   per-channel synthesis, because `MatrixSynthesisEngine.synthesize` reads
   `generalized.e_s/f_s/p_s` and ignored the `PolynomialSet.e/f` that
   `recover_filter_design` rebuilt from the perturbed S-roots. Since 0.5.0 the
   duplexer/multiplexer recovery paths write the perturbed E/F back to
   `generalized.e_s/f_s`, so the loading correction reaches the output matrices
   (`adjacent_channels_show_loading` in `tests/multiplexer.rs` asserts the
   non-degenerate effect, and
   `matrix_synthesizer_is_sensitive_to_generalized_e_s` guards the engine seam).
   **Remaining model caveat (simplified perturbation, not the full Macchiarella
   N×N p₀ spectral solve):** the correction does *not* converge to independent
   synthesis even for well-separated channels — the cross-channel term
   `D_j/S_j` approaches a non-zero constant — so far-apart channels pay a small
   passband-peak cost. This is a documented approximation, not a propagation
   gap. (The pre-0.5.0 bit-identical state is preserved in this document's
   git history.)

## Post-review status

The review findings and their resolutions are recorded in
`docs/multiplexer-code-review.md`. As of 0.5.0 (P1) the mutual-loading
perturbation **does** propagate into the returned channel matrices (deviation
item 5's no-op is resolved); the *default* simplified path remains a simplified
model (not the full N×N p₀ spectral solve) and does not reach machine
tolerance for well-separated channels. The produced channel filters are valid
passband filters, which is the verifiable contract those tests assert.

**B2 (2026-10-04): the full spectral path now reaches N channels.** The
converged 2-channel Macchiarella–Tamiazzo loop is shared with the multiplexer
via `src/synthesis/machiarella.rs` (`MultiplexerSynthesis::full_iteration(true)`
and `DuplexerSynthesis::full_iteration(true)` call the same code). For
N ≥ 3 it generalizes to the N×N p₀² solve (LU), per-channel outer-edge
return-loss selection, and contiguous N-block root assignment. The N = 2
branch is bit-identical to the historical loop, so the gate-anchored diplexer
results are unchanged. Gated by
`multiplexer_full_iteration_triplexer_converges` / `_quadruplexer_converges`
(convergence, 0 dB per-channel passbands, 22 dB antenna RL at each outer edge,
and the para-identity D·D̃ = R losslessness invariant) and
`multiplexer_full_iteration_rejects_nonideal_type_i_for_n_gt_2` (explicit
rejection of a non-ideal Type-I junction for N > 2). The full path is opt-in
and heavier (spectral factorization per iteration); the lightweight simplified
perturbation remains the default.

## Verification / Acceptance Criteria

**Status legend:** ✅ met · ⚠️ met structurally but the *quantitative* value is
not asserted (approximate model / tolerance).

1. ✅ **Diplexer parity**: `MultiplexerSynthesis` with N=2 reproduces the
   existing `DuplexerSynthesis` band placement (`tests/multiplexer.rs`).
   ("Same convergence" in the strict numeric sense is not asserted; both paths
   share the same approximate perturbation loop.) The two paths use *different*
   common-port response models (multiplexer: true (N+1)-port admittance;
   duplexer: geometric mean), so cross-implementation response tracking is not
   asserted — see the note in that test.
2. ⚠️ **Triplexer** (3 channels): synthesizes, converges to a finite bounded
   error, and the recovered channels are valid passband filters (no longer
   fully-reflective). **Caveat (deviation item 5, as of 0.5.0):** the
   mutual-loading correction now reaches the output matrices (P1), but the
   simplified perturbation model does not converge to independent synthesis even
   for well-separated channels, so adjacent channels pay a passband-peak cost.
   The dimensional N-port contract is asserted strictly. The quasi-equiripple
   RL / strong-isolation *values* are not tightly asserted (approximate model).
3. ⚠️ **Quadruplexer** (4 channels): same status as the triplexer.
4. ✅ **Response (true (N+1)-port, B1)**: `Multiplexer::response_matrix`
   returns the full complex (N+1)×(N+1) S-matrix; it is unitary and
   reciprocal for lossless channel filters and reduces exactly to the channel's
   own 2-port S-matrix at N=1 (gate tests in `src/multiplexer.rs`).
5. ✅ **Convergence (full path, B3)**: the `full_iteration(true)` path stops
   on its own criterion for well-separated, adjacent, and the vendored
   standard's N=4 case — asserted by `multiplexer_full_iteration_*` gates
   (machine-tolerance `convergence_error` for well-separated channels;
   quasi-equiripple common port for adjacent). The *default* simplified path
   still documents its bounded-residual caveat (deviation item 5). The
   vendored-standard `convergence` gate is now **PASS** on every runnable case
   (`benches/multiplexer.rs` runs the `poly` stage through the full path).
6. ✅ **No regression**: `cargo test --no-fail-fast` passes across all
   integration suites, lib tests, and doctests (the lib suite is fully green as
   of 2026-09-28; the "3 pre-existing lib failures" noted at review time are
   gone). The recovered-channel passband diagnostic is un-`#[ignore]`d and
   passes.

## Risks & Mitigations

- **Convergence for narrow-band / close channels**: multiplexers can be
  ill-conditioned. Mitigate with per-channel under-relaxation factors and a
  normalized-frequency separation guard.
- **N×N linear solve stability**: use `nalgebra`'s LU solve (already a
  dependency) with a small regularization; the matrix is diagonally dominant
  for well-separated bands.
- **Spectral factorization cost** at high N: the current crate uses a
  simplified perturbation path; keep it as the default and expose the full
  spectral path behind a flag, mirroring the diplexer's approach.
- **API churn**: keep `Duplexer*` untouched and treat N=2 as a thin wrapper, so
  existing users are unaffected.

## References

- G. Macchiarella and S. Tamiazzo, "Novel Approach to the Synthesis of
  Microwave Diplexers," IEEE TMTT, vol. 54, no. 12, 2006.
- G. Macchiarella, "Synthesis of a duplexer / multiplexer" later works
  generalizing the common-port polynomial method to N channels.
- Existing `docs/duplexer-synthesis.md` (2-channel math).