# Duplexer Synthesis: Macchiarella-Tamiazzo Method

## Reference

G. Macchiarella and S. Tamiazzo, "Novel Approach to the Synthesis of Microwave
Diplexers," IEEE Trans. Microwave Theory Tech., vol. 54, no. 12, pp. 4281-4290,
December 2006.

## Overview

This method synthesizes a diplexer by iteratively computing the characteristic
polynomials of the complete three-port network, accounting for the mutual loading
between the TX and RX filters through the junction. Unlike the naive approach of
designing filters independently, this method produces filters whose combined
response at the antenna port exhibits quasi-equiripple return loss (within ~2 dB
of ideal).

## Diplexer Architecture

```
        ┌─── TX Filter ─── Port 2 (TX)
Port 1 ─┤ Junction (N)
(Antenna)└─── RX Filter ─── Port 3 (RX)
```

Two junction types:
- **Type-I** (waveguide H-plane T): transformer ratio `n`, susceptance `jb₀`
- **Type-II** (coaxial resonating node): extra resonator with capacitance `c₀`

## Mathematical Model

### Filter Polynomials

Each filter is characterized by:

$$
S_{11}^{TX} = \frac{F_{TX}(s)}{E_{TX}(s)}, \quad S_{21}^{TX} = \frac{p_{0TX} P_{TXn}(s)}{E_{TX}(s)}
$$

$$
S_{11}^{RX} = \frac{F_{RX}(s)}{E_{RX}(s)}, \quad S_{21}^{RX} = \frac{p_{0RX} P_{RXn}(s)}{E_{RX}(s)}
$$

### Auxiliary Polynomials

$$
S_{TX} = \frac{E_{TX} + F_{TX}}{2}, \quad D_{TX} = \frac{E_{TX} - F_{TX}}{2}
$$

$$
S_{RX} = \frac{E_{RX} + F_{RX}}{2}, \quad D_{RX} = \frac{E_{RX} - F_{RX}}{2}
$$

### Diplexer Polynomials (Type-I Junction)

The three-port S-parameters are:

$$
S_{11} = \frac{n_0 N(s)}{D(s)}, \quad S_{21} = \frac{p_{0t} P_t(s)}{D(s)}, \quad S_{31} = \frac{p_{0r} P_r(s)}{D(s)}
$$

Where:

$$
N(s) = S_{TX} S_{RX} - n^2 \frac{D_{TX} S_{RX} + S_{TX} D_{RX}}{1 - jn^2 b_0}
$$

$$
D(s) = S_{TX} S_{RX} + n^2 \frac{D_{TX} S_{RX} + S_{TX} D_{RX}}{1 - jn^2 b_0}
$$

$$
P_t(s) = P_{TXn} \cdot S_{RX}, \quad P_r(s) = P_{RXn} \cdot S_{TX}
$$

### Diplexer Polynomials (Type-II Junction)

The resonating node adds one *pure-real* reflection zero (the paper's
suggested value is 1.5, which determines `c0 = 1/1.5 = 2/3`) and one extra
pole, so `N` and `D` have degree `n_tx + n_rx + 1`. In the monic/n0-factored
normalization (n0 = −1, |n0| = 1):

$$
N(s) = S_{TX} S_{RX}\left(s - \tfrac{1}{c_0}\right)
       + \tfrac{1}{c_0}\left(D_{TX} S_{RX} + S_{TX} D_{RX}\right)
$$

$$
D(s) = S_{TX} S_{RX}\left(s + \tfrac{1}{c_0}\right)
       + \tfrac{1}{c_0}\left(D_{TX} S_{RX} + S_{TX} D_{RX}\right)
$$

$$
P_t(s) = P_{TXn} \cdot S_{RX}, \quad P_r(s) = P_{RXn} \cdot S_{TX}
$$

## Iterative Synthesis Algorithm

### Step 1: Initialization

1. Synthesize TX and RX filters independently (generalized Chebyshev)
2. Get initial S_TX, S_RX from the roots of F⁰_TX and F⁰_RX

### Step 2: Compute Pt, Pr

$$
P_t(s) = \text{conv}(P_{TXn}, S_{RX}), \quad P_r(s) = \text{conv}(P_{RXn}, S_{TX})
$$

### Step 3: Evaluate |p₀r| and |p₀t|

Impose return loss at passband edges (s = ±j):

$$
|N(j)|^2 + |p_{0r}|^2 |P_r(j)|^2 + |p_{0t}|^2 |P_t(j)|^2 = 10^{-RL_{TX}/10}
$$

$$
|N(-j)|^2 + |p_{0r}|^2 |P_r(-j)|^2 + |p_{0t}|^2 |P_t(-j)|^2 = 10^{-RL_{RX}/10}
$$

Solve the 2×2 linear system for |p₀r|² and |p₀t|².
> **Diagnostic (2026-09-29, D2):** Step 3 as written is internally
> inconsistent with the unitary condition below in this step: Step 3
> imposes the constant `|N|² + |p₀|²|P|² = 10^(−RL/10)`, whereas the
> unitary condition is `|n₀|²|N|² + |p₀|²|P|² = |D|²`. The code's
> `solve_scaling_constants` implements the unitary-condition form with
> the antenna RL spec substituted in (`|D|² = |N|²/ε` at each edge) —
> OCR of the paper confirms this is the intended form (eq. 12:
> `|n₀| = 1` for both junctions; eq. 13 reduces to
> `|p₀|²|P|² = |N|²(1−ε)/ε` via the unitary condition). Measured
> (pre common-domain
> rewrite): spectral factorization is sound (with p₀ held fixed the
> iteration is a clean contraction, root displacement 0.55 → 0.02 over
> 37 iters); the 2×2 edge-RL p₀ system is singular for the standard
> order-2 probe (equal order/RL/adjacent bands) — `det → 0`, degenerate
> fallback p₀² ≈ 707; fixing p₀ = 1 removes the divergence but
> converges to a physically wrong diplexer (passband peaks −4.9 / −5.4
> dB vs spec ≈ 0, antenna ≈ −1 dB vs spec ≈ −22). After the common-
> domain rewrite: an asymmetric 2/4 probe diverges (err → 2.89); the
> order-2 gate is still RED (root displacement 4.8e-2 after 30 iters,
> antenna S11 ≈ 0 dB, channel peaks −14 / −19 dB). The paper's
> Section V worked design (see Future Work) pins the remaining
> construction mismatch.


Apply spectral factorization to the unitary condition to get D(s):

$$
D(s) D^*(-s) = |n_0|^2 N(s) N^*(-s) + |p_{0r}|^2 P_r(s) P_r^*(-s) + |p_{0t}|^2 P_t(s) P_t^*(-s)
$$

### Step 4: New S_TX and S_RX estimation

Compute roots of:

$$
\frac{a \cdot N(s) + b \cdot D(s)}{2} = S_{TX}(s) \cdot S_{RX}(s)
$$

where `a = 1 - jn²b₀` and `b = 1 + jn²b₀` (Type-I). For a Type-II junction
the two monic polynomials satisfy `D(s) - N(s) = (2/c0)·S_TX·S_RX`, so the
product is recovered as

$$
S_{TX}(s) \cdot S_{RX}(s) = \tfrac{c_0}{2}\,\bigl(D(s) - N(s)\bigr)
$$

(the degree-`(n_tx+n_rx+1)` leading terms cancel; `D` is normalized exactly
monic before the subtraction, since `spectral_factor` only fixes its leading
coefficient's *phase*).

Sort roots by imaginary part. Assign first np_RX roots to S_RX, last np_TX to S_TX.

### Step 5: Convergence Check

Compare roots of S_TX and S_RX with previous iteration.
Converges in 5-10 iterations typically.

## Filter Recovery

After convergence:

$$
F_{TX}(s) = S_{TX}(s) - D_{TX}(s), \quad E_{TX}(s) = S_{TX}(s) + D_{TX}(s)
$$

$$
F_{RX}(s) = S_{RX}(s) - D_{RX}(s), \quad E_{RX}(s) = S_{RX}(s) + D_{RX}(s)
$$

Where D_TX and D_RX are recovered from D(s) evaluated at the roots of S_TX and S_RX
via polynomial interpolation.

## Implementation in MFS

The implementation adds a `DuplexerSynthesis` struct that:

1. Takes TX and RX filter specifications (order, RL, zeros, bands)
2. Runs the iterative polynomial computation
3. Returns modified TX and RX filter polynomials accounting for mutual loading
4. Synthesizes coupling matrices for each filter independently

This produces a diplexer design where the antenna port reflection is quasi-equiripple
across both passbands, even for contiguous channels.


## Implementation Status in MFS

The MFS library implements a **simplified iterative perturbation** variant of the
Macchiarella-Tamiazzo algorithm:

1. **Independent synthesis**: Both TX and RX filters are synthesized independently
   using the generalized Chebyshev method.
2. **Mutual loading correction**: The S-polynomial roots of each filter are
   iteratively perturbed based on the cross-loading impedance of the other filter.
3. **Frequency-offset aware**: The correction accounts for the spectral separation
   between TX and RX bands by evaluating the loading polynomials at frequency-
   shifted points in the common normalized domain.
4. **Convergent iteration**: Under-relaxation with annealing ensures stability.

### Usage

```rust
use mfs::duplexer::{DuplexerSynthesis, JunctionType};

// Waveguide diplexer: TX 14.9-15.1 GHz, RX 15.15-15.35 GHz
let result = DuplexerSynthesis::new()
    .tx(7, 22.0, 14.9e9, 15.1e9)
    .rx(7, 22.0, 15.15e9, 15.35e9)
    .tx_zeros([-2.5, 2.5])
    .rx_zeros([-2.5, 2.5])
    .junction(JunctionType::TypeI { n: 1.0, b0: 0.0 })
    .synthesize()
    .unwrap();

let tx_response = result.tx.response(14.8e9, 15.2e9, 201).unwrap();
let rx_response = result.rx.response(15.1e9, 15.4e9, 201).unwrap();
```

### Full Spectral Path (as built)

The complete Macchiarella algorithm with spectral factorization is implemented
behind `DuplexerSynthesis::full_iteration(true)` (default `false`). The
infrastructure includes:

- `solve_scaling_constants`: |p₀t|² and |p₀r|² from the return-loss
  conditions (eq. 13)
- `spectral_factor` (`src/spectral.rs`): cepstrum (rustfft) with a
  **Parlett–Reinsch-balanced** faer companion-matrix fallback. The companion
  matrix is balanced by diagonal similarity (eigenvalue-preserving) before the
  solve: a companion whose coefficients span many decades (the §VI Type-II
  R(s) spans ~19) is otherwise so unbalanced that faer's eigenvalues fail the
  backward-residual gate (residual ~0.75); balancing brings it to ~1e-9. This
  replaced the earlier Wilson iteration, whose `wilson_spectral_factor` /
  `fit_polynomial_ls` helpers were removed in the 0.5.0 dead-code cleanup.
- `compute_s_product`: S_TX × S_RX product from N(s) and D(s) (eq. 11 form
  for Type-II: `(c0/2)(D − N)` with D normalized exactly monic)
- `assign_roots`: root assignment by increasing imaginary part (paper Step 4)

**Status (2026-10-03):** the full path **converges to a physically valid
quasi-equiripple 3-port at both junction types**, in the paper's exact bandpass
common-domain variable `Ω(f) = (f/f₀ − f₀/f)/B`:

- **Type-I:** orders 2–8 gated by `full_iteration_order{2,3,4,5,6,8}_quasi_equiripple`
  plus the §V printed-design golden benchmark `full_iteration_order7_section5_golden`
  (matches the paper's printed order-7 N/D/Pt/Pr/p0 to print precision).
- **Type-II** (resonating node): gated by `full_iteration_type2_order4_quasi_equiripple`
  and the §VI GSM-design anchor `full_iteration_type2_section6_anchor`
  (RX 10 / TX 9 poles, transmission zeros, real reflection zero 1.5 ⟹ c0 = 2/3,
  degree-40 R(s)). The anchor reproduces the paper's Fig. 9: |S11| peaks at the
  22 dB spec with peak-to-peak ripple < 1.5 dB, 0 dB channel passbands, and a
  lossless para-identity. This is the highest-degree design the balanced
  spectral factorization has been validated on.

Changing the *default* to `full_iteration(true)` remains a separate, deliberate
behavior change with its own changelog entry. The continuation plan is section
A of `docs/continuation-plan.md`.

### Future Work

- **DONE (2026-10-03): normalization pinned against the paper's worked
  examples.** The Macchiarella–Tamiazzo paper is vendored in the repo
  (`tests/…Diplexers.pdf`); its image-rendered equations were recovered by
  OCR and implemented exactly (eq. (8)–(11)). Eq. (12) states `|n₀| = 1` for
  both junction types; the p₀ return-loss equation (eq. 13) reduces, via the
  unitary condition plus the edge RL spec, to
  `|p₀t|²|Pt|² + |p₀r|²|Pr|² = |N|²(1−ε)/ε` — the same form
  `solve_scaling_constants` implements. Two worked designs are regression
  anchors: the §V Type-I golden (order 7, matches printed N/D/Pt/Pr/p0 to
  print precision) and the §VI Type-II GSM anchor (order 10/9, reproduces the
  Fig. 9 quasi-equiripple |S11|).
- **N=2 unification** — fold the duplexer's full-path iteration loop into
  a shared N-channel module so `MultiplexerSynthesis(N=2)` and
  `DuplexerSynthesis` maintain the Macchiarella–Tamiazzo loop once
  (continuation plan B2; depends on D3).
- Robust spectral factorization (matrix-fraction decomposition) as a
  fallback when cepstrum is ill-conditioned.
- Interval-arithmetic root isolation for high-degree polynomials.
