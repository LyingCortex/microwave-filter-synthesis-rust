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

$$
N(s) = s \cdot S_{TX} S_{RX} - \frac{1}{c_0}(S_{TX} S_{RX} - D_{TX} D_{RX})
$$

$$
D(s) = s \cdot S_{TX} S_{RX} + \frac{1}{c_0}(S_{TX} S_{RX} + D_{TX} S_{RX} + S_{TX} D_{RX})
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

Apply spectral factorization to the unitary condition to get D(s):

$$
D(s) D^*(-s) = |n_0|^2 N(s) N^*(-s) + |p_{0r}|^2 P_r(s) P_r^*(-s) + |p_{0t}|^2 P_t(s) P_t^*(-s)
$$

### Step 4: New S_TX and S_RX estimation

Compute roots of:

$$
\frac{a \cdot N(s) + b \cdot D(s)}{2} = S_{TX}(s) \cdot S_{RX}(s)
$$

where `a = 1 - jn²b₀` and `b = 1 + jn²b₀` (Type-I) or `b = -a = 1` (Type-II).

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

- `compute_n_polynomial`: diplexer reflection numerator N(s)
- `solve_scaling_constants`: |p₀t|² and |p₀r|² from return-loss conditions
- `spectral_factor` (`src/spectral.rs`): cepstrum (rustfft) with a
  balanced-companion-matrix fallback — this replaced the earlier Wilson
  iteration, whose `wilson_spectral_factor` / `fit_polynomial_ls` helpers
  were removed in the 0.5.0 dead-code cleanup
- `compute_s_product`: S_TX × S_RX product from N(s) and D(s)
- `assign_roots`: root assignment via minimum-distance matching

**Status (2026-09-28):** the full path is reachable but **not numerically
converged for order ≥ 4** (typed error; order 2 runs but exhausts its
iterations), and no test or example invokes it. The continuation plan for
converging it is section A of `docs/continuation-plan.md`.

### Future Work

- Convergence of the full spectral path at order ≥ 4 (continuation plan A)
- Robust spectral factorization (matrix-fraction decomposition) as a
  fallback when cepstrum is ill-conditioned
- Interval-arithmetic root isolation for high-degree polynomials
