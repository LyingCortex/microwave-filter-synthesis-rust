# Multi-Band Synthesis: Design & Plan

C0 of `docs/continuation-plan.md` section C ("Multi-band — true multi-passband
synthesis, scope-1"). Written the way `docs/multiplexer-synthesis.md` was:
route decision and math **before** any code, so the approximation and synthesis
work items (C2) start from a fixed target.

Status: **route committed (C-route-α), design complete, no code yet.**

---

## Purpose

Give `mfs` a **true multi-passband** synthesis: one filter structure whose
transfer function is equiripple in $K$ passbands ($K \ge 2$) with rejection in
the $K-1$ inter-band gaps and outside the outer band edges — not $K$ independent
single-band designs wired together (that composition path already exists via
`MultiplexerSynthesis` and is a *different* product).

**Scope decision (maintainer, 2026-09-29):** scope-1 confirmed, true
multi-passband synthesis. This document fixes the route, the math, the
specification surface, the transmission-zero placement law, the target
topology, and the first acceptance fixture.

---

## Current State

What already exists (2026-10-05, before any C2 code):

- `src/freq.rs:398` — **`MultiBandMapping`** (`bands: Vec<BandPassMapping>`,
  center-sorted, disjoint). Piecewise forward map
  `map_hz_to_normalized(&self, f) -> Result<(usize, NormalizedSample)>`
  (`src/freq.rs:499`) assigns each frequency to the band minimizing
  $|\Omega_b(f)|$; the inverse `map_normalized_to_hz(band, sample)`
  (`src/freq.rs:545`) is **band-local** and takes the band index — which is why
  the type deliberately does **not** implement the `FrequencyMapping` trait
  (the trait's inverse signature has no slot for the band; doc comment
  `src/freq.rs:392-396`). C1 closed 2026-09-29 with per-band round-trip and
  edge-consistency tests (`src/freq.rs:833-966`).
- The single-band approximation family (`src/approx/`): generalized Chebyshev
  producing a `PolynomialSet` (`order`, `ripple_factor`, `eps`, `eps_r`,
  `transmission_zeros_normalized`, and the `E(s)`/`F(s)`/`P(s)`
  `ComplexPolynomial`s — `src/approx/polynomial.rs:32`).
- `MatrixSynthesisEngine::synthesize(&PolynomialSet) -> Result<CouplingMatrix>`
  (`src/synthesis/engine.rs:69`) — the residue-expansion path that turns a
  `PolynomialSet` into an $N+2$ transversal coupling matrix.
- **E1 is not done.** The recurrence helper `safe_sqrt_term`
  (`src/approx/generalized_ops.rs:96-105`) rejects any zero with $|z| < 1$ with
  `MfsError::NotImplemented("...only supports |zero| >= 1 real zeros")`. The
  spec layer is *permissive*: `validated_transmission_zero`
  (`src/freq.rs:47`) only rejects NaN, and `FilterSpec` stores a passband zero
  ($0.75$) in its own test (`src/spec/types.rs:256`). So a passband TZ flows
  through `FilterSpec`/`validated_transmission_zeros` and fails **later**, at
  `approx`, as an untyped-ish `NotImplemented`.
- No multi-band approximation, no multi-band synthesis, no multi-band API.

---

## Design Goal

A $K$-band bandpass filter: $K$ disjoint passbands inside one structure, with
inter-band transmission zeros providing the rejection.

```
            ┌─ passband 1 ─┐        ┌─ passband 2 ─┐
   ─────────┘   equiripple  └── gap ─┘  equiripple  └─────────
                              ↑
                  inter-band transmission zeros
```

Deliverable (C2/C3): a `MultiBandSpec` + an approximation producing a
`PolynomialSet`-shaped output + reuse of the existing engine to get a coupling
matrix whose response is per-band equiripple and gap-rejecting.

---

## Route Decision: C-route-α (multi-band Chebyshev)

`docs/continuation-plan.md` §C names two candidates. **C0 commits to
route-α.** Rationale, and the honest reason route-β is rejected:

| | C-route-α — multi-band Chebyshev | C-route-β — step-impedance / comb-line ladder |
|---|---|---|
| Core math | Band-Chebyshev factors on a **frequency transformation**, zeros placed in every stopband | Impedance-level extraction from a classical ladder |
| Reuses | `MatrixSynthesisEngine` (output is a `PolynomialSet`) | Nothing in `synthesis` — a genuinely new path |
| Verification surface | Existing residue-expansion + response solver | A new extraction path with its own verification |
| Literature | Macchiarella–Tamiazzo 2005 (frequency-transformation prototype method, "classical multiple-coupled cavity topology") | Deep, but a different algorithmic family |

Route-α is the one consistent with the crate's model (`PolynomialSet` →
`MatrixSynthesisEngine` → `CouplingMatrix` → `ResponseSolver`) and is
**publication-anchored**: G. Macchiarella and S. Tamiazzo, *Design techniques
for dual-passband filters*, IEEE TMTT 53(11), 2005, pp. 3265–3271
(DOI `10.1109/TMTT.2005.855749`). Its abstract states the method directly:
two approaches for passbands separated by a rejection band created by
"suitably placing transmission zeros in the stopband"; both "based on defined
frequency transformations [and] employ well-established prototype techniques;
practical implementation through classical multiple-coupled cavity topology."
That is route-α.

Route-β is **not** rejected as invalid — it is deferred: it needs a new
synthesis engine and a new verification surface, so it earns its own doc if a
later decision calls for it. C0 commits to α to unblock C2/C3 now.

---

## Mathematical Model (route-α)

### 1. Frequency transformation

The single-band design maps its band to the prototype low-pass $\Omega$ via a
`BandPassMapping`. For $K$ bands, the prototype band-Chebyshev factor for band
$b$ is (the analogue of the standard $(1+s^2)$/$\Omega$ factor) evaluated on
the **band-local** transform of `MultiBandMapping`:

$$\Omega_b(f) \;=\; \frac{1}{\mathrm{FBW}_b}\left(\frac{f}{f_{0b}} - \frac{f_{0b}}{f}\right),
\qquad f_{0b}=\sqrt{f_{1b}f_{2b}},\quad \mathrm{FBW}_b=\frac{f_{2b}-f_{1b}}{f_{0b}}.$$

`MultiBandMapping` already computes exactly this and assigns each $f$ to its
nearest band (`src/freq.rs:499`). In the normalized prototype domain the $K$
passbands become the intervals $[-1,+1]$ per band; the $(K-1)$ gaps between
them and the two outer regions are the stopbands.

### 2. Characteristic function

The multi-band Chebyshev characteristic function is the **product of the
per-band band-Chebyshev factors**, with the prescribed transmission zeros
folded in as the numerator of the characteristic polynomial:

$$K(s) \;=\; \frac{F(s)}{P(s)}, \qquad
P(s) \;=\; \prod_{m} \left(s^2 + \Omega_{z,m}^2\right) \cdot \prod_{n}(s^2+\Omega_{z,n}^2)\Big|_{\text{gap+outer}},$$

where the $P(s)$ roots are the transmission zeros (see §Zeros below) and
$F(s)$ is built so that $|S_{11}|^2 = \epsilon^2 K^2/(1+\epsilon^2 K^2)$ is
equiripple in every passband simultaneously. Equivalently, in the
"frequency-transformation" reading of Macchiarella–Tamiazzo: the characteristic
variable is replaced by a function that is $\pm 1$ at every band edge, so the
ripple extrema are shared across all $K$ bands.

> **Implementation note (for C2).** This is the direct analogue of the existing
> recurrence in `generalized_chebyshev_helpers.rs`, whose band factor today is
> the single-band $(1+s^2)$-type term. The extension is to substitute the
> multi-band characteristic variable and to place $F(s)$'s extrema at the
> $2K$ band edges. The `PolynomialSet` output shape
> (`src/approx/polynomial.rs:32`) is unchanged — a new approximation function
> produces the same struct, so `MatrixSynthesisEngine::synthesize` is reused
> verbatim.

### 3. Unitary / Feldtkeller condition

$$|D(s)|^2 = |F(s)|^2 + \epsilon^2 |P(s)|^2,
\qquad E(s) = \text{Hurwitz factor of } D,$$

identical in form to the single-band case (with the multi-band $F$, $P$). The
existing spectral-factorization / `E(s)` normalization path applies.

### 4. Filter recovery → coupling matrix

$(E,F,P)$ is a `PolynomialSet`; `MatrixSynthesisEngine::synthesize` returns the
$N+2$ transversal matrix, exactly as for a single-band finite-TZ design. No new
synthesis machinery in the well-posed case — **provided E1 is solved** (see
§Zeros).

---

## Specification Surface (`MultiBandSpec`, C2)

New public type (in `spec`, alongside `FilterSpec`), *not* a widening of
`FilterSpec`:

```rust
/// One passband of a multi-band design, in the multi-band mapping's domain.
pub struct BandSpec {
    pub center_hz: f64,
    pub bandwidth_hz: f64,
    pub return_loss_db: f64,
}

/// A K-band (K >= 2) multi-passband specification.
pub struct MultiBandSpec {
    /// Disjoint, center-sorted passbands (same invariants as MultiBandMapping).
    pub bands: Vec<BandSpec>,
    /// Out-of-band rejection requirements, per region.
    pub rejection: Vec<RejectionSpec>,
}

/// A rejection requirement: a region and the minimum attenuation it must meet.
pub struct RejectionSpec {
    /// `"outer_lower"`, `"outer_upper"`, or an inter-band gap index.
    pub region: RejectionRegion,
    pub min_attenuation_db: f64,
}
```

Design constraints carried over from `MultiBandMapping`: $\ge 2$ bands, each
validated by `BandPassMapping::new`-equivalent rules, sorted by center, gap
**strictly positive** between adjacent bands (`lo_upper < hi_lower`,
`src/freq.rs:411`). Rejection regions are declared, not derived, so the
placement law below is *specified*, not guessed.

---

## Transmission-Zero Placement Law

This is the part that couples multi-band to **E1**, and the reason C0 must be
explicit.

For $K$ bands there are $K-1$ inter-band gaps and two outer regions. The
characteristic polynomial's zeros (transmission zeros) are placed:

1. **Outer regions (2 zones):** one or more finite zeros in each, or zeros at
   infinity (all-pole) if no outer rejection is required. As in single-band
   designs, a finite outer zero at $\Omega_z$ with $|\Omega_z| > 1$ is
   **stopband** relative to the nearest band prototype — already supported.
2. **Inter-band gaps ($K-1$ zones):** these zeros sit **between** two passband
   prototypes. Relative to the per-band prototype of the *nearer* band they
   fall **inside the passband** ($|\Omega_z| < 1$) — the asymmetry of the
   `BandPassMapping` on the far side of a band's edge means the same physical
   frequency is $|\Omega| > 1$ for one band and $|\Omega| < 1$ for the other.

**Consequence:** route-α *requires* passband transmission zeros. This is
exactly the plan's stated coupling: "the inter-band TZs are *passband* zeros
relative to the per-band prototypes — this is the direct, large-scale
application of E1 (passband TZ handling), which must therefore be designed
together with C, not after it" (`docs/continuation-plan.md:466-470`).

**E1 work (blocks C2):**

- `src/approx/generalized_ops.rs:96-105` — `safe_sqrt_term` must accept
  $|\Omega|<1$ (the term $1 - 1/\Omega^2$ is then *positive* if we write the
  square in the right variable; today the guard rejects it). Relaxing the guard
  is the mechanical part.
- The spec layer should stay permissive and the **approximation layer** should
  produce a typed error for genuinely unsupported cases, not
  `NotImplemented` from a helper. (Plan §E item 1.)
- **Verification:** a passband TZ must yield a real finite-frequency null in
  the response. The first test is $K=2$ with one zero in the single gap.

The $K=2$ case has exactly **one** inter-band gap (one passband zero family)
and two outer regions — the minimal case, and the first acceptance fixture.

---

## Target Matrix Topology

The output of `MatrixSynthesisEngine::synthesize` is the standard $N+2$
transversal matrix (order $N = \sum_b N_b$ across the bands). The
inter-band TZs appear as **trisection/extracted-zero structures** in that
matrix once transformed, which is exactly what the literature fixture's topology
row says: `"trisection_cascade_with_extracted_zero"`. No new topology type is
needed for C2; the existing transversal → folded/arrow/section transforms apply.
The topology decision is a *what the extracted sections are*, not a new
`MatrixTopology` variant.

---

## First Acceptance Fixture (2-band)

**Use the fixture already in the repository** rather than extracting a new PDF:

`tests/cameron2011_fixtures.json`, entry
`"case_id": "Cameron2011_10_5_dual_band_symmetric"`:

- Source: **R. J. Cameron, "Advanced Filter Synthesis," IEEE Microwave
  Magazine, 2011, Figure 25** (DOI `10.1109/MMM.2011.942007`).
- Spec: order **10**, return loss **22 dB**, double termination; lower band
  $[-1.0,-0.35]$, upper band $[0.35,1.0]$ (normalized); **5 transmission
  zeros**: 2 outer (30 dB) + 3 between bands (20 dB, one at $\Omega=0$).
- Topology: `trisection_cascade_with_extracted_zero`, `N+2` format, printed
  $12\times12$ matrix in the fixture.
- Extra: the same file's first entry
  `Cameron2011_4_2_asymmetric_trisection` (Fig. 21, order 4, 2 TZs) is the
  single-band sanity anchor.

This fixture is publication-anchored, in-repo, and carries a printed matrix —
the strongest available acceptance bar. The vendored benchmark standard has no
multi-band family (`docs/continuation-plan.md:447`), so **this fixture defines
the gate** (the fixture tool is `tests/extract_from_pdf.py`, already present).

**Note on band-mapping normalization.** The fixture's bands are given in the
normalized prototype domain ($\pm1$ edges), not in Hz. C3's test constructs a
`MultiBandMapping` whose two bands map to $[-1,-0.35]$ and $[0.35,1.0]$ and
checks the synthesized response against the fixture's $N+2$ matrix.

---

## Work Items (C2, C3 — sequenced after E1)

1. **E1 (prerequisite).** Complex/passband transmission-zero handling in
   `generalized_ops.rs` (relax `safe_sqrt_term`), a typed `approx` error for
   truly unsupported cases, and a passband-null response test.
2. **C2a — `MultiBandSpec` + builder** in `spec` (constraints of §Specification).
3. **C2b — multi-band approximation** (candidate `src/approx/multi_band.rs`)
   producing a `PolynomialSet`; $K=2$ first, $K=3$ as the extension check.
4. **C2c — synthesize** via `MatrixSynthesisEngine::synthesize` (no new engine).
5. **C3 — verification surface:** per-band response equals each band's
   standalone equiripple spec within ripple tolerance; gap rejection asserted;
   an `examples/` rendering of the multi-band response; and the Fig. 25
   fixture as the RED→GREEN anchor.

---

## Verification / Acceptance Criteria

**Status legend:** ✅ met · ⚠️ met structurally but no quantitative value
asserted · ⛔ blocked.

1. ✅ **Route committed** — route-α, with the literature anchor
   (DOI `10.1109/TMTT.2005.855749`) and the rejected-route rationale recorded.
2. ✅ **Spec structure fixed** — `MultiBandSpec`/`BandSpec`/`RejectionSpec`
   above, with `MultiBandMapping`'s invariants.
3. ✅ **TZ placement law stated** — outer zeros stopband, inter-band zeros
   passband; E1 named as a hard prerequisite.
4. ✅ **Topology fixed** — $N+2$ transversal, inter-band TZs as
   trisection/extracted-zero sections; no new `MatrixTopology` variant.
5. ✅ **First fixture identified** — `Cameron2011_10_5_dual_band_symmetric`
   (Fig. 25), in-repo, printed matrix.
6. ⛔ **2-band design reproduces the fixture** — C3, blocked on E1 + C2.

---

## Risks & Mitigations

- **E1 is a hard blocker.** Route-α is not implementable without passband TZs.
  Mitigation: sequence E1 first; C0 deliberately commits to α *because* E1 is
  already on the roadmap as a High-priority cross-cutting enabler (§E), so the
  two are co-designed rather than serialized by accident.
- **Equiripple across $K$ bands is not automatic.** A product of per-band
  Chebyshev factors gives equiripple *per band* only if the band edges are the
  ripple extrema; the $K$-band edge equations must be solved jointly.
  Mitigation: C2b solves the $2K$ edge equations for the ripple parameters; the
  $K=2$ fixture (Fig. 25) is the check.
- **Band-mapping normalization mismatch.** The fixture is in normalized
  coordinates; `MultiBandMapping` is in Hz. Mitigation: C3 builds the mapping
  so its per-band domain matches the fixture's $\pm$ edges, and asserts the
  round-trip first (C1's tests already prove per-band round-trip).
- **Overlap with the multiplexer product.** A $K$-band filter and a $K$-channel
  multiplexer are different products (one structure vs. $K$ channels on a
  junction). Mitigation: `MultiBandSpec` is a *separate* type; no change to
  `ChannelSpec`/`Multiplexer*`.
- **Route-β deferral risk.** If the fixture cannot be matched by α, the
  fallback is documented: route-β gets its own doc. Not a silent pivot.

---

## References

1. G. Macchiarella and S. Tamiazzo, "Design techniques for dual-passband
   filters," IEEE Transactions on Microwave Theory and Techniques, vol. 53,
   no. 11, pp. 3265–3271, Nov. 2005. DOI: `10.1109/TMTT.2005.855749`.
   *(Route-α anchor: frequency-transformation + prototype method, implemented
   via multiple-coupled cavity topology.)*
2. R. J. Cameron, "Advanced Filter Synthesis," IEEE Microwave Magazine, 2011.
   DOI: `10.1109/MMM.2011.942007`. *(§10.5 Fig. 25 = the first acceptance
   fixture; in-repo as `tests/cameron2011_fixtures.json`.)*
3. R. J. Cameron, "General coupling matrix synthesis methods for Chebyshev
   filtering functions," IEEE TMTT, vol. 47, no. 4, pp. 433–442, Apr. 1999.
   DOI: `10.1109/22.754877`. *(The engine this route reuses.)*
4. `docs/continuation-plan.md` section C (this plan, scope-1) and section E
   item 1 (the E1 prerequisite).
5. `docs/cameron-coupling-matrix-algorithms.md` (the transversal → folded /
   arrow / trisection workflow the multi-band matrix flows through).
