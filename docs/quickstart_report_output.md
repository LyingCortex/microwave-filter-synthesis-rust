# Quickstart Report Output

## Spec

- order: `4`
- return_loss_db: `20`
- normalized_transmission_zeros: `[-1.5]`
- unloaded_q: `2000`

## Default Prototype Convention

- `generalized_chebyshev(&spec)` runs normalized prototype synthesis without a frequency mapping
- `FilterSpec` transmission zeros must already be normalized prototype values
- physical Hz zeros should be converted ahead of time with `normalize_transmission_zeros_hz(...)`
- coefficient lists below use ascending powers
- index `0` is the constant term

## Prototype

- order: `4`
- ripple_factor: `0.10050378152592121`
- eps: `1.0524892641310224`
- eps_r: `1`

\
### Corresponding Singularities

| i | Reflection Zeros (Roots of F(s)) | Transmission Zeros (Prescribed) | Transmission/Reflection Poles (Roots of E(s)) |
| --- | --- | --- | --- |
| 1 | +j0.2257 | -j1.5 | -0.9117 + j0.2645 |
| 2 | -j0.9542 | j∞ | -0.6262 - j0.7733 |
| 3 | +j0.8983 | j∞ | -0.4397 + j1.2545 |
| 4 | -j0.5518 | j∞ | -0.1542 - j1.1276 |

## Canonical Synthesis

- approximation_kind: `GeneralizedChebyshev`
- matrix_method: `ResidueExpansion`
- topology: `Transversal`
- shape: `MatrixShape { rows: 6, cols: 6 }`

```text
┌───┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
│   │ S       │ 1       │ 2       │ 3       │ 4       │ L       │
├───┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ S │ +0      │ +0.3159 │ +0.5549 │ +0.6582 │ +0.4744 │ +0      │
│ 1 │ +0.3159 │ +1.2222 │ +0      │ +0      │ +0      │ -0.3159 │
│ 2 │ +0.5549 │ +0      │ +0.9011 │ +0      │ +0      │ +0.5549 │
│ 3 │ +0.6582 │ +0      │ +0      │ -0.3591 │ +0      │ -0.6582 │
│ 4 │ +0.4744 │ +0      │ +0      │ +0      │ -1.3822 │ +0.4744 │
│ L │ +0      │ -0.3159 │ +0.5549 │ -0.6582 │ +0.4744 │ +0      │
└───┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```

## Transform

- requested_topology: `Folded`
- source_topology: `Transversal`
- result_topology: `Folded`
- pattern_verified: `true`
- response_check_passed: `false`
- overall_passed: `false`
- notes: `["used the current folded reduction backend", "response invariance check failed on the supplied normalized grid"]`

```text
┌───┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
│   │ S       │ 1       │ 2       │ 3       │ 4       │ L       │
├───┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ S │ +0      │ +1.0324 │ +0      │ +0      │ +0      │ +0      │
│ 1 │ +1.0324 │ -0.0631 │ +0.9089 │ +0      │ +0      │ +0      │
│ 2 │ +0      │ +0.9089 │ +0.1085 │ +0.5659 │ -0.4903 │ +0      │
│ 3 │ +0      │ +0      │ +0.5659 │ +0.6166 │ +0.7653 │ +0      │
│ 4 │ +0      │ +0      │ -0.4903 │ +0.7653 │ +0.0631 │ +1.0324 │
│ L │ +0      │ +0      │ +0      │ +0      │ +1.0324 │ +0      │
└───┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```

## Response Summary

- sample_count: `21`

| sample | frequency_hz | normalized_omega | s11 | s21 | group_delay |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.5 | +0.5 | +0.042 - j0.0185 | -0.9029 + j0.4274 | +3.0778 |
| 10 | 1 | +1 | +0.3313 + j0.0045 | +0.058 + j0.9417 | +3.4731 |
| 20 | 1.5 | +1.5 | +0.2617 - j0.8482 | +0.4464 + j0.1132 | +2.5872 |
