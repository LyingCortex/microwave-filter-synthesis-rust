# Quickstart Report Output

## Spec

- order: `4`
- return_loss_db: `20`
- normalized_transmission_zeros: `[-1.5]`

## Default Prototype Convention

- `generalized_chebyshev(&spec)` uses the default normalized generalized-Chebyshev prototype domain
- coefficient lists below use ascending powers
- index `0` is the constant term

## Prototype

- order: `4`
- ripple_factor: `0.10050378152592121`
- eps: `1.0524892641310224`
- eps_r: `1`

### Transfer and Reflection Function Polynomials

| i | E(s) | F(s) | P(s) |
| --- | --- | --- | --- |
| 0 | +0.9866 + j1.034 | +0.1068 | +j1.5 |
| 1 | +2.6551 + j1.4411 | +j0.2865 | +1 |
| 2 | +3.2359 + j0.9488 | +0.9635 |  |
| 3 | +2.1318 + j0.382 | +j0.382 |  |
| 4 | +1 | +1 |  |

### Corresponding Singularities

| i | Reflection Zeros (Roots of F(s)) | Transmission Zeros (Prescribed) | Transmission/Reflection Poles (Roots of E(s)) |
| --- | --- | --- | --- |
| 1 | +j0.8983 | -j1.5 | -0.9117 + j0.2645 |
| 2 | +j0.2257 | j∞ | -0.6262 - j0.7733 |
| 3 | -j0.9542 | j∞ | -0.4397 + j1.2545 |
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
- response_check_passed: `true`
- overall_passed: `true`
- notes: `["used the current folded reduction backend", "response invariance check passed on the supplied normalized grid"]`

```text
┌───┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
│   │ S       │ 1       │ 2       │ 3       │ 4       │ L       │
├───┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ S │ +0      │ +1.0324 │ +0      │ +0      │ +0      │ +0      │
│ 1 │ +1.0324 │ -0.0631 │ +0.9089 │ +0      │ +0      │ +0      │
│ 2 │ +0      │ +0.9089 │ -0.1085 │ +0.5659 │ -0.4903 │ +0      │
│ 3 │ +0      │ +0      │ +0.5659 │ +0.6166 │ +0.7653 │ +0      │
│ 4 │ +0      │ +0      │ -0.4903 │ +0.7653 │ -0.0631 │ +1.0324 │
│ L │ +0      │ +0      │ +0      │ +0      │ +1.0324 │ +0      │
└───┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```

## Response Summary

- sample_count: `21`

| sample | frequency_hz | normalized_omega | s11 | s21 | group_delay |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.5 | +0.5 | -0.0199 - j0.0853 | -0.9702 + j0.226 | +2.9934 |
| 10 | 1 | +1 | +0.0967 + j0.0255 | -0.2536 + j0.9621 | +3.4824 |
| 20 | 1.5 | +1.5 | +0.3527 - j0.723 | +0.5338 + j0.2604 | +2.9027 |
