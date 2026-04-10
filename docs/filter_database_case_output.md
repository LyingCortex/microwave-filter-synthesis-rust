# Filter Database Case Output

case_id: `Cameron_passband_symmetry_4_2`

## Fixture

- order: `4`
- return_loss_db: `22`
- center_hz: `3500000000`
- bandwidth_hz: `200000000`
- normalized_transmission_zeros: `[1.3217, 1.8082]`

## Synthesis

- approximation_kind: `GeneralizedChebyshev`
- matrix_method: `ResidueExpansion`
- matrix_order: `4`
- response_samples: `21`
- eps: `1.1547462982469874`
- eps_r: `1`

## Coefficient Ordering

- all polynomial coefficient lists in this report use **ascending powers**
- index `0` is the constant term
- the last entry is the highest-order coefficient
- for this Cameron case, monic polynomials therefore end with `+1`

## PolynomialSet

```text
e = [-0.1268 - j2.0658, +2.4873 - j3.6256, +3.6705 - j2.1951, +2.4015 - j0.7592, +1]
f = [+0.0208, -j0.5432, +0.7869, -j0.7592, +1]
p = [-j2.3899, +3.1299, +j1]
```

## Reference Mathematical Model

```text
reference E = [-0.1268 - j2.0658, +2.4874 - j3.6255, +3.6706 - j2.195, +2.4015 - j0.7591, +1]
reference F = [+0.0208, -j0.5432, +0.7869, -j0.7591, +1]
reference P = [-j2.3899, +3.1299, +j1]
reconstructed F from reflection_zeros = [+0.0208, -j0.5432, +0.7869, -j0.7592, +1]
reconstructed E from reflection_poles = [-0.1267 - j2.0659, +2.4875 - j3.6255, +3.6707 - j2.195, +2.4016 - j0.7591, +1]
```

## Actual Generalized Helper Polynomials

```text
actual e_s = [-0.1268 - j2.0658, +2.4873 - j3.6256, +3.6705 - j2.1951, +2.4015 - j0.7592, +1]
actual f_s = [+0.0208, -j0.5432, +0.7869, -j0.7592, +1]
actual p_s = [-j2.3899, +3.1299, +j1]
```

## Coefficient Comparison

### E(s)
| idx | reference | actual | delta |
| --- | --- | --- | --- |
| 0 | -0.1268 - j2.0658 | -0.1268 - j2.0658 | +0 |
| 1 | +2.4874 - j3.6255 | +2.4873 - j3.6256 | -0.0001 - j0.0001 |
| 2 | +3.6706 - j2.195 | +3.6705 - j2.1951 | -0.0001 - j0.0001 |
| 3 | +2.4015 - j0.7591 | +2.4015 - j0.7592 | -j0.0001 |
| 4 | +1 | +1 | +0 |
### E(s) from reflection_poles vs reference E(s)(s)
| idx | reference | actual | delta |
| --- | --- | --- | --- |
| 0 | -0.1268 - j2.0658 | -0.1267 - j2.0659 | +0.0001 - j0.0001 |
| 1 | +2.4874 - j3.6255 | +2.4875 - j3.6255 | +0.0001 |
| 2 | +3.6706 - j2.195 | +3.6707 - j2.195 | +0.0001 |
| 3 | +2.4015 - j0.7591 | +2.4016 - j0.7591 | +0.0001 |
| 4 | +1 | +1 | +0 |
### F(s) from reflection_zeros vs reference F(s)(s)
| idx | reference | actual | delta |
| --- | --- | --- | --- |
| 0 | +0.0208 | +0.0208 | +0 |
| 1 | -j0.5432 | -j0.5432 | +0 |
| 2 | +0.7869 | +0.7869 | +0 |
| 3 | -j0.7591 | -j0.7592 | -j0.0001 |
| 4 | +1 | +1 | +0 |
### F(s)
| idx | reference | actual | delta |
| --- | --- | --- | --- |
| 0 | +0.0208 | +0.0208 | +0 |
| 1 | -j0.5432 | -j0.5432 | +0 |
| 2 | +0.7869 | +0.7869 | +0 |
| 3 | -j0.7591 | -j0.7592 | -j0.0001 |
| 4 | +1 | +1 | +0 |
### P(s)
| idx | reference | actual | delta |
| --- | --- | --- | --- |
| 0 | -j2.3899 | -j2.3899 | +0 |
| 1 | +3.1299 | +3.1299 | +0 |
| 2 | +j1 | +j1 | +0 |

## Root Comparison

```text
reference reflection zeros (roots of F(s)) = [-j0.8593, -j0.0365, +j0.6845, +j0.9705]
actual f_s roots = [+j0.9705, +j0.6845, -j0.8593, -j0.0365]

actual raw e_w roots = [+1.0976 + j0.0977, +0.1267 + j1.1031, -1.4178 - j0.7437, +0.9526 - j0.4571]
actual reflected e_w roots = [+1.0976 + j0.0977, +0.1267 + j1.1031, -1.4178 + j0.7437, +0.9526 + j0.4571]
actual a_s roots = [+j2.8566, +j1.4587]
reference reflection poles (roots of E(s)) = [-1.1031 + j0.1267, -0.7437 - j1.4178, -0.4571 + j0.9526, -0.0977 + j1.0976]
actual e_s roots = [-0.0977 + j1.0976, -1.1031 + j0.1267, -0.7437 - j1.4178, -0.4571 + j0.9526]
```

## Comparison Summary

- E(s): basic match, max |delta| ~= 0.0001
- E(s) reconstructed from reflection_poles: basic match, max |delta| ~= 0.0001
- F(s): basic match, max |delta| ~= 0.0001
- F(s) reconstructed from reflection_zeros: basic match, max |delta| ~= 0.0001
- P(s): basic match, max |delta| ~= 0
- reference E length: `5`
- actual e_s length: `5`
- reference F length: `5`
- actual f_s length: `5`
- reference P length: `3`
- actual p_s length: `3`

## Center Response Sample

```text
frequency_hz = 3500000000
normalized_omega = 0
s11 = +0.0006 - j0.01
s21 = +0.9981 + j0.0613
```
