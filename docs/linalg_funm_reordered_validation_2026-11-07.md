# Reordered Schur-cluster `funm` validation (2026-11-07)

This validation extends the baseline `funm` record without replacing it.
The implementation is the original MIT SciFort Schur--Parlett/Taylor method,
now with stable adjacent unitary swaps for noncontiguous clusters of repeated
or nearly repeated Schur eigenvalues. The public `funm` API is unchanged.

## Environment

- GNU Fortran 14.2.0 (Debian 13, x86-64), `-std=f2018`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.
- Checked flags: `-O0 -g -fcheck=all -fbacktrace`.
- Optimized flags: `-O3`.

## Reproducible commands

```bash
python tools/linalg_reference/run_funm_checks.py --mode both --cases 550 --max-order 9 --reorder-cases 350
python tools/linalg_reference/run_funm_checks.py --mode both --cases 180 --seed 20261108 --max-order 12 --reorder-cases 250
python tools/linalg_reference/run_matrix_functions_checks.py --mode both --cases 60
python tools/linalg_reference/run_logarithm_checks.py --mode both --cases 100 --max-order 7
python tools/linalg_reference/run_trigonometric_checks.py --mode both --cases 100 --max-order 8
python tools/linalg_reference/run_expm_condition_checks.py --mode both --cases 100 --max-order 7
```

## Results

| Test | Checked | Optimized |
|---|---:|---:|
| Native Fortran `funm` test and example | PASS | PASS |
| General randomized `funm` (seed 20261107) | 550/550 | 550/550 |
| General randomized `funm` (seed 20261108) | 180/180 | 180/180 |
| Analytic stress (two runs) | 164/164 | 164/164 |
| Reordered-cluster comparisons (seed 20261107) | 350/350 | 350/350 |
| Reordered-cluster comparisons (seed 20261108) | 250/250 | 250/250 |

The 1,200 new reordered comparisons use specialized SciPy `expm`, `logm`,
`sinm`, and `A@A` reference results (rather than SciPy's generic `funm`,
which is unreliable on nontrivial defective blocks). Their matrices span
orders 3--40, real/complex types, nonadjacent repeated and nearly repeated
eigenvalues, arbitrary upper triangular couplings, and a limited selection
of orthogonal/unitary similarities. In the reordered-cluster tests:

- Maximum relative Frobenius error for triangular matrices was approximately
  `1.12e-14` in either build, with a strict `2e-8` threshold.
- Maximum relative error for harder dense similarity-transformed defective
  matrices was `5.31e-8` in the checked build and `1.682e-7` in the optimized
  build on seed 20261107; both are below a **separately stated** `3e-7`
  threshold. The second seed reached `1.523e-9` (checked) and `1.026e-9`
  (optimized) on those transformed cases. The looser transformed-case
  threshold explicitly reflects the conditioning of an imperfect Schur QR
  on defective eigenvalues, not the accuracy expected for triangular inputs.

For the general randomized comparisons, the largest relative errors on
seed 20261107 were `2.408e-14` (checked) and `2.130e-14` (optimized).
The second seed reached `2.827e-14` in both builds. The analytic stress
set had 82 cases per run and a maximum relative error `1.148e-14`.

Existing focused `sqrtm`/`signm`, `logm`/fractional-power,
trigonometric/hyperbolic, and `expm_cond` regression runners passed with
checked and optimized builds. No other SciFort numerical module was changed.

## Limitations

- `fpm test` was not run: FPM is not installed in this environment.
- Unreliable Schur eigenvalues on sufficiently defective, nonnormal dense
  matrices can still cause sizable forward errors despite small commutator
  residuals. A successful status is not a forward-error certificate.
- Some clustered spectra are still rejected if the Taylor truncation or
  adjacent swap fails numerical screening. This remains a small-dense
  method, not a full LAPACK-backed blocked Schur--Parlett implementation.
- No third-party source was imported, and runtime dependencies are unchanged.
