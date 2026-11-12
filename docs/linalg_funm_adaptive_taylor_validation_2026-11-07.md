# Adaptive Taylor evaluation for clustered `scipy.linalg.funm` (2026-11-07)

## Change

This extension to `src/linalg_funm.f90` keeps the existing block
Schur--Parlett evaluator, globally reordered spectral clusters and the
triangular Sylvester solver. It improves the diagonal-block Taylor
calculation when a cluster has *different but close* eigenvalues:

- The default API and derivative callback contract are **unchanged**:
  derivatives are requested only to the cluster order minus one, with the
  historical near-nilpotence screen. Exact Jordan blocks are still supported.
- The opt-in `taylor_max_order=...` (1--128) permits more derivatives for
  each coupled cluster, beginning with the existing Taylor polynomial and
  continuing until three consecutive terms of order at least the block size
  are less than `256*epsilon(1.0_dp)` times the current matrix one-norm, or
  the requested order limit is reached.
- An exact-nilpotence fast path examines the **scaled Schur deviation**.
  A physical-unit Taylor power that merely underflows to zero is not treated
  as proof of nilpotence (important for derivatives that grow for tiny inputs).
- Both evaluation paths check the Schur diagonal of the computed matrix
  function against direct scalar `func` evaluations, independently of the
  existing matrix commutator check.
- A requested order less than a coupled cluster's degree `block_size-1`
  returns `linalg_status_invalid_shape`; reaching the order limit without
  convergence returns `linalg_status_no_convergence`.

For a block `T_B = cI+N` in normalized Schur coordinates, the computed series
is `f(s*c) I + sum_{k>=1} f^(k)(s*c) (s*N)^k / k!`, where `s` is the input
normalization factor. The extension is an **adaptive Taylor series**, not a
Schur--Pade algorithm, not a contour integral and not a certified numerical
error-bound method. Callbacks opting in must provide derivatives of every
requested order on the intended analytic branch. Near a pole or branch
point the series may fail even if the matrix function is mathematically
well defined; ill-conditioned Schur decompositions remain a limitation.

## Reproduction

Environment: GNU Fortran 14.2.0, Python 3.13.5, NumPy 2.3.5,
SciPy 1.17.0. The runner compiles each driver from its local SciFort sources,
with checked `-O0 -g -fcheck=all -fbacktrace` or optimized `-O3` flags.

```bash
python tools/linalg_reference/run_funm_checks.py --mode both --cases 450 --max-order 9 --reorder-cases 280 --diagonalizable-cases 280 --taylor-cases 450
python tools/linalg_reference/run_funm_checks.py --mode both --cases 250 --seed 20261108 --max-order 12 --reorder-cases 200 --diagonalizable-cases 180 --taylor-cases 260
python tools/linalg_reference/run_matrix_functions_checks.py --mode both --cases 60
python tools/linalg_reference/run_logarithm_checks.py --mode both --cases 100 --max-order 7
python tools/linalg_reference/run_trigonometric_checks.py --mode both --cases 100 --max-order 8
python tools/linalg_reference/run_expm_condition_checks.py --mode both --cases 100 --max-order 7
```

## Results

The focused `funm` suites passed **5,284/5,284 comparisons and stress checks**
across the two seeds and two build modes:

| Test group | Checked, two seeds | Optimized, two seeds |
|---|---:|---:|
| General `funm` against SciPy | 700/700 | 700/700 |
| Reordered Schur clusters | 480/480 | 480/480 |
| Semisimple repeated spectra | 460/460 | 460/460 |
| Existing Jordan/scaling stress | 164/164 | 164/164 |
| **Adaptive Taylor clustered spectra** | **838/838** | **838/838** |

The adaptive cases comprise 710 randomized matrices of orders 2--12 plus
128 larger matrices of orders 16, 24, 32 and 40, with real and complex
arithmetic, single and multiple close eigenvalue clusters, nonnormal
upper-triangular blocks, unitary similarities and inter-cluster coupling.
The scalar functions are `exp`, `log`, `sin`, `cos`, `sinh` and squaring;
reference results use SciPy's specialized matrix functions rather than its
generic `funm` for the more difficult near-repeated cases.

Maximum relative Frobenius error in the new adaptive suite was approximately
`3.34e-14` across both modes and seeds. The hardest existing reordered dense
similarity comparisons had errors around `1.43e-9`, well below their separate
`3e-7` threshold; these reflect the numerical sensitivity of the upstream
complex Schur decomposition. The checked and optimized native tests, including
new invalid-parameter and nonconvergence tests, and the runnable example passed.

Additional existing regression suites passed for `sqrtm`/`signm` (60 per
build), `logm` and fractional powers (100 per build), matrix trigonometric
functions (100 per build), and `expm_cond` (100 randomized plus 80 analytic
per build). These do not use the new Taylor option.

The monolithic `fpm test` was **not** run: FPM is unavailable in this
environment. These tests use the documented standalone GNU Fortran builds.

## Provenance

This addition is original MIT-licensed Fortran, tests and documentation
implementing standard analytic Taylor-series identities and a heuristic
stopping criterion. The general Schur--Parlett mathematical framework is
covered by Nicholas J. Higham, *Functions of Matrices: Theory and
Computation*, SIAM, 2008. SciPy is consulted only through its public
matrix-function interfaces for independent numerical comparisons. No new
external source, license or runtime dependency is introduced.
