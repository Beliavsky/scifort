# Dense `funm` validation (2026-11-07)

Environment: GNU Fortran 14.2.0 (`gfortran -std=f2018`), Python 3.13,
NumPy 2.3.5, SciPy 1.17.0, Linux x86-64.

Command:

```bash
python tools/linalg_reference/run_funm_checks.py --mode both --cases 550 --max-order 9
```

Results:

| Configuration | Native test | Runnable example | Randomized reference | Analytic stress |
|---|---|---|---|---|
| `-O0 -g -fcheck=all -fbacktrace` | PASS | PASS | 550/550 (max rel. error 2.41e-14) | 81/81 (max rel. error 1.15e-14) |
| `-O3` | PASS | PASS | 550/550 (max rel. error 2.13e-14) | 81/81 (max rel. error 1.15e-14) |

Additional checked-build randomized comparisons:

- 400/400, orders 1--7, seed 20261107: max relative Frobenius error 1.74e-14.
- 650/650, orders 1--10, seed 20261108: max relative Frobenius error 2.13e-14.
- 130/130, orders 1--20, seed 20261109: max relative Frobenius error 3.68e-14.
- 100/100, orders 1--40, seed 20261110: max relative Frobenius error 5.15e-14.

The randomized reference tests cover nine analytic functions (`exp`, `log`,
`sin`, `cos`, `sinh`, `sqrt`, `z*z`, `(1+z)**(-1)` and `exp(-z*z)`),
real and complex inputs, diagonal, triangular, orthogonal/unitary
similarities, repeated semisimple eigenvalues, and Jordan blocks with
analytic derivatives. They compare against SciPy `funm` for well-separated
cases and specialized SciPy matrix functions where SciPy's generic `funm`
has known repeated-eigenvalue limitations.

The 81 additional stress cases include repeated Jordan blocks up to order 40,
complex and real matrices, two different contiguous Jordan blocks, failure
without required derivatives, explicitly unsupported noncontiguous repeated
blocks, scaled logarithms for input magnitude 1e-200 through 1e200, expected
rejection of overflowing exponential outputs, and nonnormal triangular
matrix exponentials up to order 40.

The largest reported relative errors use the Frobenius norm normalized by
`max(1, ||reference||_F)`. Passing a relative check does not certify
conditioning on arbitrary input. Some difficult inputs are **intentionally
rejected**, e.g. noncontiguous repeated eigenvalues with coupled Schur blocks.

The full `fpm test` suite was **not** run in this environment because the
Fortran Package Manager is not installed; all listed focused GNU Fortran
builds and tests completed. No other SciFort modules were modified aside
from making the pre-existing `dense_complex_schur` helper public.
