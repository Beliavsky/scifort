# Frobenius condition number of the dense matrix exponential

SciFort provides a binary64, dependency-free, real/complex implementation of
`scipy.linalg.expm_cond`:

```fortran
use scifort_linalg_expm_condition, only : expm_cond
real(dp) :: a(3,3), condition
integer :: status, iterations
call expm_cond(a, condition, status, tol=1.0e-9_dp, iterations=iterations)
```

`a` is a square `real(dp)` or `complex(dp)` array with order 1 through 48.
`condition` is a real scalar in either case. The same generic name selects the
matching type overload. The definition is the **relative Frobenius** condition
number of the matrix exponential,

\[
  \kappa_{\mathrm{rel}}(\exp,A)
    = \frac{\lVert A\rVert_F}{\lVert\exp(A)\rVert_F}
      \max_{\lVert E\rVert_F=1}\lVert L_{\exp}(A,E)\rVert_F.
\]

`L_exp(A,E)` is the Fréchet derivative in matrix direction `E`. When `A=0`,
this *relative* condition number is zero. For `A = cI`, it equals `abs(c)`.
For real inputs the norm is over real perturbations; for complex inputs it is
over complex perturbations. Both definitions agree with SciPy for these input
types.

## Algorithm

This original SciFort implementation reuses
`scifort_linalg_expm_frechet.expm_frechet` to apply the Fréchet linear map
`L(E)` and its Frobenius adjoint `L*(E) = L_exp(A^H,E)`. A deterministic,
fully reorthogonalized Lanczos iteration applies the positive-semidefinite
operator `L* L` and estimates its largest eigenvalue using a symmetric
Jacobi eigensolve of a small tridiagonal projection. Its square root is the
induced Frobenius operator norm.

Unlike SciPy's dense `n^2`-by-`n^2` Kronecker/SVD reference algorithm,
SciFort does not construct a Kronecker matrix. Each iteration calls the
existing Fréchet kernel twice. A scalar diagonal shift is applied when
computing these calls, since multiplying both the derivative and the
exponential by the same `exp(c)` leaves the relative condition number
unchanged. For exactly diagonal matrices, SciFort evaluates the induced
norm directly from scalar exponential divided differences; this avoids
iteration and handles large diagonal spectral spreads.

## Control and failure semantics

`tol` (optional, default `1e-9`) controls the relative Ritz-residual
stopping test. `max_iterations` (optional, default 80) limits Lanczos
steps; it is capped at `n*n` and `iterations` (optional output) reports
the number completed. For diagonal or zero inputs `iterations=0`.

`status`, from `scifort_linalg`, is:

- `linalg_status_success` (0) on completed evaluation;
- `linalg_status_invalid_shape` (1) on invalid dimensions, nonfinite
  input, invalid tolerance, or nonpositive iteration cap;
- `linalg_status_no_convergence` (4) if a kernel fails, the Ritz
  tolerance is not met, or overflow/underflow prevents reliable evaluation.

On failure `condition` is zero. The input is not modified.

**Accuracy and limitations:** Lanczos produces a *numerical estimate*,
not a guaranteed upper bound and not an exact SVD. The Ritz residual
criterion does not prove the global maximum eigenvalue was reached.
Its convergence can be slow with nearly multiple dominant singular values;
there is no randomized starting state. Difficult nonnormal matrices,
extreme spreads beyond scalar-shift correction, or derivatives too close
to binary64 underflow can fail. A small iteration limit may return
`linalg_status_no_convergence`. Typical small/moderate random cases agree
with SciPy within far tighter tolerance than `1e-6`, but this is not a
universal error guarantee. There is no batching, no C ABI entrypoint,
and no external LAPACK/BLAS dependency.

## Reproducing independent checks

GNU Fortran and Python with NumPy/SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_expm_condition_checks.py \
    --mode both --cases 600 --max-order 8
```

This compiles and runs the native Fortran suite, streams independent real
and complex randomized matrices through a Fortran driver, compares
relative condition numbers with `scipy.linalg.expm_cond`, and checks
80 additional diagonal/scalar-shift analytic references per build.

See `test/test_linalg_expm_condition.f90`,
`example/linalg_expm_condition_example.f90`, and
`docs/linalg_expm_condition_validation_2026-11-07.md`.

Reference: Awad H. Al-Mohy and Nicholas J. Higham, "Computing the Fréchet
Derivative of the Matrix Exponential, with an Application to Condition
Number Estimation", *SIAM Journal on Matrix Analysis and Applications*
30(4), 1639-1657 (2009), DOI: 10.1137/080716426. No SciPy implementation
code has been imported.
