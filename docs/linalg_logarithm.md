# Dense matrix logarithm and fractional powers

`scifort_linalg_logarithm` adds small-dense-matrix analogues of
[`scipy.linalg.logm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.logm.html)
and [`scipy.linalg.fractional_matrix_power`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fractional_matrix_power.html).
The implementation is original MIT-licensed Fortran and uses **no new**
third-party numerical code or BLAS/LAPACK dependency.

## Native Fortran API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_logarithm, only : logm, logm_real, &
    fractional_matrix_power, fractional_matrix_power_real
real(dp) :: a(2,2), err, xr(2,2)
complex(dp) :: x(2,2), p(2,2)
integer :: status, roots, terms

! Allocate outputs before calling; no hidden allocation at the boundary.
call logm(a, x, status, residual=err, square_roots=roots, terms=terms)
call fractional_matrix_power(a, 0.3_dp, p, status)
call logm_real(a, xr, status)
call fractional_matrix_power_real(a, 0.3_dp, xr, status)
```

Both generics accept `real(dp)` **or** `complex(dp)` matrix inputs and return
`complex(dp)` outputs. Real input may have a non-real principal logarithm or
power (for example, `log(-1) = i*pi`). The convenience routines ending in
`_real` accept real input and return real output **only** when the principal
result has a negligible imaginary part. A non-real result is signaled by
`linalg_status_no_convergence` rather than silently discarding its imaginary
part. `x` and `p` have shape `(n,n)` and are caller-allocated.

`logm(a,x,status,residual,square_roots,terms,tol,max_sqrt,max_terms)` has optional
controls:

| Argument | Meaning | Default |
|---|---|---|
| `residual` | `||exp(log(A/alpha))-A/alpha||_1 / ||A/alpha||_1`; `alpha=maxval(abs(A))` | output, optional |
| `square_roots` | inverse-scaling steps taken | output, optional |
| `terms` | odd-power atanh terms evaluated | output, optional |
| `tol` | relative atanh-series termination tolerance | `8*epsilon(1.0_dp)` |
| `max_sqrt` | upper limit on repeated square roots | 32 |
| `max_terms` | upper limit on odd-power atanh series | 80 |

Status codes are imported from `scifort_linalg`: `success=0`,
`invalid_shape=1` (also nonfinite input or invalid option),
`no_convergence=4`, `singular=5` (detected numerically). The output matrix is
zero on failure. The backward residual is **not** an estimate of the forward
error in `log(A)`, which can be much greater for nearly defective matrices.

## Method

1. Normalize the matrix by its maximum element magnitude to avoid overflow
   or underflow in the matrix-function kernels. Record the scalar logarithm
   `log(alpha)` separately.
2. Check nonsingularity using an ordinary pivoted complex dense solve.
3. Repeatedly take the principal square root until
   `||B-I||_1 <= 0.25`. This is **inverse scaling and squaring**.
4. Solve `(B+I)*Y = B-I`, then evaluate
   `log(B) = 2*(Y + Y**3/3 + Y**5/5 + ...)`.
5. Multiply the result by `2**square_roots`, restore
   `log(alpha)*I`, and verify that exponentiating the normalized logarithm
   reproduces the normalized input.

`fractional_matrix_power(A,t)` uses `exp(t*logm(A))` for noninteger powers,
with optimized special cases for `t=0`, `t=1`, and `t=1/2` and direct binary
exponentiation for all representable integer exponents. Negative integer
powers use a pivoted linear solve. The half-power normalizes before calling
`scifort_linalg_matrix_functions.sqrtm`, enabling very small nonzero scales.

This deliberately **differs** from SciPy's Schur-Pade algorithms, which are
more robust and efficient for difficult, large, or ill-conditioned problems.
See Al-Mohy and Higham (2012), *SIAM Journal on Scientific Computing* 34,
C152-C169, and Higham and Lin (2011), *SIAM Journal on Matrix Analysis and
Applications* 32, 1056-1078, for authoritative mathematical descriptions.
No implementation source from either paper or SciPy was imported.

## Validity and known limitations

- Matrices must be square, finite, binary64, and of order `1 <= n <= 40`.
  There are no batch dimensions and no BLAS/LAPACK acceleration.
- A principal matrix logarithm is well-defined only if the spectrum avoids
  the **closed negative real axis**; for exactly negative eigenvalues the
  implementation chooses the upper complex branch where feasible (so scalar
  `log(-1)=+i*pi`), but the answer can be branch-sensitive and nonunique.
- `logm` rejects singular or numerically near-singular matrices, instead of
  regularizing them as SciPy sometimes does. Rank-deficient matrices can
  nevertheless have valid positive *integer* powers or a matrix square root;
  those cases are handled separately by `fractional_matrix_power`.
- Other fractional powers of rank-deficient matrices are not yet supported,
  even if a valid fractional power exists (e.g. a singular PSD matrix).
- The inverse-scaling series can return nonconvergence for highly nonnormal
  or poorly conditioned matrices, branch-cut proximity, or insufficient work
  limits. Each matrix square root uses SciFort's QR/Schur implementation, so
  the high-performance SciPy Schur-Pade methods remain preferable for these
  cases. A floating-point exponent giving an unrepresentable result also
  causes a nonconvergence status.
- The inverse-scaling convergence test is numerical rather than a fully
  backward-stable Schur-Pade error bound. The optional reported residual
  checks reconstruction only, not selection of the correct logarithm branch.

## Verification

```bash
python tools/linalg_reference/run_logarithm_checks.py --mode both --cases 800
```

The script compiles focused Fortran tests in checked and `-O3` modes, then
compares randomized real/complex matrices against NumPy/SciPy. Only the Python
**test harness** imports SciPy; `fpm` builds the native code without Python.
See [the validation record](linalg_logarithm_validation_2026-11-07.md) and
[the executable example](../example/linalg_logarithm_example.f90).
