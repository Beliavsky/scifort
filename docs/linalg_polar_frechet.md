# Dense polar decomposition and matrix-exponential Fréchet derivatives

SciFort adds dependency-free, real/complex binary64 versions of two
`scipy.linalg` functions. No SciPy, LAPACK, or BLAS source is copied or required
by these new routines.

## `scifort_linalg_polar`

```fortran
use scifort_linalg_polar, only : polar
call polar(a, u, p, status, side='right', residual=err, orthogonality=unit_err)
```

`a` and `u` are caller-allocated rank-two `real(dp)` or `complex(dp)` arrays
with shape `(m,n)`. `p` uses the *same arithmetic type*, with shape `(n,n)`
for the default `side='right'` or `(m,m)` for `side='left'`.
`status` is an integer defined in `scifort_linalg`.

For `side='right'`, `A = U P`; for `side='left'`, `A = P U`. `P` is symmetric
(real) or Hermitian (complex) positive semidefinite. For `m >= n`, `U^H U=I`;
for `m < n`, `U U^H=I`. Square `U` is unitary. The decomposition is valid for
rank-deficient matrices, including all-zero input. The partial isometry `U`
is not unique when rank-deficient, so do not compare `U` entry by entry against
SciPy for those cases; compare `P`, reconstruction, and unitary constraints.

`polar` uses a scaled one-sided complex Jacobi orthogonalization, with
orthonormal completion of null singular-vector directions. It does **not**
square the condition number by diagonalizing `A^H A`. `max_sweeps` optionally
sets the Jacobi sweep limit (default 120).

Limits: `1 <= min(m,n) <= 64` and `max(m,n) <= 256`. Inputs must be finite;
shapes and `side` must match. Very ill-conditioned cases or an insufficient
sweep limit can return `linalg_status_no_convergence`. The optional outputs
`residual` and `orthogonality` are relative Frobenius diagnostics. For optimal
performance on large matrices, SciPy's LAPACK-backed SVD remains preferable.

## `scifort_linalg_expm_frechet`

```fortran
use scifort_linalg_expm_frechet, only : expm_frechet
call expm_frechet(a, e, exp_a, derivative, status, terms=terms_used, scalings=scalings_used)
```

All four arrays are caller-allocated `(n,n)`, of type `real(dp)` or
`complex(dp)`. The routine computes the matrix exponential and its Fréchet
sensitivity in direction `E`:

\[
\mathrm{d}\exp_A(E)=\left.\frac{\mathrm{d}}{\mathrm{d}t}\exp(A+tE)\right|_{t=0}.
\]

The algorithm differentiates the Taylor recurrence term-by-term on a scaled
matrix and propagates the derivative through each squaring via
`D(X**2)[E] = X*D(X)[E] + D(X)[E]*X`. This correctly handles
noncommuting `A` and `E` and does **not** construct the `2n × 2n` block matrix
or require the specialized SciPy scaling-Padé-squaring kernel.

Optional `tol` (positive) and `max_terms` (positive) control Taylor
truncation; `terms` reports Taylor iterations, and `scalings` the squarings.
Default tolerance is `8*epsilon(1.0_dp)` and default maximum is 120 terms.

Limits: `1 <= n <= 48`, at most 60 scaling steps, and finite inputs.
The algorithm targets small and medium well-scaled matrices; it may return
`linalg_status_no_convergence` for extreme matrix norms, too few terms, or
floating-point overflow. Tolerance controls the Taylor stopping criterion
but is **not** a rigorous forward-error estimate. There is no batching and
no `compute_expm=False` option; both outputs are always produced.

## Status and validation

Both routines use existing SciFort status values:
`linalg_status_success` (0), `linalg_status_invalid_shape` (1), and
`linalg_status_no_convergence` (4). On failure the primary outputs are zeroed.
The original input arrays are not modified.

Build and validate with GNU Fortran and Python SciPy/NumPy:

```bash
python tools/linalg_reference/compare_polar_frechet_scipy.py --mode both --cases 1600
```

The standalone native test is `test/test_linalg_polar_frechet.f90`. An example
is `example/linalg_polar_frechet_example.f90` and a validation record is at
`docs/linalg_polar_frechet_validation_2026-11-06.md`.
