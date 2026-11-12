# Orthogonal and unitary Procrustes alignment

`scifort_linalg_procrustes` provides a SciPy-style
`orthogonal_procrustes` generic for real and complex binary64 arrays.
Given matrices `a(m,n)` and `b(m,n)`, the routine minimizes

    ||matmul(a,r) - b||_F

subject to `transpose(r)*r = I` (real) or `conjg(transpose(r))*r = I`
(complex), where the caller supplies `r(n,n)`.

## Native Fortran interface

```fortran
use scifort_kinds, only : dp
use scifort_linalg_procrustes, only : orthogonal_procrustes
real(dp) :: a(5,3), b(5,3), r(3,3), scale_factor
integer :: status
! Initialize a, b.
call orthogonal_procrustes(a, b, r, scale_factor, status)
```

Arguments are, in positional order:

| Argument | Intent | Contract |
| --- | --- | --- |
| `a` | `in` | Source matrix, `m`-by-`n` |
| `b` | `in` | Target matrix, the same shape as `a` |
| `r` | `out` | Orthogonal/unitary transformation, `n`-by-`n` |
| `scale_factor` | `out` | Sum of singular values of `A^H B` (real and nonnegative) |
| `status` | `out` | `linalg_status_*` code |
| `max_sweeps` | optional `in` | Positive SVD Jacobi sweep limit; default 120 |

The generic accepts `(real,real)` -> `real` and all three other input
combinations `(complex,complex)`, `(real,complex)` and `(complex,real)` ->
`complex`. There are no implicit allocations of caller-owned output arrays.

On success, the solution is `R = U*Vh`, where `U * diag(s) * Vh` is the
SVD of `A^H B`; `scale_factor = sum(s)`. This is the classical Procrustes
solution of Schoenemann (1966; DOI: 10.1007/BF02289451). The result is
not unique when the cross-covariance is rank-deficient. Unlike spatial
Procrustes analysis, no translation, centroid removal, or shape scaling
is performed.

All finite empty inputs are permitted. If `n=0` the result is the empty
`0`-by-`0` factor. If `m=0`, either input matrix is zero, or both, the
result is the `n`-by-`n` identity, and `scale_factor=0`.

## Error handling and numerical limits

- `linalg_status_success = 0`: valid input and successfully computed result.
- `linalg_status_invalid_shape = 1`: unequal input shapes, wrong output
  shape, or nonpositive `max_sweeps`.
- `linalg_status_no_convergence = 4`: nonfinite input, failure in the Jacobi
  SVD, or unrepresentably large returned scale.

On failure, the output factor and scale are zero. The routine checks all
matrix entries for finiteness; it does not offer SciPy's `check_finite=False`
fast path. It does not currently support prepended batch dimensions.

Each input is normalized by its maximum element magnitude before the
`A^H B` product. For complex input, maximum component magnitude is used,
avoiding `abs(z)` overflow. A binary-exponent product restores the returned
scale while avoiding premature over/underflow. If the *mathematical scale*
is too large for binary64, the routine reports failure rather than returning
infinity. If it is too small, it may underflow to zero.

For ill-conditioned cross-covariances, the result inherits SciFort's
one-sided Jacobi SVD rank cutoff and convergence limits. In rank-deficient
problems, distinct orthogonal/unitary factors can be equally optimal;
compare objective values and unitarity rather than the matrix entries.
Large problems remain slower than LAPACK-backed SciPy.

## Validation

Run the native test (requires GNU Fortran):

```sh
gfortran -std=f2018 -O0 -fcheck=all src/kinds.f90 src/linalg.f90 \
  src/linalg_svd.f90 src/linalg_procrustes.f90 \
  test/test_linalg_procrustes.f90 -o test_procrustes
./test_procrustes
```

Run independent SciPy comparisons:

```sh
python tools/linalg_reference/run_procrustes_checks.py \
  --mode both --cases 500 --max-order 32 --seed 4907
```

See [numerical validation](linalg_procrustes_validation_2026-11-07.md)
and [runnable example](../example/linalg_procrustes_example.f90).
