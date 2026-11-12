# Randomized interpolative decomposition

New module: `scifort_linalg_interpolative_random`

SciPy analogue: [`scipy.linalg.interpolative.interp_decomp(A, eps_or_k, rand=True)`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.interp_decomp.html).

The new API returns a **compact interpolative decomposition** without an SVD.
Existing deterministic `interp_decomp` and randomized interpolative SVD
remain unchanged.

## Interfaces

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_interpolative_random, only : randomized_interp_decomp

integer, allocatable :: idx(:)
real(dp), allocatable :: proj(:, :)
integer :: k, status

! Fixed rank; k is input:
k = 3
call randomized_interp_decomp(a,k,idx,proj,status, &
    seed=12345_int64,oversampling=8)

! Adaptive relative Frobenius tolerance; k is output:
call randomized_interp_decomp(a,1.0e-6_dp,k,idx,proj,status, &
    seed=12345_int64,oversampling=8)
```

There are corresponding `complex(dp)` overloads. `idx` and `proj` must be
**allocatable** variables; the routine allocates them. For `A(m,n)`, `idx`
contains a one-based permutation of `1:n`, and `proj` has shape `(k,n-k)`.
The defining relation is

```fortran
A(:,idx(k+1:n)) ≈ matmul(A(:,idx(1:k)), proj)
```

To reconstruct the full matrix, use the existing
`scifort_linalg_interpolative:reconstruct_matrix_from_id` with the skeleton
`A(:,idx(1:k))`, `idx`, and `proj`.

Optional keyword arguments:

- `seed`: `integer(int64)` deterministic local seed; default `20261110`.
  The generator is local, not the Fortran global `random_number` state.
  The same seed does not produce NumPy's randomized pivots.
- `oversampling`: nonnegative integer, default 8. For fixed-rank calls,
  sketch row count is `min(min(m,n),k+oversampling)`; for adaptive calls,
  the initial sketch size depends on this parameter and grows as needed.

## Algorithm

1. Check dimensions, rank/tolerance, and finite matrix entries. Rescale the
   input by its largest entry magnitude (or largest real/imaginary component)
   to avoid overflow and underflow during factorization.
2. Draw a seeded Gaussian row sketch `Y=G*A` (usually much shorter than `A`).
   Use SciFort's column-pivoted QR on `Y` to find candidate original columns.
3. Solve a minimum-norm least-squares system **on the original scaled A**
   to obtain interpolation coefficients. This also supports dependent columns.
4. For adaptive rank, certify the relative **Frobenius norm** of the ID
   residual on the original input. Grow the sketch if the certificate fails.
   At full sketch size a deterministic pivoted-QR fallback is available.

The fixed-rank interface selects exactly `k` columns. `k=0` uses no columns
and an empty `(0,n)` coefficient matrix. Zero arrays and zero-sized dimensions
are supported. Invalid rank/tolerance/oversampling returns
`linalg_status_invalid_shape`; nonfinite entries or failed numerical fitting
returns `linalg_status_no_convergence`.

## Differences from SciPy

- SciPy's adaptive randomized ID targets a **spectral-norm** relative error;
  this implementation checks **Frobenius-norm** residuals, with a small
  floating-point roundoff allowance. Selected ranks can differ from SciPy.
- SciFort returns **one-based** rather than zero-based column indices.
- The local Park–Miller/Box–Muller RNG does not share NumPy's RNG sequence.
- SciPy supports `LinearOperator` input; this implementation covers dense
  real and complex arrays only. No arbitrary batch dimensions are provided.
- The original-matrix coefficient refit improves numerical reliability, but
  consumes additional work. Faster results than deterministic full QR are
  workload-dependent and have not been established by a benchmark.

## Testing

Compile these sources in order (or use FPM with the existing manifest):

```text
src/kinds.f90
src/linalg.f90
src/linalg_svd.f90
src/linalg_qr.f90
src/linalg_interpolative.f90
src/linalg_interpolative_random.f90
```

Link and run `test/test_linalg_interpolative_random.f90`, then compile
`validation/reference_drivers/driver_linalg_interpolative_random.f90` and compare using

```bash
python test/compare_linalg_interpolative_random_scipy.py ./random_id_driver
```

The stand-alone example is
`example/linalg_interpolative_random_example.f90`.
See [validation report](linalg_interpolative_random_validation_2026-11-10.md).
