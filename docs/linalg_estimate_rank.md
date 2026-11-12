# Randomized rank estimation (`scipy.linalg.interpolative.estimate_rank`)

`src/linalg_estimate_rank.f90` adds **original MIT-licensed** SciFort code for
estimating the numerical rank of dense matrices and implicit linear operators.
The implementation is *inspired by the public SciPy API*, not a source-code
translation of its private randomized ID backends.

## Algorithm

For each randomized unit vector `x`, form `y = A*x` and remove its projection
onto the current orthonormal basis of the approximate column space using
twice-repeated modified Gram–Schmidt. If the residual has norm larger than
`eps` times the largest observed `||A*x||`, add it to the basis. Stop after
`oversampling` consecutive rejected probes (default 8), or at full rank.

This approach uses only matrix-vector products and a basis of at most
`min(m,n)` vectors; unlike the previous `estimate_spectral_norm` routines, it
does **not** need an adjoint callback. The matrix-free variant is useful for
sparse matrices and fast transforms where `A*x` can be applied directly.

The rank is **an estimate**. It is not a certified spectral-norm threshold,
and can differ from SVD numerical rank when singular values cluster close
to the requested threshold. In particular, SciPy's `estimate_rank` dense-array
backend typically **overestimates the numerical rank by about 7–8**, and its
`LinearOperator` backend uses a different randomized algorithm. SciFort
returns the number of accepted numerical basis directions, without artificially
adding that padding. Numerical rank results need not agree exactly with SciPy.

## Dense matrices

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_linalg_estimate_rank, only : estimate_rank
real(dp) :: a(60,30)
integer :: numerical_rank, status
! ... initialize a ...
call estimate_rank(a,1.0e-8_dp,numerical_rank,status,seed=123_int64)
```

`estimate_rank` is a generic interface accepting real(dp) or complex(dp)
rank-2 arrays. Its dense path scales by the maximum matrix-entry magnitude
to avoid unnecessary underflow and overflow. Zero-dimensional or all-zero
matrices yield zero rank. The result is always between 0 and min(m,n).

## Matrix-free interfaces

```fortran
use scifort_linalg_estimate_rank, only : estimate_rank_operator
integer :: numerical_rank, status
call estimate_rank_operator(m,n,matvec,1.0e-8_dp,numerical_rank,status)
```

`matvec(x,y,status)` has `real(dp)` arrays with `intent(in) :: x(:)` of
length `n`, `intent(out) :: y(:)` of length `m`, and integer `status` output
(0 for success). For complex(dp) callbacks, call
`estimate_rank_operator_complex`. Separate named entry points are needed
because Fortran generic resolution cannot distinguish only procedure dummy
argument types. An adjoint callback is not required.

The operator path does not scale a matrix it cannot inspect, so callers
must ensure their products remain finite; the returned status reports a
callback failure or nonfinite result. It stores an `m`-by-`min(m,n)` basis,
plus work vectors, but does not store the input matrix.

## Parameters and error handling

- `eps`: `real(dp)` relative residual cutoff; `0 <= eps < 1`.
- `rank`: output integer; initialized to zero on error.
- `status`: `linalg_status_success` (0), `linalg_status_invalid_shape` (1), or
  `linalg_status_no_convergence` (4). The last indicates nonfinite inputs,
  nonfinite products, or a failed matrix-vector callback.
- `seed`: optional `integer(int64)` local Park–Miller pseudorandom seed,
  deterministic and independent of the Fortran global RNG. The same
  integer seed is **not** bit-compatible with NumPy's `Generator`.
- `oversampling`: optional integer, default 8, must be positive; number of
  consecutive candidate vectors that must be rejected to stop. Increasing
  it reduces early termination probability at the cost of more products.

The constructor returns approximate rank rather than eigenvalues or
singular vectors. Matrix-free operation costs at most `min(m,n)+oversampling`
products, with up to `O(m min(m,n)^2)` orthogonalization work and
`O(m min(m,n) + m + n)` storage. Dense products additionally cost `O(mn)`
per probe.

## Reproducibility

From the project root, without FPM or external libraries:

```sh
mkdir -p /tmp/scifort-rank && cd /tmp/scifort-rank
root=/path/to/scifort

gfortran -std=f2018 -O0 -g -fcheck=all -c \
    "$root/src/kinds.f90" "$root/src/linalg.f90" \
    "$root/src/linalg_estimate_rank.f90"
gfortran -std=f2018 -O0 -g -fcheck=all *.o \
    "$root/test/test_linalg_estimate_rank.f90" -o rank_native
./rank_native
gfortran -std=f2018 -O0 -g -fcheck=all *.o \
    "$root/validation/reference_drivers/driver_linalg_estimate_rank.f90" -o rank_driver
python "$root/test/compare_linalg_estimate_rank_scipy.py" ./rank_driver
```

To run the example, link the same three object files with
`example/linalg_estimate_rank_example.f90`.

See [validation report](linalg_estimate_rank_validation_2026-11-09.md).
