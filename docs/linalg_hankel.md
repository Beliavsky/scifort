# Hankel matrices and FFT-based multiplication

`src/linalg_hankel.f90` provides the SciPy-compatible `hankel` constructor.
`src/linalg_matmul_hankel.f90` provides an additional `matmul_hankel`
convenience routine for efficient matrix-vector and matrix-matrix products.
**The latter is a SciFort extension, not a function exported by SciPy.**
Both modules are independently written, MIT-licensed Fortran and introduce
no dependencies beyond those already present in SciFort.

## API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_hankel, only : hankel
use scifort_linalg_matmul_hankel, only : matmul_hankel

real(dp) :: c(4), r(5), h(4,5), x(5), y(4)
integer :: status
c = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
r = [4.0_dp,7.0_dp,7.0_dp,8.0_dp,9.0_dp]
x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
call hankel(c,h,status,r)
call matmul_hankel(c,x,y,status,r)
```

* `c` is the first column, of length `n`.
* `r`, if supplied, is the last row, of length `m`. Its first value `r(1)`
  is ignored; the actual entry is always `c(n)`. If omitted, `m=n` and
  the entries to the right of `c(n)` are zero (not conjugated).
* `h` has shape `(n,m)`. `hankel` supports real and complex `c`, real and
  complex `r`, including mixed-type arguments (the output must be complex
  if either input is complex). It initializes `h` to zero on failure.
* `matmul_hankel` supports real and complex `c`, with an optional `r` of
  the same type as `c`; `x` can be a vector `(m)` or an array `(m,k)` of
  multiple right-hand sides, with real or complex values. The output `y`
  has shape `(n)` or `(n,k)` and is complex when required. It is zeroed on
  failure. No matrix is constructed in this routine.
* `check_finite` is available for `matmul_hankel` and defaults to `.true.`;
  the constructor always checks finiteness. As with the other SciFort
  structured-matrix routines, a nonfinite input returns
  `linalg_status_no_convergence` rather than producing a matrix.
* The integer `status` returns `linalg_status_success` (0),
  `linalg_status_invalid_shape`, or `linalg_status_no_convergence`.
* Zero-sized dimensions are supported, with products of the correct shape
  containing zeros.

These routines use the SciFort explicit-output/status style rather than
SciPy's dynamically allocated return values and exceptions. The constructor
currently handles one column/row pair per call (no batched constructor).

## Algorithm

A Hankel matrix has constant *anti*-diagonals:

`H(i,j) = s(i+j-1)` with `s = [c, r(2:m)]`, or `[c, zeros(m-1)]`
when the last row is omitted. With `J` the column-reversal matrix, set
`T = H*J`. Then `T` is Toeplitz with

```
T(i,1) = s(m+i-1)           i = 1, ..., n
T(1,j) = s(m-j+1)           j = 1, ..., m
H*x    = T*reverse(x)
```

`matmul_hankel` constructs only the `n` and `m` Toeplitz edge vectors,
reverses `x` by array slicing, then calls SciFort's existing
`matmul_toeplitz` FFT/circulant-embedding solver. The Toeplitz first row
is **always passed explicitly**; using the Toeplitz default would
incorrectly conjugate complex coefficients. Construction costs `O(n*m)`
and stores the resulting dense matrix; multiplication costs
`O((n+m) log(n+m))` per right-hand side with `O(n+m)` FFT workspace
(excluding input/output buffers). The inherited FFT implementation has
no parallel-worker option.

## References

* [SciPy's `hankel` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hankel.html)
* [SciPy's `matmul_toeplitz` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.matmul_toeplitz.html)
* `docs/linalg_matmul_toeplitz.md` for the FFT implementation.
* [Independent validation report](linalg_hankel_validation_2026-11-09.md).
