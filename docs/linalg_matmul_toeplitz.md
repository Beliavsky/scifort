# FFT Toeplitz multiplication (`scipy.linalg.matmul_toeplitz`)

Module: `scifort_linalg_matmul_toeplitz` (independently authored, MIT).

## Fortran API

```fortran
use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
call matmul_toeplitz(c, x, y, status, r, check_finite)
```

- `c(:)` is the first column of the Toeplitz matrix (length `n`).
- `r(:)` (optional) is its first row (length `m`). `r(1)` is **ignored**;
  `c(1)` supplies the diagonal. If `r` is absent, it is implicitly
  `conjg(c)` and `m=n`. For real `c`, this is simply `c`.
- `x(:)` or `x(:,:)` is a multiplier of length `m` or shape `(m,nrhs)`.
- `y(:)` or `y(:,:)` is a **caller-allocated** output of length `n` or
  shape `(n,nrhs)`. The routine sets `y` to zero on failure.
- `status` is an integer output: `linalg_status_success` (0),
  `linalg_status_invalid_shape` (1) or `linalg_status_no_convergence` (4).
- `check_finite` (optional) defaults to **`.true.`**, consistent with
  SciFort's other Toeplitz routines. If enabled, all entries of `c`, `r`
  (including the ignored `r(1)`), and `x` must be finite. Nonfinite or
  overflowed FFT results return `linalg_status_no_convergence`.

Overloads cover real `c` and real `x` with real `y`; complex `c` and complex
`x` with complex `y`; real `c` and complex `x` with complex `y`; and
complex `c` and real `x` with complex `y`. Each combination accepts a
rank-one or rank-two `x`. The optional `r` must have the same type as
`c`, as with the existing `scifort_linalg_toeplitz` constructor.

Empty dimensions and zero-column right-hand-side arrays return a
correctly shaped zero product when input/output dimensions agree. Unlike
`scipy.linalg.matmul_toeplitz`, this interface does not expose Python's
broadcasting over arbitrarily shaped input batches, integer dtype
promotion, the `workers` keyword, or `check_finite=False` as a default.
Fortran arrays use one-based indexing and no new FFT dependencies.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
real(dp) :: c(3), r(4), x(4), y(3)
integer :: status
c = [1.0_dp, 2.0_dp, 3.0_dp]
r = [0.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
x = [1.0_dp, 2.0_dp, -1.0_dp, 3.0_dp]
call matmul_toeplitz(c,x,y,status,r)
if (status /= linalg_status_success) error stop 'Toeplitz product failed'
! y = [22.0_dp, 15.0_dp, 18.0_dp]
```

See the [full executable example](../example/linalg_matmul_toeplitz_example.f90).

## Algorithm and performance

For `n>0`, `m>0`, embed the `n`-by-`m` Toeplitz matrix in a circulant
matrix of order `p=n+m-1`. The first column of the circulant embedding is

`g = [ c(1:n), r(m:2:-1) ]`, with `r=conjg(c)` when absent.

For each RHS column `x`, zero-pad to length `p` and calculate

`y = ifft(fft(g) * fft([x, zeros(n-1)]))(1:n)`.

Because the circulant embedding acts by cyclic convolution, the first
`n` entries are **exactly the desired Toeplitz product in exact arithmetic**.
The embedding's Fourier transform and a single FFT plan are reused across
all RHS columns; the full dense Toeplitz matrix is never materialized.

The algorithm requires `O(p log p)` arithmetic to transform the coefficient
vector, then `O(p log p)` per right-hand side, plus `O(p + n*nrhs)` work
storage. SciFort's existing `fft_plan` supports power-of-two, mixed-radix
and Bluestein transforms for arbitrary lengths. The implementation uses
complex FFTs even for all-real inputs; unlike SciPy's real-input optimized
`rfft` path, this leaves further constant-factor speedups possible.

Numerical products have normal floating-point FFT roundoff. Very large
inputs may overflow the intermediate FFT or pointwise spectrum; this is
reported as a numerical failure (the output is zeroed), not silently
accepted as a trustworthy result. Relative error can be large in components
where cancellation makes the exact product close to zero.

## Compatibility and verification

The public behavior of `scipy.linalg.matmul_toeplitz` was used to define
input conventions and output shapes:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.matmul_toeplitz.html

Both direct native Fortran tests and independent SciPy 1.17.0 numerical
comparisons (3,600 cases total, debug and optimized builds) passed.
See [validation details](linalg_matmul_toeplitz_validation_2026-11-09.md).
