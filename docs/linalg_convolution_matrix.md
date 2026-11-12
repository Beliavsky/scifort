# SciPy-style convolution matrices and FFT multiplication

The independent MIT-licensed Fortran implementation contains two modules:

- `scifort_linalg_convolution_matrix` implements **`scipy.linalg.convolution_matrix`** for integer, `real(dp)` and `complex(dp)` coefficients; includes one batch dimension.
- `scifort_linalg_matmul_convolution_matrix` adds the SciFort extension **`matmul_convolution_matrix`**, which multiplies by that matrix in `O((m+n) log(m+n))` time per right-hand side using SciFort's existing FFT code, without forming the dense matrix.

## API and dimensions

```fortran
use scifort_kinds, only : dp
use scifort_linalg_convolution_matrix, only : convolution_matrix, convolution_matrix_rows
use scifort_linalg_matmul_convolution_matrix, only : matmul_convolution_matrix

real(dp) :: a(3), matrix(5,5), x(5), y(5)
integer :: status

a = [-1.0_dp, 4.0_dp, -2.0_dp]
x = [1.0_dp, 2.0_dp, 0.0_dp, -3.0_dp, 0.5_dp]
call convolution_matrix(a, matrix, status, mode='same')
call matmul_convolution_matrix(a, x, y, status, mode='same')
! y = [2.0_dp, 6.0_dp, -1.0_dp, -12.5_dp, 8.0_dp]
```

The Fortran API uses **caller-allocated output arrays** and an explicit integer `status`. The output matrix's second dimension `n` is the length of the input to be convolved, and the first dimension must be `k`:

| Mode | `k`, output length | Starting offset in full convolution |
|---|---:|---:|
| `full` (default) | `m+n-1` | `0` |
| `same` | `max(m,n)` | `(min(m,n)-1)/2` (integer division) |
| `valid` | `max(m,n)-min(m,n)+1` | `min(m,n)-1` |

Here `m = size(a)`, and `n = size(matrix,2)` or `size(x,1)` / `size(x)`.
`convolution_matrix_rows(m,n,mode,status)` returns `k` to help allocate
output arrays, or zero for invalid arguments. A supplied `mode` is
case-sensitive and must be `full`, `same`, or `valid`.

The matrix entries satisfy `matrix(i,j) = a(i-j+shift+1)` when that
one-based coefficient index is in `1:m`, and zero otherwise. For all
modes, `matmul(matrix,x)` agrees with `numpy.convolve(a,x,mode)` up to
FFT floating-point roundoff. **The `same` crop when `m>n` uses `min(m,n)`,
not `m`,** matching SciPy's implementation; an overly literal reading of
the basic `same` formula in SciPy's notes would misplace the cropped rows.

The constructor overloads take an `a(:)` and `matrix(:,:)`, or
`a(:,:)` with batch size in the first dimension and
`matrix(:,:,:)` of shape `(batch,k,n)`. Integer coefficients produce an
integer matrix; `real(dp)` produces real and `complex(dp)` complex. No
real/imaginary conjugation is applied to coefficients. The constructor
preserves NaN/Inf values in their corresponding matrix locations, as SciPy
does; unused padded entries remain zero. An empty coefficient vector,
nonpositive `n`, unrecognized `mode`, or wrongly sized output returns
`linalg_status_invalid_shape` (1), with the output initialized to zero.

The FFT product accepts real/complex coefficients and real/complex `x`,
with all four combinations of coefficient and RHS type. `x` may be a
vector `(n)` or matrix `(n,nrhs)`; the output has shape `(k)` or
`(k,nrhs)`. The output type is complex whenever either input is complex.
`check_finite` is optional and defaults to `.true.` (the same convention
as SciFort's Toeplitz FFT multiplication). When enabled, nonfinite inputs
or FFT intermediates cause a numerical status
`linalg_status_no_convergence` (4); bad shapes or unsupported modes
return status 1. All failed products leave the output zero. Real inputs
are internally transformed using the complex FFT backend.

## Algorithm

For every mode, a convolution matrix is a Toeplitz matrix with the
mode-specific row offset. If the selected output has `k` rows, define

```
col(i) = a(i+shift)     for in-range coefficients, else 0
row(j) = a(2-j+shift)   for in-range coefficients, else 0
```

with one-based Fortran array indices. The routine explicitly supplies
**both** `col` and `row` to `matmul_toeplitz`; this is important because
`matmul_toeplitz` otherwise defaults its first row to the *conjugate*
of its first column, which would be wrong for complex convolution.

The existing Toeplitz circulant embedding performs an FFT of length
`k+n-1`, reuses it across right-hand sides, multiplies their spectra,
and inverse-transforms. The dense `k×n` convolution matrix is never
allocated by `matmul_convolution_matrix`. Execution requires
`O((m+n) log(m+n))` work per right-hand side and `O(m+n+k*nrhs)`
working storage. Memory sizes are subject to the available array index
range; oversized FFT embeddings report invalid shape.

## Scope and differences from SciPy

- SciPy accepts any numeric NumPy dtype; this implementation supports
  default Fortran integer, `real(dp)` and `complex(dp)` matrix construction,
  and `real(dp)`/`complex(dp)` multiplication.
- SciPy accepts arbitrarily many leading batch dimensions for the
  constructor; this interface currently supports no batch or one
  leading batch dimension. Multiplication supports multiple RHS columns,
  not batched coefficient vectors.
- `matmul_convolution_matrix` is a **SciFort extension**, not a named
  `scipy.linalg` routine. SciPy expresses the same operation as
  `convolution_matrix(a,n,mode) @ x`, or with `numpy.convolve`.
- Relative errors on extremely ill-conditioned/cancellation-dominated
  FFT products can be larger than those from direct convolution.

See [the executable example](../example/linalg_convolution_matrix_example.f90)
and [the validation report](linalg_convolution_matrix_validation_2026-11-09.md).
The API was checked against [SciPy's convolution_matrix reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.convolution_matrix.html).
