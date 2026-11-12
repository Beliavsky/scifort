# Hermitian Moore--Penrose pseudoinverse (`scipy.linalg.pinvh`)

`scifort_linalg_pinvh` implements the real symmetric and complex Hermitian
pseudoinverse of a square matrix, with SciPy-style eigenvalue thresholding.
Unlike a general SVD pseudoinverse, `pinvh` reuses SciFort's self-contained
symmetric/Hermitian `eigh` implementation and requires no system LAPACK.

For Hermitian `A = Q diag(w) Q**H`, the pseudoinverse is

```text
pinvh(A) = Q diag(1/w_j if abs(w_j) > cutoff, else 0) Q**H
cutoff   = atol + rtol * max_j abs(w_j)
```

Negative and positive eigenvalues count toward the numerical rank according
to their *absolute* magnitudes; `A` need not be positive semidefinite.

## Fortran usage

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_pinvh, only : pinvh
real(dp) :: a(3,3), ai(3,3)
integer :: status, rank

! Populate symmetric input A.
call pinvh(a, ai, status, rank, atol=0.0_dp, rtol=1.0e-12_dp)
if (status /= linalg_status_success) error stop 'pinvh failed'
```

The generic routine accepts:

- `real(dp) a(n,n), inverse(n,n)`;
- `complex(dp) a(n,n), inverse(n,n)`;
- and corresponding arrays `a(batch,n,n), inverse(batch,n,n)` with one
  *leading* batch dimension, batch >= 1.

The caller allocates the output. The fourth argument `rank` is optional: it
is a scalar integer for an individual matrix and an integer array of length
`batch` for the batched interface. An order-zero **single** matrix is valid
and returns rank zero; a zero-length batch is rejected.

## Options

- `atol` is the nonnegative absolute cutoff, default zero.
- `rtol` is the nonnegative relative cutoff, default `n * epsilon(1.0_dp)`
  for both `real(dp)` and `complex(dp)` inputs.
- `lower` defaults to `.true.`. The eigensolver reads only this triangle;
  the other triangle need not agree with it. For complex inputs, conjugation
  of the selected triangle reconstructs the Hermitian matrix.
- `check_finite` defaults to `.true.`, validating **both** input triangles
  (as in SciPy). With `.false.`, an unused triangle may contain nonfinite
  values, but using nonfinite selected entries remains unsupported.
- `max_iterations` controls the underlying SciFort eigensolver, default
  120. This is a SciFort-specific extension, not a SciPy argument.

No Python-style `return_rank` boolean is needed: supply `rank` to request
it. Only double-precision real and complex are currently supported, and
SciPy's arbitrarily many leading batch dimensions are not implemented.

## Status codes and accuracy

- `linalg_status_success`: computed successfully, including all-zero and
  rank-deficient inputs.
- `linalg_status_invalid_shape`: nonsquare matrix, wrong output or rank
  array shape, zero-sized batch, invalid tolerance, or nonpositive
  `max_iterations`.
- `linalg_status_no_convergence`: nonfinite input or result, or an
  eigensolver numerical failure.

On error, all outputs are zeroed; with a batch, failure of one member
zeros the entire result and rank array. The routines do not modify `a`.

Near a cutoff boundary, changing floating-point rounding, the eigensolver,
or the tolerance may change the returned numerical rank. This is also true
for numerical pseudoinverses generally. Extremely ill-conditioned problems
may warrant increasing `rtol` or comparing with the general SVD-based
`scifort_linalg_svd::pinv`.

An independent reference comparison is in
`reference/check_linalg_pinvh.py`, with the accompanying Fortran driver
`reference/driver_linalg_pinvh.f90`; native tests are in
`test/test_linalg_pinvh.f90` and an example in
`example/linalg_pinvh_example.f90`.

SciPy reference: <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinvh.html>.
