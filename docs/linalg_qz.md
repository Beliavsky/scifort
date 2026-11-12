# Generalized complex Schur (`scifort_linalg_qz.qz`)

`src/linalg_qz.f90` provides the complex-output generalized Schur factorization
of real, complex, and mixed real/complex square matrix pairs. It computes

\[
A=Q\,AA\,Z^H,\qquad B=Q\,BB\,Z^H,
\]

where `Q` and `Z` are unitary, and `AA` and `BB` are complex upper
triangular, with real nonnegative diagonal entries in `BB`. The optional
`alpha` and `beta` arrays contain their diagonals; a generalized eigenvalue
is represented homogeneously as `(alpha(i), beta(i))`. If `beta(i)=0` and
`alpha(i)≠0`, the eigenvalue is infinite.

The interface is based on `scipy.linalg.qz(A, B, output='complex')`,
with a positional output argument convention idiomatic to Fortran:

```fortran
use scifort_kinds, only : dp
use scifort_linalg_qz, only : qz
real(dp) :: a(n,n), b(n,n)
complex(dp) :: aa(n,n), bb(n,n), q(n,n), z(n,n)
complex(dp) :: alpha(n), beta(n)
integer :: status

call qz(a,b,aa,bb,q,z,status,alpha,beta,method='iterative')
```

An optional **trailing** `method` argument is supported by all four
real/complex input overloads:

- `method='iterative'`: inverse-free Hessenberg–triangular reduction followed
  by implicit single-shift two-sided complex QZ iterations. Singular and
  ill-conditioned pencils may fail to converge; the routine returns a
  nonzero status rather than silently returning a doubtful result.
- `method='inverse'`: retain the previous construction using an ordinary
  Schur decomposition of `B**(-1)A`, `A**(-1)B`, or a shifted inverse pencil.
  This can handle some singular pairs but may lose accuracy when the inverse
  matrix is poorly conditioned.
- `method='auto'` (default): first try the new iterative method; on reported
  numerical failure, try the established inverse construction.

All previously supported calls without `method` continue to compile.

## Iterative algorithm

1. Normalize each matrix independently by the maximum magnitude of its
   real/imaginary parts; retain the factors to restore physical units.
2. Compute a QR factorization of `B`, apply its left transformation to `A`,
   and reduce the pair to Hessenberg–triangular form using two-sided
   complex Givens rotations.
3. Iterate on each unreduced block using the complex roots of the
   trailing 2×2 generalized characteristic polynomial as Wilkinson-style
   shifts. Chase the QZ bulge by alternating left and right unitary
   transformations, deflating converged subdiagonal elements of `A`.
4. Rotate phases of the left basis to make `diag(BB)` real and
   nonnegative, then verify both scaled reconstruction residuals before
   returning success.

The iterative method does **not** form `B**(-1)A`, and normally avoids
amplifying errors through the condition number of `B`. The algorithm is a
compact teaching/reference implementation of the complex single-shift QZ
method, **not** a replacement for the extensive scaling, exceptional-shift,
zero-diagonal and convergence safeguards of LAPACK `ZHGEQZ`.

## Limits and status

- Supported order: `0..40`, inherited from the current project-wide dense
  matrix-function limit. The iteration itself is not fundamentally limited
  to 40, but the public contract is unchanged.
- Output: complex Schur factors only. Real quasi-triangular generalized Schur,
  generalized left/right eigenvectors, batching, and condition
  estimates are supplied by other modules or remain future work. Complex
  `ordqz` reordering is implemented in `scifort_linalg_ordqz` (see
  [ordered QZ guide](linalg_ordqz.md)).
- Poorly conditioned, singular, irregular, and nearly defective pencils
  can still cause nonconvergence or large **forward** eigenvalue differences
  despite small reconstruction and homogeneous eigenvalue backward errors.
- The user-supplied matrices are not modified; on failure, outputs should not
  be interpreted as a computed decomposition.
- `linalg_status_success` (`0`): internal reconstruction passed.
- `linalg_status_invalid_shape` (`1`): incompatible sizes, unsupported order,
  or unknown `method` value.
- `linalg_status_no_convergence` (`4`): numerical failure or nonfinite input.
- `linalg_status_singular` (`5`): singular/unsupported irregular pencil.

See [new native tests](../test/test_linalg_qz_two_sided.f90),
[existing native tests](../test/test_linalg_qz.f90),
[example](../example/linalg_qz_two_sided_example.f90),
[reference harness](../tools/linalg_reference/run_qz_checks.py), and
[validation](linalg_qz_two_sided_validation_2026-11-08.md).
