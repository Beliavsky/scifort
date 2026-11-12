# General nonsymmetric eigenproblems (`scipy.linalg.eig` and `eigvals`)

The MIT-licensed module `scifort_linalg_general_eigen` computes ordinary
real/complex **nonsymmetric** dense-matrix eigenproblems. Unlike the
`scifort_linalg_dense_eigen.eigh` module, it does *not* assume a symmetric or
Hermitian matrix, and its eigenvalues are complex even for real input.

```fortran
use scifort_kinds, only : dp
use scifort_linalg_general_eigen, only : eig, eigvals

complex(dp) :: w(n), vr(n,n), vl(n,n)
integer :: status
call eig(a, w, status, vr=vr, vl=vl) ! left and right eigenvectors
call eig(a, w, status, vr=vr)        ! only right eigenvectors
call eigvals(a, w, status)           ! eigenvalues only
```

`a` can be `real(dp) :: a(n,n)` or `complex(dp) :: a(n,n)`.
The output eigenvalues and eigenvectors are always `complex(dp)`.
The vectors are Euclidean-unit normalized, with the defining equations

\[
A v_j = w_j v_j, \qquad A^H u_j = \overline{w_j} u_j.
\]

The eigenvalue order is **unspecified**, as in SciPy. Eigenvectors of
repeated eigenvalues are not unique. Defective matrices do not have a complete
set of independent eigenvectors, so SciFort does **not** guarantee a basis.
`vr` and `vl` are independent optional arguments, and either can be omitted.

Status constants from `scifort_linalg`:

- `linalg_status_success=0` for success;
- `linalg_status_invalid_shape=1` for non-square input, mismatched output
  dimensions, or matrix order greater than 40;
- `linalg_status_no_convergence=4` for nonfinite input, Schur nonconvergence,
  unresolved eigenvectors, or nonfinite computed eigenvalues.

Zero-order and zero matrices are supported. On failure, output arrays must not
be used as eigenpairs. There is no generalized `b` matrix, homogeneous
representation, batch API or optional balancing. Those are potential future
extensions, not implemented functionality.

## Algorithm and numerical limits

The implementation reuses the pre-existing MIT-licensed
`scifort_linalg_matrix_functions.dense_complex_schur` to compute
`A = Q T Q^H` (after normalization by the largest component of `A`).
For upper triangular inputs, it uses the exact triangular form; for lower
triangular inputs it reverses the matrix order, obtaining an upper triangular
similarity without iterative QR. This avoids avoidable spectral perturbations
in highly nonnormal triangular matrices.

For each eigenvalue `T(j,j)`, triangular substitution recovers right and left
vectors, with protected near-zero pivots and intermediate scaling. Multiplying
by `Q` maps the vectors back to the input matrix. The method checks
triangular eigenpair residuals before declaring success. It has roughly
`O(n^3)` arithmetic and `O(n^2)` working storage.

**Numerical scope:** order `0..40`; unblocked, dependency-free Schur QR, not a
replacement for LAPACK's robust large-matrix nonsymmetric eigensolvers. In
particular, clustered/defective eigenvalues may be ill-conditioned, eigenpair
residuals need not imply accurate individual eigenvalues, and Schur iteration
can return explicit nonconvergence. Small and large input scales are normalized
before decomposition. Treat `status` and eigenpair residuals as essential.

## Testing

Run native GNU Fortran checks and randomized comparisons against local SciPy:

```bash
python tools/linalg_reference/run_general_eigen_checks.py \
  --mode both --cases 520 --max-order 40 --seed 20261108
```

The comparison matches eigenvalues using an optimal assignment because ordering
is arbitrary, then independently tests **left and right** normalized
residuals. It covers real and complex matrices, triangular and diagonal forms,
repeated eigenvalues, nonnormal/zero matrices, and large/small scales. The
Fortran test executable is built with runtime checks and floating-point traps
in the checked mode. See `docs/linalg_general_eigen_validation_2026-11-08.md`.

Mathematical algorithms: Golub and Van Loan, *Matrix Computations*, 4th ed.,
2013, chapters 7–8. SciPy is used only as an independent numerical reference
and for public API conventions, not as imported code:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html
