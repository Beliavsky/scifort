# Complex CSR/CSC sparse matrices and iterative solvers

SciFort supports complex binary64 compressed sparse matrices and three complex
iterative solvers without introducing an external BLAS/LAPACK dependency:

- `scifort_sparse_complex`: `csr_complex_matrix`, `csc_complex_matrix`,
  canonical COO and dense constructors, conversions between CSR/CSC and the
  existing N-D complex COO type, transpose and Hermitian transpose,
  conjugation, matrix-vector products and CSR-times-dense products.
- `scifort_sparse_linalg_complex`: `cg_complex`, `bicgstab_complex`, `gmres_complex`,
  `bicg_complex`, and `qmr_complex`, with CSR/CSC overloads.
- `scifort_sparse_linalg_complex_direct`: complex pivoted `splu_complex` and
  `spsolve_complex` (CSR/CSC), reusable `N`/`T`/`H` sparse-direct solutions;
  see [direct solver API](sparse_complex_direct.md).
- `scifort_sparse_complex_ilu`: no-fill complex `ilu0_complex` factors,
  forward and conjugate-adjoint triangular solves, and reusable factors.

Pointers (`indptr`) and coordinates (`indices`) are **one-based**. Each
CSR row / CSC column has strictly increasing unique indices; converting
unsorted COO sums duplicates and preserves explicitly stored zeros. Stored
matrix values and RHS/initial estimates must have finite real and imaginary
components. `sparse_success=0` signals a successful matrix operation;
invalid matrix input gives `sparse_invalid_input=1`.

## Minimal example

```fortran
program solve_complex
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, &
        csr_complex_from_dense, csr_complex_matvec
    use scifort_sparse_linalg_complex, only : cg_complex
    implicit none
    type(csr_complex_matrix) :: a
    complex(dp), allocatable :: b(:), x(:)
    complex(dp) :: mat(2,2)
    integer :: info
    real(dp) :: residual

    mat(1,:) = [cmplx(4.0_dp,0.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp)]
    mat(2,:) = [cmplx(1.0_dp,-1.0_dp,dp),cmplx(5.0_dp,0.0_dp,dp)]
    a = csr_complex_from_dense(mat,info)
    if (info /= 0) error stop 'invalid matrix'
    b = csr_complex_matvec(a,[cmplx(1.0_dp,2.0_dp,dp), &
        cmplx(-1.0_dp,0.5_dp,dp)])
    call cg_complex(a,b,x,info,rtol=1.0e-12_dp,jacobi=.true., &
        residual_norm=residual)
    if (info /= 0) error stop 'CG did not converge'
    print '(a,es12.4)', 'true residual norm: ', residual
    print *, 'x = ', x
end program solve_complex
```

For a general nonsymmetric system, use `bicgstab_complex`, `gmres_complex`,
`bicg_complex`, or `qmr_complex`.
`gmres_complex` accepts `restart=` (maximum Arnoldi basis size, default
`min(20,n)`) and `maxiter=` (number of restart **cycles**, not Arnoldi
steps); `iterations=` returns the total Arnoldi steps. CG, BiCG, QMR, and BiCGSTAB use
`maxiter=` for individual solver steps. All five accept `x0=`, `rtol=`,
`atol=`, `jacobi=`, `iterations=` and `residual_norm=`. Defaults are
`rtol=1e-5`, `atol=0`, `maxiter=max(1,10*n)`.

The stopping condition is checked against the actual unpreconditioned
Euclidean residual:

\[
\|b-Ax\|_2 \le \max(\mathrm{atol},\mathrm{rtol}\|b\|_2).
\]

Codes: `0` success, positive iteration/cycle limit, `-1` invalid input,
`-2` arithmetic/algorithmic breakdown, `-3` missing/zero Jacobi diagonal.
These codes are published as `sparse_complex_linalg_*` constants.
The solver calls do not change the input matrix or RHS. Solver results
are returned in an allocatable `complex(dp)` vector.

**Restrictions:** CG is for Hermitian positive-definite matrices and does not
prove that its input is HPD. `bicg_complex`, `qmr_complex`, `bicgstab_complex`, and
`gmres_complex` also accept a reusable complex ILU(0) factor through
`ilu=`, mutually exclusive with `jacobi=.true.`. CG intentionally does not
accept unpivoted ILU(0), which need not be Hermitian positive definite.
Matrices are two-dimensional, indices are default Fortran integer and values are `complex(dp)`. CSC solver
overloads convert once to CSR internally, without creating a dense matrix.
The complex solver names are suffixed `_complex` so they can be imported
alongside the original real-valued solvers without generic-name conflicts.

## Validation

```bash
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
```

The native runner incorporates `test/test_sparse_complex.f90`. The independent
randomized reference test can be compiled from
`tools/sparse_reference/complex_sparse_reference_driver.f90` and run with:

```bash
python tools/sparse_reference/compare_sparse_complex_scipy.py \
    /path/to/complex_sparse_reference_driver --per-solver 200
```

The Python comparison requires NumPy and SciPy only for **tests**, not for
the Fortran implementation. It checks both CSR and CSC, COO duplicate/zero
handling, complex matrix-vector products and all three solvers against dense
NumPy solves and `scipy.sparse.linalg`.

## Complex ILU(0), BiCG, and QMR

Use `scifort_sparse_complex_ilu` to build a zero-fill (no row/column
permutation or pivoting) factor on the existing one-based CSR pattern. The
matrix must be square, with a **stored diagonal in every row**. The optional
`diagonal_shift=` adds a nonnegative real shift to the stored diagonal before
factorization; it cannot create a missing diagonal. Explicitly stored zeros
are retained in the pattern. Factorization has no external library dependency.

```fortran
program complex_bicg_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
    use scifort_sparse_complex_ilu, only : complex_ilu0_factor, ilu0_complex
    use scifort_sparse_linalg_complex, only : qmr_complex, bicg_complex
    implicit none
    type(csr_complex_matrix) :: a
    type(complex_ilu0_factor) :: ilu
    complex(dp) :: matrix(2,2), rhs(2)
    complex(dp), allocatable :: x(:)
    integer :: info
    real(dp) :: residual

    matrix(1,:) = [cmplx(5.0_dp,1.0_dp,dp),cmplx(1.0_dp,-2.0_dp,dp)]
    matrix(2,:) = [cmplx(0.5_dp,0.5_dp,dp),cmplx(4.0_dp,0.5_dp,dp)]
    a = csr_complex_from_dense(matrix,info)
    if (info /= 0) error stop 'bad input'
    call ilu0_complex(a,ilu,info)
    if (info /= 0) error stop 'ILU0 failed'
    rhs = [cmplx(1.0_dp,2.0_dp,dp),cmplx(-1.0_dp,1.0_dp,dp)]
    call qmr_complex(a,rhs,x,info,rtol=1.0e-12_dp,ilu=ilu, &
        residual_norm=residual)
    if (info /= 0) error stop 'QMR did not converge'
    print *, 'QMR residual:', residual
    call bicg_complex(a,rhs,x,info,ilu=ilu,rtol=1.0e-12_dp)
    if (info /= 0) error stop 'BiCG did not converge'
end program complex_bicg_example
```

`ilu0_complex_solve(factor,b,x,info,adjoint=...)` solves either
`(LU)x=b` or `(LU)^H x=b` for the approximate ILU factor, not necessarily
for the original matrix. The low-level `ilu0_complex_apply` assumes
that its factor and vector shapes were checked by the caller.

Complex BiCG uses the adjoint `A^H`, not the ordinary transpose `A^T`.
The shadow recurrence and the adjoint preconditioner conjugate complex
coefficients. QMR is a two-sided Lanczos quasi-minimal-residual iteration
applied to `M^{-1} A` when preconditioning is requested; SciPy's public
QMR interface instead exposes optional separate left/right preconditioners
`M1` and `M2`. SciFort currently exposes the one-factor left-preconditioned
variant. Both use the **unpreconditioned** `||b-Ax||_2` for convergence.

Unlike SuperLU's `spilu` (`ILUTP`), SciFort complex ILU(0) does not
reorder, pivot, or create fill. It can report a bad pivot even when `A`
is nonsingular. Use a positive diagonal shift where appropriate.

The new complex ILU factor, BiCG recurrence, and wrappers are original
MIT-licensed SciFort code. Complex QMR adapts SciFort's pre-existing real
QMR recurrence and retains its SciPy BSD-3-Clause provenance; see
`LICENSES/scipy-qmr-BSD-3-Clause.txt`, `THIRD_PARTY_LICENSES.md`, and
`CODE_PROVENANCE.md`.

Run the focused complex-ILU/solver regressions in both compiler modes using:

```bash
python tools/sparse_reference/run_complex_ilu_qmr_checks.py --mode both --per-solver 180
```

For the recorded test coverage and limitations, see
`docs/sparse_complex_ilu_qmr_validation_2026-11-05.md`.

For complex rectangular least squares, see [complex sparse LSQR](sparse_complex_lsqr.md).
