# Complex sparse direct LU and solves

`scifort_sparse_linalg_complex_direct` extends SciFort's complex binary64
CSR and CSC support with pivoted sparse LU and reusable direct solves.
These are Fortran counterparts of the key `scipy.sparse.linalg.splu` and
`scipy.sparse.linalg.spsolve` operations. All sparse pointers and indices
remain **one-based**. SciFort's matrix constructors aggregate duplicate
coordinates before factorization and can retain explicitly stored zeros.

## Example

```fortran
program complex_sparse_direct_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
    use scifort_sparse_linalg_complex_direct, only : splu_complex_factor, &
        splu_complex, splu_complex_solve, spsolve_complex
    implicit none
    type(csr_complex_matrix) :: a
    type(splu_complex_factor) :: lu
    complex(dp) :: d(2,2), b(2), target(2)
    complex(dp), allocatable :: x(:)
    integer :: info

    d(1,:) = [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp)]
    d(2,:) = [cmplx(3.0_dp,-1.0_dp,dp),cmplx(4.0_dp,0.0_dp,dp)]
    target = [cmplx(1.0_dp,2.0_dp,dp),cmplx(-1.0_dp,0.5_dp,dp)]
    a = csr_complex_from_dense(d)
    b = matmul(d,target)
    call splu_complex(a,lu,info)
    if (info /= 0) error stop 'LU factorization failed'
    call splu_complex_solve(lu,b,x,info)
    if (info /= 0) error stop 'LU solve failed'
    if (maxval(abs(x-target)) > 1.0e-12_dp) error stop 'inaccurate solution'
    call spsolve_complex(a,b,x,info)
    if (info /= 0) error stop 'spsolve failed'
end program complex_sparse_direct_example
```

## Public API

| Procedure | Meaning |
|---|---|
| `splu_complex(a,factor,info[,pivot_tol])` | Reusable factorization of square canonical CSR or CSC `a` |
| `valid_splu_complex(factor)` | Validate combined factor matrix and row permutation |
| `splu_complex_solve(factor,b,x,info[,trans])` | Solve with vector `b(:)` or dense RHS block `b(:,:)` |
| `spsolve_complex(a,b,x,info[,pivot_tol])` | One-shot factorization and solve; CSR/CSC and vector/block overloads |

`pivot_tol` is an optional **absolute** minimum accepted pivot magnitude;
its default is zero. In `splu_complex_solve`, `trans='N'` solves `A*x=b`
(default), `'T'` solves `transpose(A)*x=b`, and `'H'` solves
`conjg(transpose(A))*x=b`. The factor object caches `lu` and `row_perm`
satisfying `P*A=L*U`, where `P(i,row_perm(i))=1`. Only the lower-factor
multipliers (strict lower triangle) and upper-factor entries are stored in
`lu`; the diagonal of `L` is implicitly one.

`x` is allocatable and has size zero on failure. RHS values and matrix
entries must have finite real and imaginary parts. For `n=0`, factoring and
solving the empty system is valid. RHS blocks with zero columns are valid.

| Status | Constant | Meaning |
|---|---|---|
| `0` | `complex_direct_success` | Success |
| `-1` | `complex_direct_invalid_input` | Invalid shape, data, tolerance, or transpose code |
| `-2` | `complex_direct_breakdown` | Nonfinite intermediate or factor-storage overflow |
| `-4` | `complex_direct_singular` | Pivot of magnitude <= `pivot_tol` |

## Scope and differences from SuperLU

The factorization uses sparse row storage with sorted column indices,
partial **row** pivoting, and dynamic fill-in. Neither the input matrix nor
the combined factor is converted into a dense matrix. Multiple RHS reuse
one factorization. Both CSR and CSC are supported; CSC is converted to CSR
internally. No extra numerical dependencies are required.

Unlike SciPy's SuperLU, this baseline algorithm **does not implement**
COLAMD/AMD column permutations, equilibration, iterative refinement,
supernodes or parallel factorization. It can be much slower and produce
more fill-in on large unsymmetric matrices. It is intended for small and
moderately sized sparse problems, not a replacement for optimized SuperLU
on industrial-scale systems. The pivot tolerance is absolute, rather than
SuperLU's relative diagonal pivot threshold, so the two controls are not
interchangeable.

## Testing

Native test: `test/test_sparse_complex_direct.f90`.

Independent reference tests: `tools/sparse_reference/complex_sparse_direct_driver.f90`
and `tools/sparse_reference/compare_sparse_complex_direct_scipy.py`. To
compile and test without FPM:

```bash
python tools/sparse_reference/run_complex_direct_checks.py --mode both
python tools/sparse_reference/run_native_sparse.py --mode checked
```

The randomized tests compare CSR and CSC one-shot/reusable solves, dense
right-hand sides, and `N`/`T`/`H` against SciPy's `splu` and `spsolve`, and
independently verify the matrix equations with NumPy.
