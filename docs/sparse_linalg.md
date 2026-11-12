# Sparse iterative linear solvers (`scipy.sparse.linalg` subset)

`scifort_sparse_linalg` provides real-binary64 sparse linear algebra over
SciFort's **one-based canonical** `scifort_sparse.csr_matrix`. The core does
not require BLAS, LAPACK, or an external solver library.

## Matrix-function Krylov action (`funm_multiply_krylov`)

`funm_multiply_krylov` approximates `y = f(t*A)*b` without assembling the
full dense matrix function `f(A)`. The public generic has two signatures:

- `funm_multiply_krylov(a,b,f,y,info,...)` for a **real, one-based CSR** `a`;
- `funm_multiply_krylov(n,b,matvec,f,y,info,...)` for a matrix-free square
  linear operator supplied through a user-defined procedure.

The user supplies a **dense matrix-function callback** `f(h,value,status)`:
`h(:,:)` contains a small real projected matrix, `value(:,:)` returns its
full matrix function, and `status=0` signals success. The second overload
also accepts `matvec(x,y)` for applying the large operator. Both interfaces
use `real(dp)` arrays, and all dummy arguments have explicit `INTENT`.

```fortran
program sparse_polynomial_function_example
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : funm_multiply_krylov
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: y(:)
    integer :: info, cycles, products

    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]))
    call funm_multiply_krylov(a,[1.0_dp,2.0_dp],cubic_minus_linear,y,info, &
        cycles=cycles,products=products)
    if (info /= 0) error stop 'matrix-function action failed'
    print *, y   ! f(diag(1,2))*[1,2] = [0,12]
contains
    subroutine cubic_minus_linear(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Projected real square matrix.
        real(dp), intent(out) :: fh(:,:) !! h**3-h as a dense matrix.
        integer, intent(out) :: status !! 0 on success.
        fh = matmul(h,matmul(h,h))-h
        status = 0
    end subroutine cubic_minus_linear
end program sparse_polynomial_function_example
```

All optional arguments follow `info` in this order:
`t, rtol, atol, restart_every_m, max_restarts, assume_hermitian, cycles,
products, update_norm`.

- `t` defaults to `1.0_dp`. `rtol` defaults to `1.0e-6_dp`, and `atol`
  defaults to `0.0_dp`. The stopping threshold is
  `max(atol, rtol*norm2(b))`.
- `restart_every_m` defaults to `min(20,n)` and is clamped to `n` for
  oversized requests; it must be positive when provided.
  `max_restarts` defaults to `20` and counts complete Krylov cycles.
- `assume_hermitian=.true.` selects the cheaper **Lanczos** recurrence.
  On real data this requires a symmetric operator; symmetry is *not*
  verified. The default is a two-pass modified-Gram-Schmidt **Arnoldi**
  basis suitable for nonsymmetric real matrices.
- `info=0` means the Krylov update passed the tolerance test or the
  invariant Krylov subspace closed (happy breakdown). `info>0` reports
  exhausted restart cycles. `info=-1` means invalid input; `info=-2`
  indicates an invalid callback result or numerical breakdown.
  A positive return still includes the best computed approximation.
- `cycles`, `products`, and `update_norm` are optional outputs. Unlike a
  **linear-system residual**, the difference between the most recent
  approximants is *not a certified error bound* for an arbitrary `f`.
  Always verify difficult matrix functions against independent references.
- The zero right-hand side returns zeros immediately, without calling
  the callback. An empty 0-by-0 operator returns an empty result.

Only real-to-real matrix functions are currently supported; complex matrices
and functions requiring complex output are future work. The *input matrix*
remains sparse/matrix-free, but the small projected `H` and `f(H)` are dense.
The stored projected block-Hessenberg matrix grows with the number of cycles.
SciFort honors a user-specified restart budget, whereas SciPy 1.17 sometimes
limits the number of cycles to `floor(n / m) + 1`, so their stopping behavior
can differ. Both use the same restarted Krylov mathematical construction.

See [Krylov matrix-function validation](sparse_funm_krylov_validation_2026-11-05.md)
and [SciPy's API](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.sparse.linalg.funm_multiply_krylov.html).

## N-dimensional grid Laplacians and sparse structure

The separate [`scifort_sparse_laplacian_nd`](sparse_laplacian_nd.md) module
implements `scipy.sparse.linalg.LaplacianNd`: matrix-free matvec/matmat,
canonical CSR construction, and analytic eigenvalues/eigenvectors for
Dirichlet, Neumann and periodic boundaries.

This `scifort_sparse_linalg` module additionally exports
`spbandwidth(a,below,above,info)` and
`is_sptriangular(a,lower,upper,info)` for CSR and CSC inputs. They inspect
stored structure, including explicit zeros, and support rectangular matrices.
See [validation results](sparse_laplacian_nd_validation_2026-11-05.md).

## Augmented/recycled Krylov solvers (`lgmres`, `gcrotmk`)

Two solvers for real, nonsymmetric CSR linear systems are available:

- `lgmres` augments a restarted GMRES basis with up to `outer_k` previous
  **full-cycle corrections** (Baker, Jessup, and Manteuffel, 2005).
- `gcrotmk` maintains an orthonormal image basis `C=A*U` and solves each
  Krylov correction in the space projected orthogonally away from `C`,
  retaining up to `k` recycled directions (de Sturler, 1999; Hicken and
  Zingg, 2010).

Both use right Jacobi or right ILU(0) preconditioning, with a true-residual
criterion `norm2(b-A*x) <= max(atol,rtol*norm2(b))`. They work without BLAS,
LAPACK or other new dependencies. The `outer_v` and `recycle_u` arguments
are **allocatable, optional, `intent(inout)` real(dp) matrices** holding
`n`-component vector columns. They may be unallocated on first use, are
filled with retained vectors on return, and may be passed to a later solve
with the same dimension. GCROT always recomputes `A*U`, so the matrix may
change between calls.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix
use scifort_sparse_linalg, only : lgmres, gcrotmk
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: b(:), x(:), outer(:,:), recycled(:,:)
integer :: info, steps
real(dp) :: err
! Populate canonical a and b before calling these routines.
call lgmres(a,b,x,info,rtol=1.0e-10_dp,inner_m=12,outer_k=3, &
    outer_v=outer,iterations=steps,residual_norm=err)
call gcrotmk(a,b,x,info,rtol=1.0e-10_dp,m=12,k=3, &
    recycle_u=recycled,iterations=steps,residual_norm=err)
```

The leading argument order for `lgmres` is
`(a,b,x,info,rtol,atol,maxiter,inner_m,outer_k,x0,jacobi,iterations,residual_norm,ilu,outer_v)`.
For `gcrotmk`, replace `inner_m,outer_k` with `m,k` and `outer_v` with
`recycle_u`. If present, `ilu` is an existing `ilu0_factor`, mutually
exclusive with `jacobi=.true.`.

- `inner_m` defaults to `min(30,n)`, `outer_k` to `min(3,n-1)`;
  `m` defaults to `min(20,n)`, `k` to `min(m,n-1)`. Zero recycling is
  permitted; an inner width of zero is invalid for nonempty systems.
- `maxiter` counts **outer cycles**, defaults to `max(1,10*n)` and may be
  zero. `iterations` reports all Arnoldi columns used, including LGMRES
  augmentation columns, but not matvecs needed to reimage recycled GCROT
  vectors. A successful call has `info=0`; positive means the cycle limit
  was hit; negative reports an invalid input, bad Jacobi diagonal or
  numerical breakdown.
- For either solver, previous vectors may be supplied with `outer_v` /
  `recycle_u`, which must have shape `(n,p)` and contain only finite values.
  They are kept at most `outer_k` or `k` columns, respectively. Passing
  the same arrays across related right-hand sides may accelerate solves.
- LGMRES uses **right** preconditioning and recomputes matrix images for
  cached corrections. Its retained history is a matrix rather than SciPy's
  list of `(v,Av)` tuples; the API does not yet support SciPy's optional
  `prepend_outer_v`, `store_outer_Av` or arbitrary preconditioner callbacks.
- GCROT currently implements FIFO (`truncate='oldest'`) history truncation.
  SciPy's `truncate='smallest'`, complex arithmetic, callbacks and general
  linear-operator input are not supported in this first implementation.

See [LGMRES and GCROT(m,k) validation](sparse_recycling_validation_2026-11-05.md)
for checked/optimized native regressions and independent SciPy comparisons.

## Quasi-Minimal Residual (QMR) and Transpose-Free QMR (TFQMR)

For nonsymmetric real linear systems, `qmr` and `tfqmr` provide two
quasi-minimal-residual Krylov alternatives to BiCGSTAB and GMRES. Both use
canonical, one-based `csr_matrix` storage, report the actual Euclidean
residual `norm2(b - A*x)`, and do not form a dense operator.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : qmr, tfqmr
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: x(:)
real(dp) :: mat(3,3), b(3), err
integer :: info, steps

mat = reshape([4.0_dp,1.0_dp,0.0_dp, &
               2.0_dp,5.0_dp,1.0_dp, &
               0.0_dp,2.0_dp,6.0_dp], [3,3])
a = csr_from_dense(mat)
b = [1.0_dp,2.0_dp,3.0_dp]
call qmr(a,b,x,info,rtol=1.0e-10_dp,jacobi=.true., &
    iterations=steps,residual_norm=err)
call tfqmr(a,b,x,info,rtol=1.0e-10_dp, &
    iterations=steps,residual_norm=err)
```

The common leading argument order is
`(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm)`.
`tfqmr` additionally accepts a final optional `ilu` of type `ilu0_factor`
(also returned by `spilu`), mutually exclusive with `jacobi=.true.`.

- The default convergence condition is the **true original-system residual**
  `norm2(b-A*x) <= max(1e-5*norm2(b),0)`; `rtol`, `atol` are separately
  configurable, nonnegative and finite. The result is independently checked
  on every completed iteration or half-iteration.
- `x0` defaults to zero; an already-converged guess or empty system returns
  without iterations. `maxiter` defaults to `max(1,10*n)`, accepts zero,
  and counts full QMR iterations for `qmr` but **half-iterations** for `tfqmr`.
- `info=0` denotes true-residual convergence; positive means iteration limit
  reached; negative denotes invalid input, a failed Jacobi diagonal, or
  numerical breakdown. `iterations` and `residual_norm` are optional.
- QMR needs both `A*v` and `A**T*v` and builds a canonical sparse transpose
  once. Right Jacobi preconditioning is supported; arbitrary transpose
  preconditioner callbacks or ILU-transpose solves are not implemented.
- TFQMR needs only forward `A*v`, with optional right Jacobi/ILU(0). It uses
  the mathematically consistent composed operator `A*M` and updates the
  approximation with `M*d`, unlike SciPy 1.17's preconditioned `M*A`
  recurrence. SciFort tests actual residual convergence, rather than
  accepting TFQMR's inexpensive residual bound alone.
- As with other short-recurrence nonsymmetric Krylov methods, ill-conditioned
  or defective systems can stagnate or break down; use GMRES when necessary.
  General complex matrices and arbitrary `LinearOperator` inputs are not yet
  supported.

Both recurrences retain SciPy 1.17 BSD-3-Clause attribution. QMR adapts
SciPy `qmr`; TFQMR adapts SciFort's existing BSD-3-Clause matrix-free TFQMR
kernel, which credits SciPy `tfqmr`, based on Freund (1993). Full notices
appear in `LICENSES/scipy-qmr-BSD-3-Clause.txt` and
`LICENSES/scipy-tfqmr-BSD-3-Clause.txt`.
See the [numerical comparison record](sparse_qmr_tfqmr_validation_2026-11-05.md).

## Biconjugate Gradient (BiCG) and Conjugate Gradient Squared (CGS)

The original SciFort implementations of `bicg` and `cgs` follow the established
real nonsymmetric Krylov recurrences described in Barrett et al., *Templates
for the Solution of Linear Systems* (SIAM, 1994). Both accept canonical
square real CSR matrices and validate input shapes, finite values, tolerances,
preconditioners and iteration limits.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : bicg, cgs
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: x(:)
real(dp) :: b(3), residual
integer :: info, steps

a = csr_from_dense(reshape([4.0_dp,1.0_dp,0.0_dp, &
                             2.0_dp,5.0_dp,1.0_dp, &
                             0.0_dp,2.0_dp,6.0_dp],[3,3]))
b = [1.0_dp,2.0_dp,3.0_dp]
call bicg(a,b,x,info,rtol=1.0e-10_dp,jacobi=.true., &
    iterations=steps,residual_norm=residual)
call cgs(a,b,x,info,rtol=1.0e-10_dp, &
    iterations=steps,residual_norm=residual)
```

Both routines have the signature
`(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm)`,
except that `cgs` also takes trailing optional `ilu` of type `ilu0_factor`.

- `rtol=1e-5` and `atol=0` default to the SciPy convergence condition
  `norm2(b - A*x) <= max(rtol*norm2(b), atol)`; each attempted iteration
  independently checks the **true**, unpreconditioned residual.
- The initial guess `x0` is optional (default all zero). `maxiter` defaults to
  `max(1, 10*n)`; `maxiter=0` is supported (returns positive nonconvergence
  unless the initial guess already satisfies the residual target).
- `jacobi=.true.` enables diagonal-inverse preconditioning. BiCG applies the
  same real diagonal preconditioner to forward and transposed systems;
  CGS also optionally accepts an `ilu0_factor` (mutually exclusive with Jacobi).
- Output `info=0` means actual-residual convergence, positive means the cap
  was reached, and negative indicates invalid input, bad Jacobi diagonal,
  or algebraic breakdown. `iterations` and `residual_norm` are optional output
  diagnostics. Empty systems and an already-exact initial guess terminate
  successfully without iterations.

BiCG requires both `A` and `A**T` matrix-vector products and stores a CSR
transpose. CGS avoids that transpose, but can show erratic residual behavior,
including numerical breakdown or stagnation; BiCGSTAB/GMRES may be preferable
on difficult systems. Unlike SciPy, these routines currently accept only real
CSR matrices, not arbitrary `LinearOperator` callbacks; BiCG does not accept
ILU transpose solves. General complex systems are unsupported.

The independent randomized test driver checks both routines against SciPy
1.17 on nonsymmetric, diagonally dominant systems, with and without Jacobi.
See [validation details](sparse_bicg_kronsum_validation_2026-11-05.md).

## Sparse matrix and axis norms

`norm(a,ord,info)` accepts canonical real `csr_matrix` or `csc_matrix`.
Omitting `ord` or using the string `ord='fro'` returns the Frobenius norm.
Real numeric orders may be `1`, `-1`, positive or negative infinity, or `2`;
integer `1`, `-1`, and `2` are accepted as separate Fortran overloads.
The supported scalar matrix norms are:

| `ord` | Mathematical definition | Algorithm |
|---|---|---|
| omitted / `'fro'` | sqrt(sum(abs(A)**2)) | O(nnz), scaled sum of squares |
| `1` / `-1` | largest / smallest absolute column sum | O(nnz+n_cols) |
| `+inf` / `-inf` | largest / smallest absolute row sum | O(nnz+n_rows) |
| `2` | largest singular value | Existing `svds` Golub--Kahan iteration; approximate |

Ord `2` may return nonconvergence instead of an inaccurate value; a zero
or undersized matrix may be invalid for this option. Ord `-2`, nuclear norm
and other unsupported matrix orders return a nonzero status and NaN. A
missing or zero-size reduction dimension returns invalid status for the
maximum/minimum norms; `norm(a)` of an empty matrix is zero. Fortran uses
`ieee_value(0.0_dp,ieee_positive_inf)` to specify infinity.

`norm_axis(a,axis,ord,info)` returns **vector** norms along a one-based
axis: `axis=1` reduces rows and returns one result per column, while
`axis=2` reduces columns and returns one result per row. Default order is 2;
integer and real `ord` are supported, with zero counting mathematically
nonzero entries, positive/negative orders and ±infinity. CSC operations are
adapted through canonical CSR. Nonfinite results and invalid inputs return
an empty vector and a negative status.

```fortran
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : norm, norm_axis
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: row_norms(:)
real(dp) :: norm_one, norm_infinity
a = csr_from_dense(reshape([1.0_dp,3.0_dp,0.0_dp,4.0_dp],[2,2]))
norm_one = norm(a,1)
norm_infinity = norm(a,ieee_value(0.0_dp,ieee_positive_inf))
row_norms = norm_axis(a,2,ord=2)
```

The implementation is original SciFort source (no external code copied).
See the [independent SciPy comparison driver](../tools/sparse_reference/compare_sparse_formats_norm_scipy.py).

## Sparse matrix powers and inverse

`matrix_power` and `inv` are real-binary64 numerical operations compatible
with the corresponding `scipy.sparse.linalg` functionality, adapted to
SciFort's canonical one-based CSR and CSC representations. These functions
preserve the **input storage type** (`csr_matrix` or `csc_matrix`):

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense, csr_to_dense
use scifort_sparse_linalg, only : matrix_power, inv
implicit none
type(csr_matrix) :: a, a5, inverse
integer :: info

a = csr_from_dense(reshape([2.0_dp, 0.0_dp, &
    1.0_dp, 3.0_dp], [2,2]))
a5 = matrix_power(a, 5, info)
if (info /= 0) error stop 'power failed'
inverse = inv(a, info)
if (info /= 0) error stop 'inverse failed'
```

- `matrix_power(a, power, info)` requires an integer `power >= 0` and a
  valid *square* matrix; `power=0` returns the sparse identity, `power=1`
  returns a copy (including explicit stored zeros). Higher powers use
  binary exponentiation in `O(log(power))` sparse products. The result can
  become dense as the power grows. Negative powers are intentionally
  rejected, matching SciPy's sparse function rather than NumPy's dense
  matrix-power behavior. A `0 x 0` input returns a `0 x 0` sparse matrix.
- `inv(a, info, drop_tol, block_size, pivot_tol)` computes an *explicit*
  inverse. `drop_tol` defaults to zero (no intentional sparsification);
  `block_size` defaults to 8 identity columns per solve; `pivot_tol` defaults
  to zero and is a minimum absolute pivot size. Pivoted `splu` factors are
  constructed **once** and reused for each RHS block. Diagonal matrices,
  including those with explicitly stored off-diagonal zeros, use direct
  reciprocals in `O(n+nnz)` work. `drop_tol > 0` drops inverse entries of
  absolute magnitude at most that value; this intentionally changes the
  inverse and is not an error-control guarantee.
- Status `0` means success, `-1` means invalid input or argument, `-2`
  means numerical breakdown/overflow or output-size overflow, and `-4`
  means singular or a pivot below `pivot_tol`. On failure the returned
  matrix is unallocated, and callers must check status before use.
- Computing `inv(A)` may incur major fill-in: the exact inverse of a sparse
  matrix is usually dense. If the goal is to compute `A^{-1} b`, prefer
  `spsolve(A,b,...)` or a reusable `splu` factor. These routines do not
  support complex sparse matrices, negative powers, or dense fallbacks.

Reproducible checks:

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_matrix_functions_scipy.py --mode both --cases 120
```

See [validation details](sparse_matrix_functions_validation_2026-11-05.md).

## Locally optimal block preconditioned conjugate gradients (`lobpcg`)

`lobpcg(a,x0,eigenvalues,eigenvectors,info,...)` computes a small number
of extremal eigenpairs of a real symmetric CSR matrix. Unlike the single-
vector restarted Lanczos method, LOBPCG optimizes a **block** of eigenvectors
in the local Rayleigh--Ritz space spanned by current approximations `X`,
preconditioned residual directions `W`, and previous conjugate directions `P`.
The reduced symmetric eigensystem has at most `3*k` rows; the sparse operator
is never copied into an `n`-by-`n` dense matrix.

| Optional argument | Meaning |
|---|---|
| `largest` | `true` (default) returns algebraically largest eigenvalues, descending; `false` returns smallest, ascending |
| `tol` | Absolute (projected, when constrained) residual 2-norm tolerance, default `n*sqrt(epsilon)`; zero requests roundoff scale |
| `maxiter` | Block iteration budget, default 20 |
| `b` | Symmetric positive-definite CSR mass matrix for `A*x=lambda*B*x`; absent means identity |
| `m` | Symmetric positive-definite CSR *preconditioner action* `M*r`, not a matrix to factor |
| `jacobi` | Use the positive diagonal preconditioner `1/abs(diag(A))`; cannot be combined with `m` |
| `y` | `n`-by-`ny` full-rank constraint directions; enforce `Y^T*B*X=0` |
| `iterations` | Number of completed block iterations |
| `residual_norms` | True residual norms, or Lagrange-multiplier-projected residual norms with constraints |
| `check_symmetry` | Check A/B symmetry by default; disable for previously validated inputs |

Example (one-based CSR, ordinary eigenproblem):

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : lobpcg
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: eigval(:), eigvec(:,:), errors(:)
real(dp) :: start(4,2)
integer :: info

a = csr_from_dense(reshape([8.0_dp,0.0_dp,0.0_dp,0.0_dp, &
    0.0_dp,7.0_dp,0.0_dp,0.0_dp, &
    0.0_dp,0.0_dp,2.0_dp,0.0_dp, &
    0.0_dp,0.0_dp,0.0_dp,1.0_dp], [4,4]))
start(:,1) = [1.0_dp,1.0_dp,1.0_dp,0.0_dp]
start(:,2) = [0.0_dp,1.0_dp,1.0_dp,1.0_dp]
call lobpcg(a,start,eigval,eigvec,info,largest=.true., &
    maxiter=80,tol=1.0e-9_dp,residual_norms=errors)
if (info /= 0) error stop 'LOBPCG did not converge'
```

The input `x0` must have independent columns after projecting out the
constraints; it determines `k`. `B` and any user-supplied `M` must be
symmetric positive definite (global positive definiteness is an assumed
precondition, not proved by this implementation). For constraints, the
reported residual removes the component in `range(BY)` from
`A*x-lambda*B*x`, as required by the constrained eigenproblem. Without
constraints it is the ordinary full eigenpair residual. Input validation
returns `-1`, invalid diagonal preconditioning `-3`, numerical breakdown
`-2`, and a **positive** status when the iteration budget is exhausted.
On nonconvergence the routine returns the **best-residual block seen so
far**, not an unverified success. In all cases check `info` and
`residual_norms`. Small or tightly clustered gaps may require more
iterations or preconditioning. This initial implementation is real-only;
it does not yet expose SciPy's callback operators, iteration histories, or
restart-control parameter, and does not implement SciPy's dense fallback
when the block is large relative to `n`.

Tests: `test/test_sparse_lobpcg.f90` and
`tools/sparse_linalg_reference/compare_lobpcg_scipy.py`, using a separate
Fortran reference driver. Independent NumPy/SciPy dense eigensolvers are
used to check eigenvalues; the script separately verifies each returned
residual and `B`-orthogonality, including non-diagonal B and general
constraints. See `sparse_lobpcg_validation_2026-11-05.md`.

## Threshold incomplete LU (`spilu` — ILUT subset)

`spilu` adds **incomplete LU with threshold dropping and bounded per-row fill**,
independently implemented in `scifort_sparse_linalg`. It accepts canonical
one-based **CSR or CSC** real `dp` matrices:

| Routine | Signature / behavior |
|---|---|
| `spilu(a,factor,info,drop_tol,fill_factor,diagonal_shift)` | Sparse ILUT, default `drop_tol=1e-4`, `fill_factor=10`, `diagonal_shift=0`; expects a square matrix. |
| `spilu_solve(factor,b,x,info)` | Apply the resulting incomplete-LU factors, with either vector or two-dimensional RHS. |
| `valid_ilu0(factor)` | Validate the common stored factor type for ILU(0) **and** ILUT. |
| `gmres(...,ilu=factor)` / `bicgstab(...,ilu=factor)` | Reuse the existing preconditioner argument; no change to these solvers. |

```fortran
use scifort_kinds, only: dp
use scifort_sparse, only: csr_matrix, csr_from_dense
use scifort_sparse_linalg, only: ilu0_factor, spilu, spilu_solve, gmres
implicit none
type(csr_matrix) :: a
type(ilu0_factor) :: preconditioner
real(dp), allocatable :: x(:), approximate(:)
integer :: info

a = csr_from_dense(reshape([4.0_dp,-1.0_dp,-1.0_dp,4.0_dp],[2,2]))
call spilu(a,preconditioner,info,drop_tol=1.0e-4_dp,fill_factor=3.0_dp)
if (info /= 0) error stop 'ILUT factorization failed'
call spilu_solve(preconditioner,[1.0_dp,2.0_dp],approximate,info)
if (info /= 0) error stop 'ILUT application failed'
call gmres(a,[1.0_dp,2.0_dp],x,info,rtol=1.0e-10_dp,ilu=preconditioner)
if (info /= 0) error stop 'Preconditioned GMRES failed'
```

The incomplete factorization uses the original row's largest absolute
coefficient as a relative threshold scale. Off-diagonal entries with
magnitude at most `drop_tol` times that scale are discarded. Additional
entries are retained by largest magnitude subject to separate per-row lower
and upper budgets, `int(fill_factor*max(1, original_strict_triangle_count))`.
The diagonal is always retained. **The fill parameter is a per-row bound,
not a strict upper bound on total factor nonzeros.** `drop_tol` must be in
`[0,1]`; `fill_factor` must be finite and at least one. With zero dropping
and sufficient fill (for example `fill_factor >= n`) the factorization is
exact LU **when no row pivot is needed**. A nonnegative diagonal shift can
stabilize a missing or tiny pivot; it factors `A + shift*I`, and the
preconditioned GMRES/BiCGSTAB still solves the original `A` system.

The `ilu0_factor` type is shared intentionally to preserve the existing
preconditioner APIs, although factors returned by `spilu` **may contain
fill beyond the original CSR pattern**. They can be applied by either
`ilu0_solve` or `spilu_solve`.

**Important compatibility distinction:** SciPy's `scipy.sparse.linalg.spilu`
is backed by SuperLU and supports column ordering, partial pivoting,
equilibration, and several different dropping rules. This SciFort subset is
**natural-order, unpivoted ILUT** and may report `info=-3` for a matrix
that SuperLU could factor successfully. It should not be presented as a
bitwise or factor-pattern equivalent of SciPy's routine. Use the pivoted
`scifort_sparse_linalg.splu` for difficult direct solves. `spilu` is intended
as a lightweight preconditioner for well-scaled sparse systems. No external
sparse solver is required.

Regression files: `test/test_sparse_spilu.f90`,
`tools/sparse_linalg_reference/sparse_spilu_reference_driver.f90`, and
`tools/sparse_linalg_reference/compare_spilu_scipy.py`. The randomized tests
compare exact-fill factors with dense no-pivot LU, approximate-factor
reconstruction and residuals, and GMRES on the original system. They use
SciPy SuperLU's independent implementation for numerical comparison but
do not assert identical dropping or pivoting decisions. See also
[`sparse_ilut_validation.md`](sparse_ilut_validation.md) for commands and results.

## Sparse triangular and pivoted direct solvers (`spsolve`, `splu`)

`scifort_sparse_linalg` provides **dependency-free real-binary64** direct
solvers for both canonical, one-based **CSR** and **CSC** storage:

| API | Signature / behavior |
|---|---|
| `spsolve_triangular(a,b,x,info,lower,unit_diagonal)` | Triangular solution; accepts `b(:)` or `b(:,:)`; `lower=.true.` and `unit_diagonal=.false.` by default. Entries outside the selected triangle are ignored. |
| `splu(a,factor,info,pivot_tol)` | Factor a square CSR or CSC matrix with partial **row** pivoting into `type(splu_factor)`; performs dynamic sparse fill-in. |
| `splu_solve(factor,b,x,info)` | Apply previously computed factors to a vector or RHS block, without refactoring. |
| `spsolve(a,b,x,info,pivot_tol)` | Convenience routine combining `splu` and `splu_solve`; accepts vectors or RHS blocks. |
| `valid_splu(factor)` | Check CSR factor shape, stored pivots, finite data, and permutation validity. |

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : spsolve, splu, splu_factor, splu_solve
implicit none

type(csr_matrix) :: a
type(splu_factor) :: factor
real(dp), allocatable :: x(:), x2(:)
integer :: info

a = csr_from_dense(reshape([0.0_dp,2.0_dp,1.0_dp,3.0_dp],[2,2]))
call spsolve(a,[1.0_dp,5.0_dp],x,info)
if (info /= 0) error stop 'spsolve failed'
call splu(a,factor,info)
if (info /= 0) error stop 'splu failed'
call splu_solve(factor,[2.0_dp,10.0_dp],x2,info)
if (info /= 0) error stop 'splu_solve failed'
```

The `splu_factor` stores the combined sparse factors `lu` and the final row
permutation `row_perm`, with **`P*A=L*U`**, `L` unit lower triangular, and
`U` upper triangular. `row_perm(i)` is the one-based source row for row `i`
of `P*A`. This differs from SciPy's `SuperLU` object, which has separate
`L`, `U`, and both row and column permutations. The **absolute** minimum
pivot `pivot_tol` defaults to zero; set a positive tolerance when near-singular
matrices should be rejected explicitly. `info=0` denotes success, `-1` invalid
input, `-2` non-finite arithmetic/breakdown, `-3` zero or missing triangular
diagonal, and `-4` a singular or below-threshold LU pivot. Check `info` before
using any returned result.

The implementation uses sorted, dynamically growing sparse rows during
Gaussian elimination, and packs the resulting factors into canonical CSR.
It does not form a dense `n*n` work matrix. However, *fill-in may be large*:
this is a basic serial row-pivoted implementation, not a replacement for
SuperLU/UMFPACK for very large or difficult problems. It currently has no
column reordering (AMD/COLAMD), equilibration, transpose solves, complex
arithmetic, or pivot-growth diagnostics. CSC input is converted once to CSR.
For positive-definite or well-conditioned matrices, iterative CG/GMRES may be
more efficient. To solve multiple RHS for the same matrix, prefer one `splu`
call and repeated `splu_solve` calls.

Regression tests: `test/test_sparse_direct.f90` and
`tools/sparse_linalg_reference/compare_sparse_direct_scipy.py`. The randomized
reference checks factor reconstruction `P*A=L*U`, solutions against SciPy
`scipy.sparse.linalg.spsolve`/`spsolve_triangular`, and both CSR/CSC storage.
The source code is original SciFort MIT-licensed implementation of standard
Gaussian elimination and triangular substitution; no SuperLU, UMFPACK, SciPy,
or external sparse-solver implementation code is incorporated.

## Supported methods

| Procedure | Matrix shape | Method | Preconditions |
|---|---|---|---|
| `spsolve_triangular` | n by n | Direct forward/back substitution; CSR/CSC | Triangular selected half with nonsingular diagonal or implicit unit diagonal |
| `spsolve` / `splu` / `splu_solve` | n by n | Pivoted sparse LU with elimination fill-in; CSR/CSC | Nonsingular square real matrix |
| `cg` | n by n | Preconditioned conjugate gradients | Symmetric positive definite A; optional Jacobi diagonal positive |
| `minres` | n by n | Paige--Saunders MINRES with streamed Lanczos QR | Symmetric real A, indefinite/singular allowed; optional SPD Jacobi |
| `bicgstab` | n by n | Right-preconditioned BiCGSTAB | General nonsingular A; optional Jacobi or ILU(0) |
| `gmres` | n by n | Restarted left-preconditioned GMRES, twice modified Gram-Schmidt | General nonsingular A; optional Jacobi or ILU(0) |
| `lsmr` | m by n | Golub--Kahan LSMR, with optional Tikhonov damping | General rectangular A, least-squares/minimum-norm solution |
| `lsqr` | m by n | Paige--Saunders LSQR bidiagonalization, optional damping and variance estimate | General rectangular A, including rank-deficient systems |
| `expm_multiply` / `expm_multiply_interval` | n by n | Trace-shifted, scaled Taylor exponential actions | Real CSR and real finite vector/block RHS |
| `svds` | m by n | Golub--Kahan right Krylov subspace and one-sided Jacobi projected SVD | Partial largest singular triplets; full-basis smallest triplets; real CSR |
| `ilu0` / `ilu0_solve` | n by n | ILU(0) no-fill incomplete LU and triangular application | Structurally present, nonzero diagonal pivots; no pivoting |
| `eigsh` | n by n | Symmetric block subspace iteration (default) or `method='lanczos'` or `method='lanczos_locked'`; all support MINRES-based shift-invert | Real symmetric A, k in 1..n-1; normal-mode shifted Ritz selection |
| `eigsh_lanczos` | n by n | Thick-restarted symmetric Lanczos/Rayleigh--Ritz with optional hard locking and real shift-invert | Symmetric A, reduced subspace for LA, SA, LM; full subspace for SM |
| `eigs` | n by n | Complex projected Schur QR (default) or `method='arnoldi'`/`'arnoldi_schur'`; all support complex-GMRES shift-invert | Real nonsymmetric A; complex-valued eigenpairs; full complex inverse for shifts |
| `eigs_arnoldi` | n by n | Thick-restarted Arnoldi with soft locking (default), optional Schur invariant-subspace locking, and complex shift-invert | Real nonsymmetric A; complex-valued eigenpairs; normal shifted selection |

CG, BiCGSTAB, GMRES, `eigs`, the `eigsh` subset, and LSQR are original
SciFort implementations of standard published numerical methods. MINRES
adapts the Paige--Saunders Lanczos/QR recurrence with attribution to SciPy's
BSD-3-Clause `minres.py` (see `LICENSES/scipy-minres-BSD-3-Clause.txt`). `lsmr` adapts the **existing** BSD-3-Clause
`scifort_optimize_lsmr.lsmr_solve` module through an internal CSR/CSR-transpose
callback pair; it does not duplicate or vendor upstream code. This adapter
inherits SciFort's current LSMR kernel numerical behavior and stopping rules.

## Examples

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : cg, minres, eigsh, eigs, gmres, lsmr, lsqr, ilu0, ilu0_solve, ilu0_factor
implicit none

type(csr_matrix) :: a, rectangular
type(ilu0_factor) :: preconditioner
real(dp), allocatable :: x(:)
integer :: info

a = csr_from_dense(reshape([4.0_dp,1.0_dp,1.0_dp,3.0_dp],[2,2]))
call cg(a,[1.0_dp,2.0_dp],x,info,rtol=1.0e-10_dp,jacobi=.true.)
if (info /= 0) error stop 'CG did not converge'

call minres(a,[1.0_dp,2.0_dp],x,info,rtol=1.0e-10_dp)
if (info /= 0) error stop 'MINRES did not converge'

call gmres(a,[1.0_dp,2.0_dp],x,info,restart=2)
if (info /= 0) error stop 'GMRES did not converge'

call ilu0(a,preconditioner,info)
if (info /= 0) error stop 'ILU0 factorization failed'
call gmres(a,[1.0_dp,2.0_dp],x,info,ilu=preconditioner)
if (info /= 0) error stop 'ILU-preconditioned GMRES failed'

rectangular = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
    0.0_dp,1.0_dp,2.0_dp,1.0_dp],[4,2]))
call lsmr(rectangular,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],x,info,damp=0.1_dp)
if (info < 0) error stop 'invalid LSMR arguments'

call lsqr(rectangular,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],x,info, &
    atol=1.0e-10_dp,btol=1.0e-10_dp)
if (info < 0) error stop 'invalid LSQR arguments'
```

## Sparse matrix-exponential action (`expm_multiply`)

The overloads accept a one-based canonical `csr_matrix` with **real double-
precision** values, including explicitly stored zeros. They compute an
exponential *action* directly, never constructing or storing `expm(A)`.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense
use scifort_sparse_linalg, only : expm_multiply, expm_multiply_interval
implicit none

type(csr_matrix) :: a
real(dp), allocatable :: y(:), grid(:,:), block_result(:,:)
integer :: info, products
real(dp) :: block(2,2)

a = csr_from_dense(reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2]))
call expm_multiply(a,[1.0_dp,0.0_dp],y,info,t=acos(-1.0_dp)/2.0_dp, &
    products=products)
if (info /= 0) error stop 'exponential action failed'
! y is approximately [0,1] and products counts CSR times block calls.
block = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
call expm_multiply(a,block,block_result,info)
if (info /= 0) error stop 'block exponential action failed'
call expm_multiply_interval(a,[1.0_dp,0.0_dp],-1.0_dp,1.0_dp,5, &
    grid,info,endpoint=.true.)
if (info /= 0) error stop 'time-grid exponential action failed'
! grid has shape (2,5); grid(:,i) is the action at grid time i.
```

```fortran
call expm_multiply(a,b,y,info,t,tracea,tol,scaling_steps,products)
call expm_multiply_interval(a,b,start,stop,num,y,info,endpoint,tracea,tol,products)
```

- With a vector `b(:)`, `expm_multiply` returns `y(:)`; with a block `b(:,:)`,
  the result is `y(:,:)`. `t` defaults to 1. `start`, `stop` and `num` are
  required for intervals. A vector interval returns `y(n,num)`, and a block
  interval returns `y(n,nrhs,num)`—**the time index is last**, unlike SciPy's
  Python array convention. `endpoint` defaults to true; unlike SciPy 1.17,
  a single (`num=1`) output sample is allowed and denotes time `start`.
- `tracea` is an optional finite *candidate* trace for scalar shifting;
  because an explicit CSR matrix has an inexpensive exact trace, SciFort
  compares the candidates and selects the shift with smaller operator
  infinity norm, avoiding unnecessary cancellation from inaccurate hints.
  It is not an API for altering the underlying matrix.
- Each scaled substep applies a truncated Taylor series to `A-mu*I`. An
  explicit infinity-norm upper bound on the omitted series tail controls
  local termination. A scaling substep has norm at most 2; at most 80 terms
  and 200,000 scaling steps are allowed. Default `tol` is `4e-15`; it must
  be finite and in `[0,1)`, with a four-epsilon machine floor. The bound
  controls the local Taylor truncation, **not the total relative error** in
  severely nonnormal/cancellation-dominated matrix exponentials.
- Each time grid is anchored at the sample nearest zero, then propagated
  toward the earlier and later samples. This avoids a zero that underflowed
  at a remote time propagating into a later, representable result.
- `info=0` indicates success; `-1` invalid arguments or an excessive scaling
  budget; `-2` nonfinite intermediate/overflow; `1` Taylor did not converge
  within the allowed terms. Inspect `info` before using output; failed grid
  computations can contain only partial results.
- `scaling_steps` and `products` optionally report the number of scaled
  Taylor substeps and calls to multiply the CSR matrix by a dense RHS block.
  The time-grid API offers `products` only, totaled over all time points.
- **Current scope:** real CSR matrices and real RHS blocks only; complex
  operators, matrix-free `LinearOperator`s, sparse-matrix exponential
  construction (`expm`), and full SciPy broadcasting semantics remain
  unimplemented. Storage is O(nnz+n*nrhs), excluding the requested grid
  output; time-grid storage is O(n*nrhs*num). No BLAS/LAPACK dependency.

Numerical references: [Al-Mohy and Higham (2011)](https://doi.org/10.1137/100788860);
[`sparse_expm_validation_2026-11-04.md`](sparse_expm_validation_2026-11-04.md).

## Partial singular value decomposition (`svds`)

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_coo
use scifort_sparse_linalg, only : svds
implicit none
type(csr_matrix) :: a
real(dp), allocatable :: u(:,:), singular(:), vt(:,:), errors(:)
integer :: info, steps

a = csr_from_coo(4,3,[1,2,3],[1,2,3],[5.0_dp,3.0_dp,1.0_dp])
call svds(a,2,u,singular,vt,info,which='LM',tol=1.0e-9_dp, &
    ncv=3,iterations=steps,residual_norms=errors)
if (info /= 0) error stop 'partial SVD did not converge'
! singular is [3, 5], u has shape (4,2), vt has shape (2,3).
```

```fortran
call svds(a,k,u,s,vt,info,which,tol,ncv,maxiter,v0,iterations,residual_norms)
```

- Works with real `csr_matrix` inputs of any nonempty rectangular shape,
  including singular, rank-deficient and explicitly stored-zero matrices.
  All CSR indices remain one-based. `k` may be from 1 through `min(m,n)`.
  Singular values are returned in **ascending numerical order**, even for
  `which='LM'` (largest magnitude). `u` is `(m,k)`; `vt` is `(k,n)`.
- `which='LM'` (default) finds the largest singular values with a partial
  Golub--Kahan basis. `which='SM'` finds the smallest, but **currently requires
  the full shorter dimension as `ncv`** (and `maxiter>=min(m,n)` if specified).
  This uses O((m+n)*min(m,n)) storage and is *not* a scalable smallest-SV
  algorithm; a future harmonic Ritz or shift-invert method is preferable.
- The default `ncv` is `min(min(m,n),max(40,8*k+16))` for `LM`. `ncv` caps
  memory and Golub--Kahan steps; no thick or implicit restarting is available
  yet. `maxiter` can reduce this budget, but cannot exceed `ncv`.
  The optional `v0` must have `min(m,n)` finite elements and a nonzero norm;
  it seeds the right side of the **tall orientation** (A or A-transpose).
- Two-pass reorthogonalization stabilizes the Golub--Kahan search. A small
  projected one-sided Jacobi SVD avoids explicitly forming `A^T A`. On lucky
  Krylov breakdown, a new orthogonal direction is seeded, permitting repeated
  singular values to be recovered in a full or sufficiently rich subspace.
- Success (`info=0`) requires **true two-sided** singular residuals, computed
  against the original sparse matrix and its transpose:
  `sqrt(||A*v-s*u||^2 + ||A^T*u-s*v||^2) <= tol*sigma_max`.
  `tol` defaults to `1e-8`; zero selects a near-machine-precision floor.
  `residual_norms` reports these quantities. Positive `info` is the number of
  steps used **when the requested triplets have not converged**. A negative
  status indicates invalid arguments (`-1`) or numerical breakdown (`-2`).
- Accuracy of singular values does not by itself certify that all *globally*
  largest values have been identified if a starting vector has no component
  in a relevant invariant subspace. Returned residuals certify individual
  triplets, not spectral completeness. Large clustered, smallest-singular-
  value and near-machine-rank problems may require a richer search basis.
- No external BLAS, LAPACK or PROPACK dependency is introduced. This is a
  **new SciFort partial SVD subset**, not a wrapper of SciPy's ARPACK,
  PROPACK or LOBPCG solvers. Complex sparse matrices and the full SciPy `svds`
  options (solver families, RNG policy, vector-return toggles) are not yet
  implemented.

Validation: `test/test_sparse_svds.f90`,
`tools/sparse_linalg_reference/compare_svds_scipy.py`, and
`docs/sparse_svds_validation_2026-11-04.md`.

## Signatures

All solvers are subroutines; `x(:)` is `allocatable, intent(out)` and
is allocated internally. For `cg`, `bicgstab`, and `gmres`:

```fortran
call cg(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm)
call bicgstab(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
call gmres(a,b,x,info,rtol,atol,maxiter,restart,x0,jacobi,iterations,residual_norm,ilu)
```

- `rtol` defaults to `1e-5`, `atol` to zero. Success requires the **true**
  Euclidean norm `||b-A*x|| <= max(rtol*||b||,atol)`. Both tolerances must
  be nonnegative and finite; a zero RHS succeeds without iterations when
  the starting guess is zero.
- `maxiter` defaults to `max(1,10*n)`. For CG and BiCGSTAB it counts inner
  iterations; for GMRES it counts **restart cycles** (as in SciPy's modern
  default `callback_type`). Zero is allowed and reports nonconvergence if
  the initial solution fails the tolerance. Default GMRES `restart=min(20,n)`;
  positive values greater than `n` are reduced to `n`.
- `x0` is a finite initial guess (zeros by default). The `jacobi` boolean
  enables a diagonal inverse; this is a more limited preconditioner API than
  SciPy's general `M` LinearOperator. Explicit zero diagonal entries and
  missing diagonal entries both cause a specific error. For CG, the user
  must ensure a positive-definite preconditioner. BiCGSTAB and GMRES can
  instead receive an `ilu0_factor` with keyword `ilu=preconditioner`.
  `jacobi=.true.` and `ilu=` cannot be combined. CG does **not** accept ILU(0),
  which may be nonsymmetric or indefinite.
- `info=0` means success; `info>0` means convergence not achieved within
  the budget, with the positive value representing the budget; `info=-1`
  means invalid shape/storage/tolerance; `info=-2` means numerical breakdown
  (for example, an unresolvable singular matrix); `info=-3` means missing,
  zero, or noninvertible Jacobi/ILU(0) pivot. Public constants are
  `sparse_linalg_success`, `sparse_linalg_invalid_input`,
  `sparse_linalg_breakdown`, and `sparse_linalg_bad_diagonal`.
  `iterations` returns inner iterations, independently of `info`.
- `residual_norm` is recomputed from the actual solution. For invalid input,
  it is initialized to `huge(1.0_dp)`, not an estimated residual.

`lsmr` has the separate signature:

```fortran
call lsmr(a,b,x,istop,itn,normr,damp,atol,btol,conlim,maxiter)
```

- `damp` defaults to 0; for positive `damp`, LSMR minimizes
  `||A*x-b||_2**2 + damp**2 * ||x||_2**2`.
- `atol=1e-6`, `btol=1e-6`, `conlim=1e8`, and `maxiter=min(m,n)` by default.
  For the meaning of nonnegative `istop` codes, see the SciPy LSMR reference.
  Invalid input returns `istop=-1`; zero-dimension systems return the
  trivial zero solution.
- `normr`, if requested, is the **actual data residual** `||A*x-b||_2`;
  unlike SciPy's LSMR diagnostic it does **not** include `damp*x` in the
  residual vector. This distinction matters when comparing damped solves.
- The current LSMR adapter constructs a CSR transpose once, requiring
  `O(nnz+m+n)` additional storage, and uses allocation-free matvec callbacks.


## MINRES — symmetric indefinite linear systems

```fortran
call minres(a,b,x,info,rtol,atol,maxiter,x0,shift,jacobi, &
    check_symmetry,iterations,residual_norm)
```

- Solve `(A - shift*I)x = b` with a real symmetric matrix that may be
  indefinite or singular. The symmetric Lanczos recursion uses `O(nnz+n)`
  matrix and vector storage and is more appropriate than CG for indefinite
  matrices. `shift=0`, `rtol=1e-5`, `atol=0`, `maxiter=max(1,10*n)` by default.
- `jacobi=.true.` applies inverse-diagonal preconditioning only when the
  diagonal is strictly **positive**; this guarantees a positive-definite
  diagonal preconditioner, which MINRES requires. `ilu0_factor` is *not*
  supported because incomplete LU is generally nonsymmetric.
- `check_symmetry=.true.` checks the CSR matrix against its transpose with
  relative entry tolerance `1e-12`; it is off by default, consistent with
  SciPy's optional check. True shifted residuals are evaluated at every step:
  success requires `||b-(A-shift*I)x|| <= max(atol, rtol*||b||)`.
- `info=0` indicates achieved true-residual tolerance; positive `info`
  indicates iteration exhaustion; `-1` invalid input or failed optional
  symmetry check; `-2` numerical breakdown or incompatible singular system;
  `-3` zero, missing or non-positive Jacobi diagonal. If a singular system
  has no exact solution, the solver may report breakdown instead of the
  least-squares optimum; use `lsqr`/`lsmr` for inconsistent systems.
- Unlike SciPy 1.17's current initial-residual handling with nonzero `x0`
  and nonzero `shift`, SciFort includes `shift*x0` in that initial residual.

## First `eigsh` subset — real symmetric eigenpairs

```fortran
real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:), residuals(:)
call eigsh(a,k,eigenvalues,eigenvectors,info,which='LA',tol=1.0e-8_dp, &
    maxiter=500,ncv=20,residual_norms=residuals)
```

- Supports the `which` values `LA`, `SA`, `LM` (default), and `SM`, with `k`
  in `[1,n-1]`, optional `v0(:)` starting direction, `ncv` block subspace
  dimension (`k < ncv <= n`), `iterations` and optional `check_symmetry`.
  Symmetry is checked by default; disable it only for known symmetric data.
- The implementation orthogonalizes a block of trial vectors, updates it
  by a spectral transformation and obtains Ritz eigenpairs from a small
  symmetric projected matrix. `maxiter` counts block updates, not Lanczos
  inner iterations. It uses `O(nnz+n*ncv+ncv**2)` storage, not a dense n-by-n
  copy. The `SM` transformation uses a squared operator and may converge
  **very slowly** for clustered near-zero eigenvalues; for small matrices,
  `ncv=n` ensures the full subspace. This is not an implicitly restarted
  ARPACK/Lanczos algorithm.
- Returned eigenvalues are in ascending algebraic order (SciPy's ordering
  when eigenvectors are requested). `eigenvectors(:,i)` is normalized and
  `residual_norms(i)` is the **true** `||A*v-lambda*v||` rather than an
  estimate. `tol` defaults to `1e-8` and applies relative to an upper bound
  on `||A||`; zero requests machine precision. `info=0` means all requested
  Ritz residuals passed; `info>0` means the iteration budget was exhausted
  and approximate eigenpairs are returned; negative status means invalid
  input or numerical breakdown.
- **Normal shift-invert:** supply `sigma=<real>`. Every orthogonal-subspace
  update solves `(A-sigma*I)Y=Q` using SciFort's shifted MINRES, without
  forming a dense inverse. As in SciPy, `which` selects according to the
  transformed eigenvalues `mu = 1/(lambda-sigma)`: `LM` finds the eigenvalues
  *closest to* the shift; `LA` the largest positive transformed values,
  typically just **above** the shift; `SA` the most negative transformed
  values, typically just **below** the shift; and `SM` the smallest
  transformed magnitudes, corresponding to eigenvalues *farthest from*
  the shift. **Shifted `SM` requires `ncv=n`**: reduced inverse-subspace
  iteration preferentially retains nearby modes, not the distant ones.
- `inner_rtol` (default `1e-12`) and `inner_maxiter` (default `10*n`)
  control the shifted MINRES linear solves. `maxiter` counts outer updates.
  A failed or singular inner solve reports `info=-2`; avoid choosing `sigma`
  exactly at an eigenvalue. Returned eigenvalues and true residuals refer to
  the **original** matrix, not its inverse. Transformed-mode selection may
  converge slowly when the desired modes are poorly represented in a small
  subspace; the current implementation is still not implicit-restart ARPACK.
- Not supported yet: generalized mass matrices `M`, complex Hermitian
  matrices, callback linear operators, `BE` mode, and eigenvalues-only
  output. This is a narrower numerical subset, not one-to-one ARPACK parity.

## Nonsymmetric sparse `eigs` (real CSR, complex eigenpairs)

```fortran
use scifort_sparse_linalg, only : eigs
complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
real(dp), allocatable :: residuals(:)
integer :: info
call eigs(a,2,eigenvalues,eigenvectors,info,which='LM', &
    tol=1.0e-8_dp,ncv=20,maxiter=600,residual_norms=residuals)
if (info /= 0) error stop 'eigs did not converge'
```

- The **input is real** canonical one-based CSR; values and eigenvectors are
  complex `complex(dp)` because real nonsymmetric matrices may have complex
  eigenvalues. Eigenvectors are columnwise normalized; `residual_norms(j)`
  measures the true Euclidean norm `||A*v_j-lambda_j*v_j||`.
- Supported selection criteria: `LM` (largest magnitude, default), `LR`
  (largest real part), `SR` (smallest real part), `LI` (largest imaginary
  part), `SI` (smallest imaginary part), and `SM` (smallest magnitude).
  For `SM`, a reduced subspace (`ncv<n`) now uses **implicit zero-shift
  inversion** via complex restarted GMRES, finding the eigenvalues nearest
  zero when the matrix is nonsingular. For `ncv=n`, the original full
  projected eigensystem remains available and also handles zero eigenvalues.
  Other criteria support `k < ncv <= n`; the default is
  `min(n,max(2*k+8,20))`.
- The algorithm orthonormalizes a complex block, applies a spectral
  transformation, and diagonalizes the small nonsymmetric Rayleigh--Ritz
  projection. Its internal **complex Householder/shifted-QR Schur solver**
  does not require external LAPACK. It stores `O(nnz+n*ncv+ncv**2)` entries,
  but projected QR can be expensive for large `ncv`. A full subspace
  (`ncv=n`) gives a dense-size spectral calculation and is best reserved
  for small matrices.
- `tol=1e-8` by default (zero requests machine precision). `maxiter=600`
  counts outer block updates; `v0` is an optional **complex** finite nonzero
  initial direction of length `n`. Status `0` means all requested true
  residuals satisfy the tolerance, a positive status means iteration budget
  exhausted, `-1` invalid input, and `-2` projection QR or numerical
  breakdown. The returned order is by the requested priority, *not* SciPy's
  particular final ordering. On real operators, selecting `LI` or `SI`
  strictly by the signed imaginary part may choose a different member of
  a conjugate pair than SciPy's real ARPACK wrapper.
- **Complex shift-invert:** provide `sigma=cmplx(real_part,imag_part,dp)`
  with `which='LM'` to select the eigenvalues closest to this complex point.
  The implementation applies the full complex inverse `(A-sigma*I)^(-1)`
  through restarted complex GMRES, operating directly on the real CSR storage.
  `inner_rtol` defaults to `1e-11`, `inner_maxiter` to `max(200,10*n)`
  total inner Arnoldi steps, and `inner_restart` to `min(40,n)`.
  All three can be overridden. The inner solve checks the *true* residual
  after each restart, and a failed or exactly singular shift is reported as
  `info=-2`; no inaccurate eigenpair is declared converged.
  In contrast, SciPy's **real-matrix** `eigs` defaults to the real part of a
  complex-shifted inverse (`OPpart='r'`). For a direct comparison with this
  SciFort routine, cast the SciPy input matrix to complex so it uses the
  full complex inverse.
- This is an **initial subspace-iteration solver, not implicitly restarted
  Arnoldi**. Closely clustered, highly nonnormal, or poorly separated
  eigenvalues can converge slowly or fail. True residual checks verify
  the eigenpair equations, but **small residuals do not imply accurate
  eigenvalues for highly nonnormal or defective matrices** (the familiar
  pseudospectral sensitivity of Jordan blocks). For those problems use a
  robust ARPACK/PRIMME backend with conditioning analysis when available.
  Complex-valued CSR *inputs*, generalized eigenproblems, a separate
  user-supplied inverse operator, complex-shift `OPpart` projections, and
  other shifted selection criteria remain future work. For reduced-subspace
  `SM`, the automatic zero-shift solve is inappropriate for exactly singular
  operators; use `ncv=n` or explicitly shift away from the singularity.


## Thick-restarted symmetric Lanczos eigenvalues (SciFort extension)

```fortran
use scifort_sparse_linalg, only : eigsh, eigsh_lanczos
real(dp), allocatable :: vals(:), vecs(:,:), errors(:)
integer :: info, cycles, converged, products

! Direct access to restart and locking diagnostics:
call eigsh_lanczos(a,2,vals,vecs,info,which='LA',ncv=20, &
    tol=1.0e-8_dp,maxiter=150,iterations=cycles, &
    locked_count=converged,matvecs=products,residual_norms=errors)

! Or use the method keyword on the existing symmetric API:
call eigsh(a,2,vals,vecs,info,method='lanczos',which='LA',ncv=20)
```

`eigsh_lanczos` is an opt-in, real-symmetric Krylov/Rayleigh--Ritz method.
It generates a Krylov basis through sparse matrix-vector products with
**full two-pass reorthogonalization** (a stabilized symmetric Lanczos-style
expansion). It solves a small symmetric projected eigenproblem, retains
wanted Ritz vectors across restarts, and preferentially retains already
converged vectors (**soft locking**). A Ritz-residual direction seeds the
next expansion. As in the existing Arnoldi routine, the selected Ritz pairs
must pass the *true* `||A*v-lambda*v||_2` test against the original CSR
matrix before `info=0` is returned.

- `which='LM'` (default), `'LA'`, and `'SA'` work with reduced subspaces.
  `'SM'` requires `ncv=n`; interior modes without a suitable inverse
  spectral transform are not guaranteed to converge in a small subspace.
- `ncv` defaults to `min(n,max(2*k+8,20))`; it must satisfy
  `k<ncv<=n`. `maxiter` defaults to 150 **restart cycles**. `matvecs`
  counts outer sparse products or shifted inverse applications, separately from cycles.
- `tol` defaults to `1e-8`. A zero tolerance requests near machine
  precision. `info>0` means the requested cycles were exhausted without
  satisfying the true residual bound, `-1` means invalid input, and `-2`
  means numerical breakdown. `locked_count` counts currently converged
  *wanted* Ritz vectors; soft locking does not hard-deflate them.
- Eigenvalues are returned in ascending algebraic order, together with
  matching columns of `eigenvectors` and `residual_norms`. Relative
  tolerance uses `max(1,||A||_infinity,|lambda|)` as scale.
- The optional real `v0` must be a finite, nonzero n-vector. By default,
  a deterministic starting vector is used. Symmetry checking is enabled
  unless `check_symmetry=.false.` is supplied.
- Real **shift-invert** is supported through `sigma=` directly or through
  `eigsh(...,method='lanczos',sigma=...)`. The recurrence uses shifted MINRES
  solves and ranks Ritz values of `(A-sigma*I)^(-1)`. `which` refers to the
  transformed spectrum, following SciPy normal mode. `inner_rtol` defaults
  to `1e-12`; `inner_maxiter` defaults to `10*n` for each inverse solve.
  `matvecs` counts outer inverse applications, **not** the MINRES steps.
  All returned eigenvalues and convergence residuals refer to original `A`.
  `which='SM'` still requires `ncv=n`, even with a shift.
- Memory is `O(nnz+n*ncv+ncv**2)`; the small symmetric Jacobi
  eigensolver can dominate for large `ncv`. Unlike an implicitly restarted
  Lanczos/ARPACK implementation, soft restarting does not guarantee rapid
  convergence for clustered or poorly isolated spectral edges. In
  particular, a *small* `ncv` can stagnate even if more restart cycles
  are allowed. The routine then correctly returns `info>0`.

The algorithm is independently implemented SciFort MIT-licensed code from
standard numerical methods. It does not import ARPACK, PRIMME, SciPy, or
third-party solver sources. The source, dedicated regression tests, and
randomized comparisons are described in
[`sparse_lanczos_validation_2026-11-04.md`](sparse_lanczos_validation_2026-11-04.md) and
[`sparse_restarted_shift_validation_2026-11-04.md`](sparse_restarted_shift_validation_2026-11-04.md).

### Optional hard locking and symmetric invariant-subspace deflation

Two opt-in interfaces are supported.  The original soft-locking behavior
remains the default, including for `method='lanczos'`.

```fortran
call eigsh_lanczos(a,2,values,vectors,info,which='LA',ncv=16, &
    hard_locking=.true.,locked_count=nlocked,residual_norms=residuals)
call eigsh(a,2,values,vectors,info,method='lanczos_locked',which='LA',ncv=16)
```

In hard-locking mode, a selected symmetric Ritz vector is **frozen** when
its original-operator residual meets the same public tolerance as the
original solver.  Subsequent search directions are reorthogonalized
against all frozen vectors using two-pass modified Gram--Schmidt.  The
remaining Rayleigh--Ritz problem is therefore projected into the
orthogonal complement of the accepted approximate invariant subspace.
`locked_count` now counts *persistently locked* vectors, rather than
currently converged but still active Ritz vectors.  A positive `info` can
coexist with `locked_count>0`, since some pairs may have converged while
others remain active.  Both returned eigenvectors and reported residuals
are ordered by ascending eigenvalue as in the soft solver.

This orthogonal projection is appropriate for **real symmetric** matrices;
it must not be applied to general right eigenvectors of a nonnormal matrix.
`eigs_arnoldi` remains soft-locked by default; its opt-in Schur-subspace
locking retains nonsymmetric coupling instead of discarding it.  Hard locking is not
ARPACK's implicitly restarted Lanczos algorithm, and it is not guaranteed
to accelerate every problem.  Inexact locked eigenvectors can introduce
small perturbations into the projected complement, so tolerance selection
matters.  The iteration cap still counts restart cycles, and inner MINRES
failures in shift-invert are reported as breakdown rather than success.

**Repeated eigenvalues:** A single Krylov seed cannot span every direction
of a repeated eigenspace.  The direct routine therefore also provides
`block_size=` in hard-locking mode (default `1`, valid integers `1..k`).
For example, to search for both directions of a double dominant eigenvalue:

```fortran
call eigsh_lanczos(a,2,values,vectors,info,which='LA',ncv=12, &
    hard_locking=.true.,block_size=2,maxiter=80)
```

`block_size>1` starts the first cycle with multiple deterministic,
independent seed vectors and expands them as a block Krylov space.  It can
recover multiplicities missed by a single seed, but may require more
restarts and cannot guarantee finding all eigenvalues for arbitrary
initial subspaces.  Passing `block_size` without `hard_locking=.true.`
returns `info=-1`.  The `eigsh(method='lanczos_locked')` dispatcher uses
the default single seed; use `eigsh_lanczos` directly when `block_size`
needs to be specified.  The deterministic seeds do not use global RNG.

The new routines and tests are described in
[`sparse_lanczos_locked_validation_2026-11-04.md`](sparse_lanczos_locked_validation_2026-11-04.md).
Their projection/retention design draws on the mathematical principles of
Kesheng Wu and Horst Simon, *Thick-Restart Lanczos Method for Large
Symmetric Eigenvalue Problems*, SIAM J. Matrix Anal. Appl. 22 (2000),
602--616, [DOI:10.1137/S0895479898334605](https://doi.org/10.1137/S0895479898334605).
No source from the paper or a third-party library is copied.

## Thick-restarted Arnoldi eigenvalues (SciFort extension)

```fortran
use scifort_sparse_linalg, only : eigs_arnoldi
! a is a validated one-based real CSR matrix.
call eigs_arnoldi(a,2,values,vectors,info,which='LM',tol=1.0e-8_dp, &
    ncv=12,maxiter=150,iterations=cycles,residual_norms=residuals, &
    locked_count=nlocked,matvecs=products)
```

`eigs_arnoldi` is a separate direct entry point. The same implementation
can also be selected through `eigs(...,method='arnoldi')` without changing
the default method. The opt-in `method='arnoldi_schur'` selects the new
invariant-subspace locking path. For compatibility, `method='subspace'` explicitly selects
the original algorithm. Unknown method names are rejected. In shifted
mode the recurrence expands a Krylov basis of `(A-sigma*I)^(-1)` using
complex restarted GMRES. On the first cycle it expands an orthonormal
Krylov basis using a two-pass Arnoldi recurrence. It then diagonalizes the
small projected matrix and computes the true residuals of the requested
Ritz eigenpairs. On restart it keeps leading Ritz vectors, places converged
vectors first (**soft locking**), reorthogonalizes the retained subspace,
and expands using a Ritz-residual starting direction followed by Arnoldi
steps. This is a **thick-restart / Ritz-retention method**; it is *not*
Sorensen's implicitly restarted Arnoldi or a hard-deflation algorithm.

- `which` supports `LM`, `LR`, `SR`, `LI`, and `SI` for reduced subspaces.
  `SM` is currently supported only with `ncv=n`; use the original `eigs`
  for automatically selected zero-shift smallest-magnitude problems.
- `ncv` defaults to `min(n,max(2*k+8,20))` and must satisfy `k<ncv<=n`.
  `maxiter` is a count of **restart cycles**, not sparse products; default
  `150`. `matvecs` counts outer sparse products or shifted inverse
  applications; the inner GMRES products are not included.
- `tol` defaults to `1e-8`; `tol=0` requests near-machine precision.
  Success (`info=0`) requires each selected eigenpair to satisfy
  `||A*v-lambda*v||_2 <= tol * max(1, ||A||_infinity, |lambda|)`.
  As in SciFort's other eigenvalue APIs, the output is ordered by requested
  spectral priority, which may differ from SciPy's ordering within ties.
- `locked_count` reports how many *currently requested* eigenpairs meet
  the true-residual threshold. The vectors are retained across restarts,
  but **are not removed from the projected eigenproblem**. This is soft
  locking rather than exact, hard convergence locking.
- `info>0` means the cycle cap was reached; `-1` is invalid input and `-2`
  is small projected eigensolver breakdown. `v0` accepts an optional complex
  starting vector; without one, an explicit deterministic vector is used.
- Memory is `O(n*ncv+ncv**2)` in addition to input CSR storage. No BLAS,
  LAPACK, or third-party eigenvalue sources are added.

**Shift-invert.** Supply complex `sigma=` to `eigs_arnoldi` or to
`eigs(...,method='arnoldi',sigma=...)`. The Arnoldi recurrence applies
`(A-sigma*I)^(-1)` via the pre-existing `complex_shift_gmres` kernel,
using `inner_rtol` (default `1e-11`), `inner_maxiter` (default
`max(200,10*n)`), and `inner_restart` (default `min(40,n)`).
The `which` mode selects *transformed* inverse eigenvalues; final output
uses original eigenvalues `lambda = sigma+1/mu` and certifies each with
`||A*v-lambda*v||`. Because the original CSR input is real but the inner
operator is fully complex, this differs from SciPy's real-matrix default
`OPpart='r'` for nonreal `sigma`; for matching SciPy reference results,
cast its input to complex. As for the unshifted solver, reduced `SM`
is not supported. Nearly singular shifted systems may fail the inner solve
and return `info=-2`; no claim of convergence is made.

**Limitations.** Ritz-vector retention and soft locking do not reproduce
ARPACK's implicit QR shifts or guarantee separation of tightly clustered
or defective eigenvalues. Highly nonnormal matrices may have accurate
residuals but poorly conditioned eigenvalues. Complex input CSR,
generalized operators and user-supplied inverse operators remain future work.
Schur locking is available as an opt-in mode, but it does not guarantee
convergence or complete invariant-subspace recovery for defective operators.

This algorithm follows the mathematical Arnoldi/Rayleigh--Ritz approach
and the motivation for retaining Ritz vectors discussed in Ronald B.
Morgan, *On Restarting the Arnoldi Method for Large Nonsymmetric
Eigenvalue Problems*, Mathematics of Computation 65 (1996), 1213--1230,
[DOI:10.1090/S0025-5718-96-00745-4](https://doi.org/10.1090/S0025-5718-96-00745-4).
No implementation code was taken from that paper or ARPACK.

### Optional Schur locking for nonsymmetric Arnoldi

```fortran
call eigs_arnoldi(a,3,values,vectors,info,which='LR',ncv=8, &
    maxiter=100,schur_locking=.true.,schur_locked_count=nfrozen, &
    locked_count=nconverged,residual_norms=residuals)

! The same algorithm is available through the general eigs dispatcher:
call eigs(a,3,values,vectors,info,which='LR',ncv=8, &
    method='arnoldi_schur')
```

- `schur_locking` defaults to `.false.`. The legacy `eigs_arnoldi`
  subroutine, `eigs(method='arnoldi')` and the default `eigs` algorithm
  are unchanged. `schur_locked_count` is a new optional output of the **direct**
  Arnoldi interface; it counts persistently frozen Schur directions. By
  contrast, `locked_count` continues to count the currently selected pairs
  satisfying the true residual tolerance, even if no directions were frozen.
- When one of the requested Ritz eigenpairs passes the original-matrix true
  residual test, the code tentatively adds its orthogonal component to a
  locked subspace. Before freezing it, the code checks the **invariant
  subspace residual** `||(I-Q*Q^H)*OP*Q||_F` with the active operator
  `OP=A`, or `OP=(A-sigma*I)^(-1)` for shift-invert. Nearly linearly
  dependent directions and insufficiently invariant subspaces are not locked.
  A small projected complex Schur decomposition then rotates the frozen
  block into an orthonormal Schur basis.
- Subsequent Krylov directions are orthogonalized against the locked Schur
  basis, but the **full projected Rayleigh--Ritz matrix** still includes the
  locked block and its upper-right coupling to the remaining Krylov space.
  This distinction is essential for nonsymmetric/nonnormal matrices: simply
  deleting converged right eigenvectors from the original projected
  eigenproblem generally changes the remaining right eigenvectors.
- This remains a thick-restarted method with explicit Ritz-vector retention,
  not implicitly restarted Arnoldi (IRAM), a general Schur reordering scheme,
  or a guarantee of a complete invariant subspace for repeated/defective
  eigenvalues. Freezing an approximate invariant subspace can affect the
  convergence trajectory, and highly nonnormal eigenvalues can be poorly
  conditioned despite small true residuals. For close spectral boundaries,
  approximate Krylov extraction may select a different nearly tied pair.
- `ncv`, `maxiter`, `tol`, `which`, and complex `sigma` retain their
  existing meanings. In this mode, `matvecs` also counts operator
  applications made to validate and Schur-rotate the locked block; inner
  GMRES products are still not included. Extra checks can make it slower
  than soft locking when the spectral edge is already well separated.

The new native and randomized SciPy eigenpair tests are in
`test/test_sparse_arnoldi_schur.f90` and
`tools/sparse_linalg_reference/compare_schur_locked_scipy.py`. Reproduction
and numerical limitations are recorded in
[`sparse_arnoldi_schur_validation_2026-11-04.md`](sparse_arnoldi_schur_validation_2026-11-04.md).

### Optional block Arnoldi for repeated and clustered eigenvalues

`eigs_arnoldi` also accepts `block_size` (default 1) to retain several
independent Krylov search directions instead of one. The same option is
forwarded by `eigs` when `method='arnoldi'` or `'arnoldi_schur'`:

```fortran
call eigs_arnoldi(a,2,values,vectors,info,which='LR',ncv=12, &
    block_size=2,schur_locking=.true.,residual_norms=residuals)

call eigs(a,2,values,vectors,info,method='arnoldi_schur', &
    which='LR',ncv=12,block_size=2)
```

- The accepted range is `1 <= block_size <= k`. With `block_size=b>1`,
  reduced subspaces require `ncv >= k+b+1`. A full-size `ncv=n` space is
  also accepted (no subsequent restart). The unmodified subspace iteration
  (`eigs` with no method, or `method='subspace'`) rejects `block_size`.
- The first block contains deterministic complex initial search vectors;
  subsequent block-Krylov expansions use two-pass orthogonalization. At
  restart, separate Ritz residuals (and guard seeds if needed) seed the
  next block instead of collapsing all residuals into a single vector.
  The default `block_size=1` retains the original algorithm.
- This improves recovery of *semisimple* repeated eigenvalues and closely
  spaced clusters. Eigenvalue multiplicity requires separate independent
  right eigenvectors; accurate eigenvalues or residuals alone do not
  establish multiplicity. This is especially important for nonnormal and
  defective matrices, where no completeness guarantee is provided.
- Both soft and Schur locking, and the complex shift-invert path, support
  block sizes greater than one. Complex eigenvalues and original-matrix
  true residuals are preserved. The algorithm is not ARPACK/IRAM, and
  multiseed search can be slower than single-vector Arnoldi.

Native and SciPy reference tests, including a 1,400-row sparse double
root, are documented in
[`sparse_arnoldi_block_validation_2026-11-04.md`](sparse_arnoldi_block_validation_2026-11-04.md).

## LSQR sparse least-squares API

```fortran
call lsqr(a,b,x,istop,itn,r1norm,r2norm,damp,atol,btol, &
    conlim,iter_lim,x0,anorm,acond,arnorm,xnorm,var,calc_var)
```

This is a new **independently written** Paige--Saunders LSQR implementation,
not a copy of the upstream Fortran LSQR repository. It uses Golub--Kahan
bidiagonalization and only CSR/CSR-transpose products, retaining `O(nnz+m+n)`
matrix and vector storage. `x` is allocated by the solver, even for empty
matrices. Optional arguments follow the numerical meaning of
[`scipy.sparse.linalg.lsqr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html):

- `damp=0`, `atol=btol=1e-6`, and `conlim=1e8` by default. Setting
  `conlim=0` disables the condition-limit check. The default `iter_lim` is
  `2*min(m,n)`; it can be set to zero. All parameters must be finite and
  nonnegative. There is no logging (`show`) option.
- `x0` is an optional finite initial estimate. With positive `damp`, it is
  also the **penalty center**: the objective is
  `||A*x-b||_2**2 + damp**2*||x-x0||_2**2`. The solver works on the
  correction `x-x0` rather than explicitly augmenting the matrix.
- `istop=0` indicates a zero initial residual/zero bidiagonalization,
  `1` a compatible-system residual criterion, `2` the least-squares
  normal-equation criterion, `3` the condition bound, and `7` iteration
  exhaustion. Negative codes `-1` and `-2` mean invalid arguments and
  numerical breakdown, respectively. Codes `4..6` for machine-precision
  detection in SciPy are not currently reported separately.
- `r1norm` and `r2norm` are recomputed from the final iterate. `r2norm`
  includes the damping penalty relative to `x0`; `arnorm` is the **true**
  normal-equation residual rather than a bidiagonalized estimate.
  `anorm`/`acond` are standard iteration-based estimates.
- If `calc_var=.true.` and `var` is requested, `var(:)` contains the
  *approximate* covariance-diagonal estimate accumulated by the LSQR
  recurrences; it is not an exact inverse diagonal and can vary with the
  iteration history, particularly for ill-conditioned matrices. If `var`
  is requested with `calc_var=.false.`, zeros are returned.

## ILU(0) factorization and Krylov preconditioning

```fortran
use scifort_sparse_linalg, only : ilu0_factor, ilu0, ilu0_solve, valid_ilu0

type(ilu0_factor) :: preconditioner
call ilu0(a,preconditioner,info,diagonal_shift=0.0_dp)
if (info == 0) call ilu0_solve(preconditioner,b,x,info)
```

`ilu0` is a deterministic **zero-fill** incomplete LU decomposition: it
stores strict lower-triangular `L` entries (implicit unit diagonal) and
upper-triangular `U` entries (including pivots) on the original canonical CSR
pattern. `preconditioner%lu` holds those entries, and
`preconditioner%diag_ptr` stores the one-based location of each pivot.
The factorization allocates `O(nnz+n)` storage and does not introduce
new nonzeros, even when they would improve approximation quality.
`ilu0_solve` uses forward/backward triangular substitutions.

- This is **not** SciPy SuperLU's `spilu`, which implements a more general
  thresholded and pivoted incomplete factorization. No row or column
  permutations, pivoting, drop tolerances, or extra fill are available.
- Each row must have a structurally stored diagonal entry. Zero or near-zero
  pivots fail with `sparse_linalg_bad_diagonal=-3`. A nonnegative finite
  `diagonal_shift` is optionally added to the stored diagonal *before*
  incomplete factorization to help avoid breakdown, but changes the
  preconditioning matrix, not the original system.
- `ilu0_solve` accepts any finite RHS of the right length and returns
  `x(:)` allocated by the routine. `valid_ilu0` verifies factor storage,
  pivots, and one-based CSR invariants. Callers should check `info` before
  passing factors to a solver.
- Use `ilu=preconditioner` on `bicgstab` (right preconditioning) or `gmres`
  (left preconditioning). ILU and Jacobi are mutually exclusive. The supplied
  factors must correspond to the intended coefficient matrix: shape is
  validated, but numerical matrix identity is not checked.

The factorization and solver are new SciFort MIT-licensed code implementing
standard ILU(0) recurrences. No third-party library source is embedded.

## Tests

Standalone checked/optimized native tests:

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
```

Randomized reference checks against SciPy:

```bash
python tools/sparse_linalg_reference/compare_expm_multiply_scipy.py --mode both --cases 260
python tools/sparse_linalg_reference/compare_sparse_linalg_scipy.py --mode both --cases 240
python tools/sparse_linalg_reference/compare_lsqr_ilu0_scipy.py --mode both --cases 320
python tools/sparse_linalg_reference/compare_minres_scipy.py --mode both --cases 220
python tools/sparse_linalg_reference/compare_eigsh_scipy.py --mode both --cases 120
python tools/sparse_linalg_reference/compare_eigs_scipy.py --mode both --cases 120
python tools/sparse_linalg_reference/compare_eigs_scipy.py --mode both --cases 45 --reduced
python tools/sparse_linalg_reference/compare_shiftinvert_scipy.py --mode both --cases 40
python tools/sparse_linalg_reference/compare_shiftinvert_modes_scipy.py --mode both --cases 84
python tools/sparse_linalg_reference/compare_eigs_shift_scipy.py --mode both --cases 90
python tools/sparse_linalg_reference/compare_eigs_shift_scipy.py --smallest --mode both --cases 65
python tools/sparse_linalg_reference/compare_arnoldi_scipy.py --mode both --cases 120
python tools/sparse_linalg_reference/compare_arnoldi_block_scipy.py --mode both --cases 108
python tools/sparse_linalg_reference/compare_lanczos_scipy.py --mode both --cases 72
python tools/sparse_linalg_reference/compare_lanczos_locked_scipy.py --mode both --cases 80
python tools/sparse_linalg_reference/compare_lanczos_locked_shift_scipy.py --mode both --cases 36
python tools/sparse_linalg_reference/compare_restarted_shift_scipy.py --mode both --cases 48
```

A detailed record of this update's tests is in
[`sparse_shiftinvert_validation_2026-11-04.md`](sparse_shiftinvert_validation_2026-11-04.md)
[`sparse_arnoldi_validation_2026-11-04.md`](sparse_arnoldi_validation_2026-11-04.md),
and [`sparse_lanczos_validation_2026-11-04.md`](sparse_lanczos_validation_2026-11-04.md) and
[`sparse_restarted_shift_validation_2026-11-04.md`](sparse_restarted_shift_validation_2026-11-04.md).

The tests compare solution vectors and **true residuals**, not just status
codes. They include nonsymmetric and symmetric positive definite systems,
rectangular overdetermined systems, explicit stored zeros, duplicate COO
entries, damped LSMR, invalid inputs, singular breakdowns, an 800-row
tridiagonal matrix, and new LSQR/ILU(0) tests including 1,500-row sparse preconditioning.

## Scope and future work

Not yet supported: arbitrary linear-operator callbacks in all public sparse
interfaces; complex sparse systems; AMG/SSOR preconditioners; SuperLU-style
column ordering/equilibration; or a full ARPACK-compatible implicitly
restarted or generalized sparse eigensolver. Earlier SciFort optimize modules do contain additional
matrix-free Krylov methods, but these are not yet exposed through this CSR
interface.

The [solver-library reuse survey](sparse_reuse_repositories.md) includes
[`ivan-pi/kiss`](https://github.com/ivan-pi/kiss), which is BSD-3-Clause but
links BLAS and LAPACK in its present FPM manifest. That candidate has **not**
been copied or included as a dependency. Future reuse should inspect the
specific upstream files and their license/provenance separately.

## Matrix norm estimation and sparse exponential

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense, csr_to_dense
use scifort_sparse_linalg, only : onenormest, expm

real(dp) :: d(2,2), estimate
real(dp), allocatable :: v(:), w(:), exponent(:,:)
type(csr_matrix) :: a, e
integer :: info

d = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
a = csr_from_dense(d)
estimate = onenormest(a, t=2, itmax=5, v=v, w=w, info=info)
e = expm(a, info, block_size=8)
exponent = csr_to_dense(e)
```

- `onenormest(a, t, itmax, v, w, info)` evaluates a **lower bound** on
  `||A||_1` for a real, square, canonical one-based CSR matrix. `t` defaults
  to 2 probe vectors and `itmax` to 5 iterations, both positive. The
  deterministic block sign/adjoint search avoids dependence on global RNG
  state. When `t >= n`, it evaluates the exact induced 1-norm from sparse
  column sums in O(nnz) time. Returned `v` satisfies `sum(abs(v))=1`,
  `w=matmul(A,v)`, and `sum(abs(w))=estimate` (up to floating-point error).
  If no status is supplied, inspect the returned value only for valid input.
  Empty 0-by-0 inputs return norm zero and empty witnesses. In contrast with
  SciPy, the current Fortran interface accepts only explicit CSR matrices;
  a matrix-free LinearOperator variant is not implemented.
- `expm(a, info, drop_tol, block_size, tol)` computes and **returns** a
  canonical CSR matrix approximating `exp(A)`. The input must be square.
  `block_size` defaults to 8, `drop_tol` defaults to zero, and `tol` uses
  `expm_multiply`'s default scaled-Taylor error target. Positive `drop_tol`
  omits output elements with absolute magnitude at most that value, which
  intentionally changes the mathematical exponential. Status is zero on
  success, -1 for invalid inputs, -2 for numerical overflow/breakdown, and
  positive for an unconverged Taylor series. A failure returns an unallocated
  CSR result. Large or poorly conditioned inputs may require substantial
  work. The full result can be dense even when `A` is sparse; sparse storage
  does not eliminate mathematical fill-in.

**Implementation distinction:** SciPy's sparse `expm` uses a Padé
scaling-and-squaring method; SciFort uses the existing trace-shifted,
scaled-Taylor `expm_multiply` on identity blocks, optionally retaining only
entries larger than `drop_tol`. No intermediate n-by-n dense array is
required during assembly, except when the returned matrix itself has dense
fill-in. The diagonal shortcut requires only O(n+nnz) time and O(n) storage.

**Reproduce checks:**

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_onenorm_expm_scipy.py --mode both --cases 240
```

See `docs/sparse_onenorm_expm_validation_2026-11-04.md` for the numerical
comparison methodology and known scope limitations.
