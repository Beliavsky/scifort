# Complex sparse `funm_multiply_krylov`

`scifort_sparse_linalg_complex_funm` evaluates the **action** of a matrix
function on a complex vector or block,

\[
    Y = f(tA) B,
\]

without constructing the dense full-size `f(tA)`. This extends SciFort's real
`funm_multiply_krylov` to complex binary64 arithmetic and follows the general
interface of `scipy.sparse.linalg.funm_multiply_krylov` in SciPy 1.17.
The implementation uses one-based sparse indices and Fortran `info` return codes.

## Explicit sparse interface

```fortran
use scifort_kinds, only : dp
use scifort_sparse_linalg_complex_funm, only : funm_multiply_krylov_complex

! a: canonical complex CSR or CSC matrix; b: complex vector or RHS block.
! The callback evaluates f(H) on a small *dense complex* matrix H.
call funm_multiply_krylov_complex(a,b,my_dense_function,y,info, &
    t=cmplx(0.5_dp,0.25_dp,dp), restart_every_m=12, max_restarts=10, &
    rtol=1.0e-8_dp,cycles=cycles,products=products,update_norm=change)
```

The `a` overloads accept `csr_complex_matrix` or `csc_complex_matrix`.
Vector `b(n)` produces allocated `y(n)`; block `b(n,nrhs)` produces allocated
`y(n,nrhs)`. Each block column uses its own Krylov space; the returned
`cycles` and `products` are summed across columns, and `update_norm` is the
largest final correction magnitude. CSC input is converted to CSR internally.

The optional **complex** `t` defaults to `cmplx(1.0_dp,0.0_dp,dp)`.
`restart_every_m` defaults to `min(20,n)`; `max_restarts` defaults to `20`.
`rtol` defaults to `1e-6` and `atol` defaults to zero. A complete basis of
size `n` is limited to a single cycle.

## Matrix-free operator and dense-function callback

```fortran
call funm_multiply_krylov_complex_operator(n,b,apply_a,my_dense_function,y,info, &
    restart_every_m=12, max_restarts=10)
```

The two user callbacks have these interfaces:

```fortran
subroutine apply_a(x,ax,info)
    complex(dp), intent(in)  :: x(:)     ! length n
    complex(dp), intent(out) :: ax(:)    ! A*x, length n
    integer, intent(out)     :: info     ! 0 on success
end subroutine

subroutine my_dense_function(h,fh,info)
    complex(dp), intent(in)  :: h(:,:)   ! small square projected matrix
    complex(dp), intent(out) :: fh(:,:)  ! f(h), same shape
    integer, intent(out)     :: info     ! 0 on success
end subroutine
```

A caller may use a small dense matrix function from an independent library,
or implement a polynomial, exponential, trigonometric, or other *analytic*
function in `my_dense_function`. SciFort supplies the Krylov projection, **not**
a general-purpose complex dense `funm` routine. The callback must return the
**entire** dense `f(H)`, not merely its eigenvalues or a first-column action.

The matrix-free path does not require a norm estimate or an adjoint `A^H`
callback. `assume_hermitian=.true.` selects the faster three-term complex
Lanczos recurrence (real tridiagonal projection); users must guarantee that
**A** is Hermitian. The default general path uses two-pass modified
Gram–Schmidt Arnoldi with conjugate inner products.

## Convergence and status

The restarted algorithm constructs lower-block-Hessenberg projections and
checks the norm of each correction to `f(tA)b`. The stopping threshold is
`max(atol,rtol*norm2(b))`. A numerically invariant Krylov subspace also ends
successfully. This is **an update-based heuristic, not a rigorous error bound**.
The method is intended for analytic matrix functions, and convergence can be
poor for strongly nonnormal operators or difficult functions. For the
exponential alone, `scifort_sparse_linalg_complex_expm` provides a separate
Taylor-remainder-controlled action implementation.

- `complex_funm_success = 0`: completed by full/invariant basis or small update.
- `complex_funm_no_convergence = 1`: iteration budget exhausted (partial output).
- `complex_funm_invalid_input = -1`: invalid arguments or shapes.
- `complex_funm_breakdown = -2`: nonfinite input to a numerical operation or
  nonfinite operator/projected-function output.
- `complex_funm_callback_failure = -3`: a caller-supplied callback failed.

Nonfinite matrix entries, right-hand sides, tolerances and time scales are
rejected. A zero RHS is returned as zero without invoking either callback;
rank-zero operators and zero-column blocks are handled. On invalid input,
`y` has zero extent; for numerical failure, it may contain a partial result.
Large iteration budgets cause the projected matrix storage to grow on demand,
but these matrices still cost quadratic memory and cubic-or-worse callback
work in the *accumulated* Krylov dimension. Excessive `m` or restart counts
are not appropriate for production-scale matrix-free use.

## Validation

```bash
python tools/sparse_reference/compare_sparse_complex_funm_scipy.py --mode both --cases 150
python tools/sparse_reference/run_native_sparse.py --mode checked
```

The tests compare CSR, CSC, and matrix-free implementations for vectors and
blocks with independent `scipy.linalg.expm`, `scipy.linalg.sinm`, and polynomial
reference actions. Cases include Hermitian and non-Hermitian operators,
triangular/defective matrices, explicit zeros, diagonal matrices, zero input,
complex scale factors, restarts, invalid callbacks, and a 1,024-dimensional
matrix-free problem. SciPy's public `funm_multiply_krylov` is also cross-checked
on supported cases. All source in this module is original MIT-licensed
SciFort code; no SciPy implementation was copied.
