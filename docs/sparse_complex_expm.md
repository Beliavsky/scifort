# Complex sparse matrix exponential actions

`scifort_sparse_linalg_complex_expm` computes `exp(t*A)*B` for complex binary64
CSR and CSC matrices, as well as matrix-free operators, without assembling a
dense matrix exponential. The API follows the mathematical behavior of
`scipy.sparse.linalg.expm_multiply` (SciPy 1.17), with Fortran's explicit
`info` status codes and one-based sparse indices.

## Single-time API

```fortran
use scifort_sparse_linalg_complex_expm, only : expm_multiply_complex
! a is a canonical complex CSR or CSC matrix; b is a complex vector or block.
call expm_multiply_complex(a,b,y,info,t=0.5_dp,tol=1.0e-12_dp, &
    scaling_steps=steps,products=matvec_calls)
```

`b` may have shape `(n)` or `(n,nrhs)`, with matching allocatable result `y`.
The optional real `t` defaults to `1.0_dp`. The complex `tracea` optional
argument provides a possible shift. When a matrix is stored explicitly,
SciFort computes its exact trace and chooses the supplied shift only if it
reduces the shifted infinity norm. Explicit zero entries are supported.

The algorithm uses scaling with at most 200,000 substeps and 80 Taylor terms
per substep. It applies a scalar trace shift and verifies a geometric bound
on the remaining terms in the infinity norm, separately for each RHS. The
requested `tol` is a relative *Taylor-tail* target, not a global rigorous
bound on floating-point roundoff error. By default it is `4e-15`, floored at
`4*epsilon(1.0_dp)`. The result is computed with only matrix-vector
multiplications, complex saxpy operations, and scalar exponentials.

## Matrix-free API

```fortran
call expm_multiply_complex_operator(n,b,matvec,norm_bound,y,info, &
    t=0.5_dp,tracea=trace_of_a,products=matvec_calls)
```

The callback has the interface

```fortran
subroutine matvec(x,y,info)
    complex(dp), intent(in) :: x(:)  ! length n
    complex(dp), intent(out) :: y(:) ! A*x
    integer, intent(out) :: info     ! 0 for success
end subroutine
```

**Important:** `norm_bound` must be a finite, nonnegative *conservative upper
bound* on the induced infinity norm `||A||_infinity`. An underestimate can
invalidate the Taylor-tail convergence guarantee. Matrix-free operations do
not estimate the trace automatically; pass `tracea` when known. The bound on
the shifted operator is `norm_bound + abs(tracea/n)`, which is conservative
but can require more substeps than the exact CSR/CSC norm.

## Time intervals

```fortran
call expm_multiply_complex_interval(a,b,-1.0_dp,1.0_dp,9,y,info, &
    endpoint=.true.,products=matvec_calls)
call expm_multiply_complex_operator_interval(n,b,matvec,norm_bound, &
    -1.0_dp,1.0_dp,9,y,info,tracea=trace_of_a)
```

Vector input produces `y(n,num)`; block input produces `y(n,nrhs,num)`.
The evenly spaced grid includes `stop` unless `endpoint=.false.`; with
`num=1`, only `start` is evaluated. The implementation starts at the
requested point closest to time zero and propagates forward/backward.
This avoids irreversible underflow when the first requested time is far
from zero on a decaying exponential. CSR/CSC and operator interfaces support
both vector and block RHS inputs.

## Status codes

- `complex_expm_success = 0`: action completed successfully.
- `complex_expm_no_convergence = 1`: Taylor budget exhausted.
- `complex_expm_invalid_input = -1`: invalid shape/argument or excessive scale.
- `complex_expm_breakdown = -2`: arithmetic generated a nonfinite result.
- `complex_expm_callback_failure = -3`: callback returned nonzero status.

All input values must be finite. The output is allocated with zero extent on
validation failure and may hold a partial value on computational failure.
`products` counts individual matrix-vector actions, not time steps, and
`scaling_steps` counts Taylor substeps for a single-time action.

## Validation and limitations

Run the native test (also registered in `tools/sparse_reference/run_native_sparse.py`):

```bash
python tools/sparse_reference/run_native_sparse.py --mode checked
```

Run the independent reference comparisons:

```bash
python tools/sparse_reference/compare_sparse_complex_expm_scipy.py --mode both --cases 180
```

These compare against `scipy.linalg.expm` and `scipy.sparse.linalg.expm_multiply`
on small general, triangular, diagonal, singular, Hermitian, and nonnormal
complex systems, with positive and negative times and both interval endpoint
conventions. Numerical comparisons are reference tests, not copied SciPy code.
For matrices with large exponential growth, SciPy's interval implementation
can differ from the independently evaluated dense exponential by around
`1e-8` relative in some tests; the dense exponential is used as the primary
reference. The implementation does not expose SciPy's trace estimation,
nonuniform time grids, Krylov matrix-function callbacks, or dense `expm`.
