# Cosine–sine decomposition (`scipy.linalg.cossin`)

SciFort's `scifort_linalg_cossin` module implements the **complete**
cosine–sine decomposition for real orthogonal and complex unitary matrices:

\[
X = U\,CS\,V^H,
\]

where `X` is `m × m`, `U` is block diagonal with blocks of orders `p` and
`m-p`, and `V^H` is block diagonal with blocks of orders `q` and `m-q`.
The real-valued middle factor `CS` has `r = min(p,m-p,q,m-q)` cosines and
sines. Both the default SciPy sign arrangement and `swap_sign=.true.` are
supported, including unequal partitions.

This is an independent MIT-licensed implementation based on SciFort's
existing one-sided Jacobi SVD. **No new dependency on LAPACK or BLAS** is
introduced. Unlike LAPACK `ORCSD/UNCSD`, it does not attempt to match SciPy's
specific signs, phases, or bases inside degenerate singular subspaces;
compare `CS` and `theta` directly, and compare factors through
reconstruction/orthogonality.

## Public APIs

```fortran
use scifort_kinds, only : dp
use scifort_linalg_cossin, only : cossin, cossin_separate

! Real: X = U * CS * VT
call cossin(x, p, q, u, cs, vt, status, swap_sign=.false.)

! Complex: X = U * CS * VH
call cossin(x_complex, p, q, u_complex, cs, vh, status)

! Separate real/complex subfactors and r principal angles:
call cossin_separate(x, p, q, u1, u2, theta, v1t, v2t, status)
```

The input and output types are `real(dp)` or `complex(dp)` (matching on each
call). `p` and `q` are **integers**, with `1 <= p,q < m`. The output `cs` is
always `real(dp)`, including for complex input. `u` and `vt/vh` are `m × m`.
`u1` is `p × p`, `u2` is `(m-p) × (m-p)`, `v1t/v1h` is `q × q`, `v2t/v2h`
is `(m-q) × (m-q)`, and `theta` has length `r`.

All outputs are **caller-allocated**. Status is
`linalg_status_success` (0), `linalg_status_invalid_shape` (1), or
`linalg_status_no_convergence` (4). Nonfinite and appreciably nonunitary
inputs return status 4. This unitary-input validation is deliberately stricter
than SciPy's `cossin`, which delegates directly to LAPACK. Invalid shapes
return 1. Do not interpret output factors following an unsuccessful call.

For the full-factor interface, `compute_u=.false.` and
`compute_vh=.false.` each suppress the corresponding output. The disabled
array can be either `m × m` (filled with zero) or `0 × 0`; in the latter
case it mirrors SciPy's empty return. The separate-factor interface always
computes and returns the four square subfactors. The sign switch
`swap_sign=.true.` affects only the gauge convention, not the CS angles.

### Compatibility limitations

- Operates on one matrix per call; SciPy's leading batch dimensions are
  not yet exposed in this Fortran interface.
- The SciPy form accepting four separate input blocks is not yet implemented;
  concatenate the blocks into one square matrix and supply `p` and `q`.
- `compute_u/compute_vh` avoid returning factors, but **do not currently skip**
  their internal computation. The method uses a full SVD and orthogonal
  completions; its complexity is roughly cubic in matrix order.
- This is independent SVD-based CSD rather than a literal port of LAPACK's
  more specialized `ORCSD/UNCSD`. Signs/phases of nonunique factors can differ.

## Running native tests and independent SciPy comparison

With GNU Fortran and SciPy installed, from the project root:

```sh
mkdir -p /tmp/scifort-cossin
cd /tmp/scifort-cossin
gfortran -O0 -g -fcheck=all -c \
  /path/to/scifort/src/kinds.f90 \
  /path/to/scifort/src/linalg.f90 \
  /path/to/scifort/src/linalg_svd.f90 \
  /path/to/scifort/src/linalg_cossin.f90
gfortran -O0 -g -fcheck=all -o test_cossin *.o \
  /path/to/scifort/test/test_linalg_cossin.f90
./test_cossin
gfortran -O0 -o driver_cossin *.o \
  /path/to/scifort/tools/linalg_reference/driver_cossin.f90
python /path/to/scifort/tools/linalg_reference/compare_cossin_scipy.py ./driver_cossin
```

`example/linalg_cossin_example.f90` is a minimal runnable demonstration.
See [the validation notes](linalg_cossin_validation_2026-11-09.md).
