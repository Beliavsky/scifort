# Randomized interpolative SVD

Module: `scifort_linalg_interpolative_svd_random`

SciPy analogue: [`scipy.linalg.interpolative.svd(A, eps_or_k, rand=True)`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.svd.html).

This **new**, independent module supplements the deterministic
`scifort_linalg_interpolative_svd`; no existing source module or interface
has been modified.

## Interfaces

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_interpolative_svd_random, only : randomized_interpolative_svd

! Fixed-rank; caller allocates U(m,k), s(k), V(n,k):
call randomized_interpolative_svd(a, k, u, s, v, status, &
    seed=12345_int64, oversampling=8)

! Adaptive; U, s, V must be allocatable and are allocated by the routine:
call randomized_interpolative_svd(a, eps, k, u, s, v, status, &
    seed=12345_int64, oversampling=8)
```

The same generic interfaces accept `real(dp)` and `complex(dp)` matrices.
The output `s` is always `real(dp)`. All outputs represent the **thin** SVD,

\[
A \approx U\,\operatorname{diag}(s)\,V^{H},
\]

with `U(m,k)`, `s(k)`, and `V(n,k)`; `V` is **not** the conjugate transpose.
Singular values are nonnegative and sorted from largest to smallest.
The integer `k` selects the fixed-rank overload. A `real(dp)` tolerance
`eps` selects the adaptive overload and returns `k` as an output.

Optional arguments:

- `seed`: `integer(int64)`, local reproducible seed. A deterministic default
  seed of `20261110` is used when omitted. This does not affect the global
  Fortran random-number state. It does **not** reproduce NumPy's `rng` stream.
- `oversampling`: nonnegative integer, default `8`. In a fixed-rank call,
  the sketch has at most `k + oversampling` rows. In a tolerance call, the
  algorithm starts with `max(8, 2*max(1,oversampling))` rows, capped at
  `min(m,n)`, and grows the sketch if necessary.

`status` is one of the existing `linalg_status_*` codes. The routines reject
wrong-sized fixed-rank outputs, invalid ranks or tolerances, and negative
oversampling; nonfinite input produces `linalg_status_no_convergence`.
Empty matrices, rank zero, exact rank deficiency, and all-zero inputs are
supported. For adaptive input with an empty dimension or all-zero entries,
`k=0`, `U(m,0)`, `s(0)`, `V(n,0)` are returned.

## Algorithm

1. Scale `A` by its largest entry magnitude to moderate numerical range.
2. Generate a **Gaussian row sketch**, `Y = G A`, with `l << m` rows when
   the requested rank is small. Entries of `G` come from a local seeded
   Park–Miller generator and a Box–Muller transform.
3. Apply SciFort's column-pivoted Householder QR to `Y`. The chosen column
   indices refer to the **original matrix** `A` (one-based).
4. Compute the interpolative coefficients by least-squares fitting **on the
   original scaled A**, not just on the sketch. This correction improves
   approximation accuracy and handles rank-deficient selected columns.
5. Convert the resulting interpolative decomposition to an SVD with the
   existing SciFort QR and Jacobi SVD routines.

For tolerance-based calls, the algorithm computes the **actual relative
Frobenius residual** of the interpolative approximation and grows the sketch
until the requested tolerance is achieved or the full possible sketch size
is reached. At the full sketch size it can use all `min(m,n)` columns;
small rounding errors in the returned SVD are unavoidable. The fixed-rank
interface estimates a rank-k approximation but **does not promise an error
tolerance**. With `k=0`, the approximation is zero.

The input is dense; matrix-free `LinearOperator` and sparse interfaces are
not implemented in this module. Batch dimensions are not currently supported.

## Differences from SciPy

- SciPy uses a randomized ID scheme and a spectral-norm accuracy target;
  SciFort uses a Gaussian row sketch with a verified **Frobenius** residual
  for the adaptive interface. The algorithms' pivot choices and selected
  ranks need not match.
- SciFort accepts a reproducible local `integer(int64)` seed, not a NumPy
  `Generator`. The default seed is deterministic (SciPy normally draws a
  fresh generator).
- This method can be faster than QR of the entire input for large,
  numerically low-rank matrices, but no general speedup is guaranteed; the
  original-matrix least-squares correction costs additional work.
- No code from SciPy or its bundled Fortran libraries is copied or linked.

## Testing

With GNU Fortran and the dependencies built in this order:

```text
src/kinds.f90
src/linalg.f90
src/linalg_svd.f90
src/linalg_qr.f90
src/linalg_interpolative.f90
src/linalg_interpolative_svd_random.f90
```

Compile `test/test_linalg_interpolative_svd_random.f90` against those objects.
Compile `validation/reference_drivers/driver_linalg_interpolative_svd_random.f90` the same way, then
run

```bash
python test/compare_linalg_interpolative_svd_random_scipy.py ./random_svd_driver
```

A small example is in `example/linalg_interpolative_svd_random_example.f90`.
See the accompanying [validation report](linalg_interpolative_svd_random_validation_2026-11-10.md).
