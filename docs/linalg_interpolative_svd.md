# Interpolative SVD (`scipy.linalg.interpolative.svd`)

The new `scifort_linalg_interpolative_svd` module computes a thin singular
value decomposition through a **deterministic column-interpolative
decomposition** (ID). It reuses `scifort_linalg_interpolative`, Householder
column-pivoted QR, and the Jacobi SVD in SciFort without a new dependency.

Given `A(m,n)` and rank `k`, it obtains selected columns `B=A(:,idx(1:k))`
and interpolation coefficients `P`, with `A ≈ B*P`. The existing
`id_to_svd` then computes QR(B), an SVD of the `k`-by-`n` core, and the thin
factors

```text
A ≈ U * diag(s) * transpose(V)        ! real
A ≈ U * diag(s) * transpose(conjg(V)) ! complex
U(m,k), s(k), V(n,k)
```

## API

The generic name is `interpolative_svd`. There are real and complex variants
of both calls. Every argument is explicit and has `intent`; inputs are
`real(dp)` or `complex(dp)`. Unlike `scifort_linalg_svd.svd`, the final right
factor here is **V, not V^H**, matching SciPy's interpolative SVD interface.

**Fixed rank** (output arrays allocated by the caller):

```fortran
use scifort_kinds, only : dp
use scifort_linalg_interpolative_svd, only : interpolative_svd
real(dp) :: a(50,30), u(50,5), s(5), v(30,5)
integer :: status
! ... fill a ...
call interpolative_svd(a, 5, u, s, v, status)
```

**Relative tolerance** (the call allocates the factors and returns `k`):

```fortran
complex(dp) :: a(50,30)
complex(dp), allocatable :: u(:,:), v(:,:)
real(dp), allocatable :: s(:)
integer :: k, status
! ... fill a ...
call interpolative_svd(a, 1.0e-6_dp, k, u, s, v, status)
```

- Fixed ranks must satisfy `0 <= k <= min(m,n)`.
- Tolerances must be finite and satisfy `0 <= eps < 1`.
- On success, `s` is nonnegative and nonincreasing, and the columns of U
  and V are orthonormal to floating-point accuracy.
- Valid `k=0` returns zero-column factors; empty matrices are supported.
- Invalid shapes or tolerances return `linalg_status_invalid_shape`; invalid
  matrix values or numerical failures return `linalg_status_no_convergence`.
- Fixed-rank outputs are initialized to zero before returning errors.
  Adaptive outputs have zero size if the call fails.

## SciPy differences and performance considerations

This corresponds to **`scipy.linalg.interpolative.svd(A, ..., rand=False)`**.
SciPy defaults to **randomized** sampling (`rand=True`), while the new
SciFort implementation is **deterministic only**. No `LinearOperator`,
randomized sampling, or batched-matrix interface is provided yet.

The tolerance variant inherits SciFort's relative **Frobenius norm** target
`||A - Aapprox||_F / ||A||_F <= eps` from `interp_decomp`, whereas SciPy's
precision target uses the spectral norm; selected ranks may differ.
For fixed rank, valid decompositions can select different pivot columns,
and corresponding factors may differ by signs/phases.

The ID uses the existing **full pivoted QR** rather than a truncated or
randomized QR, so a fixed small rank does not yet deliver the usual
large-matrix acceleration of randomized low-rank SVD. It is nonetheless a
convenient complete ID-to-SVD API for clients of SciFort.

## Testing

- `test/test_linalg_interpolative_svd.f90` checks real/complex factor shapes,
  orthogonality, reconstruction of the independent ID, tolerance selection,
  rank deficiencies, zero/empty matrices, coefficient scaling, invalid
  dimensions/tolerances, and NaNs.
- `validation/reference_drivers/driver_linalg_interpolative_svd.f90` and
  `test/compare_linalg_interpolative_svd_scipy.py` run independent comparisons
  against SciPy 1.17's deterministic algorithm and NumPy reconstruction.
- `example/linalg_interpolative_svd_example.f90` demonstrates fixed-rank
  and adaptive-rank calls.
- Reproduction commands and results are in
  `docs/linalg_interpolative_svd_validation_2026-11-09.md`.

SciPy public documentation:
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.svd.html>.
