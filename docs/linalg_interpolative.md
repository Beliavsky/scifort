# Interpolative matrix decomposition (`scipy.linalg.interpolative`)

`scifort_linalg_interpolative` provides deterministic real (`real(dp)`) and
complex (`complex(dp)`) interpolative decompositions and reconstruction
helpers. It uses the SciFort Householder column-pivoted QR and Jacobi SVD;
no system LAPACK or additional dependency is needed.

For an m-by-n matrix A and selected rank k, the column-pivoted ID has the form

```text
A(:,idx) ≈ B * [I_k, T]
B = A(:,idx(1:k))
T = proj(k,n-k)
A ≈ B * P,   where P(:,idx) = [I_k, T]
```

`idx` uses **one-based** column indices, unlike Python's zero-based indices.
The skeleton B is made of *actual columns of A*, not constructed basis vectors.
The interpolation coefficients T are obtained from the leading triangular
block of the pivoted QR, with a minimum-norm SVD fallback if that block is
numerically singular.

## Fortran interfaces

The real and complex variants share the same generic names. `idx` and `proj`
are **allocatable, intent(out)** and are allocated by the call.

```fortran
use scifort_kinds, only : dp
use scifort_linalg_interpolative, only : interp_decomp, &
    reconstruct_skel_matrix, reconstruct_interp_matrix, &
    reconstruct_matrix_from_id, id_to_svd

real(dp) :: a(20,12)
real(dp), allocatable :: proj(:,:), b(:,:), p(:,:), approx(:,:)
real(dp), allocatable :: u(:,:), s(:), v(:,:)
integer, allocatable :: idx(:)
integer :: status, k

! ... initialize a ...
call interp_decomp(a, 5, idx, proj, status)     ! fixed rank = 5
! Or:
call interp_decomp(a, 1.0e-7_dp, k, idx, proj, status) ! adaptive rank

k = size(proj,1)  ! for fixed-rank calls; adaptive calls return k
allocate(b(20,k),p(k,12),approx(20,12),u(20,k),s(k),v(12,k))
call reconstruct_skel_matrix(a,k,idx,b,status)
call reconstruct_interp_matrix(idx,proj,p,status)
call reconstruct_matrix_from_id(b,idx,proj,approx,status)
call id_to_svd(b,idx,proj,u,s,v,status)
```

All calls return a SciFort `linalg_status_*` code:

- `linalg_status_success`: success, including zero-sized or zero matrices.
- `linalg_status_invalid_shape`: invalid rank/tolerance or output dimensions,
  or a non-permutation index vector.
- `linalg_status_no_convergence`: nonfinite input, numerical failure in QR,
  the fallback solve, or small-core SVD, or nonfinite result.

The caller must allocate the arrays for the four construction/SVD helpers;
`interp_decomp` allocates only its `idx` and `proj` outputs. In the rare event
of a numerical failure, check `status` before using the outputs.

### `interp_decomp(a, k, idx, proj, status)`

Produces a rank-k factorization with `k` an integer satisfying
`0 <= k <= min(m,n)`. For `k=0`, the reconstructed approximation is zero.
If the selected columns are rank deficient, the interpolation uses a
minimum-norm SVD solution.

### `interp_decomp(a, eps, k, idx, proj, status)`

Chooses k using a deterministic **relative Frobenius-norm** tolerance
`0 <= eps < 1`. By construction, the trailing QR block satisfies
`||R_22||_F <= eps ||A||_F` in exact arithmetic. Because
`A - B P = Q_2 R_22` for the ID of the full pivoted QR, this is an error
bound for the reconstructed approximation (subject to floating-point error).

This differs from SciPy's random-sampling default and its relative
**spectral**-norm precision policy. Ranks need not match SciPy's; compare
errors and reconstructions rather than demanding identical rank or pivot
orders on tied or nearly tied columns.

### Reconstruction helpers

- `reconstruct_skel_matrix(a, k, idx, b, status)` requires B of shape (m,k).
- `reconstruct_interp_matrix(idx, proj, p, status)` requires P of shape (k,n).
- `reconstruct_matrix_from_id(b, idx, proj, a, status)` requires output of
  shape (m,n) and avoids explicitly constructing P.
- `id_to_svd(b, idx, proj, u, s, v, status)` returns U(m,k), s(k), V(n,k),
  with `Aapprox = U * diag(s) * V**H`. For real inputs `V**H` is `transpose(V)`;
  for complex inputs it is `transpose(conjg(V))`.

`id_to_svd` computes QR(B), followed by SVD of the k-by-n core; it never
constructs the full reconstructed m-by-n matrix. This can be cheaper than
an SVD of the original matrix when k is small.

## Scope and limitations

- Deterministic dense matrices only; randomized methods, `LinearOperator`
  callbacks, and higher-rank batch arrays are not yet implemented.
- The fixed-rank QR currently factorizes all columns even if k is small;
  a future truncated factorization could improve large-matrix performance.
- Inputs are double-precision real or complex; indices are default Fortran
  integers and one-based.
- Zero-sized matrices are valid. SciPy 1.17.0 can have different behavior for
  some degenerate or out-of-range requests; this implementation checks them.

## Tests

- `test/test_linalg_interpolative.f90`: native reconstruction, low-rank,
  tolerance, ID-to-SVD, error-path, real and complex cases.
- `validation/reference_drivers/driver_linalg_interpolative.f90` plus
  `test/compare_linalg_interpolative_scipy.py`: independent reference run
  against SciPy's **deterministic** `interp_decomp(..., rand=False)` and
  `id_to_svd`. Numerically valid decompositions can have different pivot
  orders, so the checker compares approximation errors as well as matched
  interpolation coefficients.
- `example/linalg_interpolative_example.f90`: executable usage illustration.
- `docs/linalg_interpolative_validation_2026-11-09.md`: validation commands
  and recorded results.

SciPy API reference:
<https://docs.scipy.org/doc/scipy/reference/linalg.interpolative.html>.
