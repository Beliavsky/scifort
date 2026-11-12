# Value-selected symmetric/Hermitian eigenproblems

`scifort_linalg_dense_eigen.eigh` and `eigvalsh` now accept a SciPy-compatible
`subset_by_value=[low,high]` option for both real symmetric and complex
Hermitian matrices. It is supported for ordinary eigenproblems and for
positive-definite generalized eigenproblems with `b` and `problem_type=1|2|3`.

## Half-open eigenvalue intervals

SciPy's convention is **`(low, high]`**: exclude `low`, include `high`.
Repeated eigenvalues count with multiplicity. IEEE positive and negative
infinities may be used as bounds. Equal or reversed bounds and NaN bounds
return `linalg_status_invalid_shape`. Do **not** specify both
`subset_by_value` and `subset_by_index` in one call.

The two new optional arguments are *trailing*, so existing positional calls
and the one-based, inclusive `subset_by_index` behavior are unchanged:

```fortran
call eigh(a, w, status, z, subset_by_value=[1.0_dp,3.0_dp], &
    num_selected=nfound)
call eigvalsh(a, w, status, b=b, problem_type=2, &
    subset_by_value=[1.0_dp,3.0_dp], num_selected=nfound)
```

### Caller-owned output buffers

The number of eigenvalues selected by a numerical interval is not known before
solving. Consistent with SciFort's established native Fortran API, `w` and the
optional eigenvector array `z` are **preallocated capacity buffers**:

- `w(:)` has capacity for some or all eigenvalues, usually `n`.
- `z(n,:)`, if requested, has capacity for the corresponding eigenvectors,
  usually `n` columns. `z` and `w` may have different capacities.
- `num_selected` (optional integer output) reports the true count `m`, even
  when the buffers prove too small after eigenvalue computation.
- On success, use `w(1:m)` and `z(:,1:m)` only. The unused entries and
  columns are initialized to zero. A zero-count selection succeeds, including
  when supplied with zero-length output arrays.
- If either supplied capacity is less than `m`, the routine returns
  `linalg_status_invalid_shape`. You may allocate larger buffers and retry.
- Without `subset_by_value`, outputs retain their original **exact-shape**
  requirements. `num_selected` reports the count for default and index modes
  too, if requested.

Example for a real symmetric matrix with `n=size(a,1)`:

```fortran
real(dp), allocatable :: w(:), z(:,:)
integer :: n, nfound, status

n = size(a, 1)
allocate(w(n), z(n,n))
call eigh(a, w, status, z, subset_by_value=[-2.0_dp,2.0_dp], &
    num_selected=nfound)
if (status /= linalg_status_success) error stop 'eigenproblem failed'
print *, w(1:nfound)
```

For generalized problems, `b` must be positive definite and use the same
real or complex kind as `a`. The returned eigenvectors preserve the previous
normalization: `z^H B z=I` for types 1 and 2, and `z^H B^-1 z=I` for type 3.
The `lower` option still governs which triangle of both coefficient matrices
is read. Invalid or unused triangles are not inspected.

## Algorithm and numerical qualifications

The current implementation computes the **entire** eigenvalue spectrum with
Householder reduction and tridiagonal QL iteration, then selects the entries
that satisfy the interval. If eigenvectors are requested, the full set of
intermediate eigenvectors is computed before slicing; generalized problems
still use a Cholesky reduction. Therefore `subset_by_value` does **not**
reduce the O(n^3) time or O(n^2) working storage of the original solver.

Comparisons of eigenvalues precisely on an interval endpoint can be sensitive
to rounding or to different LAPACK/QL algorithms, particularly for a repeated
or tightly clustered spectrum. Use a gap between eigenvalue clusters when a
numerical selection must be robust to such differences. This is distinct
from the exact half-open interval contract for the computed eigenvalues.

## Reproducing numerical validation

```sh
python tools/linalg_reference/run_eigen_by_value_checks.py \
  --mode both --cases 320 --max-order 48 --seed 1907
```

The script compiles native tests and a batch driver with GNU Fortran, then
compares ordinary and generalized `eigh`/`eigvalsh` against SciPy for real and
complex matrices, zero and repeated spectra, selected eigenvectors,
noncontiguous intervals, both triangles, scaled inputs, and all three
generalized types. It verifies eigenpair residuals and, where appropriate,
`B`/`B^-1` normalization. See the accompanying validation note dated
2026-11-08 for measured errors and compiler versions.

The implementation, tests, and examples are independently authored MIT-licensed
SciFort work. SciPy is used as a numerical reference and for public API
semantics, not as a source of copied implementation code.
