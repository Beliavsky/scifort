# QR rank updates, insertion and deletion (2026-11-09)

`scifort_linalg_qr_updates` adds `qr_update`, `qr_insert`, and `qr_delete`,
following the mathematical contracts of the corresponding SciPy routines.
The module depends only on existing SciFort QR, kinds and status definitions.

## API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_qr, only : qr
use scifort_linalg_qr_updates, only : qr_update, qr_insert, qr_delete

! Real symmetric update, A_new = Q R + u v^T:
call qr_update(q, r, u, v, qnew, rnew, status)

! Rank-k update, A_new = Q R + U V^T (real) or U V^H (complex):
call qr_update(q, r, U, V, qnew, rnew, status)

! Insert one row/column (rank-one u) or p rows/columns (rank-two u):
call qr_insert(q, r, u, k, qnew, rnew, status, which='row')
call qr_insert(q, r, U, k, qnew, rnew, status, which='col')

! Delete p consecutive rows/columns (p defaults to 1):
call qr_delete(q, r, k, qnew, rnew, status, p=2, which='col')
```

The overloaded entry points support `real(dp)` and `complex(dp)` inputs.
For complex input, `qr_update` uses **conjugate transpose** on `v`/`V`.
Mixed-precision and mixed-real/complex calls are not provided by these
interfaces; explicitly convert arguments to a common type if needed.

### Array shapes

For an original `m`-by-`n` matrix, `Q` is `m`-by-`h` and `R` is
`h`-by-`n`, where `h=m` for full QR or `h=min(m,n)` for economic QR.
`Q` must have orthonormal columns and `R` must be upper trapezoidal.
The routines do not independently verify these factor properties.
Caller-owned output arrays must have exactly the dimensions below:

| Operation | New matrix shape | Updated `hnew` |
|---|---|---|
| `qr_update` with `u(m,k)`, `v(n,k)` | `(m,n)` | `h` |
| `qr_insert` row with `u(p,n)` | `(m+p,n)` | `m+p` if full, else `min(m+p,n)` |
| `qr_insert` column with `u(m,p)` | `(m,n+p)` | `m` if full, else `min(m,n+p)` |
| `qr_delete` row, `p` removed | `(m-p,n)` | `m-p` if full, else `min(m-p,n)` |
| `qr_delete` column, `p` removed | `(m,n-p)` | `m` if full, else `min(m,n-p)` |

For a vector insertion, `p=1`. Output `Qnew` is `(mnew,hnew)` and
`Rnew` is `(hnew,nnew)`. When `Q` is square, a full factorization is
assumed, including for a wide matrix whose reduced and full shapes
coincide. This matches the observed SciPy behavior.

**Fortran positions are one-based**, unlike SciPy's zero-based `k`.
Insert at `k=1` to prepend or `k=m+1` (row) / `k=n+1` (column) to append.
Delete starting at `k` with `1 <= k <= dimension-p+1`.
`which='row'` is the default; `which='col'` selects columns.

### Status and edge cases

`status` is from `scifort_linalg`:

- `linalg_status_success = 0` on success.
- `linalg_status_invalid_shape = 1` for invalid dimensions, indices, `p`,
  `which`, or output shapes.
- `linalg_status_no_convergence = 4` for nonfinite inputs or numerical
  failure of a QR step.

Outputs are initialized to zero before validation. Rank-zero updates
and insertion of zero rows/columns preserve the supplied factors. Deleting
all rows or columns produces a mathematically well-defined empty QR
factorization. In SciPy 1.17.0, some such empty cases have unusual factor
shapes; the SciFort convention above is used consistently instead.

### Method and performance

`qr_update` constructs an orthonormal basis for the component of `U`
outside the current `Q` and refactors a reduced core of size at most
`(h+k)`-by-`n`. It reuses the factors and does **not** explicitly
reconstruct the original dense matrix. `qr_delete` on columns refactors
only `R` with the requested columns removed. Column insertion similarly
uses a reduced core, but falls back to a dense QR when inserted columns
are dependent and the augmented basis is too small for the requested
output shape. Row insertion and deletion use an explicit reconstructed
matrix and existing Householder QR for reliability; these operations are
not implemented as optimal-cost Givens updates or downdates.

Differences from SciPy: no `overwrite_*` flags, no `check_finite=False`,
no `rcond` option for economic column insertion, and no leading batch
axes. Inputs are checked for finiteness. Economic column insertion of a
dependent column is handled by fallback QR rather than SciPy's
condition-number-based rejection. Like SciPy, the factors may have
different valid signs/phases than those returned by a fresh QR.

## Tests and examples

- Native Fortran: `test/test_linalg_qr_updates.f90`.
- Independent reference driver: `reference/linalg_qr_updates_driver.f90`.
- Reference comparison: `validation/compare_linalg_qr_updates.py`.
- Executable example: `example/linalg_qr_updates_example.f90`.
- Validation results: `docs/linalg_qr_updates_validation_2026-11-09.md`.

To run the independent reference script (requires `gfortran`, NumPy and
SciPy), from the repository root:

```bash
python validation/compare_linalg_qr_updates.py
python validation/compare_linalg_qr_updates.py --optimized
```

These commands compile only the minimal QR dependency closure, without
fetching packages or using FPM.
