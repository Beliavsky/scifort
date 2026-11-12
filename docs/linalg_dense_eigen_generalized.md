# Generalized symmetric/Hermitian definite eigenproblems

`scifort_linalg_dense_eigen` supports the optional SciPy-like `b` matrix for
real symmetric and complex Hermitian `eigh` and `eigvalsh`:

```fortran
use scifort_linalg_dense_eigen, only : eigh, eigvalsh

call eigh(a, w, status, z, b=b, problem_type=1)
call eigvalsh(a, w, status, b=b, lower=.false., subset_by_index=[2,5])
```

The **existing positional Fortran arguments are unchanged**. `b` and
`problem_type` are optional trailing arguments; use the keywords for clarity.
Omitting `b` invokes the original standard eigenproblem solver, with no
change in its mathematical behavior. `problem_type` is ignored without `b`,
following SciPy's corresponding parameter contract.

Supported generalized problems, with real eigenvalues in ascending order:

| `problem_type` | Equation |
| --- | --- |
| `1` (default) | `A v = lambda B v` |
| `2` | `A B v = lambda v` |
| `3` | `B A v = lambda v` |

`A` must be symmetric/Hermitian and `B` symmetric/Hermitian **positive definite**.
Both matrices have shape `(n,n)` and the same real/complex double-precision
kind. `w` has length `m`; optional eigenvectors `z` have shape `(n,m)`, where
`m=n` or the number of requested eigenvalues. `subset_by_index` uses SciFort's
**one-based, inclusive** `[first,last]` convention (not Python's zero-based
indices). The same `lower` choice applies to both matrices; only that
triangle is read, and the other triangle may contain NaNs. The imaginary
Hermitian diagonal is ignored. The zero-order `(0,0)` problem is supported.

With optional `subset_by_value=[low,high]`, a generalized problem selects the
computed spectrum in the half-open `(low,high]` interval. Supply capacity
arrays `w` and `z`, usually of size `n` and `(n,n)`, and retrieve the count
through optional `num_selected`. Invalid intervals or conflicts with
`subset_by_index` return `linalg_status_invalid_shape`. See the
[separate selection guide](linalg_dense_eigen_value_selection.md) for output
sizes and failure semantics.

The returned columns are normalized by the mass matrix:

- Types 1 and 2: `transpose(conjg(z)) B z = I`.
- Type 3: `transpose(conjg(z)) inverse(B) z = I`.

When no `z` argument is provided, only eigenvalues are returned. The chosen
triangle, subset, and optional `max_iterations` have the same meanings as for
the previous standard routines.

## Algorithm

First compute the lower Cholesky factor `B = L L^H`. For types 1, 2, 3, form
an ordinary Hermitian eigenproblem for `C`:

- Type 1: `C = inverse(L) A inverse(L^H)`.
- Types 2 and 3: `C = L^H A L`.

SciFort's existing Householder/QL solver computes `C u = lambda u`. Map its
eigenvectors back to the original problem:

- Types 1 and 2: `v = inverse(L^H) u`.
- Type 3: `v = L u`.

The type-1 reduction uses triangular **solves** instead of explicitly
inverting `L`. The input triangles are copied into clean Hermitian matrices
before calling the existing Cholesky routine. This is the standard
Cholesky reduction for generalized Hermitian definite pencils (see
Golub and Van Loan, *Matrix Computations*, 4th edition, chapter 8).

No external BLAS/LAPACK routines or third-party implementation code are used.
Time complexity is `O(n^3)` and memory consumption `O(n^2)`. Selecting a
subset currently computes the full spectrum before filtering it.

## Status and numerical limitations

- `linalg_status_success`: converged solution (0).
- `linalg_status_invalid_shape`: malformed input/output dimensions,
  invalid eigenvalue interval, invalid type, or invalid iteration limit.
- `linalg_status_not_positive_definite`: Cholesky detected that `B` is not
  positive definite.
- `linalg_status_no_convergence`: nonfinite used data, overflow/underflow
  producing a nonfinite intermediate, or a numerical failure.

Like other Cholesky-based reductions, this method can lose accuracy if `B`
is extremely ill-conditioned. Both `A` and `B` should be scaled as necessary
to avoid overflow in the type-2/type-3 products, or extremely large
backtransformed eigenvectors. For large or difficult problems, a mature
LAPACK generalized eigenvalue driver may be more robust/faster. The
implementation does **not** yet support mixed real/complex `A`/`B`, batched
arrays, fast partial-spectrum algorithms, or LAPACK-specific drivers.

## Validation

```bash
python tools/linalg_reference/run_generalized_eigen_checks.py \
  --mode both --cases 500 --max-order 40 --seed 20261108
```

This compiles the native generalized-eigenproblem tests and compares both
builds against `scipy.linalg.eigvalsh`. For eigenvectors, the independent
checks measure the appropriate eigen-equation residual and `B` or
`B^-1` normalization. It exercises all three types, real and complex
inputs, selected eigenvalues, unused input triangles, repeated spectra,
non-positive-definite inputs, and scales ranging from `1e-100` to `1e100`.
See `docs/linalg_dense_eigen_generalized_validation_2026-11-08.md`.
