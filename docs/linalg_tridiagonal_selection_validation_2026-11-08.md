# Partial tridiagonal eigenvalue selection: verification (2026-11-08)

## Scope

This update adds index- and value-selected eigenvalues to the existing
`scifort_linalg_banded_eigen` public API. Selected eigenvalues without
vectors are computed by a scaled Sturm recurrence and bisection, with
Gershgorin initial bounds. An exactly diagonal input is handled by an
exact-value sorting path, preserving separated magnitudes. Selected
vectors still use the existing full-spectrum QL algorithm.

The new tests check:

- Inclusive, one-based index selections, and half-open `(low,high]`
  value selections, including zero selected eigenvalues.
- Repeated eigenvalues, diagonal/split matrices, all-zero spectra, and
  intervals with infinite endpoints.
- Eigenvalue ordering, selected-vector orthogonality and residuals,
  comparisons with SciPy's **full-spectrum** and **selected-spectrum**
  tridiagonal eigensolvers, and status reporting for invalid calls.
- Input scales from about `1e-150` to `1e150`, large/small bounds, and
  exactly diagonal inputs with coefficients ranging from `1e-200` to `1e200`.
- An order-10,000 Toeplitz tridiagonal matrix, selecting indices 5000–5004,
  checked against its known analytic eigenvalues rather than dense storage.

## Completed checked and optimized tests

Run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/linalg_reference/run_tridiagonal_selection_checks.py \
    --mode both --cases 720 --max-order 128 --seed 20261108

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/linalg_reference/run_tridiagonal_selection_checks.py \
    --mode both --cases 560 --max-order 88 --seed 31877
```

The two runs each passed in both checked GNU Fortran (`-O0 -g
-fcheck=all -ffpe-trap=invalid,zero,overflow`) and optimized (`-O3`)
modes: **2,560 completed SciPy-reference comparisons** in total.

| Run | Max normalized eigenvalue error | Max vector residual | Max orthogonality error |
|---|---:|---:|---:|
| 720 cases per mode | `5.711e-15` | `4.174e-15` | `1.998e-15` |
| 560 cases per mode | `6.059e-15` | `2.932e-15` | `2.255e-15` |

The normalized eigenvalue error uses the full-spectrum SciPy
`eigvalsh_tridiagonal` output, restricting it to the requested indices or
value interval. SciPy's independent selected-spectrum drivers were also
called to check the returned count and to track any reference discrepancy.
For one tiny-scale, closely clustered case, SciPy's selected-spectrum
`stebz` result differed from its full-spectrum result by `4.164e-10`
relative to the matrix scale, whereas SciFort agreed with the full-spectrum
reference to several parts in `1e-15`. We report that driver discrepancy
rather than increasing the strict numerical tolerance across all cases.

## Neighboring regressions

These completed tests passed after the update, in **both** compiler modes:

- `run_banded_eigen_checks.py`, 130 comparisons per mode.
- `run_dense_eigen_checks.py`, 120 comparisons per mode.
- `run_generalized_eigen_checks.py`, 70 comparisons per mode.
- `run_eigen_by_value_checks.py`, 36 comparisons per mode.

The original native tridiagonal/banded tests and the new tridiagonal
selection native tests also passed, including invalid arguments and
nonconvergence status. The full FPM suite could not be run because FPM is
not installed in this environment.

## Algorithmic limits

For k requested eigenvalues of a real tridiagonal matrix of order n, the
bisection path requires O(n) temporary storage and O(n k log(width/tol))
operations; it is most useful when k is small. It scales by the maximum
coefficient magnitude and has an automatic convergence width of roughly
machine epsilon times that scale. Clusters whose absolute separation is
far below that width cannot be individually resolved reliably. Exact
block-diagonal matrices with disconnected nontrivial blocks still use
global scaling, unlike a fully block-deflated LAPACK implementation.

Selected eigenvectors currently run full-spectrum QL before selecting
columns and therefore do not save memory or computation. Full selected
banded and dense eigenproblem solvers remain separate future work.
