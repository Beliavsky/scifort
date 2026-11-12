# Sparse integer matrix powers and inverse — validation (2026-11-05)

## Implementation and reference

- Public generic `scifort_sparse_linalg.matrix_power` accepts canonical
  one-based CSR or CSC real matrices and nonnegative integer exponents.
  Binary exponentiation uses `csr_matmul`, avoiding dense intermediates.
- Public generic `inv` uses a single pivoted `splu` factorization, followed
  by `splu_solve` for blocks of identity columns; a diagonal fast path
  handles stored off-diagonal zeros. Numerical outputs are canonical CSR/CSC.
- SciPy 1.17 `scipy.sparse.linalg.matrix_power` and `inv` are numerical
  references, alongside **independent NumPy dense** matrix powers/inverses.
  No upstream source is copied or vendored.

## Tests and commands

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_matrix_functions_scipy.py --mode both --cases 120
```

`test/test_sparse_matrix_functions.f90` exercises valid and invalid inputs,
real rectangular rejection, a zero-by-zero matrix, exact stored zeros, a
large integer power of a nilpotent matrix, pivoting, singular detection,
dropping, numerical overflow, and a 12,000-row sparse diagonal matrix.

`compare_matrix_functions_scipy.py` generates controlled sparse random matrices
including pivot-requiring permutations, sparse triangular matrices, SPD
matrices, repeated coordinates and explicitly stored zeros. It compares:

1. CSR and CSC Fortran powers against `numpy.linalg.matrix_power` and the
   two SciPy sparse class variants.
2. CSR and CSC Fortran inverses against `numpy.linalg.inv` and SciPy `inv`.
3. Residuals `A @ inv(A) - I` for both storage formats.

The test process uses GNU Fortran 14.2.0 with `-std=f2018` and either
`-O0 -g -fcheck=all -fbacktrace` or `-O3`. The full FPM test suite requires
FPM, which is not installed in the evaluation environment.

## Limitations

The inverse can be dense even when `A` is sparse. Storage and memory for the
output may therefore grow quadratically. Large powers may densify or produce
nonfinite values and are reported as numerical breakdown. Optional inverse
dropping is absolute and does not preserve an exact inverse. A zero pivot
under the requested minimum pivot tolerance is reported as singular.

## Observed verification results

GNU Fortran 14.2.0 and SciPy 1.17.0 in the test environment:

| Check | Observed result |
|---|---|
| New randomized checked run, seed 20261105 | 320 matrices, 1280 CSR/CSC numerical outputs passed |
| New randomized checked run, seed 222 | 240 matrices, 960 outputs passed |
| New randomized optimized run, seed 20261106 | 200 matrices, 800 outputs passed |
| Native sparse regression runner, checked | All 21 programs passed |
| Native sparse regression runner, optimized | All 21 programs passed |
| Large diagonal power and inverse | 12,000 rows passed |

The existing LSMR callback executable-stack linker warning persists.
`fpm test` was not run because FPM was not installed in this environment.
