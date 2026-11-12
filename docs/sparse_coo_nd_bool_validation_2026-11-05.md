# Complement-aware N-D COO Boolean mask validation (2026-11-05)

## Environment

- GNU Fortran 14.2.0 (`-std=f2018`).
- Checked: `-O0 -g -fcheck=all -fbacktrace`.
- Optimized: `-O3`.
- NumPy 2.3.5 and SciPy 1.17.0 for independent numerical reference arrays.

## Randomized independent comparisons

`tools/sparse_reference/compare_sparse_coo_nd_bool_scipy.py` compiles
`tools/sparse_reference/sparse_coo_nd_bool_reference_driver.f90` once per
mode and streams randomized cases through the Fortran implementation.
The Python side builds independent NumPy Boolean reference tensors and
numerical SciPy `coo_array` inputs. It checks canonical exception coordinates,
implicit Boolean defaults, duplicate aggregation, explicit zero preservation,
logical AND/OR/XOR/difference/NOT, scalar comparisons at implicit zeros,
sparse filtering, materialization restrictions, and exact true counts.

The cases cover ranks 1–5, empty axes, unsorted repeated points, all four
combinations of true/false defaults, duplicate cancellation, and both
implicit-false and implicit-true predicate results.

| Randomized run | Checked | Optimized |
|---|---:|---:|
| Seed 20261105, 900 trials | 900/900 | 900/900 |
| Seed 274829, 900 trials | 900/900 | Not completed (combined run reached execution limit) |

The second run is included as extra checked coverage only; the optimized
second-seed run was not completed and is not claimed.

## Native sparse regression

`python tools/sparse_reference/run_native_sparse.py --mode checked`
passed **40/40 sparse regression programs**, including the new
`test/test_sparse_coo_nd_bool.f90`.

`python tools/sparse_reference/run_native_sparse.py --mode optimized`
passed **40/40 sparse regression programs**, including the same new test
compiled with `-O3`.

The focused test verifies four Boolean operations for all four implicit
background combinations, De Morgan's law, double complement, scalar
predicates true at zero, duplicate cancellation, stored zeros, failure
statuses and output caps, and 64-bit true counting up to
3,999,999,999,999,999,999. It detects signed 64-bit overflow when a
full implicit-true background has too many elements, while still allowing
implicit-false masks on even larger shapes. A 12,000-entry sparse set
is tested on an input shape of `(2000000000,2000000000)` without dense
allocation. The documentation example compiles and runs with checking.

The full project-wide `fpm test` was not run because FPM is not installed
in this environment; these results cover the standalone sparse regression
runner and randomized independent comparisons, not every SciFort module.

## Reproduction

```console
python tools/sparse_reference/compare_sparse_coo_nd_bool_scipy.py --mode both --trials 900 --seed 20261105
python tools/sparse_reference/run_native_sparse.py --mode both
```

The original MIT-licensed implementation introduces no additional runtime
numerical dependency. NumPy and SciPy are used only in comparison tests.
