# N-dimensional sparse COO operations validation (2026-11-05)

## Environment

- GNU Fortran 14.2.0 (Fortran 2018, binary64)
- NumPy and SciPy 1.17.0
- Compiler modes: checked `-O0 -g -fcheck=all -fbacktrace`, optimized `-O3`

## Independent SciPy / NumPy comparisons

The reference harness creates SciPy `coo_array` objects from independently
randomized zero-based COO coordinates, including unsorted duplicate positions
and explicitly stored zeros. It computes the reference using SciPy's COO
dense conversion followed by NumPy reductions, broadcasting, arithmetic and
`tensordot`. SciFort performs all operations with sparse COO structures and
only converts the small *result* to dense for numerical comparison.

| Reference command (per compiler mode) | Checked | Optimized |
|---|---:|---:|
| `--trials 900 --seed 20261105` (integer values) | 900/900 PASS | 900/900 PASS |
| `--trials 600 --seed 274829` (integer values) | 600/600 PASS | 600/600 PASS |
| `--trials 750 --seed 11111 --floating` (Gaussian values) | 750/750 PASS | 750/750 PASS |
| **Total** | **2,250/2,250 PASS** | **2,250/2,250 PASS** |

**4,500 total comparisons passed.** Each group cycles equally among six
operation classes (multi-axis sum, broadcasting, addition, multiplication,
contraction, and scalar sum), covering zero and singleton dimensions, ranks
1-5, duplicate entries, stored zeros, rank expansion, broadcasted operands,
empty contractions, contracted-axis permutations, outer products and full
contractions. Exact equality was required for integer-valued cases up to
normal floating arithmetic; floating-valued results use `rtol=atol=2e-12`.

Run independently with:

```console
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 900 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 600 --seed 274829
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 750 --seed 11111 --floating
```

## Native Fortran tests

```console
python tools/sparse_reference/run_native_sparse.py --mode both
```

| Configuration | Native sparse regression executables |
|---|---:|
| Checked build | 36/36 PASS |
| Optimized build | 36/36 PASS |

The new `test_sparse_coo_nd_ops.f90` exercises duplicate and zero behavior,
noncanonical input, broadcasting rules, error paths and `max_nnz`, 2D/3D
contractions, an outer product, complete contractions, a tensor with a zero
extent, overflow handling, and 2-billion-by-2-billion sparse shapes without
allocating dense arrays. The focused test was additionally recompiled and
run after the large-shape case was added, in **both** checked and optimized
configurations.

The documentation's standalone Fortran example was also compiled and run
under the checked configuration.

The project-wide `fpm test` was not performed because `fpm` is unavailable
in this environment. The repository's `tools/check_ascii.py` currently
reports existing non-ASCII passages in many prior files; the new Fortran
implementation and new validation harness use ASCII-only source.
