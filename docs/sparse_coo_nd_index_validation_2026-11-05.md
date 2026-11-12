# N-dimensional COO indexing validation (2026-11-05)

## Environment

- GNU Fortran 14.2.0, Fortran 2018, checked `-O0 -g -fcheck=all -fbacktrace`
  and optimized `-O3` builds.
- NumPy and SciPy 1.17.0.

## Independent randomized comparisons

`compare_sparse_coo_nd_index_scipy.py` generates unsorted N-D COO inputs,
including duplicate positions, stored zeros, scalar indices, positive and
negative strided slices, empty selections, and repeated coordinate/index
arrays. It compares SciFort results against NumPy operations on SciPy COO
inputs. For basic-indexing cases, SciPy 1.17's direct COO indexing is
separately cross-checked against NumPy. For multiple indexed axes, the
explicit Fortran Cartesian-selection semantics are the reference.

| Seed / data | Checked | Optimized | SciPy basic indexing cross-checks |
|---|---:|---:|---:|
| 20261105 / integer values, 650 trials | 650 PASS | 650 PASS | 144 |
| 274829 / real Gaussian values, 450 trials | 450 PASS | 450 PASS | 110 |
| **Total** | **1,100 PASS** | **1,100 PASS** | **254 cases per mode** |

Run:

```console
python tools/sparse_reference/compare_sparse_coo_nd_index_scipy.py --mode both --trials 650 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_index_scipy.py --mode both --trials 450 --seed 274829 --floating
```

## Native validation

The focused `test/test_sparse_coo_nd_index.f90` covers negative strides,
Cartesian and repeated index selections, pointwise gathering, explicitly
stored zero retention, stored duplicate preservation, optional singleton
axes, all-scalar rank-one results, zero-sized dimensions, invalid inputs,
configurable entry caps, and 2-billion-by-2-billion sparse shapes without
dense allocation.

The integrated native regression runner compiled the new indexing module and
all existing sparse modules with GNU Fortran 14.2.0 and completed:

| Mode | Standalone native sparse tests |
|---|---:|
| Checked | 37/37 PASS |
| Optimized | 37/37 PASS |

The documentation's Fortran example was also compiled and executed with
bounds checking. The project-wide `fpm test` could not be run because FPM is
unavailable in this environment; native standalone tests compile without it.
The existing `sparse_linalg` module produces a nonfatal linker warning about
an executable stack; this predates the new indexing module.
