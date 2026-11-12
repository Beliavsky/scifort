# Sparse COO N-D predicate and coordinate-mask validation (2026-11-05)

## Environment

- GNU Fortran 14.2.0, `-std=f2018`.
- Checked numerical builds: `-O0 -g -fcheck=all -fbacktrace`.
- Optimized numerical builds: `-O3`.
- NumPy 2.3.5 and SciPy 1.17.0, used only as independent test references.

## Randomized comparisons

`tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py` builds a fresh
Fortran reference driver and generates N-D SciPy `coo_array` sources, including
unsorted entries, duplicate coordinates, explicitly stored zeros, rank 1–6,
zero-size dimensions, integer and floating values, and mask cancellation.
The script independently checks dense SciPy values against coordinate
aggregation, then compares Fortran results with NumPy/dictionary references,
including **stored-entry coordinates and values in their original order**.
Cases test sparse comparison masks, sparse-valued mask filtering, repeated
coordinate-set membership, complemented masks, value-predicate selection,
nonzero counting and correctly rejected predicates true at implicit zero.

| Seed and data | Checked `-O0 -fcheck=all` | Optimized `-O3` |
|---|---:|---:|
| 20261105, integer-valued, 300 cases | 300 / 300 | 300 / 300 |
| 274829, floating, 450 cases | 450 / 450 | 450 / 450 |
| 333471, floating, 450 cases | 450 / 450 | 450 / 450 |
| **Total** | **1,200 / 1,200** | **1,200 / 1,200** |

The reference comparison programs compile and run successfully with the
same build modes, without requiring the project-wide FPM dependencies.

## Native Fortran tests

`test/test_sparse_coo_nd_filter.f90` includes checks of every supported
comparison, duplicate cancellation, numerical zero versus stored zero,
complimentary masks, invalid operators and unsupported implicit-zero
predicates, output nnz caps, nonfinite duplicate sums, invalid indices,
zero-sized dimensions, and large sparse shapes. A separate 12,000-entry
source with 7,000 point-selector events validates sparse selection over
an input with extents `(2,000,000,000, 2,000,000,000)` without dense data.
The focused native program passes in both bounds-checked and `-O3` builds.
The Fortran usage example in `docs/sparse_coo_nd_filter.md` also compiles
and passes with bounds checking.

`python tools/sparse_reference/run_native_sparse.py --mode checked` passed
**all 39 of 39 native sparse regression programs**, including the prior
solvers, sparse formats, graph algorithms and N-D COO operations.

A complete `-O3` regression pass was attempted but reached the execution
limit before reporting results. A secondary `-O1` full-suite attempt passed
its first five native sparse programs before hitting the execution limit.
Thus **39/39 full-suite optimized verification is not claimed**; the new
module and all 1,200 per-mode randomized cases were independently run at
`-O3` successfully. The full project-wide `fpm test` was not run because
FPM is unavailable in this environment.

## Commands

```console
python tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py --mode both --trials 300 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py --mode both --trials 450 --seed 274829 --floating
python tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py --mode both --trials 450 --seed 333471 --floating
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
```
