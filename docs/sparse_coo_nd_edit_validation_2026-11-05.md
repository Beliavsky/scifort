# N-dimensional COO boolean selection and scatter validation (2026-11-05)

GNU Fortran 14.2.0 was used with checked `-O0 -g -fcheck=all -fbacktrace`
and optimized `-O3` builds, with NumPy 2.3.5 and SciPy 1.17.0.

The randomized driver independently constructs an N-D SciPy COO source
with unsorted duplicates and stored zeros, and checks its dense numerical
values. Each test then applies either full-mask selection, axis boolean
selection, repeated-coordinate overwriting, additive scatter, or scalar
masked assignment. NumPy supplies the independent numerical reference;
a Python coordinate dictionary independently checks the **exact count of
stored entries** after each operation, including cancellations and retained
explicit zeros. Cases include ranks 1–5, zero dimensions and floating data.

| Seed and input | Checked | Optimized | SciPy source comparisons |
|---|---:|---:|---:|
| 20261105 / integer values, 650 trials | 650 PASS | 650 PASS | 650 |
| 274829 / floating-point values, 450 trials | 450 PASS | 450 PASS | 450 |
| 333471 / floating-point values, 180 additional trials | 180 PASS | 180 PASS | 180 |
| **Total** | **1,280 PASS** | **1,280 PASS** | **1,280 per mode** |

The focused Fortran test covers stored zero and duplicate behavior, two axes
of boolean selection, last-write-wins and additive updates, cancellations,
transactional invalid-input handling, nonfinite arithmetic, output caps,
empty dimensions, a 2-billion-by-2-billion sparse shape, and 12,000
scatter updates without dense allocation.

The native regression runner compiles the new module with the previous sparse
modules using the same flags and exercises the complete standalone sparse
suite. **All 38/38 standalone sparse programs passed in both checked and
optimized builds.** The bounds-checked run was completed in two segments
when the combined invocation reached the execution limit; all remaining tests
passed on rerun. The documentation example was also compiled and executed
with bounds checking.

FPM was not available in this environment, so the full project-wide
`fpm test` was not run. FPM's auto-test discovery will find the newly added
native test when dependencies are available.
