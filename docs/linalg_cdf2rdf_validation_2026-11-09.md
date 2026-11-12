# Validation of SciFort `linalg_cdf2rdf` (2026-11-09)

## Implementation tested

The new `src/linalg_cdf2rdf.f90` accepts real or complex eigenvalues,
real or complex eigenvectors, and one leading batch dimension. It
constructs real eigenvalue blocks and a real invariant eigenvector basis
for eigenstructures of a real matrix. No code in the previously existing
Fortran source modules was modified.

## Native GNU Fortran results

Test environment: GNU Fortran 14.2.0, SciPy 1.17.0, NumPy 2.3.5.

| Test | `-g -O0 -fcheck=all` | `-O3` |
|---|---:|---:|
| `test_linalg_cdf2rdf` | 684 cases passed | 684 cases passed |
| `test_linalg_cdf2rdf_eig` | 2 cases passed | 2 cases passed |
| `test_linalg_general_eigen` | passed | passed |
| `test_linalg_lu_det_inv` | 2,416 cases passed | 2,416 cases passed |
| `test_linalg_dense` | passed | passed |
| `linalg_cdf2rdf_example` | passed | passed |

The new native tests cover zero-order matrices, real and complex values,
negative-imaginary-first ordering, eigenvector conjugate-phase differences,
nonadjacent pair members, real+complex eigenvector combinations, batched
matrices, malformed shapes, unmatched and interleaved complex pairs,
nonfinite values, and cleanup after failures. A separate test obtains
eigenpairs from SciFort's existing nonsymmetric `eig` routine and checks
`A * vr = vr * wr`.

A direct integration test exposed a real eigenvalue returned by SciFort
`eig` with a spurious imaginary component of `3.33e-15`. A selective
roundoff fallback restores successful conversion; the 3-by-3 example
reconstruction has max absolute residual `3.11e-15`. Tiny **valid** complex
pairs remain complex even with imaginary parts as small as `1e-20`.

## Independent comparison with SciPy

`reference/check_linalg_cdf2rdf.py` sends deterministic independent cases
to `reference/linalg_cdf2rdf_driver.f90`. Reference output is calculated
by *calling SciPy* `linalg.cdf2rdf`, not by reproducing the Fortran routine.
The generator uses SciPy eigendecompositions of random real matrices and
independent exact conjugate eigenpair constructions, including real
entries separating pair members and scales `1e-70`, `1`, `1e70`.

| Measure | Checked | Optimized |
|---|---:|---:|
| All comparison cases | 1,457 passed | 1,457 passed |
| Valid matrix comparisons | 1,291 passed | 1,291 passed |
| Expected input rejections | 166 passed | 166 passed |
| Largest scaled `wr` difference | 0 | 0 |
| Largest scaled `vr` difference | 0 | 0 |
| Largest independently checked eigen-identity residual | `8.443e-14` | `8.443e-14` |

The reference comparisons are especially sensitive to the sign of the
2-by-2 off-diagonal eigenvalue blocks and the orientation of the real
invariant eigenvector basis. Differences in eigenvector phase are not
assumed to be errors: SciFort derives the real basis from one eigenvector
of each pair, so it can also handle eigenvector phases that SciPy's
`cdf2rdf` does not align internally. Such phase-mismatch cases are
checked using the invariant-identity native tests, not direct entrywise
SciPy comparisons.

## Reproduction without FPM

From the repository root, with a clean writable build directory:

```sh
mkdir -p build-cdf2rdf
cd build-cdf2rdf
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -c \
    ../src/kinds.f90 ../src/linalg.f90 \
    ../src/linalg_cdf2rdf.f90
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace \
    ../test/test_linalg_cdf2rdf.f90 \
    kinds.o linalg.o linalg_cdf2rdf.o -o test_cdf2rdf
./test_cdf2rdf
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace \
    ../reference/linalg_cdf2rdf_driver.f90 \
    kinds.o linalg.o linalg_cdf2rdf.o -o cdf2rdf_driver
python ../reference/check_linalg_cdf2rdf.py --driver ./cdf2rdf_driver
```

For `-O3`, recompile the source and tests in a *fresh directory* with
`-O3` instead of the checked-build flags.

The project-wide `fpm test` suite could not be run in this environment
because the FPM executable is unavailable. These results are targeted
native builds and regressions, not a claim that every one of the project's
approximately 3,000 files was compiled or tested.

Public mathematical/API reference:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cdf2rdf.html
