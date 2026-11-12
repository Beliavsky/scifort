# Matrix-free randomized ID validation — 2026-11-10

## Environment

- GNU Fortran 14.2.0
- SciPy 1.17.0 / NumPy 2.3.5
- Checked (`-O0 -g -fcheck=all`) and optimized (`-O3 -fcheck=all`) builds
- New dependencies: none

## Native Fortran tests

| Configuration | New native checks | Result |
|---|---:|---|
| Checked | 42,911 | Passed |
| Optimized | 42,911 | Passed |
| Total | **85,822** | Passed |

The native test exercises real and complex callbacks, positive and zero
numerical ranks, fixed/adaptive modes, square and rectangular shapes,
zero arrays, tiny and large scales (`1e-80`, `1e80`), singleton/zero
sizes, exact selected skeleton columns, one-based permutations,
reproducible seeded pivots, and invalid/error/NaN callbacks. It confirms
that a fixed rank-three solve with `m=20,n=50` invokes exactly three
forward products and ten adjoint products (`k+oversampling+k`).

## Independent SciPy/NumPy comparison

`test/compare_linalg_interpolative_operator_scipy.py` generates real and
complex full- and low-rank matrices in 12 shapes and 4 data variants,
then compares the Fortran output against both:

- independent truncated NumPy SVD approximations; and
- SciPy 1.17.0 `interpolative.interp_decomp(aslinearoperator(A), ...)`,
  using both rank and tolerance signatures.

| Metric | Checked | Optimized |
|---|---:|---:|
| Fortran reference cases | 576 passed | 576 passed |
| SciPy operator reference calls succeeding | 576 | 576 |
| Maximum fixed-rank reconstruction error / optimal SVD tail | 1.451 | 1.451 |
| Maximum adaptive relative Frobenius error | 0.2134 | 0.2134 |
| Maximum error ratio to SciPy reconstruction | 1.043 | 1.043 |

The numerical tolerance for adaptive calls includes `eps=0.25`,
`0.001`, and `1e-7`. All Fortran adaptive cases meet their input tolerance
within the test's specified allowance. SciPy and SciFort use different
sampling methods and accuracy targets, so rank and pivot equality is
neither required nor expected. Not all SciPy 1.17.0 operator estimates meet
this implementation's Frobenius criterion, which is another reason not to
use exact rank agreement as a validation condition.

## Regression tests

In **both** checked and optimized builds:

- Existing deterministic interpolative decomposition: 14,944 passed
- Existing randomized dense interpolative decomposition: 37,449 passed
- Existing deterministic interpolative SVD: 32,770 passed
- Existing randomized interpolative SVD: 50,408 passed
- Existing QR and SVD native regression suites: passed
- True matrix-free example (rank-three real and complex operators): passed

The stand-alone full `fpm test` suite was not run because FPM is not
installed in the execution environment. The named module-level native
regressions above were compiled and run directly with GNU Fortran.

## Reproducing the new tests

From the project root, with a temporary `build` directory and GNU Fortran:

```bash
mkdir -p build
cd build
gfortran -std=f2008 -O0 -g -fcheck=all -c \
    ../src/kinds.f90 ../src/linalg.f90 ../src/linalg_svd.f90 \
    ../src/linalg_qr.f90 ../src/linalg_interpolative.f90 \
    ../src/linalg_interpolative_operator.f90

gfortran -std=f2008 -O0 -g -fcheck=all -I. \
    ../test/test_linalg_interpolative_operator.f90 *.o -o test_id_operator
./test_id_operator

gfortran -std=f2008 -O0 -g -fcheck=all -I. \
    ../reference/linalg_interpolative_operator_driver.f90 *.o -o id_operator_driver
python ../test/compare_linalg_interpolative_operator_scipy.py ./id_operator_driver
```

Note that stochastic residual probes are **estimates**, not formal error
certificates. Exact reconstruction error was computed independently using
the known dense reference matrices in the validation scripts.
