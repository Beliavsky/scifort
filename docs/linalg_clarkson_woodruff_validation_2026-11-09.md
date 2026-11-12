# CountSketch (`scipy.linalg.clarkson_woodruff_transform`) validation — 2026-11-09

## Environment

- GNU Fortran (gfortran) 14.2.0, Fortran 2018 mode
- NumPy 2.3.5, SciPy 1.17.0
- Checked build: `-O0 -g -fcheck=all -fbacktrace`
- Optimized build: `-O3 -march=native -fcheck=all -fbacktrace`
- Standalone build depends on SciFort's `kinds.f90` and `linalg.f90` only.

## Reproduce

From the project root, in a fresh temporary directory (shown as `build`):

```sh
mkdir -p build
cd build
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -c \
    ../src/kinds.f90 ../src/linalg.f90 ../src/linalg_clarkson_woodruff.f90
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace \
    -o test_cwt kinds.o linalg.o linalg_clarkson_woodruff.o \
    ../test/test_linalg_clarkson_woodruff.f90
./test_cwt
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace \
    -o driver_cwt kinds.o linalg.o linalg_clarkson_woodruff.o \
    ../validation/reference_drivers/driver_linalg_clarkson_woodruff.f90
python ../test/compare_linalg_clarkson_woodruff_scipy.py ./driver_cwt
```

Repeat with `-O3 -march=native` for the optimized configuration.

## Results

| Check | Checked (`-O0`) | Optimized (`-O3`) |
| --- | ---: | ---: |
| Independent comparisons with SciPy | 990 / 990 | 990 / 990 |
| Real dense, complex dense | 165 each | 165 each |
| Real CSR, complex CSR | 165 each | 165 each |
| Real batched, complex batched | 165 each | 165 each |
| Sketch elements compared | 25,728 | 25,728 |
| Maximum normalized sketch difference | 0.0 | 0.0 |
| Native CountSketch test | Pass | Pass |
| Mean squared Frobenius norm ratio (20 seeds) | 0.999119 | 0.999119 |
| Existing interpolative native suite | 14,944 passed | 14,944 passed |
| Existing subspace native suite | Pass | Pass |
| Executable CountSketch example | Pass | Pass |

**Total: 1,980 independent SciPy comparisons** across both
configurations, plus native tests and previous regression suites.

All reference cases use a freshly seeded NumPy generator for a SciPy
`cwt_matrix`, and the same seed is passed independently to the **public**
`scipy.linalg.clarkson_woodruff_transform`. The matrix's bucket/sign arrays
are transmitted to SciFort using one-based indices. The resulting values
match exactly, including SciPy's CSR sparse-matrix sketching and shared
batch mapping. The comparison includes empty dimensions and small and
rectangular matrices; it does not make a claim about a single-seed statistical
norm guarantee.

Native tests also exercise explicit mappings, deterministic Park–Miller
seeds, complex signs, duplicate and unsorted CSR entries, empty rows,
nonfinite inputs, malformed CSR offsets/columns, invalid plan components,
exclusive plan/seed arguments, and 20 large seeded norm-preservation checks.

The complete `fpm test` suite was not executed because FPM is not installed in
this environment. Only the targeted modules and regressions named above were
compiled and run.
