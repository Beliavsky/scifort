# Hankel matrix construction and products — validation (2026-11-09)

New files:

- `src/linalg_hankel.f90`
- `src/linalg_matmul_hankel.f90`
- `test/test_linalg_hankel.f90`
- `example/linalg_hankel_example.f90`
- `tools/linalg_reference/driver_hankel.f90`
- `tools/linalg_reference/run_hankel_checks.py`

Test environment: GNU Fortran 14.2.0, SciPy 1.17.0, NumPy 2.3.5.

## Native Fortran tests

The following source modules were compiled in order for the focused tests:

```
src/kinds.f90
src/linalg.f90
src/linalg_dense.f90
src/linalg_toeplitz.f90
src/signal_windows_extra.f90
src/signal_dpss.f90
src/signal_spectral.f90
src/linalg_matmul_toeplitz.f90
src/linalg_circulant.f90
src/linalg_hankel.f90
src/linalg_matmul_hankel.f90
```

With `-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`,
the new native `test_linalg_hankel` and existing
`test_linalg_matmul_toeplitz`, `test_linalg_toeplitz`, and
`test_linalg_circulant` all passed. All four also passed with
`-std=f2018 -O3 -march=native`.

New native tests cover the SciPy documented Hankel example, real and
complex constructors, mixed real/complex row constructors, default
zero-padded rows, real/complex/mixed right-hand sides of rank 1 and 2,
rectangular dimensions, ignored first entry of the supplied row, status
on shape errors and nonfinite coefficients, zero rows and columns,
lengths 1–130, scales from `1e-80` to `1e80`, and an order-4096 sparse
anti-diagonal case (FFT length 8191).

## Independent SciPy comparison

From the repository root:

```bash
python tools/linalg_reference/run_hankel_checks.py \
    --mode checked --cases 1800 --seed 20261109
python tools/linalg_reference/run_hankel_checks.py \
    --mode release --cases 1800 --seed 20261110
```

The script compiles the actual Fortran modules and a Fortran driver.
For each case, it compares the Fortran product with both
`scipy.linalg.hankel(c,r) @ x` and an independently constructed
anti-diagonal array multiplied by NumPy. Selected cases also compare
the Fortran `hankel` constructor elementwise against SciPy.

| Build | Products passed | Constructor comparisons | Maximum relative product difference | Constructor difference |
| --- | ---: | ---: | ---: | ---: |
| Checked `-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow` | 1,800/1,800 | 226 | 3.246e-14 | 0 |
| Optimized `-O3 -march=native` | 1,800/1,800 | 226 | 1.876e-14 | 0 |
| **Total** | **3,600/3,600** | **452** | **3.246e-14** | **0** |

The randomized references cover real and complex coefficient arrays,
real and complex multipliers, varying rectangular shapes including empty
dimensions, an optional last row, deliberately mismatched `r(1)`,
multiple right-hand sides, small and prime FFT lengths, and scales
between `1e-80` and `1e80`. Relative numerical tolerance: `3e-11`.

Full `fpm test` was not run because FPM is not installed in this environment.
Focused tests do not substitute for full cross-platform CI.
