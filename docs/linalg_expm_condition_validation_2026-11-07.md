# SciFort exponential condition-number validation — 2026-11-07

## Environment and implementation

- GNU Fortran 14.2.0, `-std=f2018`
- Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0
- `OPENBLAS_NUM_THREADS=1` for deterministic SciPy comparison timing
- `src/linalg_expm_condition.f90` reuses the original SciFort
  `src/linalg_expm_frechet.f90` without additional dependencies

## Exact commands and results

From the SciFort project root:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_expm_condition_checks.py \
    --mode both --cases 600 --max-order 8
```

Checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`) builds:

- Native test: `PASS test_linalg_expm_condition` in both modes.
- Each mode: **600/600** independent real/complex comparisons with
  `scipy.linalg.expm_cond`, no failures.
- Worst relative discrepancy per mode: **1.892e-9**.
- Mean/max Lanczos iteration counts: **7.0 / 32** per mode.
- Additional **80/80** analytical stress cases per mode: pass,
  worst relative discrepancy **1.626e-14**.
- These cases cover diagonal matrices up to **order 48**, negative and positive
  offsets up to magnitude 1000, close and repeated eigenvalues, complex
  diagonal entries, and shifted **nondiagonal** real/complex matrices.

Larger random inputs and an independent seed:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_expm_condition_checks.py \
    --mode both --cases 180 --max-order 12 --seed 20261108 --skip-stress
```

- Each mode: **180/180** SciPy comparisons pass, no failures.
- Worst relative discrepancy: **5.881e-10** per mode.
- Mean/max iterations: checked **9.6 / 63**; optimized **9.7 / 63**.

**Total final expm_cond comparisons:** 1,560 randomized SciPy cases and
160 analytical stress cases across checked/optimized builds, without a
failure. These are tests, not an exhaustive error guarantee.

## Regression of existing dense matrix functions

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_matrix_functions_checks.py \
    --mode both --cases 80
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_logarithm_checks.py \
    --mode both --cases 80
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/compare_polar_frechet_scipy.py \
    --mode both --cases 80
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_trigonometric_checks.py \
    --mode both --cases 80
```

All native tests passed. All 640 SciPy comparison cases passed across the
four harnesses and both build modes; no unexpected nonconvergences or
out-of-tolerance comparisons. These cover `sqrtm`/`signm`, `logm`/
`fractional_matrix_power`, `polar`/`expm_frechet`, and
`sinm`/`cosm`/`tanm`/`sinhm`/`coshm`/`tanhm`.

`example/linalg_expm_condition_example.f90` also compiled and ran; its
reference 3-by-3 result was 1.77878059, compared with the SciPy manual
value 1.7787805864469866.

## Limitations

The full `fpm test` suite was **not run**, because FPM is not installed.
Intel/NVIDIA Fortran and other compilers were not tested. The Lanczos
spectral estimate is not a rigorous upper bound and may fail to converge
within 80 steps for difficult clustered matrices. Neither a full SVD nor
an independent proof of accuracy for arbitrary inputs was performed.
