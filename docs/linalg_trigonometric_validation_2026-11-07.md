# Dense matrix trigonometric validation — 2026-11-07

## Environment

- GNU Fortran 14.2.0 (gfortran, standard Fortran 2018)
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0
- No FPM or external BLAS/LAPACK required for these standalone checks

## Native tests and SciPy reference comparisons

The native suite checks all six real/complex function families, matrix
identities, a nonnormal upper triangular example, a complex analytic
continuation identity, tiny arguments (1e-120), exact tangent poles, an
invalid output shape, non-square input, a matrix above order 40, NaN,
infinity and excessive scaling.

Executed from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_trigonometric_checks.py \
  --mode both --cases 1200
```

Results in both builds (`-O0 -g -fcheck=all -fbacktrace` and `-O3`):

- All native regression tests: PASS.
- 1,200/1,200 SciPy comparisons per build: PASS (200 per function).
- Worst normalized 1-norm difference: 2.277e-12 checked; 2.829e-12 optimized.
- No rejected cases, failures, or out-of-tolerance comparisons.

Larger matrices were also checked:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_trigonometric_checks.py \
  --mode optimized --cases 240 --max-order 16
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_trigonometric_checks.py \
  --mode both --cases 240 --max-order 40
```

- Up to order 16, optimized: 240/240 pass, maximum normalized difference
  3.338e-9.
- Up to order 40, both checked and optimized: 238/240 evaluated successfully.
  The two remaining cases were tangent/hyperbolic tangent matrices whose
  cosine/cosh denominator had condition numbers 1.186e12 and 7.850e10;
  SciFort correctly returned `linalg_status_no_convergence`, not an
  unreliable answer. The comparison harness checks the denominator condition
  number before counting a rejection as expected. Zero unexpected failures
  or mismatches. Maximum normalized difference among successful cases:
  3.083e-9 checked, 3.828e-9 optimized.

## Existing functionality regression checks

Ran the independent GNU Fortran/SciPy harnesses in both checked and `-O3`
builds after adding the module:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_matrix_functions_checks.py \
  --mode both --cases 120
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_logarithm_checks.py \
  --mode both --cases 120
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/compare_polar_frechet_scipy.py \
  --mode both --cases 120
```

All existing native and reference cases passed in both builds:
`sqrtm`/`signm` (120 per mode), `logm`/fractional powers (120 per mode), and
`polar`/`expm_frechet` (60+60 per mode). The largest differences were within
the established tolerances. SciPy emitted its own near-singular warning for
some logarithm reference cases; this did not affect the passing results.

## Scope and untested aspects

These are deterministic validation runs, not exhaustive proof of correctness.
FPM, Intel Fortran, nvfortran, larger orders, and optional LAPACK backends
were **not** tested here. Real/complex type dispatch and numerical behavior
are independently checked against SciPy, not based on translations of SciPy's
implementation.
