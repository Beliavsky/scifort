# Dense logarithm and fractional matrix-power validation — 2026-11-07

Environment: GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5, SciPy 1.17.0.
Two build modes, both `-std=f2018`: checked `-O0 -g -fcheck=all -fbacktrace`;
optimized `-O3`.

Reproduce from repository root:

```bash
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_logarithm_checks.py \
    --mode both --cases 800 --seed 20261107
```

| Build | Focused native test | Random SciPy comparisons | Maximum normalized 1-norm error | Median error |
|---|---|---:|---:|---:|
| Checked | PASS | 800/800 | 2.22e-13 | 2.27e-16 |
| Optimized | PASS | 800/800 | 1.25e-13 | 2.39e-16 |

**Cases**: real SPD, real nonnormal matrix exponentials, real triangular
Jordan-like matrices, complex nonnormal exponentials, real negative diagonal
spectra, normal complex matrices with phases near the negative axis,
very small/large diagonal scales, and complex matrices near identity.
Matrix orders 1 through 5 are used; cases alternate between `logm` and
`fractional_matrix_power`, including positive/negative integers and
positive/negative nonintegral exponents. All 1,600 randomized comparisons
succeeded with relative error below `2e-7`; the reported errors above are
much smaller. The focused test separately exercises reconstruction, principal
branch, singular cases, control-budget failures, invalid shapes/nonfinite
inputs, real-output convenience APIs, and tiny/large scalings.

SciPy may emit `LogmNearlySingularWarning` for numerically tiny diagonal
matrices, even though their logarithms are mathematically unambiguous.
Test generation limits extreme-scale **power** cases to those whose true
outputs fit within binary64 (e.g. `1e150 ** 3` is not representable); log
cases still include scales of `1e±150`. Fortran returns nonconvergence
instead of claiming a finite result when floating-point overflow occurs.

Additional independent checks succeeded for 160 matrices of orders 1-12 in
checked mode and 80 matrices of orders 1-20 in optimized mode. Existing
`sqrtm`/`signm` and `polar`/`expm_frechet` native regressions were also
compiled and passed in both compiler modes.

This tests ordinary small-matrix regimes, not a claim of uniform accuracy
for all matrices of order 40. This inverse-scaling/series method is **not**
SciPy's Schur-Pade method and can reject nearly singular, defective, highly
nonnormal, or branch-cut matrices. The project-wide `fpm test` could **not**
be run in this environment because `fpm` was not installed; the focused
Fortran tests and comparisons were built directly with `gfortran`.
