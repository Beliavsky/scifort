# Randomized interpolative SVD — validation (2026-11-10)

Baseline: `scifort-linalg-interpolative-svd-2026-11-09.zip` (3,054 files).
New public module: `src/linalg_interpolative_svd_random.f90`.

## Build and reference environment

- GNU Fortran 14.2.0, checked build: `-std=f2008 -O0 -g -fcheck=all`.
- GNU Fortran 14.2.0, optimized build: `-std=f2008 -O3`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.
- FPM unavailable: the complete `fpm test` suite was not run.

## Results

| Verification | Checked | `-O3` |
|---|---:|---:|
| New native checks | 50,408 passed | 50,408 passed |
| Independent NumPy/SciPy comparison cases | 730 passed | 730 passed |
| SciPy randomized reference calls succeeded | 372 | 372 |
| SciPy randomized reference calls skipped after exceptions | 6 | 6 |
| Previous deterministic interpolative SVD checks | 32,770 passed | 32,770 passed |
| Previous interpolative decomposition checks | 14,944 passed | 14,944 passed |
| Existing real/complex SVD regression | passed | passed |
| Executable example | passed | passed |

The reference script tests reconstruction, factor orthogonality, singular
value ordering, truncation error compared to the exact NumPy SVD, and
adaptive tolerances. For nontrivial fixed-rank tails, the maximum measured
ratio of SciFort's Frobenius residual to the **optimal rank-k Frobenius
residual** was 1.5078. The maximum factor orthogonality residual was
`4.80e-15`. Among fixed-rank tests, the largest relative deviation in the
leading singular values from those of the original matrix was 0.3033; an
approximation is not expected to recover all original singular values.

A raw maximum residual across **all** cases is 1.0 because some fixed-rank
cases deliberately request rank **zero** for a nonzero input. This is
expected, not a numerical failure. The adaptive comparisons passed their
requested relative Frobenius tolerances (up to the test's floating-point
allowance).

The six SciPy reference exceptions per build occurred inside the independent
SciPy randomized code, and were not counted as successful SciPy comparisons.
No SciFort comparison cases were excluded.

## Coverage and limitations

Tested zero-sized, square and rectangular, real and complex, exactly and
numerically rank-deficient, and very small/large coefficient matrices; all
reported native cases passed. A seeded repeatability check confirms
identical results when the same input and seed are supplied. Complex
conjugate orthogonality and the distinction between `V` and `V^H` are
checked explicitly.

Only seven files, including the new module, are new; there are no
changes to existing Fortran source modules or runtime dependencies.
The implementation uses a distinct seeded Gaussian sketch and a relative
**Frobenius** stopping condition. SciPy's randomized ID can use a different
algorithm, random generator and spectral-norm stopping rule, so factor
entries, ranks and random pivot sequences are not expected to be identical.
