# Real Schur decomposition: numerical validation (2026-11-08)

SciFort's new real-output `schur` generic is compared with SciPy 1.17.0,
using GNU Fortran 14.2.0 with `-std=f2018` in both checked
(`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`)
and optimized (`-O3`) builds.

## Independent randomized comparisons

| Seed | Cases per build | Both builds | Outcome |
| --- | ---: | ---: | --- |
| `39711` | 480 | 960 | Passed |
| `52971` | 540 | 1080 | Passed |
| `73019` | 360 | 720 | Passed |
| **Total** | **1380** | **2760** | **Passed** |

Maximum normalized error in the three main seeds (both modes):

| Metric | Maximum |
| --- | ---: |
| Reconstruction `Z*T*Z^T - A` | `5.1e-14` |
| Orthogonality `Z^T*Z - I` | `2.3e-16` |
| Entries forbidden by 1x1/2x2 quasi-triangular structure | `0` |
| Eigenvalue multiset vs. SciPy | `2.2e-14` |

The generator includes orders 0 to 40, dense nonsymmetric, upper/lower
triangular, real symmetric, skew-symmetric, zero and repeated eigenvalues,
scaled matrices (`1e110`/`1e-110`), defective Jordan-like matrices, and
similarity-transformed test cases. It checks all four built-in sort regions,
custom positive-imaginary selection, and a comparison with SciPy's real
Schur `sdim`, which counts whole conjugate pairs. The native tests additionally
check malformed dimensions, nonfinite inputs, order limits, empty matrices,
reconstruction, and real Schur 2x2 blocks.

The reference does **not** require entrywise equality of Schur factors,
which are not unique. Cases extremely close to a selection boundary may
be classified differently due to unavoidable roundoff; such cases are
validated without sorting. The real method uses the existing complex QR
engine and is not a replacement for LAPACK at large orders or with severely
nonnormal, nearly defective matrices.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/linalg_reference/run_real_schur_checks.py \
  --mode both --cases 480 --max-order 40 --seed 39711
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/linalg_reference/run_real_schur_checks.py \
  --mode both --cases 540 --max-order 40 --seed 52971
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/linalg_reference/run_real_schur_checks.py \
  --mode both --cases 360 --max-order 40 --seed 73019
```

The full FPM suite was not run because FPM is unavailable in the validation
environment. Existing complex Schur and neighboring eigenvalue regression
scripts are run separately and reported upon completion.
