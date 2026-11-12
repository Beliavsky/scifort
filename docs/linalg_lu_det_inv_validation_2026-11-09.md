# General LU, determinant and inverse validation — 2026-11-09

Independent comparisons used SciPy 1.17.0 / NumPy 2.3.5 with
GNU Fortran checked (`-std=f2018 -O0 -fcheck=all`) and optimized
(`-std=f2018 -O3 -march=native`) builds.

| Test | Checked | Optimized |
| --- | ---: | ---: |
| Native LU/det/inv cases | 2,416 passed | 2,416 passed |
| SciPy LU/det/inv reference cases | 1,888 passed | 1,888 passed |
| LU permutations matching SciPy | 1,595 / 1,598 | 1,595 / 1,598 |
| Determinant comparisons | 151 | 151 |
| Non-near-singular inverse comparisons | 109 | 109 |
| Matching singular inverse rejections | 28 | 28 |
| Ill-conditioned inverse cases excluded | 14 | 14 |
| Maximum relative LU reconstruction error | 8.853e-16 | 7.661e-16 |
| Maximum normalized LU factor difference* | 1.471e-14 | 1.948e-14 |
| Maximum determinant error / Hadamard bound | 7.131e-15 | 3.791e-15 |
| Maximum normalized inverse difference | 9.446e-14 | 9.545e-14 |

*Direct elementwise LU factor comparisons only make sense when pivot
sequences agree. The three differing pivot orderings still reconstruct the
input accurately and satisfy the same triangular factorization identities.

For inverse comparisons, inputs with condition number greater than `1e10`
were excluded from *inverse elementwise* comparisons because mathematically
singular but rounded complex matrices can yield different exact-pivot
classifications between SciPy/LAPACK and the independent Fortran backend.
Their LU factorizations and determinants were still tested, while the native
suite verifies singular inverse status codes for exact zero pivots.

The deterministic reference suite uses real and complex random, diagonal,
zero-column, duplicate-row, strongly scaled and forced-pivot matrices,
plus empty and rectangular shapes up to 24 in each dimension.

Also passed the previous dense triangular/LU, dense Cholesky, QR, and RQ
native regression tests in both compiler configurations. `test_linalg_rq`
reported 5,382 passing native cases in each build.

The complete `fpm test` suite was not run: FPM was not installed in the
execution environment. No new dependencies or third-party source code.
