# Validation of SciFort cosine–sine decomposition — 2026-11-09

## Tested functionality

- Real orthogonal and complex unitary full CSD, including `p != q`
- `swap_sign=.false.` and `.true.`
- `cossin_separate` principal angles and full-factor middle matrix
- Identity, row/column permutations, tiny rotations, random QR matrices
- Rejected invalid `p` and `q`, output shapes, and nonunitary inputs
- Disabled U/VH outputs (zero-sized or filled with zeros)

## Results

All tests used GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5, SciPy 1.17.0.

| Case group | Checked build | Optimized build |
|---|---:|---:|
| Native reconstruction, all partitions | 7,800 / 7,800 | 7,800 / 7,800 |
| Independent SciPy comparison | 3,070 / 3,070 | 3,070 / 3,070 |
| SciPy-reference negative cases | 10 included | 10 included |

Independent SciPy-reference results (maximum *absolute* entrywise errors):

| Numerical metric | Checked `-O0 -fcheck=all` | Optimized `-O3 -march=native` module build |
|---|---:|---:|
| `CS` middle-factor matrix | 8.88e-16 | 8.33e-16 |
| `theta` CS angles | 1.11e-15 | 9.99e-16 |
| Factor reconstruction, `X-U*CS*VH` | 1.67e-14 | 1.76e-14 |
| Factor unitarity defect | 2.00e-15 | 1.89e-15 |

Native generated-factor max reconstruction errors were 1.28e-13 in the
checked build and 1.34e-13 with an optimized SciFort module build.

The standalone reference driver (data I/O) was built without optimization
in both runs; the numerical library module was compiled at the indicated
optimization level. The full FPM suite has not been run (FPM unavailable).
Existing SVD unit tests were also rerun as individual native executables.

No new mandatory dependencies. SciPy/NumPy are **validation-only** tools.
