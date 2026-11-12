# Polar decomposition and exponential Fréchet derivative validation

Date: 2026-11-06. Source: `src/linalg_polar.f90` and
`src/linalg_expm_frechet.f90`. GNU Fortran checked `-O0 -fcheck=all` and
optimized `-O3`; SciPy 1.17.0 and NumPy 2.3.5 as numerical references.

`tools/linalg_reference/compare_polar_frechet_scipy.py` compares real and
complex results with `scipy.linalg.polar` and `scipy.linalg.expm_frechet`.
Polar validation checks `(U@P)` or `(P@U)`, Hermitian/symmetric PSD factors,
unitarity of the rectangular factor, and equality of `P` to SciPy's factor.
For rank-deficient input `U` itself need not match SciPy because it is not
unique. The Fréchet test checks both `exp(A)` and the directional derivative;
some cases independently use the upper-right block of an enlarged matrix
exponential rather than only SciPy's specialized derivative algorithm.

| Run | Polar cases | Fréchet cases | Result |
|---|---:|---:|---|
| Checked, seed 20261107, 300 cases | 150 | 150 | Passed |
| Optimized, seed 20261107, 300 cases | 150 | 150 | Passed |
| Checked, seed 20261110, 1,600 cases | 800 | 800 | Passed |
| Optimized, seed 20261110, 1,600 cases | 800 | 800 | Passed |
| Optimized, structured stress batch | 36 | 8 | Passed |

For the 1,600-case batches, maximum relative polar error was `5.329e-15`
(checked) and `5.301e-15` (`-O3`); maximum relative exponential/Fréchet
error was `1.098e-13` (checked) and `1.916e-13` (`-O3`).

The stress batch contains real/complex polar matrices up to `256×64`,
`64×256`, and complex Fréchet inputs up to order 48. It includes exactly
rank-deficient and zero-column matrices. Maximum relative polar error in
this batch was `2.597e-14`; Fréchet error was `3.501e-15`.

The native test checks zero and singular polar matrices, both polar sides,
complex rectangular inputs, invalid side, sweep-limited nonconvergence,
real divided-difference identities, complex derivative finite differences,
zero perturbations, invalid tolerances, and nonfinite inputs.

The full project-wide `fpm test` suite has **not** been run in this runtime.
This is a small-matrix, dependency-free implementation, not a replacement
for LAPACK-backed SVD or a specialized Padé Fréchet kernel on larger systems.
