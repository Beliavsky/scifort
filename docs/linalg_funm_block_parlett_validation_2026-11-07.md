# Dense `funm`: blocked Schur--Parlett and semisimple clusters (2026-11-07)

This release updates the MIT-licensed `scifort_linalg_funm` module without
changing its public interface. The previous Schur implementation reordered
repeated-eigenvalue clusters only when analytical derivatives were supplied;
inter-cluster values used scalar Parlett substitution. The new implementation
always reorders global clusters. Diagonalizable repeated spectra (including
nonadjacent eigenvalue positions) can now be evaluated with `func` alone.
For clustered blocks requiring derivatives, the existing Taylor expansion
remains. Inter-cluster blocks are solved via upper-triangular Sylvester
substitution in the **block Schur--Parlett recurrence**.

## Reproduction

Environment: GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5 and SciPy 1.17.0.
Each of the first two commands compiles in checked (`-O0 -g -fcheck=all
-fbacktrace`) and optimized (`-O3`) modes. No new BLAS/LAPACK dependency.

```bash
python tools/linalg_reference/run_funm_checks.py --mode both --cases 550 --max-order 9 --reorder-cases 350 --diagonalizable-cases 500
python tools/linalg_reference/run_funm_checks.py --mode both --cases 180 --seed 20261108 --max-order 12 --reorder-cases 250 --diagonalizable-cases 200
python tools/linalg_reference/run_matrix_functions_checks.py --mode both --cases 60
python tools/linalg_reference/run_logarithm_checks.py --mode both --cases 100 --max-order 7
python tools/linalg_reference/run_trigonometric_checks.py --mode both --cases 100 --max-order 8
python tools/linalg_reference/run_expm_condition_checks.py --mode both --cases 100 --max-order 7
```

## Results

| Focused test family | Checked | Optimized |
|---|---:|---:|
| Native `funm` test and standalone example | PASS | PASS |
| General randomized `funm` (two seeds) | 730/730 | 730/730 |
| Reordered-cluster SciPy comparisons (two seeds) | 600/600 | 600/600 |
| **Diagonalizable repeated spectra without derivatives** (two seeds) | **700/700** | **700/700** |
| Jordan, scaling, and analytic stress checks (two runs) | 164/164 | 164/164 |

Total focused numeric checks: **4,388/4,388 passed** across the two builds.
The derivative-free semisimple comparisons include real and complex matrices,
interleaved eigenvalue groups, nonnormal upper-triangular similarities,
unitary dense similarities, near-duplicates and orders up to 40. They compare
against both exact similarity formulas and independent SciPy `expm`, `logm`,
`sinm`, `cosm`, and polynomial references.

Maximum relative Frobenius error on the new semisimple tests was about
`1.63e-14`. The general `funm` tests peaked near `2.41e-14`. Triangular
reordered-cluster tests were around `1.2e-14`, while difficult dense unitary
similarities of defective Jordan matrices reached about `1.57e-7` in the
optimized build, within the separately stated `3e-7` stress threshold. These
hard cases are limited by the existing small-matrix Schur decomposition.

The existing checked and optimized regression runs also passed:

- `sqrtm`/`signm`: 60 comparisons per build
- `logm`/fractional powers: 100 per build
- Six trigonometric/hyperbolic matrix functions: 100 per build
- `expm_cond`: 100 comparisons plus 80 analytic stress cases per build

The full `fpm test` suite was not executed because FPM is not installed in
this environment; focused GNU Fortran builds and SciPy comparisons passed.

## Limitations and provenance

Coupled repeated-eigenvalue blocks still need an analytic derivative
callback. The block diagonal evaluation uses truncated Taylor series, not a
Schur--Padé or general contour-integral solver. Near branch cuts, eigenvalue
collisions, and severely nonnormal matrices can remain numerically fragile.
The commutator residual is not a forward-error bound. This release contains
original MIT-licensed Fortran and test code based on standard Schur--Parlett
and Sylvester identities; SciPy is used only as an independent numerical
reference and no SciPy source has been copied or translated.
