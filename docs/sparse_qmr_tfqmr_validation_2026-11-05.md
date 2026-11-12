# 2026-11-05 — CSR QMR and TFQMR numerical validation

## New coverage

Two nonsymmetric Krylov solvers were added to `scifort_sparse_linalg`:

- `qmr`: adapted from SciPy v1.17.0 `scipy.sparse.linalg.qmr`, with
  explicit canonical CSR transpose, an optional right Jacobi inverse,
  SciFort-style error codes, and checked original-system convergence.
  The full required BSD-3-Clause notice is in
  `LICENSES/scipy-qmr-BSD-3-Clause.txt`.
- `tfqmr`: adapts SciFort's existing BSD-3-Clause matrix-free TFQMR
  recurrence (`src/optimize/krylov_iterative.f90`), originally credited to
  SciPy's BSD-3-Clause implementation, using consistent right preconditioning
  with `A*M` and an optional Jacobi or ILU(0) inverse. The full license notice
  is retained in `LICENSES/scipy-tfqmr-BSD-3-Clause.txt`.

Both return `info=0` only after checking the true Euclidean residual
`norm2(b - A*x)`. The `tfqmr` iteration count denotes half-iterations.
Only canonical real binary64 CSR matrices are supported, not complex matrices
or arbitrary linear operators. The implementation has no new external
build or runtime dependency.

## Environment and commands

- Compiler: GNU Fortran 14.2.0
- Python: 3.13.5; SciPy 1.17.0 and NumPy installed locally
- Checked compiler options: `-std=f2018 -O0 -g -fcheck=all -fbacktrace
  -ffpe-trap=invalid,zero,overflow`
- Optimized compiler options: `-std=f2018 -O3`

From the repository root:

```bash
python tools/sparse_reference/compare_sparse_qmr_tfqmr_scipy.py --mode both
python tools/sparse_reference/run_native_sparse.py --mode both
```

The independent Python driver recompiles the four required SciFort modules
and the numerical Fortran test driver for each mode. It constructs 180
random, nonsymmetric, strictly row-diagonally dominant systems of sizes
2–40 in each build, including variable row scalings and nonzero starting
guesses. Its cases exercise QMR with/without Jacobi, TFQMR with/without
Jacobi, and TFQMR with ILU(0). Every SciFort solution is compared against
`numpy.linalg.solve`, with additional independent SciPy solver comparisons.
For TFQMR with Jacobi, SciPy is applied to the mathematically equivalent
right-preconditioned problem `A*D^-1`, rather than passing SciPy's `M`
argument (which uses a different product order). The comparison checks
the reported residual against the independently recomputed **true** residual,
the actual tolerance condition, and the recovered solution components.

## Results

| Check | Checked | Optimized |
|---|---:|---:|
| Randomized QMR/TFQMR comparisons | 180/180 | 180/180 |
| Native sparse regression executables | 29/29 | 29/29 |

The new `test/test_sparse_qmr_tfqmr.f90` covers residuals and solutions on a
nonsymmetric and indefinite test matrices, Jacobi/ILU(0), exact starting guesses, empty
systems, zero RHS, invalid dimensions/limits/tolerances, simultaneous
preconditioner rejection, zero-operator breakdown, missing Jacobi diagonals,
NaNs, and limited-iteration nonconvergence. A previous sparse BiCG/CGS and
Kronecker-sum SciPy comparison was also rerun with 30 cases of each kind in
checked mode and passed.

`fpm test` was **not run** because the `fpm` executable is not installed.
The standalone native harness builds only SciFort's sparse-module dependency
subset, so this result does not certify the full package or optional
external FPM dependencies.

## Numerical limitations

QMR may suffer a Lanczos breakdown even for nonsingular matrices; TFQMR may
encounter small denominators or stagnation on difficult nonsymmetric
operators. Such failures yield a negative status, not false convergence.
They are distinct from capped nonconvergence (positive status). No use of
ILU-transpose solves is supported by QMR; TFQMR's ILU(0) factors use the
existing SciFort inverse-application API. Complex arithmetic and general
preconditioner callbacks are not implemented.
