# CSR matrix exponential action: implementation and validation (2026-11-04)

## Scope

`scifort_sparse_linalg.expm_multiply` acts with `exp(tA)` on real vectors and
real blocks, using canonical one-based CSR input. `expm_multiply_interval`
returns a regular grid of actions and supports reversed intervals,
endpoint exclusion and a single requested time. Neither API forms a dense
matrix exponential. They are separate from SciFort's sparse iterative solvers.

## Method

The implementation centers on a scalar shift `mu` (the exact CSR trace / n,
or a supplied hint if its shifted infinity norm is smaller), then applies
repeated Taylor actions of `exp((t/s)(A-mu I))` and rescales by
`exp(t mu/s)`. The number of substeps `s` bounds the infinity norm of the
shifted substep by 2; the truncated series stops using an explicit norm bound
on the remaining tail. The time-grid routine anchors at its point nearest
zero and propagates in both directions to avoid underflow from a remote
starting time. This is original Fortran source based on the mathematical
approach of Al-Mohy and Higham (2011), not a port of SciPy's source.

## Reproducing native tests

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
```

Or compile the four standalone objects (`kinds.f90`, `sparse.f90`,
`optimize/lsmr.f90`, and `sparse_linalg.f90`) and link the program
`test/test_sparse_expm_multiply.f90` with GNU Fortran 2018 and `-fcheck=all`.
The native suite includes a rotation, upper-triangular nilpotent operator,
negative/reversed time grids, multiple RHS, independent shifts,
18,000-dimensional diagonal CSR, zero-size matrices, NaN/bad inputs,
overflow, and recovery from interval-start underflow.

## Randomized SciPy and independent dense reference

```bash
python tools/sparse_linalg_reference/compare_expm_multiply_scipy.py --mode both --cases 260 --seed 20261104
python tools/sparse_linalg_reference/compare_expm_multiply_scipy.py --mode optimized --cases 180 --seed 142857
```

The Python driver builds each randomized real matrix, sends its canonical
COO representation (including duplicate values and explicit zeros) to the
Fortran program, and checks **all entries** against `scipy.linalg.expm(t*A) @ B`
for the individual action and each grid time. SciPy 1.17
`scipy.sparse.linalg.expm_multiply` is independently checked for direct
and interval calculations, and discrepancies are reported rather than
allowed to corrupt the validation oracle. The randomized matrices include
nonnormal Jordan-like, triangular, skew-symmetric, symmetric, diagonal,
zero, sparse negative-diagonal, large-shift, and tiny-RHS cases. Test both
checked (`-O0 -fcheck=all`) and optimized (`-O3`) modes.

### Observations and limitations

- `tol` controls a local Taylor remainder bound, not rigorous overall
  forward error; very nonnormal and severely cancelled results can have
  larger relative error. This is not a replacement for a robust dense
  exponential in very small, ill-conditioned systems.
- A truncated Taylor implementation can require many sparse products for
  matrices with very large shifted norms. Work is bounded by an explicit
  scaling budget; extremely large or overflowing actions report status.
- SciPy 1.17's **interval** implementation can disagree with the individual
  actions or fail on some reversed or degenerate grids. The numerical oracle
  here is each time point's independent dense matrix exponential; the script
  prints counts of SciPy interval disagreements and exceptions separately.
  Supplying a very inaccurate `traceA` can also degrade SciPy's direct
  result; SciFort selects a lower-norm shift for explicit CSR inputs.
- The implementation currently handles real CSR matrices and RHS vectors
  and blocks; arbitrary linear operators, complex data, and sparse `expm`
  construction are out of scope.

No third-party solver code or new runtime dependencies were introduced.

## Recorded test results

- Checked GNU Fortran 2018 build: 260 randomized cases passed, with maximum
  entrywise relative-to-output-scale error **4.88e-12** for single actions
  and **1.11e-12** across all time-grid samples against dense `expm`.
- Optimized GNU Fortran build: the same 260 cases passed with the same bounds;
  an additional independent 180-case seed passed with maxima **3.95e-12**
  (single) and **2.22e-12** (grid).
- All 16 sparse native regression programs passed in each build mode,
  including this new regression and all older sparse-matrix, solver,
  Lanczos, Arnoldi, and `svds` tests.
- The 260-case SciPy 1.17 interval comparison reported 102 discrepancies with
  the independently computed dense time-point results and one exception;
  intentionally inaccurate `traceA` hints led to 11 direct-SciPy versus
  dense discrepancies. These are **reference mismatches**, not SciFort
  failures. The separate 180-case optimized seed reported 77 interval
  discrepancies and 10 direct discrepancies. This observation is specific
  to the installed SciPy 1.17.0 build and tested cases; other versions
  may behave differently.
- The complete FPM suite was not executed because FPM was not installed.
  The pre-existing LSMR-callback executable-stack warning remains.
