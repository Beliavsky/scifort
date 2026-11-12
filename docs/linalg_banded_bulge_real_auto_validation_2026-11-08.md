# Real Givens band reduction and automatic full eigenvalue selection — 2026-11-08

## Implemented changes

- The real `band_to_tridiagonal(ab,d,e,status,q,lower)` overload works entirely
  in real arithmetic. The Hermitian complex overload is unchanged. Both
  operate on the lower band plus one transient bulge diagonal and normalize
  the tridiagonal offdiagonal to nonnegative real numbers.
- The real implementation no longer allocates a complex packed-band copy or
  a complex `n`-by-`n` transformation workspace. If `q` is omitted, real
  workspace is `O(n*kd)`; if it is supplied, the caller owns the `O(n^2)`
  orthogonal output. Both overloads still use unblocked Givens rotations.
- For **full-spectrum eigenvalues only**, both `eigvals_banded` and an
  `eig_banded` call without `z` use bulge chasing automatically when
  `n >= 192` and `kd <= min(8,max(1,n/16))`. The automatic path retries
  the prior dense method on a numerical nonconvergence status.
- Full-spectrum calls requesting eigenvectors retain the old dense default.
  Selected eigenpairs and eigenvalues retain the existing inverse-iteration
  and banded-Sturm defaults. Explicit `method='dense'`, `method='bulge'`,
  `method='inverse'`, and `method='sturm'` keep their documented meanings.

This crossover rule is a conservative heuristic, not a fitted performance
model. Measurements on other hardware and compilers might favor a different
threshold. Full spectra calculated with `method='bulge'` still depend on the
existing tridiagonal QL eigensolver and its numerical limitations.

## Reproducible correctness checks

```bash
python tools/linalg_reference/run_banded_bulge_checks.py --mode both --cases 160 --max-order 128 --seed 24159
python tools/linalg_reference/run_banded_auto_checks.py --mode both --cases 28 --seed 918
```

The new auto-reference harness verifies four real/complex, upper/lower
storage combinations and compares SciFort against SciPy 1.17.0. It verifies
that `auto` returns exactly the result of its intended explicit branch on
both sides of the crossover, and that `eig_banded` without `z` follows the
same choice as `eigvals_banded`.

Fortran native tests: `test/test_linalg_banded_tridiagonal.f90`,
`test/test_linalg_banded_bulge.f90`, and
`test/test_linalg_banded_auto.f90` in both checked and `-O3` builds. In
particular, the auto test checks preservation of the full-eigenvector dense
path and selected-spectrum Sturm path, plus malformed output sizes.

All completed focused batches pass in both modes. In the 28-case independent
seed 918, the largest normalized eigenvalue difference from SciPy was
`1.69e-15`; the 160-case explicit bulge comparison seed 24159 achieved
worst normalized eigenvalue difference `1.69e-15` and eigenvector residual
`3.52e-15`. These are measured quantities, not global error guarantees.

## Reproducible single-environment benchmark

```bash
python tools/linalg_reference/run_banded_bulge_benchmark.py
```

GNU Fortran 14.2.0, `-O3`, Linux; CPU time per call (approximately):

| Order | Half-bandwidth | Real reduction | Complex reduction on same real input | Full eigvals bulge | Full eigvals dense |
|---:|---:|---:|---:|---:|---:|
| 320 | 4 | 0.0035 s | 0.0077 s | 0.0095 s | 0.1000 s |
| 320 | 8 | 0.0071 s | 0.0129 s | 0.0147 s | 0.0897 s |
| 500 | 4 | 0.0068 s | 0.0169 s | 0.0232 s | 0.2160 s |

The benchmark matrices have deterministic sinusoidal diagonal and band
entries. Runs can vary because of CPU load, thermal behavior, compiler
optimizations, and cache effects; the numbers above are not portable
speedup guarantees. The benchmark does not time eigenvector construction.

## Build limits

The full `fpm test` suite was not run because FPM is not installed in this
environment. The standalone tests compile with GNU Fortran 14.2.0 using
`-std=f2018`, with `-fcheck=all` and floating-point traps in the checked
configuration. No additional mandatory numerical libraries are introduced.
