# SciFort restarted shift-invert validation — 2026-11-04

## Algorithm

`eigsh_lanczos(...,sigma=...)` applies the real shifted inverse using
SciFort's existing MINRES; `eigs_arnoldi(...,sigma=...)` applies a complex
shifted inverse with the existing restarted complex GMRES. Both form
Rayleigh--Ritz projections of the inverse operator and retain selected
Ritz vectors across restarts, including soft convergence locking.
For each transformed Ritz eigenvalue `mu`, the reported eigenvalue is
`lambda = sigma+1/mu`. The selected `which` criterion ranks `mu`, not
`lambda`, in SciPy normal shift-invert mode. True residuals
`||A*v-lambda*v||_2` are computed against the original CSR matrix.

The operations use one-based CSR, no dense full inverse, no BLAS/LAPACK,
and no third-party source imports. The `matvecs` counter reports *outer*
operator applications, treating one inner inverse solve as one application;
it does not include the MINRES/GMRES inner sparse product counts.

## Reproduction

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_restarted_shift_scipy.py --mode both --cases 48
python tools/sparse_linalg_reference/compare_restarted_shift_scipy.py --mode both --cases 72 --seed 770077
```

The dedicated native test `test_sparse_restarted_shift` covers an interior
real spectrum, an isolated complex-conjugate pair, method dispatch,
inverse-spectrum LA/SA/SM selection, and invalid/near-singular shifts
and inner-solver parameters. The randomized script compares reported
values to NumPy dense spectra and SciPy `eigsh`/complex-input `eigs`;
it independently checks each true eigenpair residual and normalization.
It intentionally casts the nonsymmetric SciPy input to complex so its
complex shift applies the full inverse (rather than real `OPpart='r'`).

## Scope and limitations

- `which='SM'` still needs `ncv=n` in the restarted implementations.
- The inner solves may be expensive near a singular shift; failed inner
  solves produce a breakdown status, not false convergence.
- The restart retains Ritz vectors with soft locking, not full hard
  deflation or Sorensen-style implicit restart. Difficult clustered or
  highly nonnormal systems may stagnate and have poorly conditioned
  eigenvalues even with small residuals.
- Input CSR remains real-valued. Generalized problems, complex CSR,
  custom inverse operators, and exact convergence deflation are pending.

## Results obtained in this update

- GNU Fortran (`-O0 -g -fcheck=all -fbacktrace`) and `-O3`: all **11**
  standalone sparse regression programs passed, including the new
  `test_sparse_restarted_shift` and revised method-dispatch checks.
- Randomized test seed `20261104`: **48** cases per build passed;
  16 per build required more than one restart cycle.
- Additional seed `770077`: **72** cases per build passed;
  25 per build required more than one restart cycle.
- Total **120** distinct shifted cases per build mode, or **240**
  checked+optimized comparisons; comparison uses SciPy 1.17 and NumPy.
- Legacy unshifted Lanczos comparison printed 16 randomized + 2 difficult
  passing cases per build. Legacy unshifted Arnoldi comparison printed
  20 randomized + 4 clustered passing cases per build. The combined
  command reached the external execution limit after these outputs, so
  it was not counted as a clean process-level success.
- The full FPM suite was not run (FPM not installed); the prior GNU
  executable-stack linker warning from the LSMR callbacks remains.
