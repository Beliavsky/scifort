# Singular-value-only SVD and spectral diagnostics validation (2026-11-07)

## Scope

- `scifort_linalg_svd.svdvals`: separate values-only real/complex one-sided
  Jacobi rotations, avoiding singular-vector allocation and basis completion.
- `cond2`: real/complex spectral 2-norm condition numbers, including rectangular
  matrices, infinite results on singular input, invalid/nonfinite diagnostics.
- `matrix_rank`: real/complex SVD numerical rank with optional absolute cutoff,
  including empty/invalid input diagnostics.

The last two are SciFort extensions analogous to NumPy functions, rather than
separate `scipy.linalg` APIs.

## Evidence

GNU Fortran checked (`-O0 -g -fcheck=all`) and optimized (`-O3`) standalone
builds of `src/kinds.f90`, `src/linalg.f90`, and `src/linalg_svd.f90` passed
native `test_linalg_svd_condition.f90` and existing full-SVD native tests.
Python reference testing used SciPy 1.17.0 and NumPy's numerical rank behavior.

Two independent case generators were run using:

```
python tools/linalg_reference/run_svd_condition_checks.py --mode both --cases 350 --seed 903 --max-order 40
python tools/linalg_reference/run_svd_condition_checks.py --mode both --cases 280 --seed 1517 --max-order 50
```

Results: **1,260** passing SciPy/NumPy comparisons across the two compiler
modes, including zero/rank-deficient, real/complex, square/rectangular,
ill-conditioned, and rescaled inputs. Largest observed relative singular-value
spectrum error: `1.13e-15`. Largest relative condition-number error for
well-resolved matrices: `4.43e-9`. Conditions of singular and unresolved
near-singular matrices were checked qualitatively for large or infinite values
rather than compared to unstable finite reference estimates.

Existing regression comparison suites passed in both build modes:

- Full-size SVD: 200 cases per mode (400).
- Economy SVD, pseudoinverse, and rank-deficient least squares: 120 per mode (240).
- `orth`, `null_space`, `subspace_angles`: 100 per mode (200).

The standalone benchmark (`gfortran -O3`), at 120x80 with four repetitions,
reported real `svd` 0.27396 seconds versus real `svdvals` 0.033875 seconds
(8.09x), and complex `svd` 0.44581 seconds versus complex `svdvals`
0.40543 seconds (1.10x). These are single-host CPU timings and not
portable performance guarantees. Each benchmark also checks singular-value
agreement.

## Known limitations

- Intrinsic-only Jacobi SVD is unblocked and slower than optimized LAPACK for
  large matrices; convergence is sweep-limited and nonconvergence is reported.
- Singular values at or below the internal `8*epsilon*smax` floor are set to
  zero, so `cond2` may return infinity for matrices that LAPACK considers
  merely very ill-conditioned.
- Numerical rank can differ near the threshold. Custom `tol` is an *absolute*
  singular-value cutoff, not a relative rcond parameter.
- The complete FPM test suite was not run (FPM not installed in the environment).
