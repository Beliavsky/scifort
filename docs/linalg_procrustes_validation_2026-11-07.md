# Orthogonal Procrustes — validation record (2026-11-07)

## Environment

- GNU Fortran 14.2.0 (Debian build).
- Python with NumPy 2.3.5 and SciPy 1.17.0.
- Checked compilation: `gfortran -std=f2018 -O0 -g -fcheck=all
  -fbacktrace`; native failure-path tests also used
  `-ffpe-trap=invalid,zero,overflow`.
- Optimized compilation: `gfortran -std=f2018 -O3`.
- This environment does not have an `fpm` executable: the full `fpm test`
  suite was **not run**.

## Independent numerical reference

The reference harness runs SciPy's
`scipy.linalg.orthogonal_procrustes(a,b)` for every input and verifies:

1. SciFort returns success and an `n`-by-`n` factor.
2. `||R^H R-I||_F/n` is small (for `n=0`, the normalized error is zero).
3. The Frobenius residual `||A R-B||_F`, computed after common scaling,
   does not exceed the SciPy solution's residual beyond numerical tolerance.
4. The returned `scale_factor` matches SciPy's sum of singular values.
5. For selected nonsingular full-column-rank cases, the transformation
   matrices agree as well, since the polar factor is then unique.

Cases include real/real, complex/complex, real/complex, and complex/real
inputs, exact alignments, rectangular and wide matrices, singular and
rank-deficient matrices, zero matrices, empty dimensions, matrices with
zero cross-covariance, nearly duplicated columns, and scaling from
`1e-140` to `1e140` (on each factor). For rank-deficient problems, it is
important to compare the objective rather than individual matrix entries,
because many optimal transformations may exist.

### Commands and numerical results

All commands below were run from the project root.

| Run | Checked | Optimized | Matrix-order cap |
| --- | ---: | ---: | ---: |
| `--cases 140 --max-order 16 --seed 1117` | 140 passed | 140 passed | 16 |
| `--cases 500 --max-order 32 --seed 4907` | 500 passed | 500 passed | 32 |
| `--cases 500 --max-order 32 --seed 7219` | 500 passed | 500 passed | 32 |
| `--cases 140 --max-order 64 --seed 6317` | 140 passed | 140 passed | 64 |
| `--cases 120 --max-order 96 --seed 8831` | 120 passed | 120 passed | 96 |
| **Total new comparisons** | **1,400** | **1,400** | **2,800** |

Every row was run with
`python tools/linalg_reference/run_procrustes_checks.py --mode both`
followed by the specified arguments. Worst observed metrics over the runs:

- Relative scale-factor disagreement: approximately `1.84e-14`.
- Normalized unitarity error: approximately `5.49e-16`.
- Maximum positive difference in *scaled* objective: approximately
  `1.29e-14`.
- Maximum relative transformation disagreement in selected uniquely
  determined cases: approximately `6.71e-13`.

These are observed test statistics, not general numerical error bounds.

## Native failure paths

The standalone `test/test_linalg_procrustes.f90` program passed with the
checked and optimized compiler configurations. It covers known exact real,
complex and mixed alignments, analytical singular-value sums, zero and
rank-deficient systems, shape mismatch, invalid Jacobi sweep limits,
nonfinite real/complex inputs, and overflow of the true scale. The checked
native run enabled traps for invalid operations, zero division, and overflow.

A case with the source scaled by `1e200` and the target by `1e-200`
correctly returns a scale of `14` (without an intermediate overflow or
underflow). A case whose mathematical scale is `1e400` correctly returns
`linalg_status_no_convergence` instead of a misleading finite result.

## Existing SVD regressions

- `python tools/linalg_reference/run_svd_full_checks.py --mode both
  --cases 70 --max-order 20 --seed 885`: 140 cases passed.
- `python tools/linalg_reference/run_svd_condition_checks.py --mode both
  --cases 70 --max-order 20 --seed 903`: 140 cases passed.

These tests verify that the unchanged SVD infrastructure still operates
correctly. The `example/linalg_procrustes_example.f90` program also compiled
and ran with `-O3`, producing zero alignment error in its simple examples.

## Limitations

The method uses SciFort's intrinsic-only Jacobi SVD and inherits its
convergence, numerical rank and large-matrix performance limitations. There
is no implementation of SciPy's leading batch dimensions. Extremely small
scales may underflow, and unrepresentably large scale values are signaled
as a numerical failure. As with SciPy, singular cross-covariance makes the
minimizer non-unique.
