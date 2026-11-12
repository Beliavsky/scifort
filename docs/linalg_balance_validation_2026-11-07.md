# Matrix balancing validation (November 7, 2026)

## Build and numerical reference

- Build environments: GNU Fortran with `-std=f2018`, checked `-O0 -g
  -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`, and optimized
  `-O3`.
- Reference: installed SciPy's `scipy.linalg.matrix_balance`, wrapping
  LAPACK `xGEBAL`. No SciPy or LAPACK implementation source is included.
- The native suite exercises the transform matrix and optional separate
  outputs, power-of-two factors, status values, permutation-only and
  scaling-only cases, nonfinite and zero-size matrices, real/complex
  handling and scales near `1e280` and `1e-280`.

## Randomized independent comparisons

Tests deliberately do **not** compare the transformation entries one-for-one;
there can be multiple mathematically valid balances of the same input.
Instead, they verify the identity for the reported factors and permutation,
the trace invariant, and that scaling uses positive powers of two. The
independent SciPy computation is verified against its reported factors as a
separate control; off-diagonal norms are tracked as balancing diagnostics.

| Run | Cases | Maximum normalized similarity error | Maximum normalized trace error |
| --- | ---: | ---: | ---: |
| Seed 434, order at most 32, checked | 180 | 1.666e-16 | 9.053e-16 |
| Seed 434, order at most 32, `-O3` | 180 | 1.666e-16 | 9.053e-16 |
| Seed 20261108, order at most 96, checked | 650 | 1.495e-16 | 7.058e-16 |
| Seed 20261108, order at most 96, `-O3` | 650 | 1.495e-16 | 7.058e-16 |
| Seed 20261207, order at most 120, checked | 850 | 1.601e-16 | 6.164e-16 |
| Seed 20261207, order at most 120, `-O3` | 850 | 1.601e-16 | 6.164e-16 |

The maximum ratio of SciFort's balanced off-diagonal 1-norm to SciPy's
balanced off-diagonal 1-norm in the 650-case run was about `1.75`;
for the 850-case, order-120 run, it was about `2.53`.
Two SciPy reference cases with very large scaling factors produced an
internal SciPy `RuntimeWarning` while casting its LAPACK balancing metadata;
those cases still returned valid permutations, positive finite scaling
factors, and balanced matrices satisfying the independent similarity check.
Thus the method is not numerically identical to LAPACK `xGEBAL`, even though
the similarity identities and scale constraints pass. Balancing can be
unnecessary for already-balanced or triangular matrices.

Across the six checked/optimized comparison runs, all **3,360** randomized
SciPy-reference cases passed the similarity and trace-invariant checks.

The complete `fpm test` suite was not run in this environment. The new
standalone native test and direct GNU Fortran comparison harness were run in
both configurations; focused neighboring regression tests were run separately.
