# Two-sided complex QZ validation — 2026-11-08

This release adds `method='iterative'` to SciFort's generalized complex QZ.
Testing used GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5 and SciPy 1.17.0.
The reference compares (independently of basis/eigenvalue ordering):

- Both normalized reconstruction identities, `A=Q AA Zᴴ` and `B=Q BB Zᴴ`.
- Unitarity of both bases; strict triangularity and real nonnegative `diag(BB)`.
- Generalized homogeneous eigenvalues (`alpha`,`beta`) versus SciPy's QZ,
  using optimal matching in chordal/projective distance.
- For very nonnormal pencils with high forward sensitivity, a separate
  smallest-singular-value **backward error** for each computed eigenvalue.

Completed comparisons before packaging:

| Run | Build | Method | Cases | Outcome |
|---|---|---|---:|---|
| Seed 31209, order ≤40 | checked | iterative | 320 | passed |
| Seed 31209, order ≤40 | -O3 | iterative | 320 | passed |
| Seed 99181, order ≤40 | -O3 | iterative | 180 | passed |
| Seed 88173, order ≤40 | checked | iterative | 130 | passed |
| Seed 741, order ≤12 | -O3 | iterative | 90 | passed |
| Seed 45001, order ≤40 | -O3 | auto | 200 | passed |
| Seed 31209, order ≤32 | -O3 | inverse | 140 | passed |
| Seed 20261108, order ≤32 | checked | iterative | 120 | passed |
| Seed 20261108, order ≤32 | -O3 | iterative | 120 | passed |

Across the 320-case two-build run, the largest normalized reconstruction
error was 1.56×10^-14 and the largest unitarity error 2.61×10^-14.
Highly nonnormal pencils can disagree in *forward* generalized eigenvalues:
the worst projective difference in that run was about 0.002, while the
largest associated pencil eigenvalue backward error was approximately
1.35×10^-17. These forward differences are **not** described as 0.002
accuracy; the backward-error test is the correct numerical distinction.

Additional native testing exercises: complex-conjugate eigenvalue pairs,
real and complex pencils, method compatibility, very different input scales,
singular `B` with infinite eigenvalues, malformed shapes and method strings.

To reproduce the original suites:

```bash
python tools/linalg_reference/run_qz_checks.py --mode both --method iterative --cases 320 --max-order 40 --seed 31209
python tools/linalg_reference/run_qz_checks.py --mode optimized --method auto --cases 200 --max-order 40 --seed 45001
python tools/linalg_reference/run_qz_checks.py --mode optimized --method inverse --cases 140 --max-order 32 --seed 31209
```

The complete FPM suite was not run because FPM is unavailable in this
execution environment. The scope is complex-output QZ for order ≤40;
this is not the full LAPACK `ZHGEQZ` algorithm.
