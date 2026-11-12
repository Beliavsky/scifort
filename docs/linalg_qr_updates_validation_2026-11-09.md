# QR modification validation — 2026-11-09

Environment: GNU Fortran 14.2.0, Python with NumPy 2.3.5 and SciPy 1.17.0.

## Native Fortran tests

`test/test_linalg_qr_updates.f90` passed in both checked (`-O0 -fcheck=all`)
and optimized (`-O3`) builds. Each run checks **1,872** real and complex
reconstruction/orthogonality pairs for rank-one and rank-k QR updates,
row/column insertions, and row/column deletions, full and economic QR,
with matrix sizes 2 through 10 in each dimension.

- Checked build maximum normalized reconstruction/orthogonality difference:
  `1.7144e-15`.
- Optimized build maximum: `1.5543e-15` on the initial suite; the final
  test also includes additional deterministic edge and failure cases.
- Separate tests passed for zero-rank updates, zero-size insertions,
  deleting an entire axis, invalid indices/counts/shapes, nonfinite
  arguments, rank-one updates in the span of `Q`, and dependent-column
  insertion fallback.
- The existing `test/test_linalg_qr.f90` suite passed with `-O3`.

## Independent SciPy comparisons

`validation/compare_linalg_qr_updates.py` independently generates factors
using `scipy.linalg.qr`, calls the corresponding SciPy update/delete/insert
operation to determine its expected output shapes, and sends those same
factor and edit data to the native Fortran driver. Validation compares
**reconstruction**, Q orthonormality, and R upper-triangularity, not
individual Q/R entries (whose phases are nonunique).

| Build | Successful SciPy-comparable cases | Maximum normalized reconstruction error | Maximum orthogonality error |
|---|---:|---:|---:|
| Checked (`-O0 -fcheck=all`) | 1,004 | 1.5476e-15 | 1.3323e-15 |
| Optimized (`-O3`) | 1,004 | 1.5548e-15 | 1.3323e-15 |

Four generated reference calls per build were rejected by SciPy 1.17.0
with a LAPACK `?ormqr`/`?unmqr` illegal-argument error; they are skipped
rather than treated as numerical disagreements. Native Fortran tests
independently exercise rank-k updates for these shapes. Tests that would
delete every row or column are validated natively because SciPy's factor
shape conventions in some such cases are not orthonormal.

## Reproduction

From the repository root:

```bash
python validation/compare_linalg_qr_updates.py --reps 3
python validation/compare_linalg_qr_updates.py --optimized --reps 3
```

The reference script compiles the five needed Fortran sources and the
reference driver in a temporary directory. No new dependencies were
introduced into the SciFort FPM project.

Only the new QR routines and existing QR regression were compiled and
executed directly. The complete `fpm test` project suite was **not** run
in this environment because FPM is not installed.
