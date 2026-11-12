# SciPy-compatible structural matrix inspection and SVD diagonal construction

The independent MIT-licensed Fortran implementations in
`src/linalg_structure.f90` and `src/linalg_diagsvd.f90` add four SciPy-style
operations without external libraries:

| SciPy API | SciFort module | Fortran call form |
|---|---|---|
| `scipy.linalg.bandwidth` | `scifort_linalg_structure` | `call bandwidth(a, lower, upper, status)` |
| `scipy.linalg.issymmetric` | `scifort_linalg_structure` | `call issymmetric(a, result_value, status, atol, rtol)` |
| `scipy.linalg.ishermitian` | `scifort_linalg_structure` | `call ishermitian(a, result_value, status, atol, rtol)` |
| `scipy.linalg.diagsvd` | `scifort_linalg_diagsvd` | `call diagsvd(singular_values, m, n, sigma, status)` |

## Supported arrays

- `integer(int32)`, `integer(int64)`, `real(dp)`, and `complex(dp)`.
- For `bandwidth`, `issymmetric`, and `ishermitian`, core matrices are
  rank two. A rank-three overload accepts **one leading batch dimension**,
  with a separate scalar result for every member (one-dimensional output
  arrays for rank-three inputs). The batched output must already be allocated.
- For `diagsvd`, an unbatched vector of length `min(m,n)` produces an
  `(m,n)` output matrix. A batched `(batch,min(m,n))` array produces a
  `(batch,m,n)` output, also supplied by the caller.
- Core dimensions may be zero. A square `(0,0)` matrix is symmetric and
  Hermitian. A zero-length *batch* is rejected, matching SciPy 1.17.
- Every argument has explicit intent. `status = 0` is success and
  `status = 1` (`linalg_status_invalid_shape`) indicates an incompatible
  shape, invalid batch size, or invalid tolerance value. On failure,
  output arrays/flags are zero-filled/false.

## Conventions and edge cases

- `issymmetric` compares a matrix with its **transpose**, without
  conjugating complex entries. `ishermitian` uses the **conjugate transpose**
  and tests the imaginary part of complex diagonal entries.
- Default exact comparisons ignore the diagonal of real symmetric matrices,
  including diagonal NaNs; an explicitly positive tolerance invokes the
  approximate comparison and therefore also checks diagonal NaNs.
- Signed infinities match identical signed infinities, and NaNs in an
  inspected off-diagonal entry cause a false result.
- Nonzero `atol` and/or `rtol` are supported for floating and complex
  inputs with the usual `atol + rtol*abs(reference)` condition. Both matrix
  orientations are tested, as with NumPy's whole-matrix `allclose`.
- **SciPy 1.17 performs exact comparisons for integer input**, even when
  tolerances are supplied. SciFort preserves this behavior without losing
  precision when `int64` entries exceed `2**53`.
- For predictable Fortran error handling, optional tolerances must be
  finite and nonnegative. SciPy may exhibit incidental behavior for
  negative or NaN tolerances; SciFort instead returns status 1.
- **SciPy 1.17 rectangular `bandwidth` behavior:** the compiled SciPy
  implementation omits the strictly lower part of the final column and
  the strictly upper part of the final row. In a `(5,1)` column matrix it
  therefore reports `(0,0)` even if a lower entry is nonzero. SciFort
  deliberately matches this behavior. If the *mathematical* bandwidth
  of a rectangular matrix is required, rather than exact SciPy 1.17
  compatibility, this convention needs to be considered.
- `diagsvd` preserves the exact type and values of its singular-value
  vector, including negative values, integers above `2**53`, and complex
  values. It constructs the diagonal and zero-fills the rest without
  numerically computing an SVD.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg_structure, only : bandwidth, ishermitian
use scifort_linalg_diagsvd, only : diagsvd
real(dp) :: a(3,3), s(2), sigma(2,4)
integer :: lower, upper, status
logical :: hermitian

a = reshape([1.0_dp,0.0_dp,2.0_dp, &
             0.0_dp,3.0_dp,0.0_dp, &
             0.0_dp,0.0_dp,4.0_dp], [3,3])
call bandwidth(a, lower, upper, status) ! lower=2, upper=0
call ishermitian(a, hermitian, status) ! false
s = [5.0_dp,2.0_dp]
call diagsvd(s, 2, 4, sigma, status)
```

See `example/linalg_structure_diagsvd_example.f90` and
`test/test_linalg_structure_diagsvd.f90` for a runnable example and native
coverage. A seeded independent SciPy test harness is in
`reference/linalg_structure_diagsvd_driver.f90` and
`tools/compare_linalg_structure_diagsvd.py`.

### Reproducing the reference comparisons

From the project root, with `gfortran` and Python/NumPy/SciPy installed,
compile the four sources (`kinds.f90`, `linalg.f90`,
`linalg_structure.f90`, `linalg_diagsvd.f90`) in that order, placing
`.mod` and `.o` files in an arbitrary build directory. Then compile
`reference/linalg_structure_diagsvd_driver.f90` and link those four objects.
Run:

```bash
python tools/compare_linalg_structure_diagsvd.py /path/to/reference_driver
```

The Python script uses a fixed random seed and compares integer and
floating-point results, status outcomes, and the exact output shape with
SciPy 1.17. A full FPM test can also be run where project dependencies and
FPM are installed.
