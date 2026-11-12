# Complex N-dimensional COO validation — 2026-11-05

## Implemented

- `src/sparse_coo_nd_complex.f90`: complex128 N-D COO storage, canonicalization,
  explicit-zero handling, shape transformation, dense conversion, conjugation,
  Hermitian transpose and real COO conversion.
- `src/sparse_coo_nd_complex_ops.f90`: complex sums/reductions, singleton-axis
  broadcasting, addition, multiplication, sparse general contraction and
  conjugated inner products.
- `test/test_sparse_coo_nd_complex.f90`: native edge cases, invalid input,
  overflow/transactional handling and enormous sparse shapes.
- `tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py`: streaming
  randomized reference comparisons against SciPy 1.17.0 and NumPy 2.3.5.

## Environment

- GNU Fortran 14.2.0 (`-std=f2018`)
- Python 3, SciPy 1.17.0, NumPy 2.3.5
- Checked compiler mode: `-O0 -g -fcheck=all -fbacktrace`
- Optimized compiler mode: `-O3`

## Results

| Test | Checked | Optimized |
| --- | ---: | ---: |
| Native sparse regression executables (including complex COO) | 42/42 | 42/42 |
| Randomized complex integer-valued SciPy/NumPy comparisons | 520/520 | 520/520 |
| Randomized complex Gaussian-valued SciPy/NumPy comparisons | 520/520 | 520/520 |
| Documentation example | passed | not separately run |

**Total: 2,080 passing randomized comparisons** across both compiler modes,
plus the 42 native sparse regression programs run in each mode. An initial
additional 130-case checked random run also passed.

Native tests specifically cover duplicate coordinates, explicit zeros,
complex conjugation, Hermitian transpose, C/F reshape and dense conversion,
real-to-complex conversion, rejection of nonzero imaginary parts when converting
back to real, axis reductions, `vdot` conjugation versus `tensordot` semantics,
large sparse shapes, storage caps, rejected nonfinite values, and overflow.
The comparison script tests rank 1–5 tensors, arbitrary contraction axes,
broadcasting, empty dimensions and repeated coordinates against independent
SciPy COO/NumPy results.

## Reproduce

```console
python tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py --mode both --trials 520 --seed 12311
python tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py --mode both --trials 520 --seed 242311 --floating
python tools/sparse_reference/run_native_sparse.py --mode both
```

The project-wide FPM suite was **not** run because FPM is unavailable in the
test environment; these 42 executables are the project's standalone sparse
suite. The comparison driver links only SciFort source and standard Fortran,
using NumPy/SciPy externally as the reference oracle.
