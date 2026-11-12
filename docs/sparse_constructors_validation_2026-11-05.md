# Sparse constructors: implementation and numerical validation (2026-11-05)

## Scope

- New `src/sparse_constructors.f90`: `csr_diags_array`, `csr_random`,
  `csr_random_array` (alias), `csr_rand`, and `sparse_diagonal`.
- Native regression `test/test_sparse_constructors.f90`.
- Batch reference driver and randomized comparison harness in
  `tools/sparse_reference/`.
- Existing CSR, CSC, sparse indexing and linear-algebra APIs unchanged.

## Results

- GNU Fortran, `-std=f2018 -O0 -g -fcheck=all -fbacktrace`: new native test
  passed, as did affected CSC/diagonal, structural deletion/reshape,
  block/assignment and assembly native tests.
- GNU Fortran, `-std=f2018 -O3`: same five native regression programs passed.
- Full **checked-build** standalone sparse runner completed successfully on
  the final code: **26/26 regression programs passed**.
- Final-code independent SciPy 1.17 comparisons with two random seeds:
  **300 diagonal cases and 150 random-CSR cases per build**, in checked and
  optimized modes. The diagonal comparisons check complete dense values;
  random comparisons check exact nnz, unique valid coordinates, finite values,
  and uniform-value range. Earlier iterations were tested against a wider
  three-seed sample; the final-code counts above are conservative.
- Additional standalone native tests cover exact density 0 and 1; ties-to-even
  rounding; rectangular, empty, large-band and very wide matrices; explicit
  stored zeros; custom values; uniform/normal mean/variance sanity; repeatable
  seeding; unchanged RNG state on invalid input; and duplicate diagonal offsets.
- No Fortran compiler warnings from the new module under `-Wall -Wextra -Werror`.

## Reproducing

```bash
python tools/sparse_reference/compare_sparse_constructors_scipy.py --mode both --cases 150 --seed 20261105
python tools/sparse_reference/compare_sparse_constructors_scipy.py --mode both --cases 150 --seed 112233
python tools/sparse_reference/run_native_sparse.py --mode checked
```

The full optimized `run_native_sparse.py --mode optimized` invocation did not
finish within the execution limit of this environment while compiling the
pre-existing, large `sparse_linalg.f90` module at `-O3`. It is **not** reported
as passed. This does not affect the successfully completed optimized tests
listed above. `fpm test` was not run because FPM is unavailable. The existing
executable-stack linker warning from the LSMR callback implementation remains.

## Algorithm notes and differences from SciPy

- Shape inference and unpadded diagonals follow the documented real-valued
  `scipy.sparse.diags_array` subset. The Fortran output is always canonical
  one-based CSR, not a selectable output format.
- The number of randomly selected positions uses nearest-integer rounding,
  with ties to even. Floyd sampling plus an open-addressed hash set ensures
  no repeated coordinates. Bounded 62-bit draws avoid the limited address
  range of mapping a single binary64 uniform directly to a large integer.
- RNG state comes from SciFort's pre-existing MRG32k3a implementation, not
  NumPy's generator. Exact coordinates/values for matching numeric seeds are
  **not** expected to agree with Python.
- Inputs are real binary64; complex data types and callable custom data
  distributions are not implemented. Callers can supply `data_values(:)`.

## Provenance

The code is an original MIT-licensed SciFort implementation of standard
algorithms, not a source-code translation of SciPy. Existing third-party
attribution and license notices are retained. SciPy is used only as a
numerical comparison target.
