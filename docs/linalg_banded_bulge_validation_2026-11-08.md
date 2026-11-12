# Banded Hermitian tridiagonalization with Givens bulge chasing (2026-11-08)

## New public routine

`src/linalg_banded_tridiagonal.f90` exports the generic routine
`band_to_tridiagonal` (real symmetric and complex Hermitian):

```fortran
use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
call band_to_tridiagonal(ab, d, e, status, q, lower)
```

- `ab(kd+1,n)`: SciPy upper-packed (`lower=.false.`) or lower-packed
  (`lower=.true.`) symmetric/Hermitian matrix. Only the used triangle
  is read; unused padding entries are ignored.
- `d(n)`, `e(n-1)`: diagonal and nonnegative real offdiagonal of the
  symmetric real tridiagonal matrix `T`.
- `q(n,n)` (optional): orthogonal/unitary transformation satisfying
  `A = Q*T*transpose(Q)` in the real case and
  `A = Q*T*conjg(transpose(Q))` in the complex case.
- `status`: success, invalid shape, or numerical failure. Real arrays
  cannot contain nonfinite stored entries; complex diagonal imaginary
  components are ignored, matching Hermitian input semantics.
- Output arrays are caller-allocated. Optional `q` is omitted when only
  `d` and `e` are required.

The `eig_banded` and `eigvals_banded` generics now accept
`method='bulge'`, routing a full spectrum through band-to-tridiagonal
Givens rotations and the existing tridiagonal QL eigensolver. The result
is in ascending order. For subset options, `bulge` intentionally
calculates the *full* spectrum before filtering; keep using the existing
Sturm/inverse-iteration defaults for efficient partial spectra.

The full-eigenvector and selected-request defaults remain unchanged.
For full-spectrum eigenvalues-only calls, `auto` now selects `bulge` when
`n >= 192` and `kd <= min(8,max(1,n/16))`, falling back to `dense` on
nonconvergence; this is a conservative benchmark-informed rule.

## Algorithm and storage

For lower-banded work storage, annihilate the rightmost unwanted entry
in each successive row with a complex unitary Givens rotation. Each
rotation can create one out-of-band entry at distance `kd+1`. Chase this
bulge toward the bottom right using rotations `kd` positions apart
before annihilating the next unwanted entry in the current row.

Only the lower band plus **one extra bulge diagonal** is stored: working
matrix space is `O(n*kd)`. When `q` is omitted, no `n`-by-`n`
intermediate matrix is allocated. An explicit `q` necessarily uses
`O(n**2)` output storage and involves `O(n**2)` updates per sweep.
The complex subdiagonal is converted to nonnegative real values by a
unit-modulus diagonal phase similarity. Input coefficients are scaled
before elimination to avoid overflow in Givens-rotation formation.

This is an **unblocked** algorithm, not a claim of LAPACK DSBTRD/ZHBTRD
performance. The real overload now uses its own real-arithmetic Givens core with
real band work and real `q`, instead of complex work and transformation
scratch arrays. The complex Hermitian kernel is unchanged. The real
implementation and automatic dispatch were validated separately; see
`docs/linalg_banded_bulge_real_auto_validation_2026-11-08.md`.
For large matrices, clustered or nearly defective numerical problems,
use status and eigenvector residuals to assess accuracy. Full eigenvector
requests still require full `n`-by-`n` output, and `method='bulge'`
retains the existing QL eigenpair limitations.

## Reproducible tests

```bash
python tools/linalg_reference/run_banded_bulge_checks.py \
    --mode both --cases 400 --max-order 128 --seed 14159
```

The independent Python/SciPy reference harness builds under GNU Fortran
with checked (`-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow`) and
optimized (`-O3`) modes. It checks eigenvalue accuracy, eigenvector
residuals `||A Z-Z diag(w)||`, and orthogonality `||Z^H Z-I||`. Families
include real/complex, upper/lower bands, zero/repeated spectra, sparse
blocks, clustered spectra, and scaling by `1e+150` or `1e-150`.

For the 400-case seed above, **400/400 comparisons passed in each
compiler mode**. The largest normalized eigenvalue difference from
SciPy was `1.26e-15`, eigenvector residual `3.71e-15`, and
orthogonality error `4.38e-16` (rounded upward).

Native tests:

- `test/test_linalg_banded_tridiagonal.f90` reconstructs the real and
  complex original matrices from the tridiagonal and `Q`, for upper and
  lower bands, dimensions zero through 65, and bandwidths zero through 9.
  It also tests malformed shapes and unused/nonfinite input handling.
- `test/test_linalg_banded_bulge.f90` exercises full/selected eigenpairs,
  selected eigenvalues, error-path buffers, and large-scale coefficients.
- `example/linalg_banded_bulge_example.f90` demonstrates both the new
  transformation API and `method='bulge'`.

### Single-environment performance example

On an order-320 real symmetric band matrix of half-bandwidth 3,
compiled with GNU Fortran `-O3` on the current Linux test machine,
full-spectrum **eigenvalues-only** runs measured approximately
`0.009 s` for `method='bulge'` and `0.076 s` for the older
`method='dense'` path (a single timing, not a performance guarantee).
This does not imply that the bulge method is faster for all bandwidths,
all eigenvector modes, or all compiler environments.

### Build notes

The complete `fpm test` suite has not been executed in this environment
because FPM is not installed. Individual Fortran tests and the SciPy
comparison harness build directly with `gfortran`.

## References and provenance

The mathematical approach is the standard Givens bulge-chasing scheme
for symmetric band matrices, independently implemented from its
mathematical description. Useful public references include:

- Golub and Van Loan, *Matrix Computations*, 4th ed., section on
  symmetric/Hermitian band reduction.
- LAPACK DSBTRD public documentation, https://www.netlib.org/lapack/explore-html/d0/d62/dsbtrd_8f_source.html
  (algorithm description consulted, **no source code copied**).
- *Comparison of eigensolvers for symmetric band matrices* (2015),
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4617464/ .

See `CODE_PROVENANCE.md` for the licensing statement.
