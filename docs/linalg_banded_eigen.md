# Symmetric tridiagonal and Hermitian band eigenproblems

SciFort implements four `scipy.linalg`-compatible mathematical operations
in `src/linalg_banded_eigen.f90`, module `scifort_linalg_banded_eigen`.
They use original, MIT-licensed SciFort code without external LAPACK or BLAS.

## Public interfaces

```fortran
use scifort_linalg_banded_eigen, only : &
    eigh_tridiagonal, eigvalsh_tridiagonal, eig_banded, eigvals_banded

call eigh_tridiagonal(d, e, w, status, z, max_iterations, &
    subset_by_index, subset_by_value, num_selected, tol, eigenvector_method)
call eigvalsh_tridiagonal(d, e, w, status, max_iterations, &
    subset_by_index, subset_by_value, num_selected, tol)
call eig_banded(ab, w, status, z, lower, max_iterations, &
    subset_by_index, subset_by_value, num_selected, tol, method)
call eigvals_banded(ab, w, status, lower, max_iterations, &
    subset_by_index, subset_by_value, num_selected, tol, method)
```

The `z` output is optional, as are `max_iterations` and `lower`; use
keyword arguments when skipping `z`. All outputs are caller-allocated.

- `d(n)`: real symmetric tridiagonal diagonal (double precision).
- `e(n-1)`: first subdiagonal, the same as the first superdiagonal.
- `ab(kd+1,n)`: diagonal-ordered real symmetric or complex Hermitian band
  storage. `kd` is the bandwidth; rows not corresponding to matrix entries
  are ignored, even if they contain undefined values.
- `w(n)`: ascending real eigenvalues, including multiplicities, for full-spectrum
  calls. In selection mode, `w` is a caller-owned capacity buffer; the first
  `num_selected` elements contain the result, and unused entries are zero.
- `z(n,n)` for a full-spectrum call, or `z(n,capacity)` for a selected
  tridiagonal or general-banded call: corresponding orthonormal eigenvectors in columns.
  `capacity >= num_selected`, and unused columns are zeroed.
- `lower=.false.` (default): SciPy upper-packed band convention.
  `ab(kd+1+i-j,j) = A(i,j)` for `max(1,j-kd) <= i <= j`.
- `lower=.true.`: lower-packed convention.
  `ab(1+i-j,j) = A(i,j)` for `j <= i <= min(n,j+kd)`.
- `max_iterations=120` by default: maximum implicit-shift QL steps per
  eigenvalue or bisection steps per selected eigenvalue; the new banded
  Sturm path defaults to 140 steps.
- `subset_by_index=[first,last]`: optional inclusive **one-based** index limits;
  unlike SciPy's zero-based indices. Cannot be combined with value selection.
- `subset_by_value=[low,high]`: optional **(low,high]** interval, where
  `low < high`; positive and negative infinity endpoints are accepted.
- `num_selected`: optional number of matching eigenvalues. If the caller's
  output is too small, returns the required count and invalid-shape status.
- `tol`: optional nonnegative absolute bisection tolerance. Effective
  resolution is bounded below by a multiple of machine epsilon times the
  largest absolute tridiagonal or banded coefficient. An excessively loose tolerance
  can make the requested eigenvalues too inaccurate for inverse iteration.
- `method`: for `eig_banded` with eigenvectors, `auto` (default) uses the
  band-preserving inverse-iteration solver for selected eigenpairs, with a
  dense fallback for `n <= 96` if iteration does not converge. `inverse`
  forces the low-memory solver for **selected or full** eigenpairs; `dense`
  forces the original dense method and optionally filters its results. The
  default full-spectrum call with eigenvectors continues using the original
  dense algorithm. For full-spectrum *eigenvalues only* and sufficiently narrow
  bands (`n >= 192` and `kd <= min(8,max(1,n/16))`), `auto` selects the
  Givens bulge method with a dense fallback on nonconvergence.
  `method='bulge'` explicitly selects band-preserving Givens reduction,
  for both full-spectrum and selected requests. In selected mode it computes
  and filters the full spectrum; use `auto`/`inverse` for partial eigenpairs.
  With `method='inverse'` but no `z` argument, full-spectrum eigenvalues
  are computed by banded Sturm bisection only. In the real inverse-iteration
  path, real output is written directly without a complex n-by-k workspace.
- `method`: optional `eigvals_banded` string: `auto` (default) uses banded
  Sturm bisection when a subset is requested, Givens bulge reduction for
  qualifying narrow-band full spectra, and the dense solver otherwise.
  The same threshold and nonconvergence fallback apply; `sturm` always uses banded
  bisection (including full-spectrum calls), and `dense` forces the old
  dense reduction, including selection by filtering. `tol` is used by
  the Sturm method only. `method='bulge'` reduces to tridiagonal form
  without dense input expansion, then computes all eigenvalues with QL;
  it does not provide partial-spectrum speedups.
- `eigenvector_method`: optional `eigh_tridiagonal` string: `auto` (default)
  uses bisection plus inverse iteration for selected eigenvectors, with a
  full-QL fallback if the iteration fails and `n <= 256`; `inverse` forces
  the memory-efficient path, returning nonconvergence on failure; `ql` uses
  the original full-QL method. It has no effect on unselected eigenproblems.
- `status`: `linalg_status_success=0`, `linalg_status_invalid_shape=1`,
  or `linalg_status_no_convergence=4`. Wrong shapes or an iteration limit
  below one produce invalid-shape status. Nonfinite coefficients, failure
  to converge or nonrepresentable finite eigenvalues produce nonconvergence.
  No exception is raised within the library.

The complex Hermitian input uses only the stored triangle. Imaginary
parts of diagonal entries are ignored, as appropriate for a Hermitian
matrix. All finite stored off-diagonal elements are conjugate-reflected
into the nonstored triangle.

For example, a real upper-banded order-4 tridiagonal matrix with diagonal
2 and off-diagonal -1 uses:

```fortran
real(dp) :: ab(2,4), eigenvalues(4), vectors(4,4)
integer :: status
ab = 0.0_dp
ab(2,:) = 2.0_dp
ab(1,2:4) = -1.0_dp
call eig_banded(ab, eigenvalues, status, vectors)
```

See `example/linalg_banded_eigen_example.f90` for a runnable example.

## Band-preserving full-spectrum Givens reduction

The independent `scifort_linalg_banded_tridiagonal.band_to_tridiagonal`
module implements real symmetric and complex Hermitian band reduction,
optionally returning the orthogonal/unitary transformation `Q`. It uses
successive two-sided Givens rotations with bulge chasing; the working
matrix occupies `O(n*kd)` storage rather than an `n`-by-`n` dense
matrix. For full eigenpairs, transformation accumulation and the caller's
output eigenvector array still consume `O(n**2)` memory.

```fortran
use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
call band_to_tridiagonal(ab,d,e,status,q,lower=.true.)
call eig_banded(ab,w,status,z,lower=.true.,method='bulge')
call eigvals_banded(ab,w,status,lower=.true.,method='bulge')
```

The `bulge` method remains explicitly selectable for all band shapes.
The `auto` choice now uses it for full-spectrum eigenvalue-only requests
with `n >= 192` and `kd <= min(8,max(1,n/16))`. Full eigenvector calls
and selected requests preserve their earlier automatic selection behavior;
`dense`, `inverse`, and `sturm` still force their documented methods.
The real overload uses a native real-arithmetic Givens kernel rather than
allocating complex work arrays. See
`docs/linalg_banded_bulge_validation_2026-11-08.md` for algorithms,
input limits, validation, and the runnable example.

## Algorithms and numerical limits

For full-spectrum calls, `eigh_tridiagonal` first scales the diagonal and
first offdiagonal by their maximum magnitude and applies Wilkinson-shift
implicit QL iteration, accumulating rotations when eigenvectors are requested.
Results are sorted ascending. The unselected values-only form avoids the
n-by-n eigenvector matrix.

**Partial-spectrum path:** with either `subset_by_index` or `subset_by_value`,
SciFort uses a scaled Sturm recurrence, Gershgorin bounds, and bisection to
compute only the requested eigenvalues. For `k` values from an order-`n`
tridiagonal matrix, eigenvalue-only calls use **O(n) workspace** and
O(n*k*log(width/tolerance)) time. They never construct a dense matrix.

When selected eigenvectors are requested, SciFort performs inverse iteration
on tridiagonal systems with local row pivoting, a small regularizing shift,
and two-pass modified Gram-Schmidt against previously accepted vectors.
The result requires **O(n*k)** storage for the requested eigenvectors and
O(n) temporary workspace, not the O(n**2) full-QL eigenvector matrix.
Each inverse iteration is O(n); reorthogonalization can cost O(n*k**2) in
severely clustered spectra. Eigenvector signs or bases can differ from SciPy.
A zero-offdiagonal shortcut produces exact, stable-sorted coordinate vectors.

`eigenvector_method='auto'` falls back to the previous full-QL calculation
when inverse iteration does not converge and `n <= 256`. For large `n`,
or with `eigenvector_method='inverse'`, failure is returned explicitly
instead of silently allocating a dense eigenvector matrix. The full QL
method is also directly available with `eigenvector_method='ql'`.

Example:

```fortran
real(dp) :: d(10000), e(9999), w(5)
integer :: status, nfound

d = 2.0_dp
e = -1.0_dp
call eigvalsh_tridiagonal(d, e, w, status, &
    subset_by_index=[5000,5004], num_selected=nfound)
```

To request just five eigenvectors too, allocate `z(10000,5)` and use:

```fortran
real(dp) :: z(10000,5)
call eigh_tridiagonal(d, e, w, status, z, &
    subset_by_index=[5000,5004], num_selected=nfound, &
    eigenvector_method='inverse')
```

See `example/linalg_tridiagonal_vectors_example.f90` for a runnable example.

In value-selection mode, intervals are half-open `(low,high]`. Very close to
an eigenvalue, floating-point rounding can make boundary classification
ambiguous; SciPy and SciFort may disagree in those difficult cases. An
empty selection succeeds with `num_selected=0`; malformed intervals and
insufficient output capacity return `linalg_status_invalid_shape`.

### General Hermitian banded eigenvalues: memory-efficient selection

`eigvals_banded` now has a dedicated `method='sturm'` algorithm for real
symmetric and complex Hermitian band matrices, with both SciPy packed-band
layouts. It constructs an O(n*kd) **lower-band workspace**, scales the input,
computes Gershgorin eigenvalue bounds, and evaluates the negative inertia of
`A - lambda*I` via unpivoted Hermitian band LDL^H elimination. A bisection
search then locates only the requested eigenvalues, without allocating an
n-by-n dense array or constructing eigenvectors. The workspace is
**O(n*kd)** and the time is approximately
**O(n*kd^2*k_selected*log(accuracy))**, where `k_selected` is the number
of eigenvalues requested.

For example, on a very large narrow-band matrix:

```fortran
real(dp) :: ab(3,10000), w(5)
integer :: status, nfound

ab = 0.0_dp
ab(1,:) = 4.0_dp
ab(2,1:9999) = -1.0_dp
ab(3,1:9998) = -0.25_dp
call eigvals_banded(ab,w,status,lower=.true., &
    subset_by_index=[5000,5004],num_selected=nfound)
```

If no subset is requested, the default `method='auto'` keeps the historical
full-spectrum dense-backed solver for compatibility and typically better
performance on modest orders or wide bands. To force the banded algorithm for
the full spectrum, use `method='sturm'`; it can be slower than dense reduction
when requesting most eigenvalues of a wide-band matrix. `method='dense'`
always forces the old implementation, optionally selecting from its complete
computed eigenvalues. Both methods retain the previous positional arguments.

**Numerical caution:** Unpivoted indefinite LDL^H is not generally
backward-stable near singular shifts. The implementation regularizes tiny
pivots, checks for nonfinite growth and reports `no_convergence` on failed
factorizations. Rare clustered or boundary eigenvalues can lead to inaccurate
inertia counts even without overflow. As with other bisection schemes,
eigenvalue selection at an endpoint within rounding error is ambiguous.
For demanding ill-conditioned problems, compare against `method='dense'` or
an external LAPACK implementation. The `eigvals_banded` path returns **eigenvalues only**; the separately
implemented `eig_banded` selected-eigenvector path described below uses
pivoted banded LU and inverse iteration to obtain the requested vectors
without dense expansion.

### General Hermitian band eigenvectors: band-preserving selection

For `eig_banded` with `subset_by_index` or `subset_by_value`, a new selected-
eigenvector solver uses the banded Sturm eigenvalue brackets above, followed
by *partial-pivoted banded LU* and inverse iteration on the shifted matrix.
The solver supports real symmetric and complex Hermitian matrices, both
packed-band conventions, several right eigenvectors and two-pass modified
Gram-Schmidt orthogonalization. It also handles exactly diagonal and zero
matrices with coordinate eigenvectors, including repeated eigenvalues.

```fortran
real(dp) :: ab(3,10000), w(4), z(10000,4)
integer :: status, nfound
! Fill lower-packed bands in ab, with diagonal in row 1.
call eig_banded(ab,w,status,z,lower=.true., &
    subset_by_index=[101,104],num_selected=nfound,method='inverse')
```

For half-bandwidth `kd`, order `n`, and `m` selected eigenvectors, the new
method uses **O(n*kd+n*m) memory**, including pivot fill. For real symmetric
input the LU solves use real arithmetic rather than promoting the entire
factorization to complex arithmetic. Each distinct shift requires O(n*kd**2) banded LU work; each inverse
iteration after factorization takes O(n*kd) triangular-solve work.
Two-pass reorthogonalization adds O(n*m**2) total work. Eigenvalues are first
computed by bisection. The default full eigensystem still uses dense Householder/Hessenberg
reduction and full tridiagonal QL, with O(n**2) memory and O(n**3) time.
**Opt-in full-spectrum inverse iteration** also avoids a dense *input*
expansion and keeps its extra workspace O(n*kd); the required output
`z(n,n)` necessarily uses O(n**2) storage. This is not a band-to-tridiagonal
bulge-chasing reduction and can be slow for many eigenvectors.

`method='auto'` uses banded inverse iteration for selections, falling back
to the previous dense implementation only for `n <= 96` when convergence
fails; `method='inverse'` never falls back, whether requesting a subset
or the entire eigensystem. Any malformed arguments return
`invalid_shape`, not a fallback. A finite eigenvalue and orthogonalization
residual are checked before an eigenvector is accepted. Poorly conditioned
shifts, nearly degenerate eigenspaces, loose bisection tolerances, and
uncertain unpivoted LDL inertia may produce nonconvergence or reduced
accuracy. Unlike optimized LAPACK band-to-tridiagonal algorithms, this is
not a complete backwards-stability guarantee. Numerical interval-boundary
classification can be ambiguous close to an eigenvalue.

See `example/linalg_banded_full_inverse_example.f90` and
`docs/linalg_banded_full_inverse_validation_2026-11-08.md` for full-spectrum
API and validation. See `example/linalg_banded_vectors_example.f90` and
`docs/linalg_banded_vectors_validation_2026-11-08.md` for reproducible
native and independent SciPy comparisons.

See `example/linalg_banded_sturm_example.f90` for a runnable partial-spectrum
example and `docs/linalg_banded_sturm_validation_2026-11-08.md` for tests.

The input must describe a finite symmetric/Hermitian matrix. Eigenvectors
in degenerate eigenspaces are not unique and can differ from SciPy's,
even when the eigenvalues, residuals and orthogonality agree. Near multiple
or extremely ill-scaled eigenvalues the numerical accuracy can deteriorate.

## Testing

```sh
python tools/linalg_reference/run_banded_eigen_checks.py --mode both \
    --cases 600 --max-order 96 --seed 1701
```

The original banded-eigenvalue harness compiles native tests and checks ascending eigenvalues,
normalized eigen-equation residuals and orthogonality against SciPy, for
real/complex upper and lower packed bands and real tridiagonal matrices.
It uses `OPENBLAS_NUM_THREADS=1` in the reference comparison process.

The implementation follows standard implicit-shift QL, Hermitian
Householder-reduction and unitary-diagonal-phase mathematics described in
Golub and Van Loan, *Matrix Computations*, 4th ed., chapters 7 and 8.
SciPy's public reference APIs defined the intended mathematical behavior;
no SciPy/LAPACK numerical implementation code was copied.

### Partial-spectrum validation (2026-11-08)

```sh
python tools/linalg_reference/run_tridiagonal_selection_checks.py --mode both \
    --cases 720 --max-order 128 --seed 20261108
```

See `docs/linalg_tridiagonal_selection_validation_2026-11-08.md` for tests,
including independent SciPy interval and index selection, repeated spectra,
empty selections, exact endpoints, and an order-10000 narrow slice.

### Selected eigenvector validation (2026-11-08)

```sh
python tools/linalg_reference/run_tridiagonal_vectors_checks.py --mode both \
    --cases 600 --max-order 512 --seed 15803
```

The new harness checks real symmetric tridiagonal selected eigenpairs
against SciPy's reference `eigh_tridiagonal`, including reconstruction
residuals, orthogonality and unambiguous simple eigenvectors. The native
check includes an order-10000 five-eigenvector solve, exact multiplicities,
error statuses and explicit QL/inverse selection. See
`docs/linalg_tridiagonal_vectors_validation_2026-11-08.md` for results.

### Partial-spectrum Hermitian band validation (2026-11-08)

```sh
python tools/linalg_reference/run_banded_sturm_checks.py --mode both \
    --cases 240 --max-order 96 --seed 1092
```

This harness verifies real and complex input, both packed-band conventions,
full and selected spectra, diagonal/repeated/clustered cases and extreme
input scales against `scipy.linalg.eigvalsh`. Native tests independently
exercise error statuses and an order-10000 narrow-band stress case.
