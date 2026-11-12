# Full-spectrum banded inverse iteration (2026-11-08)

This update extends `scifort_linalg_banded_eigen.eig_banded` so that
`method='inverse'` can compute an entire symmetric/Hermitian eigensystem
without first expanding the band input to a dense matrix. The public routine
signature is unchanged. The original dense implementation remains the default
for full-spectrum calls; the new method is **explicit opt-in**.

For real inputs, the new `selected_banded_real_vectors` path stores real banded
LU factors and writes the vectors directly to the caller's real output. It
avoids the complex `(n,n)` temporary that the older partial-spectrum real
wrapper allocated. The existing complex inverse-iteration code writes directly
to the complex output. The eigenvalues are found through the existing banded
Sturm/LDL inertia calculations and bisection.

For order `n` and half-bandwidth `kd`, the new method uses O(n*kd) internal
band storage and O(n*n) required output storage for all the eigenvectors.
Its cost can be substantially greater than a specialized band-to-tridiagonal
solver, particularly for a full spectrum: each inverse-iteration shift is
factorized separately, and two-pass reorthogonalization against preceding
vectors adds up to O(n**3) arithmetic. Because the underlying unpivoted banded
LDL inertia is not fully backward stable near singular shifts, this is **not**
a drop-in replacement for LAPACK's bulge-chasing approach. It returns an
explicit nonconvergence status for unresolved shifts rather than silently
returning eigenvectors without passing its residual test.

## Reproduction

```bash
python tools/linalg_reference/run_banded_full_inverse_checks.py \
    --mode both --cases 300 --max-order 84 --seed 71109
```

Checks compile GNU Fortran in runtime-checked and `-O3` modes, generate
independent real/complex Hermitian band inputs with upper/lower layouts,
compare eigenvalues against SciPy, and validate `A*Z-Z*diag(w)` and `Z^H*Z`.
Cases include diagonal multiplicities, disconnected blocks, clustered
spectra, extreme matrix scales, and both narrow and moderately wide bands.
The existing selected-eigenvector harness checks the same module's prior API.

## Observed results

- 300 comparisons passed in the checked build and 300 in `-O3`, seed 71109,
  for matrices up to order 84. Largest normalized eigenvalue error:
  7.37e-15; largest normalized eigenvector residual: 3.57e-12;
  largest normalized orthogonality error: 2.51e-16.
- An earlier 180-case two-mode run, seed 20261108 and order up to 42, passed.
- The new native tests passed in both builds, including incorrect capacities,
  unused upper-packed cells, complex full spectra, and zero-sized matrices.
- The previous banded eigenvalue, selected-eigenvector and banded-Sturm native
  suites passed in both builds, including their 6000/10000-order stress tests.
- The existing selected eigenvector SciPy comparison suite passed 24 cases
  in each of the checked and `-O3` modes.

These results are environment-specific and do not imply that inverse iteration
is guaranteed to converge for every Hermitian band matrix. A larger
combined selected-eigenvector comparison run timed out and is not counted
as passed. The full FPM suite remains untested in this environment.
