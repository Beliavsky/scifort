# Real generalized Schur decomposition and ordering

Declare all four factor matrices as `real(dp)` to select native real output:

```fortran
use scifort_kinds, only: dp
use scifort_linalg_qz, only: qz
use scifort_linalg_ordqz, only: ordqz

real(dp) :: a(n,n),b(n,n),aa(n,n),bb(n,n),q(n,n),z(n,n),beta(n)
complex(dp) :: alpha(n)
integer :: status,nselected

call qz(a,b,aa,bb,q,z,status,alpha,beta)
call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sort='lhp',sdim=nselected)
```

The output satisfies `A=Q*AA*transpose(Z)` and `B=Q*BB*transpose(Z)`, with
real orthogonal `Q` and `Z`. `AA` is quasi-upper triangular with disjoint
1x1 real and 2x2 conjugate blocks. `BB` is upper triangular with nonnegative
diagonal; each conjugate block of `BB` is positive diagonal. Inputs are not
modified. The routines verify both reconstruction identities and orthogonality
on normalized data before reporting success.

Optional homogeneous outputs have different types: `alpha` is complex and
`beta` is real, both length n. Finite eigenvalues are `alpha/beta`; `beta=0`
and nonzero alpha represents infinity. Prefer homogeneous comparisons to
forming potentially overflowing quotients. An indeterminate `(0,0)` pair
returns singular status. Homogeneous pair scales and factor signs/order can
differ from SciPy; compare spectra and reconstruction, not unique factors.

## Sorting

`ordqz` supports the existing `lhp` (default), `rhp`, `iuc`, `ouc` and `none`
rules. Half planes and the unit-circle interior/exterior are strict; zero
belongs to neither open half plane, infinity belongs to the exterior, and
unit-modulus eigenvalues belong to neither strict circle region.

An optional `procedure(ordqz_selector)` callback receives complex homogeneous
`alpha,beta`, using the existing interface even though the real denominator
has zero imaginary part. Selecting either member of a conjugate pair selects
the entire block; `sdim` counts both. Custom predicates should be deterministic
and stable near their selection boundaries. The real path tests the returned
homogeneous values directly, without the complex path's additional beta floor.

`ordqz_factors(aa,bb,q,z,status,sort=...,sdim=...,selector=...)` also accepts
real factors. It checks quasi-triangular structure and orthogonal bases,
reorders whole blocks, and verifies the result. On failure all four input
matrices remain unchanged and `sdim=0`. The original complex specific procedure
and existing complex overloads remain available.

## Implementation and limits

The real path adapts Reference-LAPACK v3.12.1 `DGGHRD`, `DHGEQZ`, `DTGEXC`
and supporting routines under its retained BSD license. Existing native
Householder QR first triangularizes `B`. The real double-shift QZ engine and
block reordering work without forming an inverse. No linked LAPACK/BLAS package
is required. See the [exact source catalog](linalg_real_qz_sources.md).

Order **0..40** is supported, including empty inputs and singular `B` for
regular pencils. Both `method='auto'` (default) and `'iterative'` select the
real QZ engine. `'inverse'` is unsupported for real output and returns invalid
shape/status 1; it remains available for existing complex output.

Statuses retain existing meanings: 0 verified success; 1 incompatible shape,
unsupported size/method/sort; 4 nonfinite input, representability failure,
nonconvergence or rejected swap/check; 5 an indeterminate pencil detected by
the solver. Factors and homogeneous outputs on decomposition failure must not
be interpreted as valid. Not every irregular pencil is guaranteed to be
recognized as such. Extremely clustered/defective eigenvalues can have poor
forward accuracy despite small reconstruction errors; sorting near a boundary
is sensitive to rounding. Underflow/overflow guards can reject outputs that
cannot be represented reliably after rescaling.

Batching, condition estimates, larger public sizes and C ABI expansion remain
separate work. See the [validation record](../validation/real_qz_2026-11-15.md).
