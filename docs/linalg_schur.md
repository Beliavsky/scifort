# Real and complex Schur decomposition (`scipy.linalg.schur`)

`src/linalg_schur.f90` provides a public `schur` generic accepting real or
complex double-precision square matrices. With **complex output matrices**,
it constructs the complex Schur form

\[
A=ZTZ^H,
\]

where `T` is upper triangular and `Z` is unitary. The routine uses the existing
MIT-licensed `dense_complex_schur` implementation; a fast exact triangular
path avoids unnecessary QR perturbations of triangular inputs.

```fortran
use scifort_kinds, only : dp
use scifort_linalg_schur, only : schur

! Real and complex coefficient matrices are both supported.
complex(dp) :: t(n,n), z(n,n)
integer :: status, nselected

call schur(a, t, z, status)
call schur(a, t, z, status, sort='lhp', sdim=nselected)
call schur(a, t, z, status, sort='iuc', sdim=nselected)
call schur(a, t, z, status, selector=positive_imag, sdim=nselected)
```

The user-defined `selector` must be a procedure matching the public
`schur_selector` abstract interface: a logical function with one
`complex(dp), intent(in)` argument for an eigenvalue in the **original input
units**. It should be deterministic and remain stable under tiny eigenvalue
rounding errors. Do not pass a non-`none` `sort` together with `selector`.

Supported `sort` strings:

| Value | Selected eigenvalues moved to the leading Schur block |
| --- | --- |
| `none` (default) | No reordering; `sdim=0` |
| `lhp` | Real part strictly less than zero |
| `rhp` | Real part greater than or equal to zero |
| `iuc` | Modulus less than or equal to one |
| `ouc` | Modulus strictly greater than one |

`sdim` is optional and reports the number of reordered eigenvalues. The
selected invariant subspace is spanned by the first `sdim` columns of `Z`.
The eigenvalue order **within** each group is unspecified.

**Real Schur output for real matrices** is also available. Declare `t` and
`z` as `real(dp)` instead of `complex(dp)`; the generic selects the real
output overload without changing the existing complex overloads:

```fortran
real(dp) :: a(n,n), t(n,n), z(n,n)
integer :: status, nselected

call schur(a, t, z, status)
call schur(a, t, z, status, sort='lhp', sdim=nselected)
```

It computes `A = Z*T*transpose(Z)` where `Z` is real orthogonal and `T` is
real **quasi-upper triangular**. Its diagonal has 1x1 blocks for real
eigenvalues and 2x2 blocks for conjugate complex pairs. The strictly lower
part of `T`, apart from the allowed 2x2 subdiagonal entries, is zero. A
returned pair always stays within the same selected/unselected Schur region.

The same `sort` values are available for real output. The optional Fortran
`selector` uses the existing `schur_selector` interface, receiving **one
complex scalar**, unlike SciPy's Python callback for real output (which takes
separate real and imaginary arguments). If either member of a conjugate pair
is selected, both are moved to the leading Schur block, with `sdim` counting
two eigenvalues. For numerically real singletons, imaginary roundoff is
removed before evaluating a custom predicate.

The real-output implementation first calculates the existing complex Schur
form on a normalized input, groups conjugate eigenvalues by stable adjacent
unitary Schur swaps, orders the selected 1x1/2x2 groups, and constructs nested
real invariant subspaces by twice-orthogonalized real and imaginary Schur
vectors. A final real orthogonal similarity constructs `T`. It checks
unitarity, the quasi-triangular structure, and the reconstruction residual,
returning nonconvergence status if a reliable real form cannot be obtained.

Both output styles support matrix order 0 through 40 (as set by
`matrix_function_max_order`). Caller-allocated output matrices must both
have shape `(n,n)`. SciPy's batched arrays, overwrite controls, and some other
options are not implemented.

Status values from `scifort_linalg`:

- `linalg_status_success` (0): decomposition and requested reordering succeed.
- `linalg_status_invalid_shape` (1): invalid dimensions, order, sort code, or
  mutually exclusive selection options.
- `linalg_status_no_convergence` (4): nonfinite input, QR failure, unstable
  adjacent Schur swaps, a selection rule changed under reordering roundoff,
  or an unacceptable reconstruction residual.

On failure, the outputs are not a valid Schur decomposition. For near-boundary
or nearly defective matrices, the selected invariant subspace can be
ill-conditioned; a failure status is preferable to a false success. Values
extremely close to the sorting boundary may have different membership from a
LAPACK-backed SciPy calculation because of floating-point rounding.

## Numerical implementation

The input is normalized by its largest absolute real/imaginary component
before the complex QR-Schur process, and the final `T` is rescaled. Sorted
Schur forms use successive two-dimensional **unitary similarities** to move
chosen eigenvalues to the upper-left, preserving `A=Z*T*Z^H` and upper
triangularity. Reordering is verified by checking region membership and a
normalized reconstruction residual. A zero matrix is handled directly, with
correct `sdim` for predicates that select the zero eigenvalue.

Unlike LAPACK's blocked Schur QR and reordering routines, this implementation
is intended for small matrices; large and highly nonnormal cases are still
challenging. Computation is roughly cubic (plus Schur swaps), with quadratic
workspace. No new external numerical dependency is required.

## Validation

```bash
python tools/linalg_reference/run_schur_checks.py \
  --mode both --cases 480 --max-order 40 --seed 20261108
python tools/linalg_reference/run_schur_checks.py \
  --mode both --cases 384 --max-order 40 --seed 37271
```

The Fortran tests use GNU Fortran with `-fcheck=all` and floating-point traps in
checked mode and `-O3` in optimized mode. The independent Python reference
uses `scipy.linalg.schur(a, output='complex')` and its sorting rules, checking
reconstruction, unitarity, triangularity, spectra, and selected dimensions.
Comparisons do not require equality of Schur vectors or triangular factors:
they are not unique. Additional native tests cover a user-defined predicate,
zero and empty matrices, invalid sizes, and nonfinite data. Full FPM tests were
not run in this environment.

Algorithm reference: Golub and Van Loan, *Matrix Computations*, 4th ed.
(2013), chapters 7–8. Public SciPy API reference:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.schur.html

## Real-output validation

```bash
python tools/linalg_reference/run_real_schur_checks.py --mode both \
  --cases 480 --max-order 40 --seed 39711
python tools/linalg_reference/run_real_schur_checks.py --mode both \
  --cases 540 --max-order 40 --seed 52971
python tools/linalg_reference/run_real_schur_checks.py --mode both \
  --cases 360 --max-order 40 --seed 73019
```

The SciPy comparison checks real Schur reconstruction, orthogonality,
non-overlapping 2x2 blocks, eigenvalue multisets, and selected invariant
subspaces. Eigenvectors and exact 2x2 entries are not uniquely defined and
are not compared elementwise. See
`docs/linalg_real_schur_validation_2026-11-08.md`.

The real form is obtained through the existing complex QR-Schur engine rather
than native real Francis QR; it may fail for severely nonnormal or tightly
clustered spectra and is slower than an optimized LAPACK-backed approach.

## Convert real Schur factors without rerunning Schur QR (`rsf2csf`)

Use this routine when the real quasi-upper-triangular factors `(Tr,Zr)`
have already been computed and a **complex upper-triangular** form is wanted.

```fortran
use scifort_linalg_schur, only : rsf2csf

real(dp) :: tr(n,n), zr(n,n)
complex(dp) :: tc(n,n), zc(n,n)
integer :: status
call rsf2csf(tr, zr, tc, zc, status)
```

The factors satisfy `Zr*Tr*transpose(Zr) = A` and
`Zc*Tc*transpose(conjg(Zc)) = A`. `Tc` is upper triangular and `Zc` is
unitary. Complex-conjugate eigenvalues in 2x2 blocks are diagonalized
using local complex unitary similarities; the eigenvalue order and phases
are not uniquely prescribed. `Tr` and `Zr` are unchanged.

All arguments are required; the real inputs and complex outputs must be
square and of the same shape, including order zero. The numerical status
codes are `linalg_status_success` (0), `linalg_status_invalid_shape` (1)
for malformed dimensions, and `linalg_status_no_convergence` (4) for
nonfinite inputs, invalid real Schur structure, or numerical difficulties.
The outputs must not be used after a nonzero status.

Unlike the public `schur` implementation, this conversion has **no
order-40 algorithmic limit**: it operates on any already-computed real
Schur factors that fit in memory. It does not compute a real Schur
factorization itself and expects disjoint 1x1 / 2x2 diagonal blocks.
No additional numerical dependency is required.

See [native and SciPy-reference validation](linalg_rsf2csf_validation_2026-11-08.md)
and [the runnable example](../example/linalg_rsf2csf_example.f90).
