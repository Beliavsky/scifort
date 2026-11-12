# Ordered generalized complex Schur decomposition (`ordqz`)

`scifort_linalg_ordqz` implements the **complex-output** form of
[`scipy.linalg.ordqz`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ordqz.html),
with an idiomatic native Fortran calling convention. For real, complex, or
mixed-type square matrices `A` and `B`, it constructs

\[
 A=Q\,AA\,Z^H,\qquad B=Q\,BB\,Z^H,
\]

where `Q` and `Z` are unitary and `AA` and `BB` are upper triangular. The
returned homogeneous generalized eigenvalues are `(alpha(i),beta(i)) =
(AA(i,i),BB(i,i))`. Selected eigenvalues are moved to the top-left; `sdim`
reports how many were selected. Infinite eigenvalues are represented by a
zero or numerically negligible `beta` with nonzero `alpha`.

```fortran
use scifort_kinds, only : dp
use scifort_linalg_ordqz, only : ordqz

real(dp) :: a(n,n), b(n,n)
complex(dp) :: aa(n,n), bb(n,n), q(n,n), z(n,n), alpha(n), beta(n)
integer :: status, nselected

call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sort='lhp',sdim=nselected)
```

`ordqz` supports the existing `qz` methods, passed as an optional final
`method='auto'|'iterative'|'inverse'`. Its factorization order is 0..40,
inherited from `scifort_linalg_qz`. All four real/complex input overloads
return complex Schur factors, including for real inputs. Real
quasi-triangular generalized Schur output is **not** implemented.

## Sorting

The optional `sort` argument defaults to `'lhp'`. The accepted values are:

| Selector | Condition on finite `lambda=alpha/beta` |
|---|---|
| `'lhp'` | `real(lambda) < 0` |
| `'rhp'` | `real(lambda) > 0` |
| `'iuc'` | `abs(lambda) < 1` |
| `'ouc'` | `abs(lambda) > 1` |
| `'none'` | No reordering (SciFort extension) |

These agree with SciPy's **ordqz** selectors, not necessarily the boundary
conventions of ordinary `scipy.linalg.schur` sorting. An infinite eigenvalue
is in neither half-plane, is outside the unit circle, and is not inside it.
The irregular pair `(alpha,beta)=(0,0)` is never selected by a built-in rule.

A user-defined selector takes the **homogeneous** complex pair `(alpha,beta)`:

```fortran
call ordqz(a,b,aa,bb,q,z,status,sdim=nselected,selector=select_custom)
```

`select_custom` must have the public `ordqz_selector` interface. If supplied,
the callback overrides the built-in `sort` rule. The caller should evaluate
homogeneous coordinates without dividing by a potentially zero `beta`.

## Reordering precomputed QZ factors

`ordqz_factors` is a separate in-place operation on complex Schur factors:

```fortran
call ordqz_factors(aa,bb,q,z,status,sort='ouc',sdim=nselected)
```

It accepts **already triangular** factors from SciFort or another solver and
does not rerun QZ. It is not subject to the `ordqz` order-40 cap, but it uses
quadratic storage and may be expensive on large problems. The matrices are
left unchanged on a failed reorder. This can be useful when investigating a
precomputed pencil or testing QZ reordering independently of the QZ engine.

The reordering algorithm uses local two-sided unitary rotations to exchange
adjacent complex 1x1 generalized Schur blocks. The first column of the right
rotation is a generalized eigenvector of the lower block; the left rotation
aligns both resulting first columns. This preserves both matrix identities
without forming an inverse. The implementation normalizes `AA` and `BB`
independently during reordering and verifies reconstruction before returning.

## Failure modes and accuracy

- `linalg_status_success` (`0`): completed and verified numerically.
- `linalg_status_invalid_shape` (`1`): incompatible array shapes, unknown
  selector code, or unsupported factorization size.
- `linalg_status_no_convergence` (`4`): nonfinite inputs, unresolved adjacent
  swap, or a reconstruction check that exceeds the allowed tolerance.
- Other statuses from the underlying `qz` routine are propagated.

QZ can return tiny roundoff-sized `beta` values for theoretically infinite
eigenvalues. The built-in selectors therefore classify numerically zero
`beta` relative to the full second Schur factor. This is an intentional
robustness choice; eigenvalues extremely close to a selection boundary may
be classified differently from SciPy depending on rounding.

Reordering can become ill-conditioned when eigenvalues are tightly clustered
or nearly defective. A successful small reconstruction residual and small
triangularity error do **not** guarantee small forward eigenvalue errors.
Batch operations, real generalized Schur output, LAPACK-grade reordering
condition estimates, and robust sorting of nearly singular irregular pencils
remain unimplemented. No new external numerical dependencies were added.

See [native tests](../test/test_linalg_ordqz.f90),
[Python reference comparisons](../tools/linalg_reference/run_ordqz_checks.py),
[Fortran example](../example/linalg_ordqz_example.f90), and
[validation notes](linalg_ordqz_validation_2026-11-09.md).
