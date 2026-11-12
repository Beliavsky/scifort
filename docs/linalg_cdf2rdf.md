# Complex eigendecomposition to real block form (`scipy.linalg.cdf2rdf`)

The `scifort_linalg_cdf2rdf` module converts the complex eigenvalues and
right eigenvectors of a **real** nonsymmetric matrix to real factors
`wr` and `vr`. It does not itself compute eigenvectors. For a real matrix
`A` with right eigenvectors `v` and eigenvalues `w`, the identity is

```text
A * vr = vr * wr
```

This translation supports real or complex double-precision eigenvalues,
complex eigenvectors (or real eigenvectors for real eigenvalues), and a
single leading batch dimension. The output arrays are always `real(dp)`.

## Fortran interface

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_general_eigen, only : eig
use scifort_linalg_cdf2rdf, only : cdf2rdf

real(dp) :: a(3,3), wr(3,3), vr(3,3)
complex(dp) :: w(3), v(3,3)
integer :: status

call eig(a, w, status, vr=v)
if (status == linalg_status_success) then
    call cdf2rdf(w, v, wr, vr, status)
end if
```

All `cdf2rdf` interfaces use the same five positional arguments, with
caller-allocated output and an explicit integer status. The generic accepts:

- `complex(dp) w(n), complex(dp) v(n,n)`;
- `real(dp) w(n), real(dp) v(n,n)`;
- `real(dp) w(n), complex(dp) v(n,n)`;
- the corresponding one-batch-dimension variants `w(batch,n)`,
  `v(batch,n,n)`, `wr(batch,n,n)`, and `vr(batch,n,n)`.

A zero-order single matrix is valid. A zero-length *batch* is rejected.
Unlike SciPy, arbitrary additional batch dimensions are not supported.

## Pair ordering, phases, and blocks

For each conjugate eigenvalue pair at positions `j` and `k`, ordered as
`w(j) = alpha + i*beta`, `w(k) = alpha - i*beta`, the real block has

```text
wr(j,j) = alpha     wr(j,k) = beta
wr(k,j) = -beta     wr(k,k) = alpha
```

The real basis at those positions is `vr(:,j)=-aimag(v(:,j))` and
`vr(:,k)=real(v(:,j))`. The real/imaginary parts of the *first*
eigenvector suffice, so a different complex phase of `v(:,k)` cannot
break the reconstruction. All other positions in `wr` are zero except
the diagonals for real eigenvalues. Any real eigenvalues between the two
members of a complex pair remain in their original positions.

This uses the elementary identity obtained by separating real and
imaginary parts of `A*v = w*v`; no dense eigensolver is run again.
Nonreal eigenvalues must appear in conjugate pairs **in occurrence
order**. For example, a real eigenvalue may separate two members, but
two different complex pairs must not be interleaved. Unlike SciPy, which
may silently return an incorrect result for interleaved pairs, SciFort
checks pair values and returns an error for mismatches.

Fortran `eig` may occasionally give a real eigenvalue a tiny imaginary
roundoff component. SciFort first tries the exact complex-pair layout.
Only if that fails, imaginary parts of size no more than 64 machine
precisions relative to the input eigenvalue scale are treated as
roundoff. Valid tiny complex pairs are *not* discarded when their exact
pairing is valid. Eigenvectors associated with effectively real
eigenvalues are converted by taking their real part.

## Status and limitations

- `linalg_status_success`: conversion completed.
- `linalg_status_invalid_shape`: output/input shape mismatch, zero batch,
  unmatched conjugate pair, or interleaved/nonconjugate pair.
- `linalg_status_no_convergence`: NaN or infinity in either input.

On error all caller-allocated outputs are filled with zero. For multiple
batch members, an error in any member zeros the entire batch output.

The caller must supply a valid eigendecomposition of a **real** matrix;
complex-valued input matrices cannot in general be reduced to a real
eigenbasis this way. Repeated and defective eigenvalues may not have a
complete eigenvector basis; the method does not repair ill-conditioning.
The routine is not a substitute for a real Schur decomposition.

See `test/test_linalg_cdf2rdf.f90`,
`test/test_linalg_cdf2rdf_eig.f90`,
`reference/check_linalg_cdf2rdf.py`, and
`example/linalg_cdf2rdf_example.f90` for checks and usage.

Public SciPy reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cdf2rdf.html
