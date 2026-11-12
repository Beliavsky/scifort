# Generalized nonsymmetric eigenproblems (`eig`, `eigvals`)

The MIT-licensed `scifort_linalg_generalized_eigen` module solves the regular
nonsymmetric matrix pencil

\[
A v_j = \lambda_j B v_j, \qquad A^H u_j = \overline{\lambda_j}B^H u_j,
\]

for dense real, complex, and mixed real/complex matrices. It computes
complex left/right eigenvectors and represents each eigenvalue homogeneously
as `(alpha(j), beta(j))`, with `lambda = alpha/beta`. This is particularly
important when `B` is singular or when the finite ratio overflows. The
underlying SciFort QZ solver supports order `0..40` only.

## API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_generalized_eigen, only : eig, eigvals

real(dp) :: a(n,n), b(n,n)
complex(dp) :: w(n), vr(n,n), vl(n,n), alpha(n), beta(n)
integer :: status

call eig(a,b,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
call eig(a,b,w,status,vr=vr,method='iterative')
call eigvals(a,b,w,status,alpha=alpha,beta=beta)
```

The generic interfaces accept all four combinations of `real(dp)` and
`complex(dp)` for A and B. **The B argument is required**: ordinary
`eig(a,w,status,...)` and `eigvals(a,w,status)` remain available from the
unchanged `scifort_linalg_general_eigen` module. Import the two modules with
renamed generic names if both are needed in one scope.

The optional `vr` and `vl` are `complex(dp)` `(n,n)` arrays, each column a
Euclidean-unit eigenvector; phases and bases in repeated eigenspaces are
unspecified. The optional `(alpha,beta)` arrays have length n. An infinite
generalized eigenvalue is returned as `w(j)=(+Inf,0)` when `beta(j)=0` and
`alpha(j) ≠ 0`. Even when finite `w` overflows or underflows, the homogeneous
values retain the original QZ factor diagonals and should be preferred.
For left eigenvectors the defining finite/infinite identity is

\[
\overline{\beta_j}A^H u_j=\overline{\alpha_j}B^H u_j,
\]

and the analogous right identity is `beta*A*v = alpha*B*v`.

`method='auto'` (default) uses the two-sided QZ iteration with inverse-based
fallback, and `method='iterative'`/`'inverse'` select the underlying method.
Already jointly upper- or lower-triangular matrix pairs use an exact diagonal
fast path, protecting ill-conditioned triangular spectra from Schur roundoff.

## Numerical method and limitations

After QZ produces `A=Q*AA*Z^H` and `B=Q*BB*Z^H`, each upper-triangular
homogeneous pencil `beta*AA-alpha*BB` is solved by protected back or forward
substitution. A right eigenvector is `Z*y`; a left eigenvector is `Q*y`.
Independent coefficient scaling prevents overflow when the two matrices have
very different magnitudes. The eigenvectors are normalized and the triangular
pencil residual is checked, with a failure status rather than an unchecked
answer for unresolvable defective cases.

- `linalg_status_success=0`: computed requested eigenpairs.
- `linalg_status_invalid_shape=1`: incompatible arrays, order over 40, or
  unknown QZ method.
- `linalg_status_no_convergence=4`: nonfinite input, failed QZ, or poorly
  resolved eigenvectors.
- `linalg_status_singular=5`: zero/zero homogeneous eigenvalue, a singular
  or unsupported irregular pencil.

This is **not** LAPACK's generalized `xGGEV` eigensolver. Strongly nonnormal
pencils may have large forward eigenvalue errors despite small backward
residuals; defective and repeated eigenvalues may not yield a full independent
basis. No batching, generalized balancing, eigenvalue condition estimates,
or `ordqz` sorting is included. When status is nonzero, output arrays must not
be interpreted as valid eigenpairs.

## Validation

```bash
python tools/linalg_reference/run_generalized_eigen_checks.py \
  --mode both --cases 400 --max-order 40 --seed 91283
```

The native tests and reference harness include finite and infinite roots,
triangular and dense pencils, mixed input kinds, both-singular but regular
pairs, extreme scale ratios, repeated eigenvalues, and invalid inputs.
Eigenvalues are compared to SciPy's homogeneous `eig` using a projective
(chordal) distance and optimal matching; eigenvectors are checked by both
left and right residuals, not by entrywise agreement. Compile modes include
`-fcheck=all` with floating-point traps, and `-O3`. See
[`docs/linalg_generalized_eigen_validation_2026-11-08.md`](linalg_generalized_eigen_validation_2026-11-08.md)
for completed test results.

Mathematical basis: generalized Schur decomposition and triangular
homogeneous eigenvector substitution; Golub and Van Loan, *Matrix
Computations*, 4th ed. (2013), chapters 7–8. SciPy supplies only the public
interface conventions and independent numerical reference:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html
