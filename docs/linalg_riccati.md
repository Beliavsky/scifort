# Dense algebraic Riccati equations (CARE and DARE)

The module `scifort_linalg_riccati` implements two `scipy.linalg` counterparts
with **real and complex binary64** generic interfaces:

- `solve_continuous_are(a,b,q,r,x,status,...)`
- `solve_discrete_are(a,b,q,r,x,status,...)`

Their generalized equations are, respectively,

\[
E^HXA + A^HXE - (E^HXB+S)R^{-1}(B^HXE+S^H)+Q=0,
\]

\[
A^HXA-E^HXE-(A^HXB+S)(R+B^HXB)^{-1}(B^HXA+S^H)+Q=0.
\]

Here \(A,Q,E,X\) are \(n\times n\), \(B,S\) are \(n\times m\),
\(R\) is \(m\times m\), and \(H\) denotes conjugate transpose.
If omitted, `e` is the identity matrix and `s` is zero.

For the original standard equations, only `a,b,q,r,x,status` are required:

```fortran
use scifort_kinds, only : dp
use scifort_linalg_riccati, only : solve_continuous_are
real(dp) :: a(2,2), b(2,1), q(2,2), r(1,1), x(2,2), res
integer :: status, iters
! Initialize the input matrices.
call solve_continuous_are(a,b,q,r,x,status,iterations=iters,residual=res)
```

Both generic interfaces accept these optional named arguments **after** `status`:

| Argument | Description |
|---|---|
| `e` | Nonsingular descriptor matrix `(n,n)` |
| `s` | State-control cross cost `(n,m)` |
| `max_iter` | Matrix-sign or doubling iteration limit (default 100) |
| `tol` | Positive relative convergence tolerance (default `1.0e-11_dp`) |
| `iterations` | Number of principal algorithm iterations used (output) |
| `residual` | Normalized residual of the *original* generalized equation (output) |

The result `x` must have shape `(n,n)`; all arrays are caller allocated.
The main inputs have `intent(in)` and are never altered. The result is
Hermitian for complex inputs and symmetric for real inputs. Invalid or failed
solves leave `x` zero. Integer statuses come from `scifort_linalg`:

- `linalg_status_success` (0): converged with a small true residual and a
  stabilizing-branch selection under the stated control assumptions.
- `linalg_status_invalid_shape` (1): input shape, non-Hermitian input cost,
  iteration/tolerance options, or dimension exceeds limit.
- `linalg_status_not_positive_definite` (2): `R` is not numerically Hermitian
  positive definite.
- `linalg_status_no_convergence` (4): nonfinite data, iteration limit,
  unstable/undetectable case, or residual did not meet the acceptance threshold.
- `linalg_status_singular` (5): descriptor `E` was singular to pivot tolerance.

## Algorithms, assumptions, and limitations

**Continuous ARE** is solved with scaled Newton iteration for the sign of the
\(2n\times2n\) Hamiltonian matrix, extracting its stable invariant subspace.
A fallback recovers the graph when the lower-right sign block is singular.
When needed, Newton refinements solve linear Sylvester equations.

**Discrete ARE** uses the structure-preserving doubling iteration on
\(A\), \(G=BR^{-1}B^H\), and \(H=Q\); when needed, discrete Lyapunov
corrections refine its solution. Both forms use the ordinary cross-cost and
descriptor transformations before their main iteration. The returned result
is checked in the **original generalized equation**, not just the transformed
one. The Hamiltonian sign isolates the stable continuous-time graph; discrete
doubling selects the stabilizing solution under the usual detectability
assumptions. A numerical repeated-squaring screen additionally rejects
nonstabilizing zero-state-cost discrete solutions. We deliberately do not
use a general power-norm screen for nonnormal closed-loop matrices because
it can reject genuinely stable systems.

The routines require stabilizability/detectability and well-separated stable
subspaces for dependable convergence. They assume Hermitian `Q` (ideally
positive semidefinite) and Hermitian positive-definite `R`; generalized
cross-cost inputs should preserve the expected nonnegative cost structure.
No positivity verification of `Q` or `Q-S R^{-1} S^H` is currently performed.

These algorithms are **not** SciPy's balanced Schur/QZ decomposition algorithms
and do not handle every indefinite, marginally stable, ill-conditioned, or
nearly singular problem that SciPy handles. The solution is limited to
`riccati_max_states=20` and `matrix_equation_max_unknowns=400`, primarily
because the optional Sylvester/Lyapunov refinements use Kronecker-based
linear solves. The matrix-sign approach can lose stability near eigenvalues
on the imaginary axis; doubling can fail near the discrete unit circle.
The algorithms assume stabilizability and detectability and do not provide
formal certificates of closed-loop stability.

No BLAS/LAPACK dependency is introduced. Implementation and test code are
original MIT-licensed SciFort work, not translated SciPy source.

## Validation

Tests are in `test/test_linalg_riccati.f90`. A runnable example is
`example/linalg_riccati_example.f90`. The independently compiled Fortran
reference driver and the randomized SciPy/NumPy harness are in
`tools/linalg_reference/`. To reproduce:

```shell
python tools/linalg_reference/run_native_linalg.py --mode both
python tools/linalg_reference/compare_riccati_scipy.py --mode both --cases 320
```

The SciPy 1.17 randomized harness uses `balanced=False` because in some
well-conditioned generalized cases SciPy's default balancing has returned
solutions with **larger true equation residuals** than unbalanced QZ.
This is a reference-accuracy safeguard, not an argument that the Fortran
method is superior. The script verifies the SciPy result's true equation
residual before comparing the two solutions, and reports any rejected reference
cases. See `validation/linalg_riccati_2026-11-06.md`.
