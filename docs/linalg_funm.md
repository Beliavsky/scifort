# General analytic dense matrix functions: `funm`

`scifort_linalg_funm` implements a subset of the mathematical functionality of
`scipy.linalg.funm` for small real or complex square matrices. It is a
self-contained MIT-licensed Schur--Parlett implementation: no new BLAS/LAPACK,
SciPy, Python, or other runtime dependency is added.

## API and callbacks

```fortran
use scifort_kinds, only : dp
use scifort_linalg_funm, only : funm, funm_real

call funm(a, func, x, status, commutator_residual, separation_tol, &
    derivative, taylor_max_order)
call funm_real(real_a, func, real_x, status, commutator_residual, &
    separation_tol, derivative, taylor_max_order)
```

- The generic `funm` accepts a `real(dp)` or `complex(dp)` square input array
  `a`, and writes a **complex** matrix `x` of the same order.
- `func(z)` is a user-supplied function with one `complex(dp), intent(in)`
  argument, returning `complex(dp)`. The function must be analytic in a
  neighborhood of the spectrum on the selected branch.
- Optional `derivative(z, order)` returns the `order`-th complex derivative
  of `func`. **By default**, the order is a positive integer up to the
  size of a repeated eigenvalue block minus one, preserving the existing
  callback contract. It enables Taylor evaluation of Schur blocks
  with repeated (or sufficiently close) eigenvalues even when those
  eigenvalues were originally **noncontiguous** on the Schur diagonal.
  When a repeated spectral cluster is diagonalizable and its reordered
  Schur block has negligible off-diagonal entries, derivatives are **not**
  needed; only the scalar callback is used.
- Optional `taylor_max_order` (1--128) enables **extended, adaptive Taylor
  evaluation** inside clustered Schur blocks when its value exceeds the
  block size minus one. The user must supply accurate derivatives for **all
  orders requested**, potentially up to this limit. The extension checks
  for three consecutive negligible tail terms after order `block_size-1`
  and reports `linalg_status_no_convergence` if the limit is reached first.
  A limit below `block_size-1` is invalid for a coupled cluster.
  It is best suited to scalar functions such as `exp`, `sin`, and `log`
  whose arbitrary-order derivatives are available. The stopping criterion
  is heuristic, **not** a certified truncation-error estimate.
- `status` uses the existing `scifort_linalg` integer status constants:
  `linalg_status_success` (0), `linalg_status_invalid_shape` (1) and
  `linalg_status_no_convergence` (4). On failure `x` is zero.
- Optional `separation_tol` defaults to `1e-10` and must be positive and
  below `0.1`. It controls when the recurrence regards eigenvalues as too
  close to divide safely. This is a *screening* parameter, not a requested
  error tolerance.
- Optional `commutator_residual` returns the normalized Schur-coordinate
  one-norm `||TF-FT||_1/(||T||_1 ||F||_1)`. It is **not an estimate or bound
  of the numerical error**: many incorrect matrices can commute with `A`.
- `funm_real` takes a real matrix, returns a real matrix only if the complex
  result has negligible imaginary part, and otherwise reports failure.
- Matrices of orders 1--40 are supported. Array storage is caller-owned.
  Batched SciPy inputs and non-square matrices are not supported.

## Example

The file `example/linalg_funm_example.f90` illustrates a non-diagonalizable
Jordan block, using the fact that every derivative of `exp(z)` is `exp(z)`.
Its central call is:

```fortran
call funm(a, scalar_exp, x, status, derivative=exp_derivatives)
```

Callbacks are ordinary Fortran procedures, including internal procedures.
For example, `scalar_exp(z)` calls the intrinsic scalar `exp(z)`.

## Numerical method

1. Scale the input by its largest element magnitude and compute the complex
   Schur form `A/s = Z*T*Z^H` with SciFort's existing internal QR kernel.
2. Find *all* eigenvalue clusters on the Schur diagonal, **whether or not**
   derivatives are supplied, including pairs separated by other eigenvalues.
   Components are formed by the relative `separation_tol` proximity test, with transitive
   closure. Stably group them using adjacent **complex unitary** Schur swaps,
   updating both the triangular factor `T` and the Schur basis `Z`.
3. Set the diagonal of `F` to `func(s*T(i,i))`. For coupled clustered diagonal
   blocks, the default method uses a degree-`block_size-1` Taylor polynomial
   with a near-nilpotence screen. With `taylor_max_order` greater than that
   degree, it sums additional Taylor terms until three consecutive tail
   terms have negligible one-norm relative to the partial result, or
   reports nonconvergence at the limit. In both modes, the resulting
   diagonal is checked against direct scalar function values. Diagonalizable
   repeated clusters with negligible off-diagonal coupling need no derivatives.
4. For each pair of distinct spectral blocks, construct the **block Parlett**
   right-hand side from previously computed blocks and solve
   `Tii*Fij - Fij*Tjj = RHS` via triangular Sylvester back-substitution.
   Reject unresolved small spectral separations or nonfinite intermediates.
5. Return `Z*F*Z^H`, subject to finite-number and commutator checks.

The Schur swaps are unitary similarities; no eigenvector inversion or new
LAPACK dependency is required. The solver **does** use a block Schur--Parlett
recurrence, but clustered diagonal blocks still rely on Taylor series;
this is not a complete robust Schur--Padé method. No-derivative
calls still reject genuinely coupled Jordan blocks. The test for negligible
couplings is numerical, and does not guarantee forward accuracy near a
branch point or a pole.

The scalar scaling handles very small and very large matrices without using
such magnitudes directly in the Schur QR iteration. No third-party source
has been translated or vendored. The method is based on the standard
Schur--Parlett framework described in Nicholas J. Higham, *Functions of
Matrices: Theory and Computation*, SIAM (2008), and Golub & Van Loan,
*Matrix Computations*, 4th ed. (2013), section on matrix functions.

## Known numerical and compatibility limits

- This is **not** SciPy's complete `funm` implementation; it does not return
  SciPy's heuristic `errest` or accept vectorized scalar callbacks.
- Without derivatives, a non-diagonalizable matrix with repeated eigenvalues
  may have an undefined divided difference in this algorithm and will return
  `linalg_status_no_convergence`. **Diagonalizable** repeated-eigenvalue
  matrices are now supported without derivative callbacks, including
  noncontiguous repeated eigenvalues. With derivatives, Jordan blocks are
  supported by Schur reordering and Taylor evaluation, subject to screening.
  Very tight eigenvalue clusters and badly conditioned swaps may still be rejected, and no returned status can certify small forward error.
- A small commutator residual does **not** certify accuracy. Severely
  nonnormal matrices, close eigenvalues, and near poles/branch cuts can be
  difficult. Even a successful return can have significant forward error
  when eigenvalue conditioning is poor. For example, dense unitary
  similarities of large defective Jordan structures can amplify errors
  in the upstream small-matrix Schur decomposition. The implementation
  detects and rejects some, but not all, numerically unresolved cases.
- User callbacks must select a consistent analytic branch. NaN and infinity
  in the input or scalar results are rejected. In the adaptive mode the
  Taylor series at the cluster center must converge throughout the cluster;
  near a pole, branch point, or cut it may fail even when the matrix function
  is mathematically defined. Derivative callbacks can themselves overflow
  for large orders or extreme input scales; such cases report failure.
- For the usual specialized functions (`expm`, `logm`, `sqrtm`, `sinm`, etc.),
  use the dedicated SciFort algorithms when numerical robustness is important.
- For large or highly nonnormal matrices, a LAPACK-backed Schur decomposition
  and more advanced diagonal-block evaluation would be more reliable.

## Test commands

```bash
python tools/linalg_reference/run_funm_checks.py --mode both --cases 550 --max-order 9 --reorder-cases 350 --diagonalizable-cases 350 --taylor-cases 400
```

The stand-alone checked/optimized runner builds the native tests and example,
compares against `scipy.linalg.funm` for well-separated spectra, uses
specialized SciPy matrix functions for Jordan matrices, and independently
checks repeated spectra, nonnormal couplings, extreme scaling, and globally reordered repeated eigenvalues, including diagonalizable
noncontiguous repeated spectra **without derivatives**. See `docs/linalg_funm_validation_2026-11-07.md`
for the initial implementation and
`docs/linalg_funm_reordered_validation_2026-11-07.md` for this extension.

The extension to block Sylvester solves and derivative-free semisimple
clusters is validated separately in
[`linalg_funm_block_parlett_validation_2026-11-07.md`](linalg_funm_block_parlett_validation_2026-11-07.md).

The opt-in higher-order Taylor mode and clustered-spectrum SciPy comparisons
are documented in
[`linalg_funm_adaptive_taylor_validation_2026-11-07.md`](linalg_funm_adaptive_taylor_validation_2026-11-07.md).
