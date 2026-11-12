# Dense symmetric and Hermitian eigenproblems (`scipy.linalg`)

The new MIT-licensed `scifort_linalg_dense_eigen` module provides real and
complex overloads of `eigh` and `eigvalsh` for **standard** symmetric and
Hermitian dense eigenproblems, now also accepting a positive-definite `b`
for generalized problems. See [Generalized eigenproblems](linalg_dense_eigen_generalized.md).

```fortran
use scifort_linalg_dense_eigen, only : eigh, eigvalsh

call eigh(a, w, status, z, lower, subset_by_index, max_iterations)
call eigvalsh(a, w, status, lower, subset_by_index, max_iterations)
```

Inputs and outputs are caller-allocated:

- `a(n,n)`: finite real symmetric or complex Hermitian matrix, double
  precision. **Only the selected input triangle is read.**
- `w(m)`: ascending real eigenvalues, where `m=n` by default, including
  multiplicity.
- `z(n,m)`: optional real orthogonal or complex unitary eigenvectors,
  with `A z(:,i) = w(i) z(:,i)`.
- `lower=.true.` (default): use the input lower triangle; if false, use
  the upper triangle. The unused triangle can contain undefined data.
- `subset_by_index=[first,last]`: optional *one-based, inclusive* interval
  of eigenvalues in ascending order. For example, `[2,4]` selects the
  second through fourth eigenvalues; `w` then has three entries and `z`
  has three columns. **SciPy/Python uses zero-based indices; SciFort uses
  Fortran's one-based convention.**
- `subset_by_value=[low,high]`: optional *half-open* `(low,high]` interval,
  with an additional `num_selected` count output and capacity-buffer semantics.
  See [value selection](linalg_dense_eigen_value_selection.md).
- `max_iterations=120`: maximum implicit-shift QL iterations per eigenvalue.
- `status`: `linalg_status_success=0`,
  `linalg_status_invalid_shape=1` for malformed shapes/indices/iteration
  limit, or `linalg_status_no_convergence=4` for nonfinite used data,
  numerical failure, or nonrepresentable eigenvalues.

The complex Hermitian diagonal is interpreted as real; imaginary diagonal
parts are ignored. The ignored triangular region is *not* checked for
finite entries, unlike SciPy's default `check_finite=True` when applied to
an entire array. Arrays of size zero are supported.

Usage example:

```fortran
real(dp) :: a(3,3), w(2), z(3,2)
integer :: status
! a has been populated with a symmetric matrix.
call eigh(a,w,status,z,subset_by_index=[2,3])
if (status /= linalg_status_success) error stop 'eigh failed'
```

See `example/linalg_dense_eigen_example.f90` for an executable example.

## Algorithm, complexity and limitations

`eigh` constructs the symmetric/Hermitian matrix from the selected
triangle, then uses SciFort's two-sided Householder Hessenberg reduction.
For Hermitian matrices the resulting Hessenberg form is tridiagonal to
roundoff. Its complex subdiagonal is made real and nonnegative by a diagonal
unitary transformation. The resulting real symmetric tridiagonal eigenproblem
is solved by SciFort's scaled implicit-shift QL iteration, accumulating
orthogonal transformations when eigenvectors are requested. Eigenvalues
are returned in ascending order. This is the standard Householder/QL method
in Golub and Van Loan, *Matrix Computations*, 4th edition, chapters 7–8.

- Runtime is O(n^3); storage is O(n^2). This implementation is unblocked,
  dependency-free and slower than SciPy's LAPACK implementations for large n.
- `subset_by_index` computes the full spectrum first and only then extracts
  the chosen range; it **does not** accelerate partial-spectrum problems.
- Generalized positive-definite `b` problems are supported for types 1, 2,
  and 3 using a Cholesky reduction. See the linked generalized API notes.
  Value-selected eigenpairs via `subset_by_value` are now supported (and
  compute the full spectrum). See [value selection](linalg_dense_eigen_value_selection.md).
  Batched inputs and LAPACK `driver` selections remain unimplemented.
- For repeated eigenvalues, eigenvectors in the same eigenspace are not
  unique; comparing individual eigenvectors against SciPy is inappropriate.
- Strongly nonuniform or nearly multiple eigenvalues can lose accuracy;
  the QL method can return a numerical-failure status rather than converge.

The implementation reuses existing SciFort MIT-licensed routines and does
not import or translate SciPy/LAPACK implementation code.

## Testing

```sh
python tools/linalg_reference/run_dense_eigen_checks.py \
  --mode both --cases 500 --max-order 96 --seed 20261108
```

The independent SciPy comparison harness checks eigenvalue differences,
normalized eigen-equation residuals and orthogonality. It exercises real
and complex inputs, selected indices, both triangles, unused-triangle NaNs,
zero/repeated/clustered eigenvalues, and magnitudes from `1e-150` to
`1e150`. See `docs/linalg_dense_eigen_validation_2026-11-07.md`.
