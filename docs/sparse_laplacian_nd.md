# Multidimensional sparse grid Laplacians and structure queries

`scifort_sparse_laplacian_nd` implements the real, negative-semidefinite
finite-difference operator corresponding to `scipy.sparse.linalg.LaplacianNd`
for **Dirichlet**, **Neumann**, and **periodic** boundaries. The operator
uses one-based Fortran vector indexing, but the **last grid dimension varies
fastest**, matching NumPy C-order flattening and SciPy's 2-D/3-D stencil.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix
use scifort_sparse_laplacian_nd, only : laplacian_nd_operator, laplacian_nd, &
    laplacian_nd_matvec, laplacian_nd_tocsr, &
    laplacian_nd_eigenvalues, laplacian_nd_eigenvectors
implicit none
type(laplacian_nd_operator) :: op
type(csr_matrix) :: a
real(dp), allocatable :: y(:), values(:), vectors(:,:)
integer :: info, i

op = laplacian_nd([3, 4], 'neumann', info)
call laplacian_nd_matvec(op, [(real(i,dp), i=1,12)], y, info)
a = laplacian_nd_tocsr(op, info)
call laplacian_nd_eigenvalues(op, values, m=3)
call laplacian_nd_eigenvectors(op, vectors, m=3)
```

The constructor is `laplacian_nd(grid_shape, boundary_conditions, status)`.
`boundary_conditions` defaults to `neumann`. Every dimension must be positive,
and its total grid size must fit a default Fortran integer. Invalid input
returns an uninitialized operator and a nonzero optional status. Use
`valid_laplacian_nd(op)` or `laplacian_nd_shape(op, shape, n, status)`
for inspection.

The module exposes:

- `laplacian_nd_matvec(op, x, y, status)` computes `L*x` in **O(d*N)**
  work with only the output vector allocated (`N=product(grid_shape)`).
- `laplacian_nd_matmat(op, x, y, status)` applies the operator to each
  column of a `(N, nrhs)` dense matrix.
- `laplacian_nd_tocsr(op, status)` generates a canonical, one-based CSR
  matrix using SciFort's existing COO constructor. A periodic dimension
  of length two contributes **two edges to the same neighboring vertex**.
- `laplacian_nd_eigenvalues(op, eig, status, m)` returns the `m` largest
  (closest to zero) eigenvalues, **sorted ascending**; without `m`, returns
  the entire spectrum.
- `laplacian_nd_eigenvectors(op, vec, status, m)` returns the corresponding
  orthonormal columns in the same eigenvalue ordering. The full basis is
  necessarily dense; requesting `m` reduces output size, but the current
  spectrum sorter still enumerates all `N` eigenvalues.

The 1D eigenvalues used in the tensor-product spectrum are

\[
\lambda_j^{D}=-4\sin^2\!\left(\frac{\pi(j+1)}{2(n+1)}\right),\quad
\lambda_j^{N}=-4\sin^2\!\left(\frac{\pi j}{2n}\right),\quad
\lambda_j^{P}=-4\sin^2\!\left(\frac{\pi\lfloor(j+1)/2\rfloor}{n}\right),
\]

where `j=0,...,n-1` and higher-dimensional eigenvalues are sums across
axes. The associated bases are discrete sine (Dirichlet), cosine (Neumann)
and real Fourier (periodic). Values are real binary64.

## SciPy 1.17 edge cases

SciPy 1.17's `LaplacianNd.tosparse()` and `matvec()` can disagree for a
**dimension of length one** under Neumann or periodic boundaries. SciFort
follows the mathematically consistent, matrix-free operator: such a dimension
contributes zero for Neumann/periodic boundaries and -2 for Dirichlet.
This agrees with SciPy's matrix-free action and analytic spectrum, not with
its affected sparse-conversion output.

SciPy 1.17's `eigenvalues(m)` optimization can also yield spurious modes on
some grid shapes because it selects a tuple-wide reduced shape using
lexicographic tuple ordering. SciFort returns the true largest `m`
eigenvalues of the full problem, matching `LaplacianNd.eigenvalues()[-m:]`.

## Sparse structural queries

`scifort_sparse_linalg` now exposes `spbandwidth(a, below, above, info)`
and `is_sptriangular(a, lower, upper, info)` for **CSR and CSC** matrices.
The lower bandwidth is `max(0, row-column)` over stored entries; the upper
bandwidth is `max(0, column-row)`. Triangularity is a structural property:
even an explicitly stored numeric zero counts as occupying its position.
Rectangular and empty sparse matrices are allowed; empty matrices have
bandwidth `(0,0)` and are both structurally lower and upper triangular.
An invalid sparse representation returns `sparse_linalg_invalid_input`.

## Testing

```sh
python tools/sparse_reference/compare_sparse_laplacian_nd_scipy.py --mode both
python tools/sparse_reference/run_native_sparse.py --mode both
```

See [recorded checked and optimized validation](sparse_laplacian_nd_validation_2026-11-05.md).

The first script compiles two reference drivers and checks grid matvecs,
canonical sparse construction, eigenvalues, eigenpair residuals, basis
orthogonality, and CSR/CSC structure against SciPy. The second executes
standalone native sparse regression tests. No new runtime dependencies are
required in SciFort (Python, NumPy and SciPy are test-only dependencies).
