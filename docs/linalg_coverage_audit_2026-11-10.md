# SciFort dense linear algebra coverage audit — 2026-11-10

## Scope and approach

The SciFort snapshot immediately preceding this update contained approximately
3,082 files, including 60+ `src/linalg*.f90` files. I compared the public names
exported by the installed SciPy 1.17.0 `scipy.linalg` module with identifiers in
these source files, then inspected the most relevant candidates manually.
This is an **API name screen**, not a formal numerical or signature compliance
audit: a textual match may mean a mention, an internal procedure, or a partial
implementation. Conversely, a different Fortran name may implement similar
functionality. Passing native tests for a routine is also not proof of
arbitrary-input equivalence to SciPy.

The installed SciPy reports 112 public `scipy.linalg.__all__` names (including
submodules, errors, aliases and dispatch helpers) and 98 callable names.
The pre-update source already included independent modules implementing many
principal dense algorithms:

| Area | Examples inspected |
|---|---|
| Decompositions | `qr`, `rq`, `lu`, `ldl`, `svd`, `schur`, `qz`, `ordqz`, `cossin` |
| Linear equations | `solve`, `solve_triangular`, `solve_banded`, `solve_toeplitz`, `solve_circulant`, Sylvester/Lyapunov/Riccati solvers |
| Eigenproblems | `eig`, `eigh`, generalized eigenproblem, banded/tridiagonal eigensolvers |
| Matrix functions | `sqrtm`, `logm`, `signm`, `fractional_matrix_power`, `funm`, `expm_frechet`, `expm_cond`, trigonometric matrix functions |
| Structure and low rank | Toeplitz/Hankel/circulant, QR updates, interpolative decomposition, randomized SVD, sketching |

### Highest-priority gap found: standalone `expm`

`src/linalg_expm_frechet.f90` already offered `expm_frechet`, computing both
`exp(A)` and its Fréchet derivative, but the source had no public `expm`
generic. `expm` is a foundational matrix function useful beyond derivative
calculations. This update adds `src/linalg_expm.f90` with real/complex and
single-batch interfaces, an independent Pade algorithm, native tests, SciPy
reference comparisons and documentation.

### Remaining differences and follow-on work

- The public `solve_lyapunov` name is a SciPy alias of the already implemented
  `solve_continuous_lyapunov`; no duplicate numerical solver is needed.
- `find_best_blas_type`, `get_blas_funcs`, and `get_lapack_funcs` are SciPy
  backend dispatch helpers, not missing mathematical algorithms. SciFort uses
  its own Fortran APIs and FPM build configuration.
- Module namespace identifiers such as `decomp`, `matfuncs`, `blas`, and
  `special_matrices` need not map to separate Fortran modules.
- Complete SciPy type support, arbitrary batch dimensions, error details,
  performance and randomized-algorithm parity remain substantial cross-cutting
  compatibility tasks. They have **not** been verified by this name screen.
- For future new-algorithm priorities, first check whether the routine is
  already present under an alternate name, as `orth`, `subspace_angles`, and
  the matrix-equation solvers were in this snapshot.

### Provenance and reproducibility

The audit used local source inspection and SciPy 1.17.0 API names; numerical
algorithms in the new code are independently implemented from standard
mathematical identities. Reference tests and source code are in the archive.
For the next substantive SciPy addition, rerun the comparison against the
then-current public API and inspect matching source modules to avoid duplication.
