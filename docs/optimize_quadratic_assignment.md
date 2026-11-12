# Quadratic assignment (FAQ and 2-opt)

SciFort offers `quadratic_assignment` in
`scifort_optimize_quadratic_assignment` and the `scifort_optimize` umbrella.
This is an original, deterministic Fortran implementation of two approximate
algorithms also exposed by `scipy.optimize.quadratic_assignment`.

The objective for a permutation `p` is

`f(p) = sum(a(i,j) * b(p(i),p(j)), i=1..n, j=1..n)`.

The mapping `res%col_ind(i)` assigns row/node `i` of `a` to the corresponding
row/node of `b`. Unlike SciPy's zero-based Python indices, **all SciFort
permutation indices, including the input pairs, are one-based**.

## API

```fortran
use scifort_optimize, only : dp, qap_result, quadratic_assignment

type(qap_result) :: res
res = quadratic_assignment(a, b, method='faq', maximize=.false.)
res = quadratic_assignment(a, b, method='2opt', partial_guess=guess)
```

The function is `pure` and uses no global random state. The input matrices
must be equally sized, finite and square. An empty problem succeeds with
an empty permutation and objective zero. A single-node problem is trivial.

Optional arguments (all keyword-compatible):

| Argument | Meaning |
|---|---|
| `method` | `'faq'` (default) or `'2opt'`, case-insensitive |
| `maximize` | Minimize by default; maximize when true |
| `partial_match(m,2)` | Fixed `(a_row,b_row)` one-based pairs; may be unsorted |
| `partial_guess(g,2)` | Initial `(a_row,b_row)` pairs for 2-opt; not fixed |
| `p0(n-m,n-m)` | Initial doubly stochastic free block for FAQ (default barycenter) |
| `maxiter` | Max accepted FW steps (default 30) or max accepted swaps (default 10000) |
| `tol` | FAQ step-change norm tolerance (default `0.03_dp`) |

A `partial_guess` is allowed only with `method='2opt'`; a `p0` is allowed only
with `method='faq'`. `p0` must be a finite, approximately doubly stochastic
matrix indexed by the **unfixed rows and columns in ascending order**. 2-opt
completes unspecified guessed entries by assigning remaining free columns in
ascending order; it does not generate a random initial permutation.

`type(qap_result)` contains:

- `integer, allocatable :: col_ind(:)` -- permutation, or empty on error;
- `real(dp) :: fun` -- true objective at that permutation (NaN on error);
- `integer :: nit` -- accepted FW steps or accepted 2-opt swaps;
- `integer :: status` and `logical :: success` -- result code and success flag.

Status codes are `qap_success=0`, `qap_invalid_input=1`,
`qap_assignment_failure=2`, and `qap_numerical_failure=3`. An unsuccessful
result is not a certified optimum. `quadratic_assignment_score(a,b,p)` is a
separate pure objective evaluator that returns NaN for invalid shapes or
permutations.

## Methods

**FAQ** starts from a doubly stochastic matrix on the unfixed nodes and
uses Frank-Wolfe steps. Each step solves a linear assignment problem for the
relaxation's gradient via SciFort's existing `linear_sum_assignment` and
minimizes the exact quadratic along that step. It rounds the relaxed block
via a final linear assignment; fixed matches are always retained. `p0`
can alter the local solution. This implementation has no SciPy-style RNG or
`shuffle_input` option. A barycentric start is deterministic and may become
trapped at a stationary point; varying `p0` can help.

**2-opt** starts from a feasible permutation and repeatedly accepts the
best improving free-row swap. It evaluates swap *deltas* in O(n) work rather
than recomputing the O(n^2) objective for every pair; accepted swaps are
checked by evaluating the complete objective. It stops at a pairwise local
optimum unless `maxiter` is reached. Its default initialization is
index-ordered, not SciPy's randomized initialization. The result and its
`nit` therefore need not match SciPy.

For a stronger heuristic, use FAQ's output as a 2-opt initial guess:

```fortran
faq = quadratic_assignment(a, b, method='faq')
allocate(guess(size(faq%col_ind), 2))
do i = 1, size(faq%col_ind)
    guess(i, :) = [i, faq%col_ind(i)]
end do
refined = quadratic_assignment(a, b, method='2opt', partial_guess=guess)
```

This preserves or improves the FAQ *permutation objective*, but it does not
provide a bound on the global optimum. `maxiter=0` skips local iterations.

## Validation

The native Fortran test is `test/test_optimize_quadratic_assignment.f90`.
The independent SciPy check is
`validation/compare_optimize_quadratic_assignment_scipy.py`, which compiles
and runs a batch Fortran driver without needing FPM:

```text
python validation/compare_optimize_quadratic_assignment_scipy.py
python validation/compare_optimize_quadratic_assignment_scipy.py --optimize
```

The driver compares 144 FAQ and 144 2-opt solves against SciPy 1.17.0,
including symmetric/asymmetric, seeded, maximization and user-supplied `p0`
cases. It independently recomputes every returned objective and verifies
that 2-opt has no improving permitted pairwise swap. Comparisons with SciPy
2-opt assess quality rather than requiring identical permutations because
its initialization and swap order are randomized.

## Limitations and provenance

These are approximate algorithms for an NP-hard problem. No global-optimum
certificate is produced. FAQ rounds the relaxation and can be far from the
optimal discrete permutation; 2-opt can become stuck at a local minimum.
Numerical overflow in objective evaluation is treated as a solver failure.
There is no general `rng` option, randomized multistart, or solver-supplied
`shuffle_input` in this version. No new third-party code or dependency was
added. See `CODE_PROVENANCE.md`.

References:

- Vogelstein et al. (2015), *Fast Approximate Quadratic Programming for Graph
  Matching*, PLOS ONE 10(4), e0121002, DOI: 10.1371/journal.pone.0121002.
- SciPy documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.quadratic_assignment.html
- SciPy FAQ method documentation: https://docs.scipy.org/doc/scipy/reference/optimize.qap-faq.html
- SciPy 2-opt method documentation: https://docs.scipy.org/doc/scipy/reference/optimize.qap-2opt.html
