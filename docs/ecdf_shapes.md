# Empirical distribution array contracts

Source: [empirical_distribution.f90](../src/empirical_distribution.f90).
Bridge: [ecdf_bridge.f90](../src/c_api/ecdf_bridge.f90), with
[generated contracts](ecdf_bridge_shapes.md).

| Procedure or storage | Array relationship | Empty behavior |
| --- | --- | --- |
| `ecdf_uncensored` | `sample(n)`, `n = size(sample)`; finite values | `n = 0` allowed |
| `ecdf_right_censored` | events `(n)` and right-censoring times `(r)`; independent lengths, finite values | either or both may be empty |
| EDF storage | quantiles and probabilities `(m)`, sorted distinct knots; `0 <= m <= n+r` | `m = 0` for no observations |
| EDF variance storage | at-risk and event counts `(m)` for original estimates; `(0)` for confidence bounds | original empty estimates have empty variance arrays |
| `linear_greenwood`, `loglog_greenwood` | `low(m)`, `high(m)`; `m = size(self%probabilities)`; risk/event arrays also `(m)` | empty bounds permitted |
| `set_edf` | matching quantile/probability lengths; risk/event lengths depend on variance flag as above | all empty permitted |
| `stable_sort_indices` | values `(n)` and explicit output `order(n)` | no sorting for `n <= 1` |
| `merge_sort_real` | values `(n)`, workspace capacity at least `n` | immediate return for `left >= right` |
| `merge_sort_indices` | values and order `(n)`, workspace capacity at least `n` | immediate return for `left >= right` |

For nonempty merge work, indices satisfy `1 <= left <= right <= n`;
index-sort order entries are a permutation of `1:n`. Unrestricted helper vectors
have no invented equality constraints. Public and private declaration comments
now record the constrained arrays used by this ECDF path. Numerical statements
and native public signatures are unchanged.

The bridge uses exact-size retrieval: query `ecdf_info`, then pass two distinct
caller-owned buffers of length `count` to `ecdf_values`. Incorrect nonnegative
lengths return status 1 and NaN-filled buffers. Negative C extents are rejected
before native entry, preserving outputs. Empty C buffers still require valid
addresses, as documented in the generated header. Array inputs/outputs must be
disjoint; the Python constructor copies overlapping read-only event/censor inputs.

Original estimates, their CDF/SF partners and all confidence bounds own separate
storage. A bound survives closing its source; closing one function does not close
its partner. Bounds do not contain variance data, so nested intervals return
status 2. All owned outputs need fresh slots and must be destroyed exactly once;
null destruction is idempotent. Raw non-null handles must be valid and live.
