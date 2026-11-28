# Logrank and shared interval contracts

Sources: [empirical_distribution.f90](../src/empirical_distribution.f90),
[logrank_bridge.f90](../src/c_api/logrank_bridge.f90), and the
[generated bridge contracts](logrank_bridge_shapes.md).

| Procedure | Arrays and relationships | Empty cases |
| --- | --- | --- |
| `logrank_arrays` | finite `x(nx)` and `y(ny)`; independent lengths | either or both may be empty |
| `logrank_censored_array` | finite uncensored `y(ny)`; independent event/censor lengths in X | any length may be zero |
| `logrank_array_censored` | finite uncensored `x(nx)`; independent event/censor lengths in Y | any length may be zero |
| `logrank_core` / `survival_logrank` | finite events `x(nx)`, `y(ny)` and right-censor times `x_right(rx)`, `y_right(ry)`; four independent lengths | all lengths >= 0; empty groups and all-censored samples permitted |
| `unpack_censored_data` | separately allocated copies of event and censor components | an absent component becomes length zero |
| `count_equal_value`, `count_at_risk` | unrestricted vector lengths; no equality constraint between event and censor vectors | empty vectors contribute zero |

The bridge requires the total observation count to fit the native default integer
kind; overflow returns status 1 and NaN scalars. Negative C extents return 1 before
native entry and leave scalar outputs unchanged. All pointers, including empty
input buffers, must be valid; nonempty array buffers are disjoint. Python copies
overlapping read-only inputs, so comparisons of a sample with itself are allowed.

Events tied with censoring times remain in the risk set at that event time.
The statistic is signed standard-normal: positive means more observed deaths in
X than expected under the pooled risk sets. Alternatives 0/1/2 mean two-sided,
less and greater; their probabilities are respectively `2*SF(abs(z))`, `CDF(z)`
and `SF(z)`, capped at one where appropriate. The native statistic and formulas
are unchanged. Zero variance, including empty/all-censored comparisons, produces
NaN statistic and p-value with successful completion status 0. Invalid alternatives
or nonfinite input times produce status 2 and NaN scalars before native calls that
can terminate. Left and interval censoring are not accepted by this interface.

The shared `confidence_interval` value type now lives in
[confidence_interval.f90](../src/confidence_interval.f90). Its `low` and `high`
fields, order, precision and initialization behavior are unchanged. Existing
`use scifort_contingency_meta, only : confidence_interval` and the
`scifort_stats` import path re-export the same type; this is not a duplicate type.
Native consumers should rebuild their Fortran modules after the move. Existing
C symbols are unchanged, and compiler module files are not part of the C ABI.
The empirical module imports the lightweight type directly, avoiding the previous
contingency/hypothesis/descriptive dependency chain in the binding build.
