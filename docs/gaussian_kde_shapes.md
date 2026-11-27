# Gaussian KDE native array contracts

Source: [gaussian_kde.f90](../src/stats/gaussian_kde.f90).
The native implementation and public Fortran signatures are unchanged. Required
shapes now appear next to dummy declarations and in procedure documentation.
These readable contracts apply to initialized estimators unless noted below.
The [generated bridge contracts](gaussian_kde_bridge_shapes.md) describe the
separate reviewed C adapter dimensions.

| Native procedure | Required arrays and dimensions | Empty cases and constraints |
| --- | --- | --- |
| `gaussian_kde_init_1d` | `dataset(n)`, optional `weights(n)`; `n = size(dataset)` | `n >= 2`; finite observations; positive weighted variance |
| `gaussian_kde_init_nd` | `dataset(d,n)`, optional `weights(n)`; dimensions from dataset | `1 <= d <= n`, `n >= 2`; finite observations and positive-definite weighted covariance |
| `gaussian_kde_evaluate_point`, `gaussian_kde_logpdf_point` | `point(d)`; `d = self%d >= 1` | finite evaluation coordinates |
| `gaussian_kde_evaluate_points`, `gaussian_kde_logpdf_points` | `points(d,m)`; `d = self%d >= 1`, `m = size(points,2)`; allocated result `(m)` | `m >= 0`; empty results allowed; finite coordinates |
| `gaussian_kde_integrate_gaussian` | `mean(d)`, `cov(d,d)`; `d = self%d >= 1` | finite arrays; sum of kernel covariance and `cov` must be positive definite |
| `gaussian_kde_integrate_box` | `low(d)`, `high(d)`; `d = self%d >= 1` | bound vectors must match the estimator dimension |
| `gaussian_kde_resample` | `samples(d,m)`; `d = self%d >= 1`, `m = size(samples,2)` | `m >= 0`; empty batches preserve RNG state |
| `gaussian_kde_marginal` | `dimensions(k)`; `k = size(dimensions)` | `1 <= k <= self%d`; distinct one-based indices in `[1,self%d]` |
| `gaussian_kde_inv_cov` | `inv_cov(d,d)`; `d = self%d >= 1` | caller supplies matching storage |
| private `evaluate_many` | `points(d,m)`, allocated result `(m)`; dimensions as above | empty batches allowed; nonfinite points return invalid-shape status |
| private `gaussian_log_norm_from_cholesky` | square lower Cholesky factor `l(d,d)`; `d = size(l,1)` | empty factor returns zero; nonempty diagonal is positive |
| private `outer_product` | independent vectors `x`, `y`; result `(size(x),size(y))` | either length may be zero; no equality constraint |

Weights are finite, nonnegative and have a positive finite sum. They are
normalized internally; zero-weight observations are permitted. Invalid shape or
parameter statuses and singular-covariance handling remain native API behavior.
Native evaluation currently returns zero-initialized results on shape/nonfinite
point errors. The C bridge deliberately returns NaNs on native evaluation
failure and Python raises an exception; it does not treat these zeros as valid
densities. C validation failures before native entry preserve output buffers.

The estimator fields use dataset `(d,n)`, normalized weights `(n)` and covariance
and Cholesky caches `(d,d)`. Changing fields directly can invalidate the native
object; the C/Python interface keeps them private and copies input data.

Integration, resampling and marginalization are now exposed through the reviewed
C/Python bridge. Box integration/resampling require explicit RNG state; native
resampling consumes d+1 uniforms per column, and empty batches consume none.
C marginals use distinct one-based coordinates; Python accepts zero-based and
negative coordinates and preserves their order. The child copies storage,
retains the parent bandwidth factor and survives parent destruction.

Covariance and inverse getters now copy into caller-owned `(d,d)` arrays, where
`d` must equal the estimator dimension. Bandwidth updates accept Scott, Silverman
or a positive finite constant. The C/Python bridge commits updates only when
covariance/Cholesky caches remain finite with positive diagonals, preserving the
old object on failure; this guard is a SciFort extension to the wrapper contract.
The native setter is unchanged. Getter outputs never alias estimator storage.
General type-bound-method generation remains a later extension.
