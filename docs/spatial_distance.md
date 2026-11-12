# `scipy.spatial.distance`-style pairwise metrics (updated 2026-11-02)

SciFort's new `scifort_spatial` public module exposes `point_distance`,
`pdist`, `cdist`, and a rank-overloaded `squareform` for **binary64 real**
observations. The implementation is self-contained and introduces no new
Fortran dependencies or third-party source code. `pdist`, `cdist`, and
`squareform` reproduce SciPy's layout conventions. As with other SciFort
APIs, invalid arguments return an integer status rather than raising a Python
exception.

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : point_distance, pdist, cdist, squareform
implicit none
real(dp) :: x(3,2), y(2,2)
real(dp), allocatable :: condensed(:), pairwise(:,:), cross(:,:), weighted(:,:)
integer :: status

x(1,:) = [0.0_dp, 0.0_dp]
x(2,:) = [3.0_dp, 4.0_dp]
x(3,:) = [2.0_dp, 1.0_dp]
y(1,:) = [1.0_dp, 0.0_dp]
y(2,:) = [2.0_dp, 2.0_dp]

condensed = pdist(x, metric='euclidean', status=status)
! condensed = [5, sqrt(5), sqrt(10)], SciPy's (1,2), (1,3), (2,3) order
pairwise = squareform(condensed)
cross = cdist(x, y, metric='cityblock')  ! shape (3,2)
weighted = cdist(x, y, metric='minkowski', p=3.0_dp, w=[2.0_dp,1.0_dp])
print *, point_distance(x(1,:), x(2,:))  ! 5.0
```

## Supported metrics

The original 11 metrics and the five Boolean metrics below accept nonnegative
finite coordinate weights `w(:)` of length `n_features`, with uniform
weights as the default. A weight of zero
excludes its coordinate from distance computations. **The new covariance and
Jensen--Shannon metrics do not accept `w=`**, consistent with their SciPy
formulas. Passing `p=` is supported **only** for `minkowski`; its default value
is 2. A value `0 < p < 1` produces a quasi-distance rather than a metric,
as with SciPy.

| `metric=` | Computation |
|---|---|
| `euclidean` | `sqrt(sum(w * (u-v)**2))` |
| `sqeuclidean` | `sum(w * (u-v)**2)` |
| `cityblock` | `sum(w * abs(u-v))` |
| `chebyshev` | `max(abs(u-v))` over coordinates with positive weight |
| `minkowski` | `(sum(w * abs(u-v)**p))**(1/p)` |
| `cosine` | One minus weighted cosine similarity |
| `correlation` | One minus weighted centered correlation |
| `hamming` | Fraction of unequal coordinates, weighted when requested |
| `jaccard` | Weighted Boolean-set dissimilarity, using nonzero as true |
| `canberra` | `sum(w * abs(u-v)/(abs(u)+abs(v)))`, omitting `0/0` terms |
| `braycurtis` | `sum(w * abs(u-v))/sum(w * abs(u+v))` |
| `seuclidean` | `sqrt(sum((u-v)**2 / variance))` |
| `mahalanobis` | `sqrt((u-v)^T * inv_cov * (u-v))` |
| `jensenshannon` | Square root of the Jensen--Shannon divergence between normalized nonnegative arrays |
| `dice` | Weighted Dice dissimilarity of Boolean vectors |
| `rogerstanimoto` | Weighted Rogers--Tanimoto dissimilarity of Boolean vectors |
| `russellrao` | Weighted Russell--Rao dissimilarity of Boolean vectors |
| `sokalsneath` | Weighted Sokal--Sneath dissimilarity of Boolean vectors |
| `yule` | Weighted Yule dissimilarity of Boolean vectors |

In the cosine/correlation calculations, vector norms are rescaled to reduce
overflow and underflow. The Euclidean/Minkowski kernels also rescale the
coordinate differences. Results are matched against SciPy for representative
finite inputs. Extremely large differences can still overflow before
rescaling; very large exponent/weight combinations may overflow the final
returned distance. Underflow/overflow of input vector magnitudes and
undefined distances are distinct from invalid parameter errors.

For `cosine` with a zero-norm vector, `correlation` with zero variance, and
`braycurtis` with both numerator and denominator zero, the result is NaN
with **success status**. If Bray--Curtis has zero denominator but positive
numerator, the result is positive infinity. Jaccard of two empty support
sets is zero. Input NaNs and infinities are rejected; SciPy may instead
propagate them to individual output cells.

## Boolean dissimilarities (added 2026-11-02)

`dice`, `rogerstanimoto`, `russellrao`, `sokalsneath`, and `yule` are
supported by **all three** pairwise APIs (`point_distance`, `pdist`, `cdist`).
Inputs are converted to Boolean indicators (`0.0` = false, every other finite
value = true), consistently including negative and non-unit positive values.
Optional finite, nonnegative per-feature weights `w=` replace Boolean counts
with sums of weights. A zero-weight coordinate is ignored.

Let `TT`, `FF`, `TF`, and `FT` denote the (possibly weighted) numbers of
coordinates for which the two vectors are both true, both false, only the
first true, or only the second true. Set `D = TF + FT` and
`N = TT + FF + TF + FT`. The formulas are:

| Metric | Dissimilarity |
|---|---|
| Dice | `D / (2*TT + D)` |
| Rogers--Tanimoto | `2*D / (N + D)` |
| Russell--Rao | `1 - TT/N` |
| Sokal--Sneath | `2*D / (TT + 2*D)` |
| Yule | `2*TF*FT / (TT*FF + TF*FT)` |

When both inputs are entirely false, Dice and Sokal--Sneath return NaN with
**success status**, because their denominators vanish. Rogers--Tanimoto
returns 0 and Russell--Rao returns 1 for those inputs (if total weight is
positive). Yule returns 0 when its denominator is zero, matching SciPy's
current bulk implementation. All-zero coordinate weights give NaN for the
other four metrics and zero for Yule. Finite positive weights are normalized
internally, avoiding overflow of summed weights even at huge magnitudes.

SciPy's 1.17 scalar `sokalsneath` raises an exception for two entirely false
vectors, but its `pdist`/`cdist` methods return NaN; SciFort uses the latter
behavior consistently. SciPy's raw *nonbinary numeric* inputs have
inconsistent behavior among scalar and bulk Boolean metrics, including
negative Dice values. SciFort deliberately interprets any nonzero numeric
value as Boolean true; independent reference tests explicitly convert the
inputs to `bool` before invoking SciPy. Exact 0/1 inputs have the same
interpretation on both sides.

```fortran
real(dp) :: u(4), v(4), d
u = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
v = [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
d = point_distance(u,v,metric='dice')                ! 0.5
d = point_distance(u,v,metric='yule')                ! 1.0
d = point_distance(u,v,metric='russellrao', &
    w=[2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp])            ! 15/17
```

## Covariance and probability metrics (new, 2026-11-02)

The new optional keyword arguments **follow** the existing `status` argument;
existing calls, including positional calls, remain valid:

- `variance(:)` -- positive finite per-feature sample variances for
  `metric='seuclidean'`, equivalent to SciPy's `V` argument.
- `inv_cov(:,:)` -- finite inverse covariance for `metric='mahalanobis'`,
  equivalent to SciPy's `VI` argument. Its shape must be `(n_features,n_features)`.
- `base` -- logarithm base **greater than 1** for `metric='jensenshannon'`;
  omitted uses natural logarithms. Unlike SciPy's current compiled `cdist`
  implementation, the SciFort matrix APIs also accept `base=`.

For **`point_distance`**, standardized Euclidean and Mahalanobis require
`variance=` or `inv_cov=`, respectively. For **`pdist`**, they are optional:
SciFort estimates the sample variances or covariance from the rows of `x`
using `ddof=1`. For **`cdist`**, SciFort estimates them from the combined
rows of `xa` and `xb`, following SciPy's default. The covariance is factored
using a self-contained Cholesky algorithm applied to its correlation matrix,
then rescaled into an inverse covariance; no BLAS/LAPACK is required.
Singular/near-singular covariance, insufficient observations for automatic
Mahalanobis (total observations <= feature dimension), and zero sample
variance give `spatial_distance_invalid_input`. For explicit `inv_cov=`, the
caller is responsible for providing an appropriate positive-definite matrix;
SciFort validates its shape and finiteness but does not enforce symmetry.

Jensen--Shannon accepts nonnegative finite weights and normalizes each vector
separately (vectors need not sum to one). A zero-sum vector is invalid; zeros
in otherwise nonempty probability vectors are supported. The implementation
rescales large/small finite weights before normalization to avoid overflow.

```fortran
real(dp), allocatable :: dis(:), cross(:,:)
real(dp) :: inverse(2,2)

! x contains (observations, features); use sample variances estimated from x.
dis = pdist(x, metric='seuclidean')

! Supplying inv_cov skips covariance estimation entirely.
inverse(1,:) = [1.0_dp, 0.2_dp]
inverse(2,:) = [0.2_dp, 2.0_dp]
cross = cdist(x, y, metric='mahalanobis', inv_cov=inverse)

! Jensen--Shannon distance with the base-2 logarithm.
print *, point_distance([1.0_dp,0.0_dp], [0.0_dp,1.0_dp], &
    metric='jensenshannon',base=2.0_dp)  ! 1.0
```

## Shapes, status and errors

- `point_distance(u(:),v(:),metric,p,w,status,variance,inv_cov,base)` returns one real value.
- `pdist(x(:,:),...)` expects samples **in rows**, returns the condensed
  vector in SciPy's order `(1,2), (1,3), ..., (n-1,n)`.
- `cdist(xa(:,:),xb(:,:),...)` expects `(n_samples,n_features)` and returns
  a `(n_samples_A,n_samples_B)` matrix.
- `squareform(v(:),status)` infers the square matrix dimension from the
  triangular vector length; the empty vector converts to a `1 x 1` zero
  matrix. `squareform(m(:,:),checks,status)` returns its condensed vector.
  By default it requires **exact** symmetry and zero diagonal; pass
  `checks=.false.` to bypass these checks.

`spatial_distance_success=0` indicates that the call was accepted.
`spatial_distance_invalid_input=1` indicates invalid dimensions, nonfinite
observations, unsupported metric, invalid weights, invalid exponent, or
invalid squareform input. An invalid scalar call returns NaN, invalid
`pdist` returns an empty vector, and invalid `cdist` or vector-to-matrix
`squareform` returns `(0,0)`.

This version supports 19 distinct binary64 distance metrics on real arrays
with nonempty feature axes. Callback metrics, reusable precomputed
covariance factors, and `out=` buffers remain future work. `squareform`
currently handles real arrays only.

## Numerical checks

```sh
# Use the SciPy 1.17 reference implementation as an external oracle.
python tools/spatial_reference/compare_distance_scipy.py
python tools/spatial_reference/compare_distance_scipy.py --seed 8675309
```

Each seed runs 495 independent cases for each of two GNU Fortran builds,
covering 11 metrics, 1/2/3/5/8 features, weighted and unweighted vectors,
negative and zero-valued samples, and ordinary/small/large magnitudes.
The native regression `test/test_spatial_distance.f90` additionally covers
condensed index ordering, round trips, invalid parameters, degenerate
metrics, and empty-row matrices.

SciPy's cosine/correlation implementation can overflow or underflow its
intermediate norm products for input magnitudes near `1e+90` or `1e-90`.
Accordingly those two reference metrics use moderate magnitudes, while the
SciFort kernel's scale-invariance is checked independently at wider scales.

For the new metrics run `python tools/spatial_reference/compare_covariance_distance_scipy.py`
(and `--seed 84633`). The two seeds cover 600 scenarios per compiler mode,
including explicit and automatically estimated covariance/variance, positive
and partly-zero probability vectors, non-unit probability totals, and
Jensen--Shannon with natural/base-2 logarithms. SciPy's fast bulk JS path can
report spurious NaNs for proportional vectors; the independent comparison
instead uses its public scalar `jensenshannon` reference for that metric.

For the five Boolean metrics run
`python tools/spatial_reference/compare_boolean_distance_scipy.py`
(and `--seed 86123`). Two independent seeds compare 900 randomized scenarios
in each checked and optimized GNU Fortran build. Scenarios include weighted
and unweighted Boolean vectors, sparse/dense sets, nonbinary numeric inputs,
zero supports, zero-weight coordinates, all-zero weights, and very
large/small finite weight magnitudes. The reference uses explicit NumPy
Boolean inputs to avoid SciPy's nonbinary numeric inconsistencies.
`test/test_spatial_distance_boolean.f90` additionally verifies independent
hand-computed values, invalid inputs, and very large weights.

Other `scipy.spatial` algorithms: [native KD-tree pair and cross-tree searches](spatial_kdtree.md).

New independent point-cloud geometry APIs: [Hausdorff distances and geometric SLERP](spatial_geometry.md).
