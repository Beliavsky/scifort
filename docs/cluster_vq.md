# Vector quantization and K-means clustering

`scifort_cluster_vq` independently implements a subset of
[`scipy.cluster.vq`](https://docs.scipy.org/doc/scipy/reference/cluster.vq.html).
Its only SciFort dependencies are `scifort_kinds` and the existing
`scifort_random` generator. The API uses one-based labels and explicit
integer status codes rather than Python exceptions or warnings.

## Example

```fortran
program example_vq
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_seed
    use scifort_cluster_vq, only : whiten, vq, kmeans, kmeans2, &
        cluster_vq_success
    implicit none

    real(dp) :: x(6,2)
    real(dp), allocatable :: scaled(:,:), centers(:,:), optimized(:,:), distances(:)
    integer, allocatable :: labels(:)
    type(rng_state) :: state
    real(dp) :: distortion
    integer :: status

    x = transpose(reshape([0.0_dp,0.0_dp, 0.0_dp,1.0_dp, 1.0_dp,0.0_dp, &
        10.0_dp,10.0_dp, 10.0_dp,11.0_dp, 11.0_dp,10.0_dp], [2,6]))

    scaled = whiten(x, status)
    if (status /= cluster_vq_success) error stop 'invalid input'

    call rng_seed(state, 20261103)
    call kmeans2(scaled, 2, centers, labels, status, &
        state=state, init='++', iterations=20)
    if (status /= cluster_vq_success) error stop 'kmeans2 failed'

    call vq(scaled, centers, labels, distances, status)
    if (status /= cluster_vq_success) error stop 'vq failed'

    call kmeans(scaled, centers, optimized, distortion, status)
end program example_vq
```

## Public procedures

- `whiten(obs, status)` accepts a vector `(n)` or a matrix `(n, p)` and
  returns the same shape, divided feature-by-feature by the **population**
  standard deviation (no centering). It leaves exactly constant features
  unchanged, consistent with modern SciPy, and rejects NaN or infinity.
- `vq(obs, codebook, labels, distances, status)` assigns each row of an
  observation matrix to its nearest Euclidean centroid. Vector overloads
  support 1D observations. Indices are **one-based**; equal distances select
  the lowest centroid index. Distances are ordinary Euclidean distances,
  not squared distances.
- `kmeans2(obs, k, centers, labels, status, state, iterations, init, missing)`
  performs exactly `iterations` Lloyd updates (default 10). `k` is a positive
  integer no larger than the number of observations. The count overload
  requires explicit mutable `rng_state` for random initialization; `init` may
  be `'points'` (default), `'++'` (K-means++), or `'random'` (Gaussian initialization).
- `kmeans2(obs, initial, centers, labels, status, iterations, missing)` accepts
  a supplied `(k, p)` initial centroid matrix instead of a cluster count.
  Its behavior corresponds to SciPy's `minit='matrix'`; no RNG is used.
  Both K-means generics also accept one-dimensional observation vectors and
  either an integer `k` or one-dimensional initial centroids, returning a
  one-dimensional centroid vector.
- `kmeans(obs, k, centers, distortion, status, state, restarts, threshold, &
  max_iterations)` performs `restarts` independent starts (default 20),
  selecting the codebook with the lowest distortion. The RNG state is explicit.
- `kmeans(obs, initial, centers, distortion, status, threshold, &
  max_iterations)` runs a single Lloyd optimization from a supplied codebook.
  `threshold` is the absolute change in **mean unsquared** Euclidean distance
  used for stopping (default `1e-5`); `max_iterations` is 300 by default.
  Empty centroids are removed, so the returned codebook can contain fewer
  than the originally requested `k` rows.

### Status codes and empty clusters

| Code | Name | Meaning |
|---|---|---|
| 0 | `cluster_vq_success` | Normal completion |
| 1 | `cluster_vq_invalid_input` | Invalid shape, nonfinite input, or invalid option |
| 2 | `cluster_vq_empty_cluster` | `kmeans2` encountered an empty cluster |

For `kmeans2`'s default `missing='warn'`, status 2 signals a **nonfatal**
empty-cluster incident: the unused centroid is retained, and final arrays
are returned. With `missing='raise'`, status 2 is an error and the result
arrays are empty. Library routines do not print warnings or abort execution.

## SciPy differences and numerical limits

- Gaussian `init='random'` draws one codebook from a multivariate normal
  with the **sample mean and unbiased sample covariance** of the observations.
  Rather than calculate a matrix square root, it projects independent normal
  variates through the centered observation matrix, which produces the same
  covariance and works with rank-deficient data. The sampler uses the existing
  SciFort RNG, with Box-Muller normal deviates. It requires at least two
  observation rows; insufficient observations return invalid-input status.
  The time cost of sampling `k` centers is O(k * n * features), and no
  BLAS/LAPACK is required. SciPy can use Cholesky or SVD instead.
  The native default remains `'points'`, **not** SciPy's `'random'`, to
  preserve existing SciFort behavior. SciPy's unused `thresh` argument is
  not implemented. `init='matrix'` is selected by passing initial centers
  through the separate generic overload, not by a string flag.
- RNG streams differ from NumPy and must be supplied through SciFort's
  `rng_state`. Using the same SciFort seed reproduces the same SciFort run;
  it does not reproduce NumPy's pseudo-random centroids.
- `kmeans2` returns labels reassigned to the **final** centroids, whereas
  SciPy versions can return labels from the immediately preceding update.
- `whiten` does not emit SciPy's `RuntimeWarning` for constant features.
- Distances are computed with a scale-normalized Euclidean formula to avoid
  overflow from squaring large coordinates; sufficiently large true distances
  can nevertheless exceed binary64's representable range.
- This implementation emphasizes correctness and portability. Lloyd's method
  only finds a local optimum; unlike SciPy's compiled vector-quantization
  kernels, the current Fortran code does not yet have vectorized/SIMD
  performance tuning.

## Independent numerical validation

Run the standalone native regression (including Gaussian initialization) and
independent SciPy checks:

```sh
python tools/cluster_reference/run_native_vq.py
python tools/cluster_reference/compare_vq_scipy.py --cases 100
```

Each script builds both checked (`-O0 -g -fcheck=all -fbacktrace`) and
optimized (`-O3`) executables in temporary directories. The reference script
uses externally supplied initial centroids to make the deterministic
algorithms comparable across distinct random-number generators.

### Gaussian initialization example

```fortran
use scifort_kinds, only : dp
use scifort_random, only : rng_state, rng_seed
use scifort_cluster_vq, only : kmeans2, cluster_vq_success
implicit none
real(dp) :: x(4,2)
real(dp), allocatable :: centers(:,:)
integer, allocatable :: labels(:)
type(rng_state) :: state
integer :: status
x = reshape([0.0_dp, 1.0_dp, 5.0_dp, 6.0_dp, &
             0.0_dp, 1.0_dp, 5.0_dp, 6.0_dp],[4,2])
call rng_seed(state, 20261103)
call kmeans2(x, 2, centers, labels, status, state=state, &
    init='random', iterations=20)
if (status /= cluster_vq_success) error stop 'kmeans2 Gaussian failed'
```

The additional `test/test_cluster_vq_gaussian.f90` independently reproduces
Gaussian initialization from the supplied RNG and checks its agreement with
explicit-centroid Lloyd updates, reproducibility, rank-deficient matrices,
wide data, zero-variance columns and invalid inputs.
