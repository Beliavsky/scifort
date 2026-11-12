# Frequency-continuous canonical spectral modes

**SciFort extension (2026-10-30).** `scifort_signal_canonical_tracking`
provides an optional continuity-based relabeling of the canonical modes from
`scifort_signal_canonical_coherence`. This is **not** part of the public
`scipy.signal` API. No external Fortran BLAS or LAPACK library is required.

Canonical correlation modes are normally sorted from strongest to weakest at
each spectral frequency. Sorting by strength can swap the identity of two
physically continuous modes whenever their squared coherences cross. A related
issue is the arbitrary common complex phase of their paired coefficient
vectors. `track_canonical_modes` addresses these visualization and
interpretation problems **without changing** the original canonical estimate.

## Usage

```fortran
use scifort_signal, only : dp, canonical_coherence_result, &
    canonical_mode_tracks_result, multitaper_canonical_coherence, &
    track_canonical_modes, canonical_tracking_success
implicit none
real(dp) :: x(1024,5)
type(canonical_coherence_result) :: original
type(canonical_mode_tracks_result) :: tracks

! Populate x with the observations, indexed (sample,channel).
original = multitaper_canonical_coherence(x,3.5_dp,[1,2],[3,4], &
    conditioned_on=[5],kmax=6,weighting='eigen', &
    ridge=0.02_dp,return_vectors=.true.)
tracks = track_canonical_modes(original,min_overlap=0.6_dp)
if (tracks%status == canonical_tracking_success) then
    ! tracks%squared_modes(f,track) follows continuity rather than rank.
    ! tracks%source_mode(f,track) is the original descending mode index.
end if
```

The same generic procedure accepts a `canonical_coherence_tf_result` from
`multitaper_tf_canonical_coherence` or the time-frequency matrix/band APIs;
it returns `canonical_mode_tracks_tf_result`. Tracking proceeds over **frequency
only**, independently within each time window. Results produced by band
integration also work, but adjacent ordered bands should be scientifically
comparable; unrelated bands do not establish a continuous trajectory.

`return_vectors=.true.` is required for the canonical estimator. If the
coefficients or stability mask were not requested, or the result array shapes
are inconsistent, tracking returns `canonical_tracking_invalid_input = 1`.
Otherwise `status = canonical_tracking_success = 0` indicates that the
*tracking operation* was valid; it does not mean every source spectral cell
was invertible or every direction was successfully matched. The source
estimator's status is retained in `%source_status`, and its `%valid` mask is
copied. `min_overlap` defaults to `0.5` and must be finite and in `[0,1]`.

## Matching and phase convention

Let `a_i(f)` and `b_i(f)` denote the two groups' coefficient vectors of mode
`i`. The complex similarity for a previously tracked mode `i` and a current
strength-sorted mode `j` is

\[
 h_{ij}=\frac12\left(
 \frac{a_i(f-1)^H a_j(f)}{\|a_i(f-1)\|_2\|a_j(f)\|_2}
 +\frac{b_i(f-1)^H b_j(f)}{\|b_i(f-1)\|_2\|b_j(f)\|_2}
 \right),\qquad q_{ij}=|h_{ij}|.
\]

Normalize the A and B overlaps separately before averaging. The complex
averaging matters: the *same* complex phase must rotate both members of a
paired mode, preserving its nonnegative real canonical cross-spectrum.
The algorithm computes a **globally optimal one-to-one assignment**
maximizing `sum(q[i,permutation(i)])`, using an O(m^3) Hungarian algorithm
for `m = min(size(A),size(B))` modes. After a successful match it multiplies
**both** coefficient vectors by `conjg(h)/abs(h)`, making the average
cross-frequency overlap real and positive.

This is a Euclidean coefficient-coordinate similarity, not a covariance-
weighted, subspace-principal-angle, or statistically optimal correspondence.
Because the conditional spectral covariance can change with frequency,
coefficient magnitudes and orientations can vary. Scaling individual input
channels changes the similarity metric. The results should be interpreted as
heuristic continuity diagnostics rather than mode-identification guarantees.
A permutation can change mode ordering: `%squared_modes` in the **tracking
result is not generally descending**. Use the original `%squared_modes` when
strength ordering is needed.

## Outputs and interruption rules

For the static result the shapes are:

| Field | Meaning |
|---|---|
| `%squared_modes(f,track)` | squared coherence in continuity order |
| `%coefficients_a(f,i,track)` | phase-aligned A coefficients in input group order |
| `%coefficients_b(f,j,track)` | phase-aligned B coefficients in input group order |
| `%directions_stable(f,track)` | original numerical nondegeneracy flag, reordered |
| `%source_mode(f,track)` | original 1-based descending mode index; zero at invalid cells |
| `%overlap(f,track)` | `q` of the assigned pair, or NaN when not linked |
| `%matched(f,track)` | true only if accepted from the immediately preceding bin |
| `%segment_start(f)` | valid bin at which tracking begins or restarts |
| `%valid(f)` | original spectral validity mask |
| `%frequencies`, `%source_status`, `%status` | unchanged coordinates and statuses |

The TF result inserts the time index `t` after `f` in all the per-mode arrays,
with `%segment_start(f,t)` and `%valid(f,t)`. The `%times` array is preserved.
Original aggregate group coherence, information and shrinkage metadata remain
in the **original** canonical result, not the tracking result.

The algorithm starts a new segment (leaves the original strength ordering and
phase gauge for that cell, with `matched = .false.`) whenever:

- a current or preceding cell is invalid or lacks a usable coefficient vector;
- any root in the current or preceding cell is marked numerically unstable
  (zero or tied/near-tied directions);
- frequency coordinates are non-finite, repeated, decreasing, or FFT-wrapped
  between adjacent cells;
- any globally assigned similarity is below `min_overlap` or numerically zero.

For tied or nearly equal roots, use the complementary [`track_canonical_subspaces`](signal_canonical_subspace_tracking.md): it follows joint block spans rather than asserting unique individual eigenvectors.

An unstable root interrupts tracking of **the entire group**, avoiding
unwarranted claims of individual direction identity in a degenerate subspace.
A stable, valid bin after a reset serves as a new starting point and can be
matched to the next valid ascending bin. The algorithm does not interpolate
over missing bins or track modes across time windows. Successful matching
means numerical similarity passed a chosen threshold; near-equal assignment
scores can still make labels ambiguous, and no statistical confidence level
is implied. Tracking phase alignment does not change the conditional
variance normalization or the paired spectral covariance.

## Verification

```sh
python tools/run_signal_canonical_tracking.py
python tools/check_signal_canonical_tracking_public.py
```

The native tests cover genuine canonical singular-vector crossings, complex
cross-spectrum phase conservation, an explicitly non-greedy 3-mode global
assignment, FFT-wrap boundaries, unstable/invalid bins, insufficient overlaps,
parameter validation, and two independent time-frequency windows. The
reference script generates 80 reproducible complex-unitary spectral trajectories
and compares Fortran matchings, phases and continuity values to NumPy complex
matrix products and `scipy.optimize.linear_sum_assignment` for 2,307 modes.
Both checked and optimized builds are exercised by the first script.
