# Continuity tracking of nearly degenerate canonical spectral subspaces

**SciFort extension (2026-10-30).** The `scifort_signal_canonical_subspace_tracking`
module provides `track_canonical_subspaces` for static or time-frequency
results returned by `canonical_coherence_matrix`, its band variants, or the
corresponding real/complex multitaper interfaces. This operation is not part
of `scipy.signal`. It does not need external BLAS or LAPACK libraries.

The existing [`track_canonical_modes`](signal_canonical_tracking.md) matches
**isolated, nonzero individual canonical directions** and deliberately stops
when roots coincide. The new routine instead groups sufficiently close roots
and follows the **subspace spanned by the paired canonical coefficients**.
Individual canonical directions in such a group are not intrinsically unique.

For coordinate-invariant matching under fixed changes of channel units, the
new optional `covariance_metric=.true.` uses the conditioned, regularized
spectral covariance blocks returned by the canonical estimator with
`return_covariances=.true.`. See
[covariance-metric canonical subspace tracking](signal_canonical_covmetric_tracking.md)
for the midpoint inner product, mathematical invariance, edge cases, and
memory implications. The default Euclidean method documented below is
unchanged.

## Usage

```fortran
use scifort_signal, only : dp, canonical_coherence_result, &
    canonical_subspace_tracks_result, multitaper_canonical_coherence, &
    track_canonical_subspaces, canonical_subspace_success
implicit none
real(dp) :: x(1024,5)
type(canonical_coherence_result) :: canonical
type(canonical_subspace_tracks_result) :: tracks

! Populate x(sample,channel) before calling.
canonical = multitaper_canonical_coherence(x,3.5_dp,[1,2],[3,4], &
    conditioned_on=[5],kmax=6,weighting='eigen',ridge=0.02_dp, &
    return_vectors=.true.)
tracks = track_canonical_subspaces(canonical, &
    gap_tolerance=0.001_dp,min_similarity=0.6_dp)
if (tracks%status == canonical_subspace_success) then
    ! tracks%track_id(f,rank) labels the rank's near-degenerate block.
    ! tracks%basis_a and %basis_b are continuity-aligned BLOCK BASES.
end if
```

The generic also accepts a `canonical_coherence_tf_result` and returns a
`canonical_subspace_tracks_tf_result`. Each time window is followed
independently along frequency. Existing canonical-mode and mode-tracking
APIs remain unchanged.

`return_vectors=.true.` is required in the canonical estimator. The source
`directions_stable` mask is **not** required: its false values typically
identify precisely the near ties that this procedure is designed to handle.
The output `%status` reports successful tracking *input validation*, not that
all bins were valid or matched. `%source_status` and `%valid` preserve the
canonical estimator's diagnostics.

## Grouping and matching

At a frequency, adjacent, descending squared canonical strengths `s_i` and
`s_(i+1)` belong to the same block when

\[
 |s_i-s_{i+1}| \leq
 \texttt{gap\_tolerance}\,\max(s_i,s_{i+1},0.05).
\]

The blocks form a consecutive partition of canonical ranks. The default
`gap_tolerance` is `0.001` and must lie in `[0,1]`. This tolerance is a
**heuristic squared-coherence tolerance**, not a statistical confidence
interval. For a block of dimension `d`, define the paired coefficient
matrix, with one scale **per whole A and B block**,

\[
 X=\frac1{\sqrt{2}}\begin{bmatrix}
 A/\sqrt{\|A\|_F^2/d}\\
 B/\sqrt{\|B\|_F^2/d}
 \end{bmatrix}.
\]

The single scale per group is unchanged by unitary column rotations and
avoids assigning a preferred direction within a tied subspace. Its columns
span the concatenated, paired coefficient subspace. Orthonormalize these
columns separately at two frequencies to obtain `Qprev` and `Qcur`. The
singular values of `Qprev^H Qcur` are the **cosines of the principal angles**.
The smallest cosine is the block similarity, with 1 indicating identical
subspaces and 0 indicating an orthogonal direction.

Blocks are eligible to match only when their dimensions agree, their
coefficients are finite and full rank, every squared coherence is positive,
and their smallest principal cosine meets `min_similarity` (default `0.5`,
finite, within `[0,1]`). A rectangular block-pair score table is padded with
zero-score dummy blocks and globally assigned using the Hungarian algorithm.
The objective prioritizes the **number of accepted block matches** and then
the total similarity. Thus a genuine crossing of two independent blocks can
retain their IDs even if their descending strength ranks change.

For an accepted block, let `H = Xprev^H Xcur = U Sigma V^H`. Apply the
right unitary Procrustes rotation

\[
 R=VU^H
\]

to **both** `Acur` and `Bcur`. This minimizes their shared, group-balanced
Euclidean Frobenius distance from the previous aligned basis while
preserving within-group conditional spectral variance normalization. The
whole block rotates as one unit: it does **not** assign a physical meaning
to any particular vector within a repeated-root group.

**Important:** For *nearly* equal but distinct roots, a general unitary block
rotation need not leave paired canonical cross covariance diagonal. Therefore
`%basis_a` and `%basis_b` are **aligned subspace bases**, not individually
canonical singular vectors. `%squared_modes` deliberately retains the
*original* descending singular strengths (and must not be interpreted as the
strength of each rotated basis column). For exact ties the canonical
cross-covariance within a block is proportional to the identity.

## Result fields and interruptions

The static result contains:

| Field | Meaning |
|---|---|
| `%frequencies(f)` | Original frequency grid |
| `%squared_modes(f,k)` | Original descending squared strengths |
| `%basis_a(f,channel,k)`, `%basis_b(f,channel,k)` | Jointly rotated block bases |
| `%block_first(f,k)`, `%block_size(f,k)` | Starting original rank and dimension of the rank's block |
| `%track_id(f,k)` | Persistent subspace ID (shared by all ranks of one block; 0 if invalid) |
| `%principal_cosines(f,k)` | Descending principal cosines within matched block; NaN if unmatched |
| `%similarity(f,k)` | Smallest cosine in the corresponding matched block; NaN if unmatched |
| `%matched(f,k)` | True if block was linked to the immediately preceding bin |
| `%segment_start(f,k)` | True when a new block ID starts at this bin |
| `%valid(f)`, `%source_status`, `%status` | Input spectral validity and two statuses |

Time-frequency results insert `time` immediately after `frequency` in
per-bin arrays, with `%times(:)` preserved. Track IDs are numbered
independently for each time window.

An invalid frequency bin has `track_id=0`; a subsequent valid bin starts new
labels. A nonfinite, repeated, decreasing, or FFT-wrapped frequency breaks
connections. **Splits and merges** of near-degenerate blocks do not invent
one-to-many identities: affected blocks start new labels while unrelated
blocks may continue. Zero-strength and rank-deficient blocks cannot be
matched. A block that fails the overlap threshold starts a new ID.

This is a descriptive, coordinate-dependent tracking diagnostic: Euclidean
coefficient-subspace comparisons are not invariant to rescaling individual
input channels, and do not provide calibrated identification confidence.
Unequal frequency-bin uncertainties or smoothing-induced dependence are not
modeled. Band-integrated canonical results can be tracked only when adjacent
band coordinates really form an interpretable sequence.

## Reproducible verification

```sh
python tools/run_signal_canonical_subspace_tracking.py
python tools/check_signal_canonical_subspace_tracking_public.py
```

The first script exercises checked and optimized Fortran builds, native
repeated-root/split/merge/FFT-wrap/zero-root tests, preceding canonical
regressions, and 100 seeded complex-matrix scenarios against NumPy SVD and
`scipy.optimize.linear_sum_assignment` (2,215 coefficient-mode comparisons
per build). The public check builds the full `scifort_signal` dependency
closure and runs the new and preceding public-interface regression tests.
