# Banded large-input spline selection

Native `generate_knots` and automatic positive-target `make_splrep` /
`make_splprep` support up to **10,000 sites and 256 channels**, degrees 1-5,
for nonperiodic and periodic data. Both constructors consume the same final
selected knot vector and retain their existing derivative-jump smoothing solve.
The [soft `nest` contract](interpolate_modern_selection.md) applies throughout
these limits, including feasible overshoots at thousands of sites.

Larger residual solves reuse existing nonperiodic B-spline banded QR and the
attributed periodic FITPACK cyclic QR. The latter rejects a deficient triangular
factor before back substitution. Small problems (at most 256 sites) and short
periodic coefficient cycles retain their existing dense minimum-norm path.
A failed banded solve can use dense QR/SVD only with at most 256 independent
coefficients and at most one million active sample-by-coefficient entries.
Beyond that budget it returns status 3. This is a bounded fallback, not a
promise to solve every rank-deficient or ill-conditioned large problem.

Banded least-squares storage grows with sites and coefficients times degree
and channels, rather than their product. At 10,000 sites, degree 3 and one
channel, the cyclic kernel's declared buffers occupy about **2.481 MiB**;
a full dense design alone would occupy **763.092 MiB**. These are buffer-size
estimates with binary64 and four-byte default integers, not measured process
resident memory. They exclude caller state, expression temporaries, allocator
and runtime overhead. The end-to-end pipeline also retains input, output and
smoothing workspace. Time still depends on the number of refinement steps;
this addition does not give a fixed execution-time guarantee.

Above 10,000 sites or 256 channels, constructors retain their prior banded
legacy automatic path and hard `nest` capacity. The standalone generator
returns status 3 outside its size limits. Supplied-knot paths, zero-target
constructor interpolation and legacy tuple APIs keep their existing contracts.

Normalization and range guards remain those in the
[generator guide](interpolate_generate_knots.md). A successful selection does
not guarantee that a near-interpolation or ill-conditioned smoothing target
is attainable. Explicit failures never silently return an unmet residual.
Canonical periodic exteriors intentionally differ from stale exteriors in
SciPy v1.18.0 selection. Root stopping within the permitted 0.001 residual
band can still produce different coefficients and higher derivatives.

The [validation record](../validation/large_selection_2026-11-14.md) includes
301-10,000-site scalar/shared fits, complete independent knot sequences,
soft-count overshoots, cyclic rank failures and a periodic solve with 9,999
independent coefficients. Modern selection beyond the size limits, repeated
interiors, broader conditioning guarantees and exact periodic/root-iteration
parity remain gaps. C ABI and wrappers remain deferred.
