# Remaining SciPy functionality

Initial inventory: 2026-11-11; coverage entries updated 2026-11-13.
This preliminary table summarizes gaps identified from
SciFort's current source, documentation, roadmap, and an initial API-name
comparison. It is not an exhaustive compatibility audit. A matching name or
passing native test does not establish complete SciPy behavior, option,
data-type, dimensionality, accuracy, or performance compatibility.

Some older roadmap entries describe work that has since been completed,
including dependency-backed FPM tests. Verify candidates against current
source and topic guides before starting an implementation.

For the complete list of remaining censored-distribution candidates, see
the [censored distribution coverage inventory](censored_distribution_gaps.md).

| Area | Remaining functionality or compatibility |
| --- | --- |
| FFT | Standalone 1D, 2D, and rank 1-15 multidimensional FFT pairs, DCT/DST types I-IV, selectable axes, frequency bins, shifts, backend-specific fast-length helpers, and [FFTLog fast Hankel transforms](hankel.md) are available in [`scifort_fft`](fft.md). Remaining gaps include workers, public reusable transform plans beyond the existing 1D FFT plan, and backend/option parity. Hankel singular/range handling and wide-bias accuracy limits are documented. |
| Integration | Some QMC integration options. All 32 [Lebedev sphere rules](lebedev.md) are now available. Native real/complex [`solve_bvp`](bvp.md) now covers adaptive two-point systems, unknown parameters, singular terms, analytic Jacobians, and dense output; large-system backend and iteration-path parity remain limited. |
| Interpolation | Native fixed-knot [spherical least-squares splines](interpolate_sphere.md), [scattered spherical smoothing](interpolate_smooth_sphere.md), and [rectangular spherical splines](interpolate_rect_sphere.md) (interpolation and adaptive smoothing, with pole options) are available within documented limits. Additional spline-construction APIs, larger-problem limits, and broader rank-deficient compatibility audits remain. Native real/complex [AAA rational approximation](interpolate_aaa.md) is available; its poles, zeros, residues, and cleanup now support up to 512 retained support points, with documented conditioning limits. Native nonperiodic scalar `splrep(s>0)` and parametric `splprep(s>0)` now provide automatic-knot weighted smoothing for degrees 1-5; periodic and other automatic-knot smoothing interfaces still rely on the optional FITPACK backend. Modern constructors and smoothing-iteration/conditioning parity remain gaps. |
| Optimization | Stored sparse Jacobians and worker parallelism in least squares; stronger sparse LP/MILP capabilities; SciPy-equivalent SHGO homology stopping and sampling options. |
| Special functions | Native [complex Riemann zeta](complex_zeta.md) and [complex Bessel Y zeros](complex_bessel_zeros.md) are available within documented limits. Higher-imaginary-argument zeta coverage and a broader function-by-function complex-domain/option audit remain; the three Y-zero helpers support up to 1000 complex roots with explicit companion-sign and ordering details. |
| Statistics | Censored likelihoods and bounded fits are available for normal, exponential, Weibull, lognormal, gamma, beta, logistic, and both Gumbel distributions (lognormal/gamma fits require fixed location, beta fits require fixed location and scale, and gamma/beta fits exclude exact support endpoints); extend censoring support to other families. Observation weights and additional method/option compatibility remain. See [censored fitting](censored_fitting.md). Distribution-name coverage alone does not establish complete compatibility. |
| Sparse graphs | Native sparse weighted full matching is available as `min_weight_full_bipartite_matching` and the original `minimum_weight_full_bipartite_matching`; see [graph compatibility details](sparse_csgraph.md). Additional graph conversion helpers and input-format/option parity remain. |
| Linear algebra | Broader batching, larger-problem backends, stronger conditioning/error handling, and real generalized Schur output. Most major algorithm families already exist. |

C ABI coverage and the planned Python, R, MATLAB, and Octave wrappers also
remain incomplete. These are interoperability tasks distinct from implementing
numerical algorithms. Current work prioritizes native Fortran SciPy coverage;
C ABI and wrapper expansion are deferred.

For the project's statistics-first objective, prioritize distribution fitting,
special-function gaps, and a current function-by-function compatibility inventory
before expanding into additional SciPy subpackages. New scope remains subject
to the approval requirements in [repository rules](../AGENTS.md).

## References

- [SciPy public API](https://docs.scipy.org/doc/scipy/reference/index.html)
- [SciFort roadmap](../ROADMAP.md)
- [SciFort API reference](api.md)
- [Dense linear algebra coverage audit](linalg_coverage_audit_2026-11-10.md)
- [SHGO guide](optimize_shgo.md)
- [MILP guide](optimize_milp.md)
- [FITPACK compatibility and optional backend](interpolate_fitpack_tck_backend.md)
