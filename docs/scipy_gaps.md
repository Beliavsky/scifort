# Remaining SciPy functionality

Initial inventory: 2026-11-11; statistics entry updated 2026-11-12.
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
| FFT | A complete standalone [`scipy.fft` API](https://docs.scipy.org/doc/scipy/reference/fft.html), including cosine transforms, real/Hermitian multidimensional transforms, and fast Hankel transforms. SciFort already uses FFTs internally. |
| Integration | [`solve_bvp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html) for boundary-value ODEs, Lebedev quadrature, and some QMC integration options. |
| Interpolation | AAA rational approximation, spherical spline interfaces, and additional spline-construction APIs. Automatic-knot smoothing currently relies on the optional FITPACK backend. |
| Optimization | Stored sparse Jacobians and worker parallelism in least squares; stronger sparse LP/MILP capabilities; SciPy-equivalent SHGO homology stopping and sampling options. |
| Special functions | Additional complex-valued cases, including complex zeta and complex zeros of selected Bessel Y functions. |
| Statistics | Censored likelihoods and bounded fits are available for normal, exponential, Weibull, lognormal, gamma, beta, logistic, and both Gumbel distributions (lognormal/gamma fits require fixed location, beta fits require fixed location and scale, and gamma/beta fits exclude exact support endpoints); extend censoring support to other families. Observation weights and additional method/option compatibility remain. See [censored fitting](censored_fitting.md). Distribution-name coverage alone does not establish complete compatibility. |
| Sparse graphs | Functions such as `min_weight_full_bipartite_matching` and additional graph conversion helpers. |
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
