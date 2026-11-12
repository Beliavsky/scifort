# Root-finding upstream reuse review (2026-10-31)

Source catalog: <https://beliavsky.github.io/Fortran-code-on-GitHub/#root-finding>.
Reviewed to identify prospective *attributed* reuse in SciFort without mandatory
system BLAS/LAPACK libraries, especially for FPM builds on Windows.

| Package | Licensing finding | SciFort fit | Decision |
|---|---|---|---|
| [polyroots-fortran](https://github.com/jacobwilliams/polyroots-fortran) (Jacob Williams) | BSD-style; separate provenance/attribution conditions exist for incorporated legacy algorithms | All real/complex roots of general polynomials (`rpoly`, `cpoly`, Aberth-family methods), potentially a high-degree `PPoly.solve` backend | **First choice for polynomial-root backend**, after selecting a precise implementation and checking the notices for that algorithm |
| [roots-fortran](https://github.com/jacobwilliams/roots-fortran) (Jacob Williams) | BSD-style | Scalar bracketing methods (Brent, bisection, regula falsi); good candidate for `scipy.optimize.root_scalar` and for the bracketing stage of `PPoly.solve` | **First choice for scalar root finding** |
| [gf_solver](https://github.com/gfell/gf_solver) | MIT-labeled repository; inspect its attribution/citation language before redistribution | Hybrid Newton and bisection, somewhat overlapping `roots-fortran` | Optional; no strong immediate gain over the above |
| [FORTRAN-Roots](https://github.com/bryand4140/FORTRAN-Root-Finding-Toolbox) | A clearly applicable source-redistribution license was not established | Scalar and polynomial roots | **Do not copy without explicit permissions** |
| [GetRootsF](https://github.com/physicsfan/GetRootsF) | A clearly applicable source-redistribution license was not established | Root finding with/without derivatives | **Do not copy without explicit permissions** |
| [gssFortran](https://github.com/huijunchen9260/gssFortran) | Not a standalone SciPy root-solving library (golden-section optimization demonstration) | None here | No integration |

## Preferred integration strategy

1. Keep SciFort's current MIT-licensed **independent** piecewise-polynomial
   algorithm as the default. It recognizes roots of even multiplicity and
   discontinuity crossings, emits NaN sentinels for identically-zero intervals,
   handles descending knots, and preserves SciPy's interval-order semantics.
   These are *additional* tasks beyond finding all roots of one polynomial.
2. Evaluate `polyroots-fortran` as an **optional, version-pinned FPM backend**
   for polynomial degrees above a measured threshold. Filter roots with a
   numerical-imaginary-part tolerance, clip to the interval domain, recheck
   residuals, and retain the existing zero/discontinuity handling and ordering.
3. Evaluate `roots-fortran` as an optional backend for safeguarded scalar
   root isolation and eventually `scipy.optimize.root_scalar`. Its FPM package
   supports `real64`, matching SciFort's `dp` convention.
4. Add tests for badly scaled and clustered polynomial roots, repeated roots,
   interval endpoints, extrapolation, and both real and complex coefficients.
   Benchmark the upstream solver against the independent SciFort solver.
5. Before incorporating **any** source, record the exact upstream revision,
   original author, algorithm-level provenance, copyright, license conditions,
   and required notices. A top-level BSD or MIT license is not sufficient to
   override restrictions or notices carried by an embedded third-party work.

No upstream source was copied, no dependency was added to `fpm.toml`, and no
upstream implementation was claimed to be part of SciFort in this update.

Reference documentation:
- <https://github.com/jacobwilliams/polyroots-fortran>
- <https://github.com/jacobwilliams/roots-fortran>

## Update: opt-in adapter (2026-10-31)

A pinned `polyroots-fortran` v1.2.0 FPM subpackage is now supplied at
`optional/polyroots_backend/` with an `rpoly` callback adapter. It is not
a main SciFort dependency. The adapter was interface-compiled using a stub,
but it has **not** been built or benchmarked against the genuine upstream
package in this environment; see `docs/interpolate_polyroots_optional_backend.md`.
