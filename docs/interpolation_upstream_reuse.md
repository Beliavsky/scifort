# Interpolation upstream reuse and license review (2026-10-31)

This review covers the projects listed at
<https://beliavsky.github.io/Fortran-code-on-GitHub/#interpolation>.
The goal is to expand SciPy-compatible interpolation **without vendoring whole
libraries** or imposing new system BLAS/LAPACK dependencies on FPM Windows builds.
Links below point to the upstream repositories and their license files.

| Upstream project | License reviewed | Candidate SciFort capability | Decision |
|---|---|---|---|
| [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran) (Jacob Williams) | [BSD-style with incorporated NIST/SLATEC/LAPACK notices](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) | B-spline setup, evaluation, derivatives, integrals, 1D–6D tensor-product splines | **First integration candidate.** Prefer an FPM dependency and adapt the SciFort API. Its bundled BLAS is a better fit than a system BLAS dependency. |
| [finterp](https://github.com/jacobwilliams/finterp) (Jacob Williams) | [BSD-style](https://github.com/jacobwilliams/finterp/blob/master/LICENSE) | 1D–6D nearest/linear grid interpolation, `RegularGridInterpolator`, `interpn` | **Optional back-end / cross-validation candidate.** SciFort now implements 1–6D nearest/linear grids independently; retain finterp's license if integrating its source later. |
| [regridpack](https://github.com/jacobwilliams/regridpack) (Jacob Williams; originally John C. Adams, NCAR) | [BSD-style](https://github.com/jacobwilliams/regridpack/blob/master/LICENSE.txt) | 1D–4D linear/cubic regridding on orthogonal grids | Good alternative/complement to finterp and bspline; verify exact SciPy boundary/coordinate semantics. |
| [PCHIP](https://github.com/jacobwilliams/PCHIP) (Jacob Williams) | [BSD-style, also identifies public-domain SLATEC code](https://github.com/jacobwilliams/PCHIP/blob/master/LICENSE.txt) | PCHIP interpolation/inversion and error handling | Optional alternative, but `scifort_interpolate` already has independent SciPy-compatible PCHIP. Avoid duplicate core implementations. |
| [NdSpline](https://github.com/Takayuki-Miyagi/NdSpline) (Takayuki Miyagi) | [BSD 2-Clause](https://github.com/Takayuki-Miyagi/NdSpline/blob/master/LICENSE) | Multidimensional B-splines | Possible alternative; assess build systems, maintained precision/kind conventions, and API scope. |
| [interpolation](https://github.com/furstj/interpolation) (Jiří Fürst) | [MIT](https://github.com/furstj/interpolation/blob/main/LICENSE) | `interp1` conveniences | Compatible, but much overlaps with currently supported PCHIP/linear interpolation. |
| [akima86](https://github.com/NRLMMD-GEOIPS/akima86) (H. Akima / D. Ryglicki, NRL packaging) | [NRLMMD license and accompanying original permission notice](https://github.com/NRLMMD-GEOIPS/akima86/blob/main/LICENSE) | The improved 1986 Akima algorithm, a distinct extension from SciPy's 1970 Akima and Makima | Review per-file rights before incorporation; avoid confusing the algorithms. |

## Candidates requiring caution

- [interpolation2d3d](https://github.com/arjenmarkus/interpolation2d3d):
  its **wrapper** has an MIT license, but the author says the underlying Renka
  numerical routines remain under the ACM Software License Agreement. Do not
  assume that the MIT wrapper licenses the bundled ACM code for copying.
- [PPR](https://github.com/dengwirda/PPR): its own conditions restrict
  commercial redistribution and how modified code may be distributed. Those
  conditions do **not** fit a simple MIT-licensed SciFort source copy.
- [WenOOF](https://github.com/Fortran-FOSS-Programmers/WenOOF) and
  [FOLLIA](https://github.com/giacrossi/FOLLIA) describe a special
  GPLv3/commercial multi-license scheme. Clarify which terms legally apply
  to SciFort before linking or copying instead of treating all options as
  automatically interchangeable.
- [fortran-rbf](https://github.com/dongli/fortran-rbf) and
  [grid_interpolation](https://github.com/ChristianSteger/grid_interpolation):
  no explicit reusable license was established in this review. A public
  repository without a license does **not** imply permission to reuse code.

## Practical integration sequence

1. SciFort already provides an independently tested, self-contained
   `RegularGridInterpolator`-style nearest/multilinear implementation (1–6D,
   real/complex and trailing channels) in `src/interpolate_regular_grid.f90`.
   Consider `finterp` as an **optional revision-pinned FPM back-end** if its
   performance and behavior prove advantageous. Do not require it for the
   existing tested API. Any future source adaptation must preserve BSD notices.
2. Evaluate `bspline-fortran` as an **optional pinned performance backend** for
   SciFort's now-independent degree-1/3/5 tensor-product interpolation.
   Benchmark per-axis banded solvers, check its built-in BLAS and its FPM
   build on Windows, and preserve notices if copying any implementation code.
3. Add `regridpack` only for capabilities not already supplied by the above
   two packages, with exact validation against `scipy.interpolate`.
4. Retain full upstream copyright, permission and disclaimer text for any
   **copied** source file, and add an entry to `CODE_PROVENANCE.md` naming
   the origin, revision, modifications, and license. For a normal FPM dependency,
   the upstream dependency keeps its own license; document the linkage in
   the README/provenance and verify distribution obligations before shipping.

No upstream interpolation source code was copied in the 2026-10-31
`CubicSpline` increment. Its native Fortran tridiagonal and periodic cyclic
solvers were written independently from the standard spline continuity and
endpoint equations. SciPy was used only for independent numerical regression.

## 2026-10-31 update: native grid interpolation

`scifort_interpolate_regular_grid` provides the nearest/linear 1–6D behavior
listed above, independently implemented and validated against SciPy. The
upstream `finterp` package remains a possible **optional FPM dependency**,
not a required dependency for this new feature. The newly implemented public
API and compatibility limitations are documented in
[`interpolate_regular_grid.md`](interpolate_regular_grid.md). No upstream
interpolation code has been copied as part of this update.

## 2026-10-31 update: SciPy-style 1D B-splines

`scifort_interpolate_bspline` now provides direct `(knots, coefficients,
degree)` evaluation, analytic derivatives/integrals and independent
`make_interp_spline` collocation, plus a dense design matrix and knot
insertion. No upstream code or additional runtime dependency was used.
The richer 1D–6D tensor-product interpolation capabilities of Jacob
Williams' BSD-style `bspline-fortran` remain valuable future targets,
but its degree convention (`order = degree + 1`) and sampled-grid setup
need a separate SciPy-facing adapter. Adopt a version-pinned FPM dependency
only after validating its public API, built-in BLAS portability and test
behavior on Windows and GNU/Intel Fortran. Upstream documents an example
FPM version tag `7.3.0` but that is **not** a pin in SciFort's current manifest.

## 2026-10-31 update: tensor-product cubic and quintic grids

`src/interpolate_regular_grid.f90` now includes 1–6D not-a-knot B-splines of
degrees 1, 3 and 5 with partial derivatives. The fit reuses SciFort's existing
B-spline basis and independently written per-axis collocation code; it does
**not** incorporate source from any of the above projects. `bspline-fortran`
remains a potential optional pinned performance backend. The FPM manifest
remains unchanged, and no third-party license files are newly required.

## 2026-10-31 update: recursive grid PCHIP and fast banded collocation

SciFort now also provides real-valued `RegularGridInterpolator(method='pchip')`
with the SciPy-compatible last-to-first recursive evaluation. The existing
`'slinear'`/`'cubic'`/`'quintic'` fit was refactored to use banded B-spline
collocation and an original no-pivot LU solver, eliminating the earlier
quadratic-factor-storage and cubic-factor-time bottleneck. This reduces the
immediate performance motivation for an optional `bspline-fortran` backend;
its other numerical features remain of interest. No external implementation
source was copied or new FPM dependency introduced.
