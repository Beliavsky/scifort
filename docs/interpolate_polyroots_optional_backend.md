# Optional `polyroots-fortran` backend for SciFort `PPoly.roots` and `.solve`

The MIT-licensed SciFort `PPoly` module now accepts an optional procedure
callback for finding all real zeros of one local polynomial. Its existing
internal derivative-isolation algorithm remains the default; nothing new
is required for normal `fpm` builds.

```fortran
use scifort_interpolate, only : ppoly, ppoly_type, piecewise_native_real_roots
! The following module exists ONLY in the opt-in package below:
use scifort_optional_polyroots, only : polyroots_segment_solver

type(ppoly_type) :: pp
real(dp), allocatable :: roots(:)
pp = ppoly(reshape([1.0_dp,-3.0_dp,3.0_dp,-3.0_dp,2.0_dp,0.0_dp], &
    [6,1]),[0.0_dp,3.0_dp],extrapolate=.false.)
roots = pp%roots(root_solver=polyroots_segment_solver)
! Likewise: roots = pp%solve(y=0.5_dp,root_solver=polyroots_segment_solver)
```

The core `%roots(root_solver=...)` and `%solve(root_solver=...)` preserve
SciFort's SciPy-compatible interval semantics: local coordinate shifts,
zero-polynomial NaN sentinels, discontinuity sign crossings, ascending and
descending knots, channel selection, extrapolation policy, and segment
ordering. Callback outputs are sorted, interval-clipped, and validated to be
finite (invalid outputs trigger a native-solver fallback). The generic
callback signature is:

```fortran
function my_solver(coeff,lower,upper) result(roots)
    real(dp), intent(in) :: coeff(:)   ! descending power coefficients
    real(dp), intent(in) :: lower      ! finite inclusive lower local x
    real(dp), intent(in) :: upper      ! finite inclusive upper local x
    real(dp), allocatable :: roots(:)
end function
```

`piecewise_native_real_roots(coeff,lower,upper)` is also public, providing
an explicit fallback for adapters.

## Optional FPM project and licensing

The self-contained SciFort project at `optional/polyroots_backend/` declares:

```toml
[dependencies]
scifort = { path = "../.." }
polyroots-fortran = { git = "https://github.com/jacobwilliams/polyroots-fortran.git", tag = "1.2.0" }
```

Build it explicitly only when desired:

```bash
cd optional/polyroots_backend
fpm test --profile release
```

It uses Jacob Williams's `polyroots_module:rpoly` (real-coefficient
Jenkins–Traub algorithm), with normalized filtering of candidates. Degree
below four and upstream solver failure fall back to SciFort's native
routine. This optional package is intentionally *not* part of the main
`fpm.toml`; no source from `polyroots-fortran` is redistributed in SciFort.

**Important Windows limitation:** the pinned upstream package's `fpm.toml`
declares system `lapack` and `blas` link dependencies for the *whole*
package, even though `rpoly` itself does not require either numerical
library. A standard upstream FPM build may consequently be unsuitable on
Windows systems without those libraries. Upstream package availability,
linking and actual run-time performance must be verified on the target
system before replacing SciFort's native root solver in production.

Original upstream project and license: Jacob Williams,
<https://github.com/jacobwilliams/polyroots-fortran>, BSD-style license
plus applicable algorithm-specific attribution notices. Keep all upstream
notices as required when FPM retrieves the package. No third-party code
has been copied into this repository.

## Verified offline build using the genuine uploaded upstream source (2026-10-31)

The genuine `polyroots_module.F90` supplied in
`polyroots-fortran-master.zip` has now been compiled and linked with the
optional SciFort adapter on Linux using GNU Fortran 14.2.0. **This is not a
verification of the separately pinned upstream release tag `1.2.0`: the
uploaded source ZIP is from `master`.** The upstream source is not copied into
SciFort. See `docs/polyroots_real_upstream_validation_2026-10-31.md`.

To build and execute the actual upstream solver, the SciFort callback test,
120 additional constructed-polynomial cases for degrees 5--25, and the
benchmark without FPM or a network connection, run:

```bash
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip
```

The default `--link-mode sections` uses GNU Fortran's
`-ffunction-sections -fdata-sections` and the GNU linker option
`-Wl,--gc-sections`, so unused companion-matrix routines referring to
`dgeev`/`zgeev` are discarded and the `rpoly` adapter links without LAPACK
or BLAS. This was tested successfully. A second, conventional link mode is
available and was tested too:

```bash
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip --link-mode lapack
```

The latter requires installed LAPACK and BLAS libraries. The former requires
a linker capable of garbage-collecting unused function sections. **The
original upstream FPM manifest still requests BLAS/LAPACK**, so an FPM build
using that manifest is not automatically dependency-free. Neither Windows
MinGW nor the actual FPM dependency-resolution path was tested.

Upstream cannot be compiled unchanged with `gfortran -std=f2018`, because
of one `COMPLEX*16` declaration inside the LAPACK `zgeev` interface. This
does not affect the normal GNU Fortran build above. A strictly conforming
local patch can replace this with `complex(real64)` and add
`use iso_fortran_env, only: real64` inside the `zgeev` procedure interface;
a temporary copy with those changes compiled under `-std=f2018`.

## Tests and validation limits

`test/test_interpolate_root_callback.f90` checks the real core callback
independent of polyroots-fortran. `tools/interpolate_reference/
check_optional_polyroots_adapter.py` compiles the optional adapter against a
**small signature-equivalent test double**, in both checked and optimized
builds, and checks basic root filtering and invocation. This remains only
an interface test. The separate
`tools/interpolate_reference/check_optional_polyroots_local.py` now validates
the **genuine uploaded upstream** source in both builds, runs deterministic
root-location tests, and executes the benchmark. FPM itself and Windows
compatibility are not yet validated.

When upstream prerequisites are available, run the included benchmark:

```bash
cd optional/polyroots_backend
fpm run --profile release --target benchmark_polyroots
```

It prints CPU times for 100 repetitions for each of four odd degrees,
plus the maximum difference in root locations. It was successfully run
using the genuine uploaded upstream source on Linux; the small-sample
timings should not be interpreted as general performance rankings.
