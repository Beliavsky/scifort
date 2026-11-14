# General supplied-knot smoothing: validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and results

The native modern scalar and parametric constructors now accept ordered
nonclamped/partially repeated exterior knots and rank-deficient regularized
fits for positive `s`. Interiors remain simple. Periodic fits enforce
endpoint derivatives explicitly on the supplied basis. Conventional banded
paths are retained; a bounded dense QR/SVD path handles other forms and
banded numerical failures. Automatic knot selection is unchanged.

Both `-O0` and `-O3` builds with bounds/runtime checks and floating-point
invalid/zero/overflow traps passed:

- `test_interpolate_fixed_general`: **206 checks**, degrees 1-5, noncyclic
  exterior knots, partial clamping, scalar/shared residuals, unwrapped seam
  derivatives, retained knots/sites, constant limits, rank-deficient clamped
  and cyclic observation matrices, infeasible targets, tiny/huge physical
  intervals, common coordinate/weight scaling, an ignored large final periodic
  weight, enclosing ordinary bounds, and clustered degree-five constant limits.
- Existing periodic fixed-knot (1,548), periodic constructor (206), scalar
  constructor (117), parametric constructor (184), closed tuple (134), periodic
  scalar tuple (115), scalar smoothing (108), parametric smoothing (125),
  parametric interpolation, TCK, and LSQ regressions.
- Existing independent scalar/parametric constructor comparisons (165 fits
  each) and periodic constructor comparisons (200 fits), with prior printed
  metrics preserved for their existing conventional/default cases.
- The refreshed periodic fixed-knot checker: **130 cases**, **60 active
  constraints**, and **six rank-deficient observation fits**, including three
  active cases previously rejected. Other printed metrics remain unchanged.
- The new guide's nonclamped periodic example compiled and ran successfully.

`check_fixed_general.py`, seed **20261114**, checks **160 cases**, including
**80 active constraints** and **28 rank-deficient observation cases**. It
covers scalar and three-coordinate curves, ordinary and periodic boundaries,
degrees 1-5, nonclamped, left/right partially clamped, and clustered clamped
vectors. Constant/polynomial and active smoothing limits are checked.

The independent reference uses SciPy `BSpline` identity coefficients,
`scipy.linalg.null_space` for periodic endpoint constraints, NumPy augmented
least-squares solves, and a scalar root solve matching the achieved residual.
Polynomial limits are independently constructed by collocating global
monomials on a grid in every knot interval; production instead uses Marsden's
identity. No FITPACK smoothing routine serves as the new coefficient oracle.

| Check | Maximum in the two checked builds | Acceptance criterion |
|---|---:|---|
| Active relative residual-target error | 0.000962695 | 0.0011 |
| Normalized unwrapped seam mismatch | 6.3105e-16 | 2e-12, coefficient/knot-gap derivative scaling |
| Independent tuple representation | 1.28838e-10 | 2e-8, including explicit derivative cancellation allowance |
| Scaled normal-equation backward error | 1.40418e-12 | 2e-8 |
| Coefficient error divided by conditioned tolerance | 0.579586 | 1 |
| Raw normalized coefficient difference | 3.79321e-5 | Condition-dependent; see below |
| Augmented-system condition estimate | 2.54248e12 | Diagnostic |

Coefficient comparisons use `2e-7 + 16*epsilon*cond(augmented_system)`.
Polynomial limits are compared directly with `2e-7`. The extreme degree-five
clustered active case has a condition number above `2e12`; its raw coefficient
difference must not be presented as a uniform high-accuracy result. The
worst coefficient error/tolerance ratio is 0.579586 in debug and 0.351056 in
the optimized checked build. Backward errors are 1.40418e-12 and 1.09685e-12,
respectively. Other table maxima cover both builds.

Representation checks allow `2e-8*(1+abs(reference))` plus an absolute
rounding allowance `64*epsilon*(abs(basis_derivative) @ abs(coefficients))`.
This addresses cancellation in clustered-knot derivatives rather than
loosening all tolerances. The normal-equation error is divided by
`norm(A)*norm(A*c-b) + lambda*norm(D)**2*norm(c)` in the endpoint-constraint
null basis. The multiplier is independently recovered by matching the
achieved residual. Endpoint tests never use periodic query wrapping.

The nonperiodic banded fixed-knot residual lookup was corrected to search
across multiple knots between samples. Active-target assertions in the new
independent cases caught that defect. Existing automatic-knot comparisons
retain their prior values and knots.

## Limits and provenance

The dense fallback is capped at 256 coefficients, 10,000 sites, and 256
channels. Normalized exterior magnitude is capped at `1e4`. Existing gap,
relative-weight, scaling, finite-range, and physical-residual checks apply.
These range limits do not guarantee good conditioning; rank handling cannot
restore precision lost to ill-conditioned binary64 problems. Larger dense
problems, repeated interiors, legacy periodic `task=-1`, exact modern automatic
selection, and soft capacities remain gaps. C ABI and wrappers remain deferred.

The new module, tests, and validation code are original MIT SciFort code
composing existing native numerical routines. The polynomial mathematical
reference is Carl de Boor, [*B(asic)-Spline Basics*, Section 4](https://ftp.cs.wisc.edu/Approx/bsplbasic.pdf),
Marsden's identity and polynomial coefficient functionals. The smoothing
objective follows Dierckx (1993) and
[SciPy's public specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html).
No implementation is copied from those mathematical references.

The interval-search correction extends the already attributed perazz/FITPACK
kernels at pinned commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`, retaining
all notices and the full BSD-3-Clause license. No new third-party implementation
or dependency is imported, and no new upstream full continuation-suite run is
claimed. See [provenance](../CODE_PROVENANCE.md) and the
[guide](../docs/interpolate_fixed_general.md).

## Commands

Full repository suite:

```bat
fpm test
```

The full run uses Python `subprocess.run`, redirecting stdout/stderr directly
to `build/fixed-general/full-test.log` and recording the native return code
in `full-test-exit.txt`. The full suite **passed with native exit code 0**,
including all 206 new general fixed-knot checks and the existing constructor
regressions. Existing array-temporary runtime warnings remain informational.

The checked builds compile the changed core, new fixed-smoothing module,
compatibility module, and public interpolation module freshly, linking them
before the existing FPM support archive. These CMD commands reproduce the
new debug checks in this workspace; substitute `build\general-release` and
`-O3` for the optimized checked build. FPM directory hashes can differ.

```bat
set "out=build\general-debug"
if not exist "%out%" mkdir "%out%"
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -I%out% -Ibuild\gfortran_2806A042DAA8BBF7 -J%out%"
gfortran %fitflags% -c src\interpolate_curfit_core.f90 -o %out%\interpolate_curfit_core.o
gfortran %fitflags% -c src\interpolate_fixed_smoothing.f90 -o %out%\interpolate_fixed_smoothing.o
gfortran %fitflags% -c src\interpolate_fitpack_compat.f90 -o %out%\interpolate_fitpack_compat.o
gfortran %fitflags% -c src\interpolate.f90 -o %out%\interpolate.o
set "fitobjects=%out%\interpolate_curfit_core.o %out%\interpolate_fixed_smoothing.o %out%\interpolate_fitpack_compat.o %out%\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a"
gfortran %fitflags% test\test_interpolate_fixed_general.f90 %fitobjects% -o %out%\test_interpolate_fixed_general.exe
%out%\test_interpolate_fixed_general.exe
gfortran %fitflags% validation\fixed_general_driver.f90 %fitobjects% -o %out%\fixed_general_driver.exe
C:\Programs\Python313\python.exe validation\check_fixed_general.py %out%\fixed_general_driver.exe
```

Dependent regressions and existing independent checkers use the same
compile/link/run pattern. Cached sources, orchestration scripts, logs,
executables, and object files remain ignored build artifacts. Later edits
limited to documentation/comments do not require repeating numerical tests.
