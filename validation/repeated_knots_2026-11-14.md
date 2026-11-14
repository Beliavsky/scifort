# Repeated supplied knots validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and upstream audit

The exact SciPy v1.18.0 `fpcheck`, `_make_splrep_impl` and `disc` source and
root BSD license were inspected. Supplied interior knots must be distinct
in SciPy's smoothing constructors, and `disc` assumes simple interiors.
This feature is a **SciFort extension**, correcting the earlier gaps entry;
no repeated-knot SciPy smoothing parity is claimed. Installed SciPy's public
`make_lsq_spline` is separately checked on full-rank ordinary repeated-knot
cases; that interface minimizes residual without the smoothing penalty.

The extension retains SciFort's kth-piecewise-derivative-jump objective,
counting each distinct knot once. The zero-jump space can contain lower-
order corners and discontinuities. Tests therefore check its full numerical
null space instead of assuming polynomial/constant limits. See the
[contract](../docs/interpolate_repeated_knots.md) and
[provenance](../CODE_PROVENANCE.md).

## Native coverage

Checked debug `-O0` and optimized `-O3` builds passed **365 checks** in
`test_interpolate_repeated_knots`:

- Degrees 1-5 and every interior multiplicity from 2 through k+1, with
  clamped and nonclamped exterior knots, scalar and shared-channel fits.
- Zero-jump reproduction of truncated lower-degree powers and steps;
  repeated knots and parameter arrays retained exactly.
- Feasible positive-target smoothing and comparison with existing ordinary
  repeated-knot `make_lsq_spline` least-squares residuals.
- Periodic piecewise-constant recovery with discontinuous interiors,
  right-continuous knot values and endpoint derivatives through k-1.
- Excess multiplicity, extra base endpoint copies, NaN knots, unattainable
  targets, invalid objects and parametric-output cleanup.

Both profiles also passed general fixed-knot **206**, supplied periodic
**1,548**, large modern selection **271**, modern selection **123**, scalar
constructor **112**, parametric constructor **179**, and periodic constructor
**196** checks. The existing simple-knot general checker passed **160 cases**,
including 80 active constraints and 28 rank-deficient cases. Its existing
conditioning-aware bounds remained unchanged; its maximum normalized
coefficient difference was 3.79321e-5 on clustered high-degree cases.

## Independent repeated-knot checker

`check_repeated_knots.py` uses seed **20261114**, 61 observation sites,
nonuniform weights, one/three channels, ordinary/periodic boundaries,
clamped/nonclamped exterior knots, uniform and mixed interior multiplicities,
and clustered high-multiplicity knots. SciPy `PPoly.from_spline` leading
coefficients give kth derivatives on nonempty polynomial intervals; their
differences define one jump per distinct interior knot. This avoids mirroring
the native midpoint-derivative construction.

SciPy null spaces impose periodic endpoint constraints and construct the
zero-jump limit. Independent NumPy augmented least-squares systems and a
Brent multiplier search reconstruct each active solution at its actual native
residual. Comparing at that residual avoids attributing permitted root-stop
differences to coefficient error. Values, piecewise derivatives through
`min(3,k)`, and integrals are independently checked against SciPy B-splines.
Queries include each interior break, its immediately adjacent representable
values, and extrapolation. Periodic derivatives through k-1 must match at
the endpoints. Common-multiplier stationarity and minimum-norm components
in joint unobserved/unpenalized null directions are checked separately.

| Check | Debug | Optimized |
| --- | ---: | ---: |
| Successful fits | 120 | 120 |
| Active residual constraints | 60 | 60 |
| Rank-deficient observation cases | 40 | 40 |
| Maximum normalized coefficient difference | 0.000236308 | 0.000236308 |
| Maximum normalized representation difference | 1.39234e-10 | 1.41684e-10 |
| Maximum scaled periodic seam difference | 7.57693e-16 | 5.84317e-16 |
| Maximum relative target error | 0.000992259 | 0.000992259 |
| Maximum scaled stationarity error | 2.98087e-15 | 2.98087e-15 |
| Maximum relative joint-null component | 2.53304e-11 | 2.53303e-11 |
| Maximum retained-subspace condition | 5.55122e13 | 5.55099e13 |
| Coefficient error / conditioned tolerance | 0.00289291 | 0.00289291 |

Active residuals use a `0.0011` relative check, stationarity `2e-8`, and
joint-null components `2e-7`. Coefficient comparisons use normalized error
and `2e-7 + 16*epsilon*condition`, with the condition of the numerically
retained singular subspace, not a formally infinite full-matrix condition.
The retained cutoff is `max(rows,columns)*epsilon*sigma_max`, matching the
rank convention of the existing native solvers. Representation comparisons
use the earlier rounding-aware `2e-8` check. The reported large condition
numbers and coefficient differences are deliberate limits: small backward
errors do not imply high forward accuracy for clustered repeated knots.

## Commands

Local ignored orchestration scripts:

```bat
python build\repeated-knots\checks.py debug
python build\repeated-knots\checks.py release
```

They rebuilt these sources in dependency order:

```text
src/interpolate_curfit_core.f90
src/interpolate_fixed_smoothing.f90
src/interpolate_periodic_lsq.f90
src/interpolate_generate_knots.f90
src/interpolate_fitpack_compat.f90
src/interpolate.f90
```

The debug source compile, for example, used:

```bat
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\repeated-knots-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\repeated-knots-debug -c src\interpolate_fixed_smoothing.f90 -o build\repeated-knots-debug\interpolate_fixed_smoothing.o
```

Tests and drivers linked the fresh objects followed by
`build\gfortran_3C159EB802C22FB9\scifort\libscifort.a`, using the same flags.
The optimized build used `-O3` and `build\repeated-knots-release`, retaining
checks and traps. Supporting archive objects were not independently rebuilt
at both optimization settings. The final expanded 365-check native test was
rebuilt/run separately in debug after the initial 355-check run.

Final independent repeated-knot commands:

```bat
python validation\check_repeated_knots.py build\repeated-knots-debug\fixed_general_driver.exe
python validation\check_repeated_knots.py build\repeated-knots-release\repeated_knots_driver.exe
```

The unchanged existing checker ran against each build's `fixed_general_driver`:

```bat
python validation\check_fixed_general.py build\repeated-knots-debug\fixed_general_driver.exe
python validation\check_fixed_general.py build\repeated-knots-release\fixed_general_driver.exe
```

Full-suite command on the final numerical source/test state:

```bat
fpm test
```

Full-suite result: **passed, native exit code 0**, in **2325.10 s** (about
38 minutes 45 seconds), including all 365 repeated-knot checks. Existing
array-temporary runtime warnings were present; no test failed. Native
stdout/stderr and exit status were captured with Python `subprocess.run`
in ignored `build/repeated-knots/full-test.log` and `full-test-exit.txt`.

## Limits

Repeated supplied knots always use the existing dense limit: 256 physical
coefficients, 10,000 sites and 256 channels. Finite-range, normalized distinct-
knot gaps and final physical residual checks remain. Simple-interior banded
paths and automatic generation retain their prior contracts. No uniform
coefficient/derivative forward-accuracy bound is claimed, especially with
clustered knots. The objective does not penalize lower-order jumps. Legacy
tuple knot rules, C ABI and wrappers remain unchanged.
