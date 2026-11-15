# Periodic scalar short-basis correction, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python 3.13, SciPy **1.15.3**, binary64.

## Failure and correction

The 61-site reproducer uses `u=((i-1)/60)**1.1`,
`y=cos(2*pi*u)+0.08*sin(14*pi*u)` with the closing value equal to the first,
`w=0.8+u`, and `s=0.03*sum((w(:60)*y(:60))**2)`.
Before correction, degrees 3 and 5 return status 3; degrees 1, 2 and 4 succeed.
An instrumented copy of the pre-change scalar kernel reproduced this failure.

The failing fits have two independent cyclic coefficients (`n7=2<k`).
Both scalar solves extended them with `c(n7+1:n7+k)=c(1:k)`. Fortran evaluates
that overlapping right-hand side before assignment, reading stale entries
when more than one wrap is necessary. This destroys cyclic identities and
can make the residual appear inconsistent with the root bracket. Ordered
`c(n7+i)=c(i)` copying propagates each previously written coefficient. Existing
fixed-knot scalar and parametric paths already used ordered loops.

The correction changes neither QR factors, rank tests, knot selection rules,
root iteration nor tolerances. Selected knots can nevertheless change because
correct coefficients now feed interval residuals and continuation history.
No public API, new dependency or C ABI change is introduced.

## Native regression

`test_interpolate_periodic_short_basis` passes **542 checks** in checked
`-O0` and `-O3` builds. It tests degrees 1-5, fresh and continued fits at targets
`energy*[.03,.2,.001,10,0,.03]`, interpolation restart, independent one-channel
parametric fits, cyclic identities and periodic derivatives through k-1.
Cubic/quintic first fits explicitly assert the short two-coefficient basis.

Both profiles also pass the existing continuation (521), scalar smoothing
(108), parametric smoothing (125), periodic tuple (115), closed tuple (134),
repeated-knot (365), and large-selection (271) checks.

## Independent validation

`check_periodic_short_basis.py` uses seed 20261122, 31/61/97 nonuniform sites,
degrees 1-5, and five initial targets (.01, .03, .1, .2, .6 times energy), each
followed by targets .75 and 1.25 times its initial value. The 31/97 fixtures
include small random perturbations. Separate SciPy B-spline evaluation checks
residuals and endpoint derivatives. All 75 fresh fits compare knots and
coefficients against SciPy `splrep` on independently normalized inputs.

Both checked optimization profiles produced identical results: **225 fits**, including **81 short cyclic bases**; maximum
normalized coefficient difference **5.66091e-10**, active relative residual
error **0.000997233**, and scaled derivative seam error **5.40213e-16**.
Acceptance thresholds are 2e-13 for knots, 2e-7 for normalized coefficients,
0.0011 for active residuals, and 4096 binary64 epsilons for scaled seams.

The earlier 320-fit continuation comparison originally used exact upstream
`percur`, whose scalar kernel contains the same overlapping-copy defect.
`continuation_driver` now additionally supports explicit `reference-cyclic`
mode: periodic scalar references call unchanged upstream `clocur` with one
channel. Ordinary scalar, ordinary parametric and periodic parametric
references retain their appropriate unchanged upstream routines. Original
`reference` mode remains available to reproduce the defective scalar path.
No upstream source is patched. The closed interpolation-restart capacity bug
is still recorded before independently restarting unchanged upstream task 0.

Both checked optimization profiles produced identical results: **40 sequences / 320 fits**, **20 recorded upstream restart
failures** recovered independently; maximum normalized coefficient difference
**1.29856e-15**, value/derivative difference **2.2431e-12**, active relative
residual error **0.000989183**, and scaled seam error **6.47482e-15**. Thirty
scalar fits have knots differing from fresh SciPy fits, consistent with history.

The exact cached upstream source SHA-256 is
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`, commit
`1907c7bad9f5c4557970e906bdd4909eb52b2d4c`; BSD terms and attribution remain
in the repository. See [provenance](../CODE_PROVENANCE.md).

## Commands

Local ignored orchestration:

```bat
python build\periodic-reliability\checks.py debug
python build\periodic-reliability\checks.py release
python build\continuation\upstream.py debug
python build\continuation\upstream.py release
```

Fresh dependency order: `interpolate_curfit_core`, `interpolate_fixed_smoothing`,
`interpolate_periodic_lsq`, `interpolate_generate_knots`,
`interpolate_fitpack_compat`, `interpolate`. Both profiles used `-g -fbacktrace
-fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`, with
`-O0` or `-O3`, existing supporting modules and the FPM library archive. The
separately compiled exact upstream module/comparison driver omitted floating
point traps. Supporting archive modules were not rebuilt independently in
both focused profiles. Public validation commands:

```bat
python validation\check_periodic_short_basis.py build\continuation-debug\continuation_driver.exe
python validation\check_periodic_short_basis.py build\continuation-release\continuation_driver.exe
python validation\check_continuation.py build\continuation-debug\continuation_driver.exe
python validation\check_continuation.py build\continuation-release\continuation_driver.exe
fpm test
```

The full `fpm test` run passed with native exit code **0** in **1660.17 s**
(27 minutes 40 seconds), including all 542 new regression checks. Existing
array-temporary warnings appeared, with no failing test. Python subprocess
capture retained the native combined log, return code and elapsed time in
ignored `build/periodic-reliability/full-test.log`, `full-test-exit.txt` and
`full-test-seconds.txt`. Logs and compiled diagnostics stay in ignored `build`
storage.
This targeted correction does not establish general conditioning guarantees
or resolve every periodic fitting failure; broader audits remain work.
