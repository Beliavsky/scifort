# Adaptive knot generator validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Focused results

Both checked `-O0` and `-O3` builds passed **758 native checks** covering
degrees 1-5, exact interpolation recovery, interpolation knot conventions,
minimal sample counts, weighted residual monotonicity and stopping,
independent copied states, finite exhaustion, periodic shared-channel knot
cycles, constant limits, soft counts, ignored finite closing weights, invalid
shapes/parameters/nonfinite inputs, range/size guards and collapsed midpoint
representation failures. Default and failed generators return empty outputs.

The independent checker, seed **20261116**, passed in both builds:

| Check | Result in each build |
| --- | ---: |
| Base scalar/shared-channel sequence cases | 140 |
| Balanced data/weight scaling cases | 140 |
| Site scaling cases (`1e-100`, `1e100`) | 40 |
| Direct public SciPy nonperiodic sequence comparisons | 60 |
| Periodic base cases | 70 |
| Individual knots compared | 18,078 |
| Maximum absolute knot discrepancy | 2.22045e-16 |
| Independent residual checks | 837 |
| Stale upstream capacity exteriors canonicalized | 9 |

Base cases cover degrees 1-5, one/three channels, nonuniform random sites and
weights, zero targets, loose constant/polynomial targets, intermediate targets,
near-interpolation targets, capacity stops and enlarged nonperiodic bounds.
Knot comparisons use `rtol=2e-14, atol=2e-14`; balanced scaling requires exact
vector equality. Residual checks use independent NumPy least squares,
allowing `1e-10` absolute roundoff for monotonicity and checking uncapped final
residuals against `1.001*s+1e-20`.

The public nonperiodic oracle is installed SciPy 1.15.3. The periodic oracle
executes exact hash-checked SciPy v1.18.0 Python selection control flow with
NumPy cyclic least squares substituted for SciPy's internal QR. Its residual
evaluator canonicalizes periodic exteriors before every solve, and the final
capacity-limited vector is canonicalized before comparison. This deliberately
models SciFort's corrected cyclic contract; it is not an unmodified SciPy
1.18 periodic result. The upstream loop can leave those exteriors stale
within batches and at a capacity stop. No upstream SciPy test-suite run is
claimed. Source versions, full hashes and license are in
[provenance](../CODE_PROVENANCE.md).

## Commands

Focused tests rebuilt the new generator and the public interpolation facade,
then linked the preexisting native support archive. The existing support
modules were not themselves rebuilt with both optimization settings.
The debug commands issued by the local runner were equivalent to these
CMD-compatible commands (cache paths record this checkout):

```bat
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\generate-knots-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\generate-knots-debug -c src\interpolate_generate_knots.f90 -o build\generate-knots-debug\interpolate_generate_knots.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\generate-knots-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\generate-knots-debug -c src\interpolate.f90 -o build\generate-knots-debug\interpolate.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\generate-knots-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\generate-knots-debug test\test_interpolate_generate_knots.f90 build\generate-knots-debug\interpolate_generate_knots.o build\generate-knots-debug\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\generate-knots-debug\test_interpolate_generate_knots.exe
build\generate-knots-debug\test_interpolate_generate_knots.exe
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\generate-knots-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\generate-knots-debug validation\generate_knots_driver.f90 build\generate-knots-debug\interpolate_generate_knots.o build\generate-knots-debug\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\generate-knots-debug\generate_knots_driver.exe
python validation\check_generate_knots.py build\generate-knots-debug\generate_knots_driver.exe
```

The optimized run used the same commands with `-O3` and
`build\generate-knots-release` in place of `-O0` and
`build\generate-knots-debug`, retaining all runtime and floating-point checks.
The new guide example was also compiled and run with the checked debug build.

Full-suite command:

```bat
fpm test
```

Full-suite result: **passed, native exit code 0**, including the 758 new native
checks and the existing interpolation, distribution, optimization, signal,
sparse and special-function tests. The invocation captured native stdout/stderr
with Python `subprocess.run` and recorded FPM's own exit code; shell stderr
formatting did not decide test success. Existing array-temporary runtime
warnings were present; no test failure was reported.

## Remaining scope

This addition exports a generator; it does not replace the automatic paths or
hard capacities of `make_splrep` / `make_splprep`. Periodic path parity,
zero weights, repeated sites, degrees outside 1-5 and larger dense problems
are not established. Range guards do not guarantee good conditioning. C ABI
and wrappers remain deferred. See the [API guide](../docs/interpolate_generate_knots.md).
