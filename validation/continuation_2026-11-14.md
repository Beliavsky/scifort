# Explicit spline continuation validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Native tests

Checked `-O0` and `-O3` builds passed **521 checks** in
`test_interpolate_continuation`:

- Degrees 1-5, ordinary/periodic scalar and three-channel fits, with separate
  interleaved state objects. Sequences go from polynomial/constant limits to
  active targets, looser targets, interpolation and positive-target restart.
- Weighted shared residuals, retained parameters, periodic endpoint
  derivatives, scalar padding, independent deterministic deep copies.
- Default/reset state lifecycle, missing state, fixed-knot/state rejection,
  changed active data, weights, sites, degree, periodic mode and fitting family.
- NaN targets and unrepresentably small normalized targets return invalid
  tuples without losing successful history; subsequent valid reuse succeeds.
- Ignored periodic closing data/weights may change to NaN. Generated chord
  parameters, replacement by a fresh fit of a different family, enlarged
  ordinary bounds, mismatched bounds and rejected fresh initialization.

Both builds also passed existing scalar smoothing **108**, parametric
smoothing **125**, periodic tuple **115**, closed parametric tuple **134**,
repeated supplied knots **365**, and large modern selection **271** checks.
The final seven lifecycle assertions were rebuilt/run separately in each
profile after the initial 514-check run. The numerical source was unchanged
for that test expansion.

## Exact upstream continuation comparison

`validation/continuation_driver.f90` links both native SciFort and the exact
unchanged cached perazz/FITPACK module. Source:

- Commit: `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.
- [Source URL](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90).
- SHA-256: `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
- BSD-3-Clause license retained in `LICENSES/perazz-fitpack-BSD-3-Clause.txt`.

The upstream module was compiled directly from that source, not generated
from native kernels. It is cached only under ignored build storage. Native
and reference runs use separate processes. Each reference sequence retains
its full upstream real/integer work arrays, knot vector and knot count and
calls genuine `curfit`, `parcur`, `percur` or `clocur` task 1. Inputs and
residual targets are independently normalized to the same binary64 scale
used by the native wrapper. No upstream source patch is made.

`check_continuation.py` uses seed **20261121**, 47/91 nonuniform sites,
random positive weights, degrees 1-5, one/three independent channels and
ordinary/periodic boundaries. Each of 40 sequences has eight targets:
`energy * [10, .0005, .001, .002, .0002, 10, 0, .0005]`, where energy is
the total active weighted squared data norm. Selected knots and coefficients
are compared with the stateful upstream reference; SciPy B-splines
independently check residuals, derivatives through k, seams and extrapolation.

Both checked optimization profiles produced:

| Check | Result |
| --- | ---: |
| Complete sequences | 40 |
| Successful fits | 320 |
| Maximum normalized coefficient difference | 0 |
| Maximum normalized value/derivative difference | 0 |
| Maximum active relative residual error | 0.000989183 |
| Maximum scaled periodic seam difference | 6.47482e-15 |
| Scalar fits with knots differing from a fresh SciPy fit | 34 |
| Recorded upstream closed restart failures independently recovered | 10 |

Knot comparisons use `rtol=atol=2e-13`, normalized coefficient error `2e-7`,
value/derivative error `2e-6`, active residual error `0.0011`, and seam error
`4096*epsilon`. Interpolation residuals have explicit absolute roundoff
allowances. Zero printed coefficient/derivative differences are exact for
these normalized fixtures; they are not a universal parity guarantee.

The 34 fresh-fit knot differences demonstrate why a fresh SciPy task-0 fit
alone cannot validate continuation: retained knots depend on history.
No general speedup or QR-factor reuse claim is made.

## Closed-curve restart correction

All ten periodic parametric sequences encounter an upstream positive error
when restarting after interpolation. The exact upstream `fpclos` bootstrap
compares the prior `n` (already equal to full capacity) with `nest` before
replacing it by the minimal one-interior seed. Native code compares `nmin`
with `nest` at that point. The validation driver records the original error,
then runs a fresh unchanged upstream task-0 fit for that transition. The
corrected native knots and coefficients match that independent fit exactly.
Subsequent upstream task-1 calls continue from the recovered reference state.
This is a recorded intentional correction, not hidden upstream success.

## Subsequent correction

The later [short cyclic basis correction](periodic_short_basis_2026-11-14.md)
fixes the scalar failure below. The historical results in this record describe
the pre-correction code and original upstream scalar reference. Current
`check_continuation.py` explicitly uses unchanged one-channel `clocur` for
periodic scalar references because upstream `percur` contains the same defect.

## Existing numerical limits

The extension preserves existing normalization, rank and final physical-
residual guards. A diagnostic against the pre-change native archive confirmed
an existing periodic scalar status-3 case, independent of continuation:
61 sites `u=((i-1)/60)**1.1`,
`y=cos(2*pi*u)+0.08*sin(14*pi*u)` with the closing value set to the first,
`w=0.8+u`, and target `0.03*sum((w(:60)*y(:60))**2)` fails for degrees 3/5
in that build and the changed wrapper's fresh-fit path. The broader periodic
numerical reliability audit remains a gap; these failures are not silently
replaced by a different algorithm. The successful sequence fixtures above
exercise richer active fits and exact upstream state transitions.

## Commands and flags

Local ignored focused orchestration scripts:

```bat
python build\continuation\checks.py debug
python build\continuation\checks.py release
python build\continuation\one.py debug
python build\continuation\one.py release
python build\continuation\upstream.py debug
python build\continuation\upstream.py release
```

The fresh native source dependency order was:

```text
src/interpolate_curfit_core.f90
src/interpolate_fixed_smoothing.f90
src/interpolate_periodic_lsq.f90
src/interpolate_generate_knots.f90
src/interpolate_fitpack_compat.f90
src/interpolate.f90
```

For example:

```bat
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\continuation-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\continuation-debug -c src\interpolate_curfit_core.f90 -o build\continuation-debug\interpolate_curfit_core.o
```

Tests linked those fresh objects followed by
`build\gfortran_3C159EB802C22FB9\scifort\libscifort.a`. Optimized runs used
`-O3` and `build\continuation-release`, retaining bounds checks and traps.
Supporting archive modules were not independently rebuilt in both profiles.
The exact upstream module and comparison driver used the same optimization,
bounds/backtrace/interface checks, without the floating-point trap flag.
The `.F90` upstream file uses its own existing preprocessor directives.

Final independent commands:

```bat
python validation\check_continuation.py build\continuation-debug\continuation_driver.exe
python validation\check_continuation.py build\continuation-release\continuation_driver.exe
```

The full-suite command covers the final numerical source/test state:

```bat
fpm test
```

Full-suite result: **passed, native exit code 0**, in **1594.22 s** (about
26 minutes 34 seconds), including all 521 continuation checks. Existing
array-temporary runtime warnings appeared; no test failed. Python
`subprocess.run` captured native stdout, stderr and exit status under ignored
`build/continuation/full-test.log` and `full-test-exit.txt`, avoiding PowerShell
stderr-formatting artifacts.

The complete example in the continuation guide also compiled and ran with
checked `-O0` flags. Changed-source ASCII, free-form line widths, Python
syntax, local documentation links and `git diff --check` passed.

## API and portability limits

Explicit state is required for task 1; no global Python-style cache is
introduced. Matching uses exact active-input comparisons, not hashes or
approximate equality. State copies and call transactions allocate/copy input
and knot history; state is not serialized or C interoperable. QR factors and
smoothing iterations are recomputed. Ordinary zero-target bounds must equal
the sites' endpoints. General ill-conditioned forward accuracy, broader
periodic numerical reliability, C ABI and wrappers remain separate work.
