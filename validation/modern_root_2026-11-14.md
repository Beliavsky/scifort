# Modern supplied-knot smoothing stop validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python 3.13, SciPy **1.15.3**, binary64.

## Finding and correction

With identical supplied knots, native simple-knot fits previously returned
unpenalized least-squares coefficients whenever their residual was inside
`0.001*s`, even if a feasible target was strictly above that residual.
The inspected modern reference instead evaluates finite penalties before
accepting a residual. These valid but different stops can change coefficients
and high derivatives. Automatic modern fits also use the supplied-knot solve
after knot selection.

The four existing attributed scalar/parametric ordinary/periodic kernels now
skip that early exit only when `fixed` and `fp-s<0`. They enter their existing
rational search. A least-squares residual at or immediately above the target
can still be accepted inside the existing tolerance. This preserves the
near-feasibility allowance and avoids claiming a finite exact solution when
the target is below the minimum residual. Polynomial/constant limits, rank
checks, QR, root safeguarding, tolerances and legacy tuple tasks are unchanged.
General nonclamped/repeated-knot/rank-deficient fallback searches are unchanged.

No new numerical solver, public argument, dependency or C ABI is introduced.
The preliminary audit also encountered warnings and unmet targets in broader
reference cases; these outputs are not adopted to manufacture parity. In
particular, the native rational search retains its existing finite-upper-
bracket safeguard rather than copying different reference branch behavior.

## Regressions

`test_interpolate_modern_root` passes **180 checks** in both checked `-O0` and `-O3` builds.
Deterministic 61-point fixtures cover degrees 1-5, scalar/three-channel fits,
ordinary/canonical periodic knots and targets just above the least-squares
residual. Fixed numerical references check selected coefficients and the
highest derivative. Nearby targets just below the least-squares residual
verify that the existing tolerance allowance still succeeds. Shared fits
also verify retained supplied parameters.

The pre-change static archive fails the new `finite-penalty coefficients`
assertion. With the correction all new assertions pass. This negative check
was performed before the full FPM run updated the archive; its executable
and log remain only in ignored build storage.

Existing tests also pass in both checked optimization profiles: modern selection (123), scalar
constructor (112), parametric constructor (179), periodic constructor (196),
general fixed knots (206), supplied periodic knots (1,548), repeated knots
(365), large selection (271), continuation (521), and short cyclic bases (542).

## Independent root and linear-solve comparison

`validation/check_modern_root.py` uses seed **20261123**, 31/61/97 nonuniform
sites, degrees 1-5, ordinary/periodic boundaries, scalar and one/three-channel
parametric calls. Deterministic trigonometric inputs gain small random
perturbations in the 31/97-site cases. All use identical full supplied knots;
targets are the independently computed least-squares residual times
`[1.00001,1.0005,1.002]`.

Reference root control comes from the exact cached SciPy **v1.18.0**
`root_rati`, `fprati` and `Bunch`; packed jumps come from its `disc`/`prodd`.
The source is hash verified before extracting definitions, and is not patched.
Every penalty evaluation uses independently constructed dense NumPy augmented
least-squares equations. Canonical periodic coefficients are folded by cyclic
column identities. Initial scales use positive dense QR diagonals of the
weighted design; the polynomial/constant and least-squares residual brackets
are computed independently. The helper checks reference convergence, exact
returned knots/parameters, coefficients, all derivatives through degree k,
native evaluation/integration and physical residuals. Known final reference
penalties also supply an independent stationarity check.

The 135 ordinary cases additionally compare directly against public installed
SciPy 1.15.3 `make_splrep` and `make_splprep`, avoiding dependence on the custom
linear-solve oracle alone. Periodic references are the stated independent
linear solves with pinned control flow, not an installed SciPy 1.18 binary or
an unmodified upstream C++ solver suite.

Both checked optimization profiles produced the same printed results:

| Check | Result |
| --- | ---: |
| Supplied-knot fits | 270 |
| Public SciPy comparisons | 135 |
| Independent penalty evaluations | 1,057 |
| Maximum normalized coefficient error | 1.91138e-14 |
| Maximum normalized value/derivative error | 4.73461e-11 |
| Maximum public SciPy coefficient error | 1.89028e-14 |
| Maximum relative residual error | 0.000995333 |
| Maximum scaled stationarity error | 6.50057e-16 |

Acceptance bounds: coefficients `3e-11`, value/derivative comparisons `3e-8`,
active residual `0.0011*s`, stationarity `2e-11`; representation and integration
checks have separate `3e-10`/`3e-9` absolute and relative allowances. Knot and
supplied-parameter retention are exact. These are fixture-specific tolerances,
not universal forward-error guarantees.

## Existing independent constructor checks

Both checked optimization profiles also pass with the fresh numerical source:

| Checker | Result |
| --- | --- |
| Scalar constructor | 165 fits; 26,737 strict SciPy comparisons; 55 residual-aware automatic fits |
| Scalar maximum reported normalized difference | 1.28873e-11, down from 1.48599e-7 before this correction |
| Parametric constructor | 165 fits; 51,779 strict comparisons; 80 residual-aware automatic fits |
| Periodic constructor | 200 fits; 23,020 strict interpolation comparisons; 100 active constraints |
| Modern selection | 100 successful fits; 20 infeasible targets; 20 feasible soft-count overshoots; 2,685 knots compared |

Selected knots are unchanged in these fixtures. The parametric checker still
reports a raw overall difference of **1.37358** across derivatives and query
regimes; the periodic legacy diagnostic still includes different selected
knots. Those broader differences are not solved or hidden by the new
near-least-squares correction. See the earlier
[selection record](modern_selection_2026-11-14.md).

## Source and commands

Inspected cached source: SciPy v1.18.0 `_fitpack_repro.py`, SHA-256
`82449ae292bb361675439b422cec6c084f9fbfa6519eabec06b1009c51a818bd`.
The full BSD license remains in `LICENSES/scipy-generate-knots-BSD-3-Clause.txt`.
The existing perazz-derived kernels retain their notices and exact source
attribution. See [provenance](../CODE_PROVENANCE.md).

Ignored local orchestration:

```bat
python build\modern-root\build.py debug
python build\modern-root\tests.py debug
python build\modern-root\independent.py debug
python build\modern-root\build.py release
python build\modern-root\release.py
```

Sources were freshly compiled in order: `interpolate_curfit_core`,
`interpolate_fixed_smoothing`, `interpolate_periodic_lsq`,
`interpolate_generate_knots`, `interpolate_fitpack_compat`, `interpolate`.
Debug flags were `-O0 -g -fbacktrace -fcheck=all
-ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`; optimized flags
replace `-O0` with `-O3`, keeping checks/traps. Includes/module outputs use
`build/modern-root-debug` or `build/modern-root-release` and existing
`build/gfortran_2806A042DAA8BBF7` support modules. Fresh objects precede
`build/gfortran_3C159EB802C22FB9/scifort/libscifort.a` when linking. Supporting
archive modules were not independently rebuilt at both optimization levels.
The driver is the existing `validation/fixed_general_driver.f90`.

```bat
python validation\check_modern_root.py build\modern-root-debug\fixed_general_driver.exe
python validation\check_modern_root.py build\modern-root-release\fixed_general_driver.exe
fpm test
```

The full `fpm test` run passed with native exit code **0** in **2756.72 s**
(45 minutes 57 seconds), including all 180 new regression checks. Existing
array-temporary warnings appeared; no test failed. Python `subprocess.run`
captured combined native output, exit code and elapsed time in ignored
`build/modern-root/full-test.log`, `full-test-exit.txt` and
`full-test-seconds.txt`. This full run covers the final numerical source and
tests; subsequent edits only record the documentation results.

Source ASCII/free-form line widths, Python syntax, local documentation links
and `git diff --check` also passed.
General conditioning, fallback root parity, exact periodic knot selection,
larger modern limits and C ABI/wrappers remain outside this correction.
