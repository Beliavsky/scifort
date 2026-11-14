# Banded large-input selection validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Native coverage

Checked `-O0` and `-O3` builds both passed:

- **271** large-selection checks: 301, 1,001, 3,001 and 10,000 sites;
  degrees 1-5 at the first two sizes and degree 3 at the larger sizes;
  scalar/shared residuals, nonperiodic soft-count overshoots, nonperiodic
  and periodic refinement, canonical periodic knots, seams and returned sites.
- Balanced data/weight scaling by `1e100` / `1e-100` in a 1,001-site cyclic fit.
- Direct periodic banded recovery from known coefficients for degrees 1-5,
  plus explicit missing-column rank failures with no returned coefficients.
  A 10,000-site cubic problem recovered **9,999 independent coefficients**
  within `2e-11`, without a dense sample-by-coefficient matrix.
- Modern constructor **123**, scalar **112**, parametric **179**, periodic
  constructor **196**, general supplied-knot **206**, periodic supplied-knot
  **1,548**, periodic tuple LSQ **1,450**, generator **758**, closed tuple
  **134**, periodic tuple **115**, scalar smoothing **108**, parametric
  smoothing **125**, plus parametric interpolation, TCK and LSQ regressions.

The generator's oversize rejection fixture is now 10,001 sites. The legacy
constructor fallback fixture also uses 10,001 sites and loose targets giving
polynomial/constant fits; it verifies preserved dispatch, not arbitrary large
legacy smoothing accuracy. Active large modern refinement is tested separately.

## Independent large checks

`check_large_selection.py` uses seed **20261120**, nonuniform sites and
weights, one or three channels, and periodic/nonperiodic cases. It compares
complete positive-target sequences for 301 and 1,001 sites in degrees 1-5,
plus cubic 10,000-site cases. Ten additional zero-target 10,000-site
sequences cover both boundary conventions in every degree.

Nonperiodic sequences use SciPy 1.15.3's public generator. Periodic selection
uses the same hash-verified SciPy v1.18.0 control flow and canonical exterior
correction documented in [generator validation](generate_knots_2026-11-14.md).
Its independent residual evaluator uses SciPy banded QR for nonperiodic fits
and folded sparse B-spline matrices with LSQR for periodic fits. This is not
an unmodified SciPy 1.18 binary or an upstream suite run. Sparse normal
matrices are not formed; sparse design and derivative-jump operators check
the native shared-multiplier smoothing equations.

Both optimization builds produced these metrics:

| Check | Result |
| --- | --- |
| Complete large sequences | 74 |
| Individual selected-knot comparisons | 104,986 |
| Successful constructor fits | 64 |
| Feasible soft-count overshoots | 20 |
| Maximum relative residual-target error | 0.000890578 |
| Maximum scaled common-multiplier stationarity error | 2.91527e-13 |
| Maximum scaled unwrapped periodic seam difference | 3.00656e-14 |

Knots use `rtol=2e-14, atol=2e-14`; returned sites must match exactly.
Active residual errors must be below `0.0011`; stationarity below `2e-8`;
seams below `4096*epsilon`. Values and all derivatives through degree `k`
are checked against the returned B-spline representation, including queries
outside the base interval. Soft fixtures are generated from first-refinement
basis functions so their improvement over the minimal polynomial/constant
is substantial rather than hidden within the residual stopping tolerance.

The earlier scalar, parametric, periodic and combined independent constructor
checkers also passed both builds with the printed results in
[modern selection validation](modern_selection_2026-11-14.md), including
zero modern nonperiodic knot differences. Small-input dense behavior remains
unchanged. Root stopping still permits coefficient/high-derivative differences
from SciPy; matching knots alone does not establish roundoff-level parity.

## Storage and elapsed time

For 10,000 sites, cubic degree and one channel, the periodic QR helper's
explicit arrays plus returned coefficients total **2.481 MiB**, assuming
binary64 and four-byte default integers. An equivalent full dense design
alone would be **763.092 MiB**. The checker derives the buffer total as

```text
8*(m*(channels+k+2) + nt*(5*k+8+2*channels)
   + (nt-k-1)*channels) + 4*nt
```

with `nt=m+2*k`. This excludes caller data, expression temporaries, allocator
and runtime overhead and is **not a measurement of process RSS**. The
9,999-coefficient native test exercises the banded implementation at this
scale. Dense fallback is capped separately at 256 independent coefficients
and one million active independent-design entries.

Observed generator subprocess time over the 64 positive-target cases was
**4.969 s** in the checked debug run and **6.159 s** in the checked optimized
run. The batched constructor subprocess took **1.796 s** and **1.808 s**,
respectively. These include process launch and text I/O, exclude Python
oracle work, and are observations from single runs rather than a benchmark
or guaranteed speedup. The ten zero-target vectors are outside those timings.

## Commands

The local focused runners used:

```bat
python build\large-selection\checks.py debug
python build\large-selection\dependencies.py debug
python build\large-selection\checks.py release
python build\large-selection\dependencies.py release
```

The first debug runner stopped at an overstrict soft-fixture assertion in
the new Python checker after all its native and existing independent checks
passed. After correcting that fixture, the checker was rerun directly and
passed. No production numerical change was needed for that validation fix.
The final optimized runner completed in full.

These local ignored orchestration scripts compiled the following source
files in dependency order:

```text
src/interpolate_curfit_core.f90
src/interpolate_fixed_smoothing.f90
src/interpolate_periodic_lsq.f90
src/interpolate_generate_knots.f90
src/interpolate_fitpack_compat.f90
src/interpolate.f90
```

For example, the debug kernel compile was:

```bat
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\large-selection-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\large-selection-debug -c src\interpolate_curfit_core.f90 -o build\large-selection-debug\interpolate_curfit_core.o
```

Tests and drivers used those flags and fresh objects followed by
`build\gfortran_3C159EB802C22FB9\scifort\libscifort.a`. The optimized build
used `-O3` and `build\large-selection-release`, keeping checks/traps enabled.
Supporting archive objects were not independently rebuilt in both profiles.
Independent commands for each build included:

```bat
python validation\check_make_splrep.py build\large-selection-debug\make_splrep_driver.exe
python validation\check_make_splprep.py build\large-selection-debug\make_splprep_driver.exe
python validation\check_periodic_constructors.py build\large-selection-debug\periodic_constructors_driver.exe
python validation\check_modern_selection.py build\large-selection-debug\modern_selection_driver.exe
python validation\check_large_selection.py build\large-selection-debug\generate_knots_driver.exe build\large-selection-debug\modern_selection_driver.exe
```

The optimized commands used corresponding `large-selection-release` paths.
Native executables were run directly after linking.

Full-suite command (shared kernel and final numerical state):

```bat
fpm test
```

Full-suite result: **passed, native exit code 0**, in **903.34 s**. The run
included all 271 large-selection checks and all 123 modern-constructor checks.
Native stdout/stderr and exit status were captured with Python
`subprocess.run`, avoiding shell stderr-formatting artifacts. Existing
array-temporary runtime warnings were present; no test failed. The raw log
and exit status are in ignored `build/large-selection/full-test.log` and
`full-test-exit.txt`.

## Remaining limits

Modern generation and soft counts are bounded at 10,000 sites and 256 channels.
Larger constructor inputs retain legacy selection and hard capacities. Dense
rank fallback is bounded; arbitrary larger deficient fits are not supported.
Repeated interiors, broader conditioning guarantees, exact periodic selection
and root-iteration parity, C ABI and wrappers remain separate gaps.
