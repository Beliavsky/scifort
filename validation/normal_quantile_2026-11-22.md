# Direct normal quantile validation and timing

Recorded on 2026-11-22 (America/New_York). The change replaces the private
standard-normal bisection with the attributed SciPy v1.15.3 / Cephes rational
kernel. Public Fortran signatures, C symbols, validation and endpoint behavior
are preserved. The CDF, SF and logarithmic-tail implementations are unchanged.
See [provenance](../CODE_PROVENANCE.md) and the
[retained license](../LICENSES/scipy-cephes-ndtri-BSD-3-Clause.txt).

## Source review

The exact upstream commit is `e29dcb65a2040f04819b426a04b60d44a8f69c04`.
All 47 coefficient literals were checked against the inspected source, with
unchanged ordering, branch cutoffs and normalization. Source SHA-256 values:

| Upstream file under `scipy/special/xsf/cephes/` | SHA-256 |
| --- | --- |
| `ndtri.h` | `4ec2945a1a6213d094d9f07b379be8186e45754b9bac87c7466ef1c81182f4e3` |
| `polevl.h` | `eff5938ec806f56284c593b445253c7b473f8faaaf9d6bc33d811aeb72ead099` |
| `const.h` | `15f2bb7181b75bc7c2cc14deecce98f048df77fe8e2f59334acc02f4ac93364e` |

The source files and BSD-3-Clause license were inspected before adaptation.
No library or preprocessor dependency was added. The Fortran catalog was
consulted; stdlib and fstats normal modules do not provide a direct inverse CDF.

## Focused accuracy validation

Commands from the repository root:

```text
python validation/check_normal_quantile.py --profile checked
python validation/check_normal_quantile.py --profile optimized
```

Both commands passed with GNU Fortran (GCC) 17.0.0, build 20260510 (experimental),
Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3 and mpmath 1.3.0 on Windows.
The checked profile uses `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`; the optimized profile uses `-O3`.
Both also use `-std=f2018 -fimplicit-none -Werror=implicit-interface`.

Each profile passed:

- 20,184 native checks, including 32 independent 400-digit reference quantiles,
  support/probability endpoints, invalid parameters, NaN/infinity handling,
  location/scale transforms, CDF/SF and extreme log-CDF round trips, strict
  dense-grid monotonicity and scalar/array agreement.
- 48,044 probability comparisons with `scipy.special.ndtri`, including a dense
  grid, logarithmic tails, random binary64 probability bit patterns, subnormals,
  the minimum positive binary64 value, neighbors of 0.5 and approximation
  switches, and the largest representable probability below one. ISF reflection
  and endpoint infinity signs were also checked.
- 242 independent 400-digit mpmath comparisons, evaluating
  `sqrt(2)*erfinv(2*p-1)` with the actual binary64 input probability.

The maximum difference from SciPy was **1 ulp**, and the maximum relative
error against the independent reference was **4.78727e-16**, in both profiles.
The inspected upstream `TestNdtri` median, endpoint and invalid-domain
contracts are covered by equivalent native assertions. No tolerances in
existing tests were loosened.

## Normal quantile timing

The same 16,384-probability workload and seven-sample benchmark protocol were
used before and after the change. Each native kernel was compiled with `-O3`.
Compilation and input setup are excluded. These focused measurements were made
before starting the full-suite rebuild. They are illustrative local timings,
not a portable performance guarantee or a CI threshold.

| Kernel | SciFort median (ms per array) | SciPy median (ms per array) | SciFort / SciPy |
| --- | ---: | ---: | ---: |
| Original bisection | 23.207450 | 0.399301 | 58.120 |
| Direct rational kernel | 0.290316 | 0.402771 | 0.721 |

The native time decreased by about **79.9 times**.
The measured SciPy time also varied between runs; neither ratio is an overall
library or language ranking. The original bisection averaged 58.8164 CDF
computations per quantile on the two benchmark input sets; the new kernel has
no iterative CDF evaluations.

To repeat the current benchmark:

```text
python benchmarks/compare_scipy.py --case normal_ppf_array --output build/normal-quantile-investigation/current-timing.json
```

The saved local before/after JSON reports are
`build/normal-quantile-investigation/timing.json` and
`build/normal-quantile-investigation/cephes-timing.json`. These build files are
ignored by Git. Their timing samples in measurement order, in ms per array:

| Kernel | Backend | Seven samples |
| --- | --- | --- |
| Original bisection | SciFort | 23.137875, 23.378825, 23.20745, 24.9092, 23.611575, 21.845625, 21.681475 |
| Original bisection | SciPy | 0.425091406, 0.399301172, 0.405751953, 0.408902344, 0.367463672, 0.356065625, 0.363612109 |
| Direct rational kernel | SciFort | 0.292274023, 0.306252539, 0.311169922, 0.278358008, 0.290316016, 0.281588672, 0.252345117 |
| Direct rational kernel | SciPy | 0.398312109, 0.393510547, 0.402770703, 0.403930469, 0.424424219, 0.408385156, 0.366069531 |

## Dependent regression suite

The normal quantile is shared by other distributions, random sampling,
special-function inverses and native/C ABI paths. The full `fpm test` command
was started after the numerical implementation and new native test were final.
It completed in 1322.81 seconds (22 minutes 2.81 seconds). Its sole failing
executable was the hypothesis regression described below; every other
executable passed, including the new native quantile regression and the C ABI,
random-sampling, distribution and special-function paths. The initial full
command returned exit code 1, so it is not recorded as a successful full run.

The run exposed a roundoff-sensitive Fligner statistic reference in
`test_hypothesis_extended.f90`: SciFort returned `10.803687663522233`, while
the stored SciPy result was `10.803687663522238`, outside the test's absolute
tolerance of `3e-15`. Independent 100-digit inverse-normal scores and statistic
arithmetic, using the actual binary64 rank probabilities, give
`10.8036876635222341424557164318327493235371542`.
The stored reference was corrected to this independently rounded value,
retaining the original tolerance. The calculation is reproduced by
`validation/check_normal_quantile.py`. The complete hypothesis test passed
with checked `-O0` and `-O3` test builds against the rebuilt library.

After correcting only this test reference, the dependency-backed FPM rerun
also passed with exit code 0:

```text
fpm test --target test_hypothesis_extended
```

The numerical implementation did not change after the full run started.
Unchanged suites were not rerun after the isolated reference correction.
The FPM version was `0.12.0, alpha`; FPM used the local GNU Fortran 17 compiler.

## Remaining limits

Validation here used GNU Fortran 17 on Windows. Other compiler/platform
combinations remain subject to CI validation. Binary64 rounding can give tied
or few-ulp differences for adjacent probabilities at approximation boundaries;
the dense-grid and branch-neighbor checks distinguish these regimes. This
change does not extend existing logarithmic inverse-CDF APIs or affine
range guarantees beyond their documented contracts.
