# Censored finite-sample two-sided Kolmogorov-Smirnov validation

Date: 2026-11-18.

The new censored kstwo kernel is independent of the existing scalar finite-
difference density. Its contract is integer `1<=n<=32`, normal positive
scale<=huge/2, abs(loc)/scale<=100, and finite noncollapsed support. See the
[API guide](../docs/censored_fitting.md#finite-sample-two-sided-kolmogorov-smirnov-censoring)
and [provenance](../CODE_PROVENANCE.md).

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Python 3.13.3, SciPy 1.15.3, mpmath 1.3.0.

## Native checks

Both `-O0` and `-O3` builds passed 718 checks with `-g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`. Each build compiled:

```text
gfortran <optimization> -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J<out> -I<out> -Ibuild/gfortran_2806A042DAA8BBF7 src/empirical_distribution.f90 src/censored.f90 src/fit.f90 src/stats.f90 test/test_censored_kstwo.f90 build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o <out>/test.exe
<out>/test.exe
```

Here `<optimization>` is `-O0` / `-O3` and `<out>` is
`build/kstwo-native` / `build/kstwo-optimized`, respectively.
Checks cover scalar CDF/SF and quantile consistency, affine scales 1e-200 and
1e200, first-knot one-sided limits and affine crossing intervals, support and
probability endpoints, adjacent intervals, invalid/NaN/infinite inputs,
NNLF sign, fixed and free-scale fits, enumeration of integer n, and preserved
real-array fit calls.

## Independent reference checks

```cmd
C:\Programs\Python313\python.exe validation/check_censored_kstwo.py
```

The script constructs independent 100-decimal nested polynomial integrals over
ordered uniforms and obtains central density derivatives with mpmath. It does
not repeat the Fortran count-convolution or density-flux algorithms. Reference
cases include n=1,2,3,4,10,32, support-adjacent tails, density knots, central and
adjacent intervals, and affine scales through 1e307. SciPy's finite-difference
central PDF is deliberately not treated as a tight density oracle.

Both checked and optimized builds passed all 241 independent reference cases.
The maximum scaled log error was 2.19e-14; including relative errors for tiny
complementary logs, the maximum comparison error was 5.62e-14, below 3e-11.
The script prints the exact compiler commands it uses. Optional `--reuse`
reuses the ignored local reference fixtures while rebuilding both profiles.

For the mixed sample in the script and native test, SciPy evaluated the NNLF
at scale 1.2 as 9.696556989999241; SciFort's fit finds scale 1.2 and NNLF
9.696556989998978. This point puts an exact observation at the first density
knot. SciPy's smooth bounded scalar optimizer instead finds a local solution
near scale 1.1318853116319016, NNLF 9.840680042272329. That local solution is
not used as the global fit oracle: the knot has a strictly better objective.
The native fit checks compare against the independently evaluated knot value.

## Full suite

`fpm test` passed with exit code 0 in 1593.967 seconds (26 minutes 34
seconds), including all 718 new kstwo checks and the existing scalar kstwo
regressions. The run emitted existing array-temporary runtime warnings from
scattered interpolation. No full debug-suite run was performed; the focused
native and independent-reference builds above ran with bounds/runtime checks
and invalid/zero/overflow traps at both optimization levels.

The final source/docs audit passed: ASCII and standard line lengths in changed
Fortran, Python syntax, changed relative documentation links, 107 distinct
supported families plus 2 remaining, matching documented likelihood names,
and `git diff --check`. No generated reference fixtures or build artifacts
are included in the change.
