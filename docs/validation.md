# Starter validation record

Validation date: 2026-08-20

## Completed checks

The starter source was checked with GNU Fortran 14.2.0 on Linux using Fortran
2018 mode, warnings, implicit-interface errors, bounds and runtime checks, and
a backtrace-enabled debug build.

The following checks passed:

- strict compilation of every Fortran module;
- execution of `test/test_scifort.f90`;
- optimized compilation and execution of the distribution example;
- construction of a shared library;
- compilation and execution of a C program against `include/scifort.h`;
- exported C-symbol inspection;
- TOML parsing of `fpm.toml`;
- repository ASCII-policy check;
- selected comparison of 948 results with SciPy 1.17.0.

The sampled SciPy comparison covered PDF, log-PDF, CDF, survival,
log-CDF, log-survival, PPF, and inverse-survival values. The largest sampled
absolute difference for a quantile was approximately `7.11e-15`. The largest
sampled relative difference for a normal probability or log-probability was
approximately `1.03e-13` in a far tail.

These checks are smoke and regression checks, not a complete accuracy proof or
a final accuracy specification.

## FPM validation status

The execution environment used to create this archive did not contain an FPM
binary and did not permit downloading one. Therefore, the exact `fpm build`
and `fpm test` commands were not executed locally. The source dependency order
and equivalent GNU Fortran builds were exercised directly.

The repository includes a GitHub Actions workflow that installs FPM 0.13.0 and
runs debug build, tests, release build, and the example on Linux, macOS, and
Windows. Run the following after unpacking to complete local FPM validation:

```text
fpm build --profile debug
fpm test --profile debug
fpm build --profile release
fpm run --example distributions_demo
```

## Incomplete gamma, gamma, and chi-square validation

Validation date: 2026-09-19

Environment: Windows 11, GNU Fortran 17.0.0 20260510 (experimental), FPM
0.12.0 alpha, Python 3.13 with mpmath 1.3.0 and SciPy 1.15.3.

Commands run, all passing (`test_gamma` and `test_scifort`):

```text
fpm test
fpm test --profile debug
fpm test --profile release
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
fpm run --example distributions_demo
python tools/check_ascii.py
```

The trapping build produced no compiler warnings and no division-by-zero or
overflow traps. It found one division by an underflowed Newton slope in the
inverse solver, which was fixed before this record.

Additional sweeps compiled the modules directly with `gfortran -O2` and
compared them with mpmath at 50 digits. The sweep scripts were not committed.

- `log1p_safe` and `expm1_safe`, 5000 arguments on `[-0.99, 5]` and
  magnitudes down to `1e-20`: maximum relative errors `5.6e-16` and
  `2.2e-16`. The 0.1.0 versions had `5.0e-13` and `3.0e-12`.
- `gammainc`, `gammaincc`, and their logarithms, 3000 random `(a, x)` with
  `1e-8 <= a <= 1e6`, concentrated near `x = a` and spread over
  `1e-10 <= x <= 1e3`: with `V` the tail and `x f(x) / V` its condition
  number, the maximum of
  `|error| / (V (1 + |log V| + x f(x) / V))` was `6.4e-16`.
- `gammaincinv` and `gammainccinv`, 1500 random `(a, p)` with
  `1e-6 <= a <= 1e5` and `p` down to `1e-300`: the relative error in `x`
  implied by the mpmath residual, divided by `1 + V / (x f(x))`, was at most
  `4.1e-14`. The worst case was `p = 7.7e-300`, where the logarithmic
  residual limits accuracy to about `eps |log p| / a`. All 740 results equal
  to zero were confirmed to have exact roots below the smallest positive
  normal number.
- Against SciPy 1.15.3 on 20000 random `(a, x)`, the median relative
  difference was `2.3e-16` for `P` and `1.1e-16` for `Q`. The maximum was
  `2.0e-12`. In each of the eight largest differences checked with mpmath,
  SciPy had the larger error (up to `2.0e-12`) and SciFort's error was at
  most `3.7e-14`.

Known limitations:

- Evaluation cost near `x = a` grows like `sqrt(a)`. No uniform asymptotic
  expansion for large `a` is implemented yet.
- Tail probabilities below about `1e-300` carry relative error of order
  `eps |log V|` from forming the prefactor in logarithmic form.
- Inverse values in the subnormal range are returned with reduced precision,
  and roots below them are returned as zero.
- Only GNU Fortran was tested. Intel `ifx` and LLVM Flang were not run.
- The C ABI does not yet expose the new functions.
