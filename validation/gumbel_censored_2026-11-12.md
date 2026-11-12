# Censored Gumbel validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 150 new native Gumbel checks passed at both `-O0` and `-O3`, with
bounds/runtime checks and traps for invalid arithmetic, division by zero,
and overflow. All 945 high-precision reference comparisons passed in each
reference build; the maximum scaled log-likelihood error was `1.08e-12`
(an existing gamma stress case).

The Gumbel regression covers both skew directions, mixed censoring fits
compared with SciPy, fixed location/scale, objective agreement, preservation
of array likelihood/fit calls, reflection, ordinary and nested interval
masses, subnormal widths, underflowed standardized widths, adjacent intervals
in both tails, overflowing widths and displacements, the rounded subtraction
overflow threshold, the smallest scale at the location, unrepresentable log
densities/tails, infinite/equal endpoints, empty/invalid observations and
parameters, and fit error statuses. The prior private standardization helper
is renamed and shared without changing its numerical algorithm.

Checked native command (create `build\gumbel-native` first):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\gumbel-native -Ibuild\gumbel-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_gumbel.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\gumbel-native\test_censored_gumbel.exe
build\gumbel-native\test_censored_gumbel.exe
```

The optimized command replaces `-O0` with `-O3` and `gumbel-native` with
`gumbel-optimized`, retaining runtime checks and arithmetic traps.
FPM directory hashes are local artifacts; supply matching module and archive
paths for other builds. The modules listed above are compiled explicitly;
the archive supplies unchanged dependencies.

Reference commands:

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The harness uses seed 20261112 and 100 mixed cases for each of nine families,
plus one Weibull, three lognormal, nine gamma, nine beta, seven logistic, and
sixteen Gumbel stress cases (945 total). Gumbel stress cases comprise eight
right-skewed intervals and their left-skewed reflections: subnormal width,
scaled underflow, adjacent light and steep tails, a very steep tail near
standardized -700, overflowing width with finite scaled endpoints,
overflowing displacement with finite scaled value, and tiny scale.
Random references use 100 digits; beta, logistic, and Gumbel stress references
use 700 digits. Gumbel reference intervals use high-precision tail subtraction
independently of the exponential-increment identity used in Fortran.
Reference builds use `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fbacktrace` and `-O3 -fcheck=bounds`.

The mixed sample has exact values `[0.5,1.2,2.5]`, left limit `[0.3]`, right
limit `[3.5]`, and intervals `[1.1,1.4]` and `[4.0,4.4]`. At location 0.4 and
scale 1.3 the independent log-likelihoods are `-15.063330674681756` for
`gumbel_r` and `-36.48621545554803` for `gumbel_l`. SciPy fits give location
`1.187243219112296`, scale `1.4928122886183504` for `gumbel_r`, and location
`2.7284240787521883`, scale `1.7359414219851461` for `gumbel_l`.
Parameter comparisons allow `4e-4` for optimizer differences. The scaled
interval identities allow `1e-13` absolute error to account for subtraction
of log widths and log scales near magnitude 700; the observed tiny-scale
identity error was about `4e-15`.

## Full-suite verification

The full `fpm test` run passed with exit code 0, including the existing 373
censored checks and all 150 new Gumbel checks. The exact invocation used to
capture compiler/runtime output and preserve the FPM exit status was:

```cmd
cmd /d /c "fpm test > build\gumbel-full-test.log 2>&1"
```

Existing array-temporary runtime warnings were emitted; there were no test
failures. `python tools\check_test_inputs.py` passed. The final API diff,
whitespace, changed Fortran ASCII/line lengths, selected documentation links,
and file list were checked. No generated build artifact or unrelated file
is included in the change.

## Limits

Fits use the existing bounded local search. Convergence does not establish
a global optimum, identifiability, or existence of an unconstrained finite
MLE. No weights, new optimizer, dependency, or C ABI symbols are added.
Unrepresentable Gumbel log densities/tails return negative infinity; the
existing likelihood accumulation range limits still apply.
