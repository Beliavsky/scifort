# Optimized-build regressions, 2026-11-19

## Modified Fresnel phase accuracy

The previous full release run failed at the exact binary64 input
`x=100000000.25`. A standalone probe isolated the problem to the combined
`sincos` runtime call emitted by GNU Fortran at `-O3`. For the rounded square
`10000000050000000`, separate `sin`/`cos` calls at `-O0` and the combined call
at `-O3` disagreed by about `1.28e-5` in the sine. The split-square correction
was `0.0625` in both builds; the compensated multiplication was not the cause.
Assembly inspection confirmed the optimized probe called `sincos`.

`square_phase` now evaluates cosine using `1-2*sin(t/2)**2` for both the rounded
square and its correction. Distinct sine arguments prevent this combined-call
substitution while retaining small absolute phase errors. The correction also
needs this treatment: for larger inputs, its magnitude can itself be large.
No compiler flags, public interfaces, or reference tolerances were relaxed.

The mathematical references use the rotated integral
`F_+(x)=sqrt(pi)/2*exp(i*pi/4)*erfc(x*exp(-i*pi/4))`, with the documented
scaling for `K_+` and conjugation for the negative pair. These follow the
[public integral definitions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modfresnelp.html)
and [DLMF error/Fresnel interrelations](https://dlmf.nist.gov/7.5).
No upstream implementation code was copied or translated.

The new checker uses mpmath 1.3.0 at 400 decimal digits, converting binary64
inputs exactly. Its 159 inputs cover zero, ordinary arguments, adjacent large
arguments, deterministic logarithmically spaced random magnitudes, both signs,
and magnitudes up to half the square-overflow threshold. It compares all four
complex outputs using relative error, including tiny tails.

## Circular-statistics reference

The former negative-period standard-deviation test used antipodal observations.
Their exact mean resultant is zero, so the finite standard deviation reference
depended on rounding of trigonometric functions. The replacement uses a
quarter-period separation, with exact resultant magnitude `sqrt(1/2)`, mean
`1.25`, and signed standard deviation `-sqrt(log(2))/pi`. An antipodal variance
check still verifies the well-conditioned value one. Library circular-statistics
behavior is unchanged.

## Error-path calls

Ten child-process test programs assigned pure-function results to unused local
variables. Optimized executables could return normally without executing those
calls. Each tested result is now printed if the call unexpectedly returns.
The parent still fails when the child exits normally; no unconditional error
termination was added after the calls. Histogram tests consume a scalar metadata
field rather than dereferencing possibly unallocated result arrays.

The affected programs are `test_descriptive_circular_errors`,
`test_descriptive_entropy_errors`, `test_descriptive_frequency_errors`,
`test_descriptive_kstat_errors`, `test_descriptive_lmoment_errors`,
`test_descriptive_nan_raise`, `test_descriptive_proportion_trim_errors`,
`test_power_transform_errors`, `test_sample_distance_errors`, and
`test_sample_distance_nd_errors`.

## Validation

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Python: 3.13.3; mpmath: 1.3.0.

From the repository root:

```bat
python validation/check_modified_fresnel_phase.py --modules build/gfortran_100AF6AD77A3D894 --library build/gfortran_F2AB82903052A58A/scifort/libscifort.a
fpm test --profile release
```

The module/archive paths identify the existing release cache on this machine;
substitute the corresponding paths after a fresh FPM build. The checker prints
the compiler version and rebuilds the changed Faddeeva module in both profiles:

- Checked: `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
- Release: `-O3 -funroll-loops`.

Both profiles passed four native suites: modified Fresnel, Faddeeva, complex
error functions, and complex Fresnel. Both passed 636 independent complex
comparisons; maximum relative error was `8.015731e-16` in each profile.
The four additional native large-phase references, generated from the same
400-digit formula, scale their tails before comparison to avoid underflow.

The eleven circular/error-path programs also passed when compiled and linked
individually against the existing release archive. Each was compiled once with
`-O0 -g -fcheck=all -fbacktrace` and once with `-O3 -funroll-loops`, using
`-I build/gfortran_100AF6AD77A3D894` and
`build/gfortran_F2AB82903052A58A/scifort/libscifort.a`. Runtime checks here apply
to the test callers; the unchanged library dependencies retain release flags.

Full `fpm test --profile release`: **passed**, exit status zero. Elapsed time
was 2298.938 seconds (38 minutes 19 seconds), including recompilation and
relinking. All twelve previously failing executables passed in this run.
The final diff check, changed-source ASCII/line-length checks, and Python
syntax check also passed. No numerical tests were rerun solely for the final
documentation update.

These checks establish the observed fix for this toolchain and the sampled
arguments; they are not an exhaustive accuracy proof for every binary64 input
or a test of other compilers' mathematical runtimes.
