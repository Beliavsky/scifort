# Levy-stable censored skewness extension

Date: 2026-11-19.

## Scope and derivation

Extend the general native S1 stable censored likelihood, NNLF and bounded-fit
shape rectangle to 1.25<=alpha<=1.75 and abs(beta)<=0.75. The special cases,
affine restrictions, public names and existing scalar distribution API remain
unchanged. All finite physical coordinates are supported in the censored
rectangle, including overflow of the standardized coordinate.

The formulas, quadrature panels, join at abs(z)=6 and exprel polynomial are
unchanged from the [tail derivation](stable_censored_tails_2026-11-19.md).
For a=alpha, c=1-i*beta*tan(pi*a/2), theta=pi/(8*a), d=c*exp(-i*a*theta),
the extended rectangle gives

```text
Re(d) >= cos(pi/8)-0.75*abs(tan(5*pi/8))*sin(pi/8)
      = 0.25*cos(pi/8) = 0.2309698831278216... > 0.
```

The real coefficient remains positive throughout the rotation as well.
The additional damping exp(-sin(theta)*y) is still bounded above by
exp(-sin(pi/14)*y). Since exprel(w)=integral(exp(t*w),t=0..1), Re(w)<=0
implies abs(exprel(w))<=1. Thus the same polynomial-times-exponential
absolute quadrature-tail bound applies, with the larger abs(c) prefactor.
Neither the integration cutoff nor the cancellation-free tail construction
requires a change. With abs(c)<=sqrt(1+(0.75*tan(5*pi/8))**2),
the factored PDF remainder beyond y=256 is bounded by 9.11e-21 and the
factored SF remainder by 3.50e-23 (integrals of y**1.75 and y**0.75
times that exponential, respectively, divided by pi).
This argument establishes damping; independent references
below check the actual binary64 quadrature accuracy.

The mathematical specifications remain [SciPy's S1 characteristic function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levy_stable.html)
and, for validation only, equations (3) and (5) of [Saenko (2023)](https://arxiv.org/abs/2303.12488).
No upstream implementation was imported.

## Checks

`validation/check_censored_stable_skew.py` reuses the independent reference
formulas in `check_censored_stable_tails.py`: 90-digit real Fourier integration
near the center and join, and independent 32-/48-term asymptotic sums remotely.
The finite truncations must agree relatively within 1e-35 before use.
Cases cover both signs of beta=0.75 at alpha=1.25,1.5,1.75, interior shapes,
reflection, zero, central coordinates, the central/tail join, x=8,16,64,
1e100,1e300, adjacent intervals, central/join/remote intervals, affine scales
1e-200 and 1e200, and overflowed standardized coordinates. The log comparison
tolerance is 3e-10*max(1,abs(reference)); nonzero tiny complementary logs use
relative error instead.

The native regression retains the old beta=0 and +/-0.5 tail cases and adds
+/-0.25 and +/-0.75. It checks finite densities, monotonicity, complements,
reflection, join continuity, adjacent intervals, rejection just outside both
skew boundaries even for empty observations, and fixed/free mixed censored
fits. The independent SciPy scale fit uses alpha=1.5, beta=0.75, loc=0,
exact=[0.5,1,2,8], left=[-10], right=[20], interval=[30,32], scale bounds
[0.5,20]. Its reference scale is 9.026265604293837 and NNLF is
24.262217457254415. The reflected beta=-0.75 sample must give the same fit.

## Results

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: Python 3.13.3, SciPy 1.15.3 and mpmath 1.3.0.
The checked profile uses -O0 and the optimized profile uses -O3; both add
-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace. The affected
empirical_distribution, censored, fit and stats modules are compiled locally
and linked before the existing SciFort archive.

Commands run from the repository root (CMD syntax):

```text
python validation/check_censored_final.py --reuse
python -c "import sys; sys.path.insert(0,'validation'); from check_censored_stable_skew import compile_drivers; compile_drivers()"
python validation/check_censored_stable_skew.py --profile checked
python validation/check_censored_stable_skew.py --reuse --profile optimized
fpm test --profile release
```

The initial broad checker passed 953 native censored checks and the existing
81 larger-k native checks, then could not launch its reference driver because
cached executables were missing during the concurrent full FPM build. The explicit rebuild
above restored the drivers and passed all 953 censored native checks in each
profile. It did not rerun the unaffected larger-k executable.

All 560 new independent reference comparisons passed in each profile. Maximum
normalized error was 6.067496584174028e-14 in both. No reference tolerance was
relaxed. The independent SciPy fit returned the scale and NNLF recorded above;
both positive and reflected native fits passed their regression tolerances.

For a subsequent complete focused run, the rebuilt reference tool supports
`python validation/check_censored_stable_skew.py --compile` (generate fixtures)
or `--compile --reuse` (reuse fixtures). This invokes the same compilation and
native-test function used above, then the independent checks. It uses the
existing FPM module/archive cache paths, as does check_censored_final.py.

The full `fpm test --profile release` completed with exit code 1 in
9973.818 seconds (2 hours 46 minutes 14 seconds, including the first uncached
release compilation). The updated stable regression passed all 953 checks;
the existing larger-k regression passed all 81 checks. All 12 failed
executables are unchanged in this batch and exclude the changed
`scifort_censored` module, confirmed with `nm --defined-only` on each binary.

Failed executables:

- `test_descriptive_circular.exe`
- `test_descriptive_circular_errors.exe`
- `test_descriptive_entropy_errors.exe`
- `test_descriptive_frequency_errors.exe`
- `test_descriptive_kstat_errors.exe`
- `test_descriptive_lmoment_errors.exe`
- `test_descriptive_nan_raise.exe`
- `test_descriptive_proportion_trim_errors.exe`
- `test_power_transform_errors.exe`
- `test_sample_distance_errors.exe`
- `test_sample_distance_nd_errors.exe`
- `test_special_modified_fresnel.exe`

Ten failures are error-path child tests whose scalar pure-function results
are unused. In a release diagnostic, the original frequency child
`--error-case percentile_nan_raise` returned exit code 0. An otherwise
identical call to `percentileofscore([1,qnan],2,nan_policy='raise')` with its
result printed, compiled with -O3 -funroll-loops against the same release
archive, correctly returned exit code 1. Consuming the result restores the
expected termination in this configuration. These test calls do not execute
the changed censored kernel.

The other failures are numerical checks in unchanged modules:

- `circstd([1,2],high=0,low=2)` returns -2.7504585811707036, while its test
  expects -2.7504573645368673. The sample is antipodal on that circular branch,
  making its near-zero resultant sensitive to rounding.
- `modfresnelp(100000000.25)` returns F=(6.96195740127326031e-10,
  -4.95129390830604579e-9), while the test expects
  (6.96131996384738757e-10,-4.95130287082192080e-9), with tolerance 2e-23.

The full release suite is therefore not green under this compiler. These
separate release failures were diagnosed and recorded; their source and
tests were not changed as part of the stable skewness extension. The focused
checked and optimized native/reference commands above completed successfully.

## Remaining limits

General alpha outside [1.25,1.75], general skew alpha=1, and abs(beta)>0.75
remain unsupported by these censored methods. The supported normal, Cauchy,
and Levy special cases remain available. The original affine restrictions
remain. This extension does not establish broader scalar stable API accuracy.
