# Levy-stable censored skewness through 0.9

Date: 2026-11-19.

## Scope and derivation

Extend the general native S1 stable censored likelihood, NNLF and bounded-fit
rectangle from abs(beta)<=0.75 to abs(beta)<=0.9, retaining
1.25<=alpha<=1.75. Normal/Cauchy/Levy special cases, affine restrictions,
public names and the scalar distribution API are unchanged. The censored
rectangle retains all finite physical coordinates, including overflowed
standardized coordinates.

The formulas and integration rules remain those in the
[original tail derivation](stable_censored_tails_2026-11-19.md).
With a=alpha, c=1-i*beta*tan(pi*a/2), theta=pi/(8*a) and
d=c*exp(-i*a*theta), let B denote the binary64 value of 0.9. Then

```text
Re(d) >= cos(pi/8)-B*abs(tan(5*pi/8))*sin(pi/8) > 0.09238.
```

The coefficient is positive throughout the contour rotation as well. The
additional damping is at least exp(-sin(pi/14)*y). The identity
exprel(w)=integral(exp(t*w),t=0..1) gives abs(exprel(w))<=1 for Re(w)<=0.
Using abs(c)<=sqrt(1+(B*tan(5*pi/8))**2), integration of the polynomial-times-
exponential envelopes beyond y=256 bounds the factored PDF remainder by
1.06e-20 and the factored SF remainder by 4.04e-23. The powers are y**1.75
and y**0.75 respectively, and both bounds include division by pi. These
rounded-up bounds also cover the rounding of the decimal skew boundary.

Thus neither the contour angle nor the cutoff, panel width, GL16 rule,
central/tail join or exprel polynomial needs changing. The bounds establish
damping and truncation control; independent references check the actual
binary64 quadrature accuracy.

The mathematical sources remain
[SciPy's S1 characteristic function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levy_stable.html)
and, for validation only, equations (3) and (5) of
[Saenko (2023)](https://arxiv.org/abs/2303.12488).
No upstream implementation or new numerical coefficient table was imported.

## Checks

The extended `validation/check_censored_stable_skew.py` retains its previous
eight shape pairs and adds beta=+/-0.9 at alpha=1.25,1.5,1.75, plus interior
pairs (1.37,-0.825) and (1.68,0.85). It also retains the previous affine
beta=+/-0.75 cases and adds beta=+/-0.9. All 1120 cases use the independent
reference formulas in `check_censored_stable_tails.py`: 90-digit working
precision real Fourier integration near the center/join and finite asymptotic
sums remotely. The 32- and 48-term sums must agree relatively within 1e-35.

Cases cover reflection, zero, central coordinates, the abs(z)=6 join,
x=8,16,64,1e100,1e300, adjacent intervals, central/join/remote intervals,
affine scales 1e-200 and 1e200, and overflowed standardized coordinates.
The unchanged comparison tolerance is 3e-10*max(1,abs(reference)) for logs;
nonzero tiny complementary logs use relative error.

The native regression retains beta=0,+/-0.25,+/-0.5,+/-0.75 remote cases,
adds +/-0.825 and +/-0.9, and tests the new central boundary. It checks
monotonicity, complements, reflection, join continuity, adjacent intervals,
empty observations at the boundaries, and rejection at the next binary64
value outside each new boundary for both empty and nonempty samples.
Both old and new boundaries have fixed/free-scale reflected censored-fit
checks and overflowed-coordinate/adjacent-interval checks.

For alpha=1.5, beta=0.9 and loc=0, the independent SciPy scale fit uses
exact=[0.5,1,2,8], left=[-10], right=[20], interval=[30,32], and scale bounds
[0.5,20]. It returned scale 8.532460150892556 and NNLF 24.488664184760914.
The reflected beta=-0.9 sample must produce the same fit. The validation
script additionally checks the fitted objective with independent Fourier
integration for both beta=0.75 and beta=0.9.

## Validation

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: Python 3.13.3, SciPy 1.15.3 and mpmath 1.3.0.

Commands from the repository root (CMD):

```bat
python validation/check_censored_stable_skew.py --compile --profile checked
python validation/check_censored_stable_skew.py --reuse --profile optimized
fpm test --profile release
```

The compile step rebuilds empirical_distribution, censored, fit and stats
locally, links the reference driver and native regression ahead of the
existing FPM archive, and runs the native regression at both -O0 and -O3.
Both profiles add -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace.
The script uses the existing FPM module/archive cache paths documented in
the [previous extension](stable_censored_skew_2026-11-19.md). New fixtures
use the stable-skew09 prefix so old 0.75 fixtures cannot accidentally satisfy
the new --reuse command.

Both native profiles passed all 1343 checks. Both passed all 1120 independent
reference comparisons, with maximum normalized error
9.921398641842514e-14 in each profile. No tolerance was relaxed. The beta=0.75
and beta=0.9 fitted objectives also agreed with independent Fourier
integration, and both reflected native fits passed their reference checks.

Full `fpm test --profile release`: **passed**, exit status zero. Elapsed time
was 3030.294 seconds (50 minutes 30 seconds), including recompilation and
relinking. The updated stable regression passed all 1343 checks within that
run, and the unchanged large-k studentized-range regression passed all 81.
The final diff, source ASCII/line-length, Python syntax and local-link checks
passed. No numerical tests were rerun solely for the final documentation
update. Other compilers were not tested.

The current rectangle still excludes general skew alpha=1, other alpha
ranges, and abs(beta)>0.9. The boundary abs(beta)=1 is not claimed by this
general method. Existing normal/Cauchy/Levy delegated cases remain available.
