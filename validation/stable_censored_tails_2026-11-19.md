# Levy-stable censored tails: validation and derivation

Date: 2026-11-19.

## Scope

The native censored likelihood, NNLF and bounded fit retain the supported
S1 shapes 1.25<=alpha<=1.75 and abs(beta)<=0.5. The former abs(z)<=6
restriction is removed. Existing normal/Cauchy/Levy special cases and affine
parameter restrictions are unchanged. No scalar distribution API or C ABI
is changed. All finite physical observations are supported within the stated
parameter restrictions, including overflow of the standardized distance.

## Independent contour derivation

Start from the S1 characteristic function in the
[SciPy mathematical specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levy_stable.html).
For positive standardized x, let a=alpha, c=1-i*beta*tan(pi*a/2),
theta=pi/(8*a), r=exp(-i*theta), d=c*r**a, q=x**(-a).
Rotate the positive Fourier ray to t=r*y/x. Throughout this rotation,
Re(c*exp(-i*a*angle)) is positive for the supported shapes. In particular,
Re(d)>=cos(pi/8)-0.5*abs(tan(5*pi/8))*sin(pi/8)>0.46.
The exponential exp(-i*r*y) additionally damps at least as fast as
exp(-sin(pi/14)*y), independently of x.

Define w=-d*q*y**a and

```text
H(y) = -d*y**a * exprel(w) * exp(-i*r*y)
exprel(w) = (exp(w)-1)/w, exprel(0)=1
pdf(x) = x**(-a-1)/pi * Re(r * integral(H(y), y=0..infinity))
sf(x)  = x**(-a)/pi * Im(integral(H(y)/y, y=0..infinity))
```

The subtraction of 1 inside exprel removes a baseline Fourier term with
zero real density integral and imaginary survival integral -pi/2. It thus
removes cancellation of the one-half term before numerical integration.
Factoring q outside the integral preserves log density and log SF even
when q underflows. Physical log distance minus log scale avoids explicitly
forming an overflowing standardized coordinate. Reflection gives the
negative tail; stable logarithmic complementation gives the large tail.

For abs(z)<=6 the original real Fourier integral is retained. Beyond this
join, set y=v**4 and integrate 0<=v<=4 using composite GL16 with maximum
panel width 0.025. The discarded y>256 integrand is exponentially damped;
its absolute density bound after factoring q is below 1e-19 over this
shape rectangle. For abs(w)<0.5, a 21-term Taylor polynomial evaluates
exprel without subtraction. Its next-term bound is below 1e-26. These
choices are numerical implementation decisions checked against independent
references, not imported constants or coefficient tables.

Narrow intervals retain positive density integration and the physical
log width. Only the central segment needs absolute-width subdivision;
relatively narrow tail intervals use GL16. Well-separated endpoints use
log differences of the directly computed small tails. The shared positive
studentized-range interval path is unchanged.

## References and checks

`validation/check_censored_stable_tails.py` uses 90-digit mpmath real-axis
Fourier integration at the join, x=8 and x=16. At x>=64 it uses finite
asymptotic sums from equations (3) and (5) in
[Saenko (2023)](https://arxiv.org/abs/2303.12488). Convert that paper's C
parameterization by magnitude abs(1+i*beta*tan(pi*a/2)) and phase
pi*a/2+arg(1+i*beta*tan(pi*a/2)). Independent 32- and 48-term truncations
must agree relatively within 1e-35 before a fixture is accepted. The
asymptotic expansion is used only in validation, not the runtime kernel.

Fixtures cover the shape boundaries and interior nonsymmetric shapes,
both tails, coordinates through 1e300, physical scales from 1e-200 to
1e200, overflowed standardized coordinates, adjacent binary64 intervals,
and ordinary tail intervals. Nonzero near-zero complementary logs are
checked with relative error; other logs use 3e-10*max(1,abs(reference)).
Native tests cover monotonicity, reflection, complements, continuity at
abs(z)=6, adjacent endpoints and fixed/free fits. A SciPy bounded scalar
minimization independently checks the mixed-sample fit: alpha=1.5,
beta=0.5, loc=0, exact=[0.5,1,2,8], left=[-10], right=[20], interval=[30,32],
scale bounds [0.5,20]. Reference scale=9.92118608975177 and
NNLF=23.868504660682664.

## Reproduction

GNU Fortran 17.0.0 20260510 (experimental), Windows.

```cmd
python validation\check_censored_final.py --reuse
python validation\check_censored_stable_tails.py
fpm test
```

The first command reuses the previously generated central/studentized-range
fixtures; omit --reuse to regenerate them. It compiles both -O0 and -O3
profiles with -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace.
The reference tooling uses existing fpm module/archive cache paths listed
in check_censored_final.py; adjust them for another build/compiler cache.
SciPy and mpmath are optional validation dependencies only. After generating
the tail fixtures once, append --reuse to the second command to skip their
high-precision regeneration. That cached-fixture path was also checked.

Python 3.13.3, SciPy 1.15.3 and mpmath 1.3.0 supplied validation tools.
Both builds passed 500 independent tail comparisons; the maximum normalized
log error was 6.04243045010842e-14. Both builds also passed 563 native checks and 395 existing
central/studentized-range reference comparisons. Their maximum normalized
error was 3.038311982575117e-9 in the ordinary SciPy range references
(tolerance 3e-8); independent stable references use the tighter 3e-10 target.
The full `fpm test` run passed with exit code 0 in 1703.076 seconds
(28 minutes 23 seconds), including all 563 checks in the changed native test.
Existing array-temporary runtime warnings occurred in unchanged interpolation,
optimization and spatial code; no test failures were reported.
A supplementary 525-value SciPy grid through abs(z)=32 agreed within
2.65e-14 in scaled log error. Adjacent intervals at the largest positive
and negative finite binary64 endpoints remained finite under runtime traps.

## Remaining limitations

Wider alpha/beta ranges, general skew alpha=1, and the existing affine
parameter restrictions remain. The original scalar stable routine is
unchanged and may differ in remote-tail or near-mode behavior. This
extension establishes neither full SciPy parity nor additional C ABI coverage.
