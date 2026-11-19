# Studentized-range censored likelihoods through k=100

Date: 2026-11-19.

## Scope and numerical changes

The native censored likelihood, NNLF and bounded-fit overload now support
real 2<=k<=100, extending the former k<=10 limit. Finite degrees of freedom
remain 1<=df<=100. The existing df>=100000 infinite-df approximation is
unchanged, including its positive standardized-coordinate cap of 40.
The scalar distribution routines and C ABI are unchanged.

The independent normal-range/chi-scale mixture remains based on the
[SciPy mathematical specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.studentized_range.html).
No upstream numerical implementation or coefficients were imported.

For the normal range, D=Phi(z+r)-Phi(z). The positive CDF integrand is
k*phi(z)*D**(k-1), and its r derivative gives the density. The direct SF
integrand is k*phi(z)*(sf(z)**(k-1)-D**(k-1)). Composite GL16 integration
retains the original unit-width panels for k<=10 and uses width 0.5 for
larger k. This resolves the concentrated lower-tail peak, whose limiting
width is proportional to 1/sqrt(k).

At small q, weighting the scale density by S**(k-1), where
S=sqrt(chi-square(df)/df), gives a tilted chi-squared law with df+k-1
degrees of freedom. The scale upper bound now satisfies

```text
df*Smax**2 = df+k-1+200
```

rather than df+200. At df=k=100, the former bound produced an absolute log error of about
4.78e-6 in the limiting lower tail. The revised bound's tilted
chi-squared upper probability is approximately 1.67e-15 in this worst
shape corner. The existing lower bound and log-scale panel width remain.

There is also a cancellation correction in direct SF integration. With
R=sf(z+r)/sf(z) and m=k-1, the hazard increment is -m*log(1-R).
For small R, evaluate it with log1p(-exp(log(R))); use the complementary
expm1 form near R=1 and the existing linear log limit remotely. The
former expression lost accuracy while recovering R from a value near one.
At q=12, k=100 and infinite df, the former log-SF error was about 1.93e-5.
That point is now a native regression test, with independent adaptive
reference log-SF=-29.87048776184788.

## Independent references

`validation/check_censored_studentized_range_large_k.py` uses independent
positive adaptive quadrature over the normal minimum and the physical
chi scale. Each integral is normalized near its numerical peak before
integration, so small absolute probabilities do not disappear under a
blanket absolute error threshold. Small Gaussian interval masses are
integrated directly. CDF/SF references resolve the smaller probability and
use stable logarithmic complementation for the larger probability.

80-digit mpmath erf integrals independently check small infinite-df
normal ranges. Analytic lower-end coefficients check finite and infinite
degrees of freedom:

```text
log C = log(k)/2 - (k-1)*log(2*pi)/2
        + (k-1)*log(2/df)/2
        + log_gamma((df+k-1)/2) - log_gamma(df/2)
CDF(q) ~ C*q**(k-1)
PDF(q) ~ C*(k-1)*q**(k-2)
```

For infinite df, retain only log(k)/2-(k-1)*log(2*pi)/2 in log C.
For remote finite-df q, independent normal-range moments M=E[R**df]
provide the leading asymptotes:

```text
c = log(2) + df*log(df/2)/2 - log_gamma(df/2)
log SF(q)  ~ c + log(M) - log(df) - df*log(q)
log PDF(q) ~ c + log(M) - (df+1)*log(q)
```

At q>=1e100 the next scale-density terms are negligible at the test target.
Moment quadrature through R=40 captures the Gaussian tail far beyond the
moment's peak for df<=100. These formulas also check overflowing
standardized coordinates and adjacent physical interval widths.

The grid includes k just above 10, 11, 20, 37.5, 50, 75 and 100; finite
df=1,10,100 and the infinite-df limit; ordinary CDF/SF/PDF values; tiny
lower-end probabilities; remote tails; and physical scales 1e-200 and
1e200. Nonzero complementary logs near zero use relative error. Other
logs use the existing 3e-10*max(1,abs(reference_log_value)) target.

Native tests use representative shape/df pairs, endpoint and invalid-input
checks, monotonicity, complementary tails, adjacent intervals, affine
identities, overflowing standardized coordinates, and fixed/free fits.
The Python reference grid covers the broader parameter cross-product.

A separate SciPy bounded scalar minimization checks a mixed sample at
k=df=100, loc=0: exact=[3,4,5], left=[2.5], right=[7], interval=[4.5,4.75],
scale bounds [0.5,2]. Reference scale=0.800935178253828 and
NNLF=25.123232284243887. Independent scaled adaptive integration at that
scale gives NNLF=25.12323228367478, a difference below 6e-10.

## Reproduction and results

GNU Fortran 17.0.0 20260510 (experimental), Windows. Python 3.13.3,
SciPy 1.15.3 and mpmath 1.3.0 are optional validation tools.

```cmd
python validation\check_censored_final.py --reuse
python validation\check_censored_studentized_range_large_k.py
fpm test
```

The first command compiles both -O0 and -O3 with -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace, runs both native programs,
and checks the existing stable/studentized-range references. Omit --reuse
to regenerate those fixtures. The second command generates and compares
the larger-k references; add --reuse after the initial generation to skip
that work. The saved fixtures can also be checked concurrently in separate
processes with --reuse --profile checked and --reuse --profile optimized.
The existing fpm module/archive cache paths in the first script
must be adjusted when reproducing with a different build cache/compiler.

Both checked and optimized builds passed 563 existing native checks and
81 new larger-k checks, plus 395 existing independent reference comparisons.
The latter's maximum normalized error was 3.0366159424440845e-9 in the
ordinary SciPy range comparisons, within their unchanged 3e-8 tolerance.
Independent larger-k targets retain 3e-10, with relative checks on nonzero
near-zero complementary logs. Stable-complement reference generation was
also spot-checked against the saved fixtures.

Both profiles also passed all 239 new larger-k reference comparisons. Their
maximum normalized error was 8.548717289613705e-14, including relative checks
on nonzero complementary logs near zero. In total each profile passed
644 native checks and 634 reference comparisons. Both profile-specific
--reuse commands completed successfully. No comparison tolerance was relaxed.

The full `fpm test` run passed with exit code 0 in 3324.084 seconds
(55 minutes 24 seconds), including all 81 new and 563 existing checks in
the two focused native programs. Existing array-temporary runtime warnings
occurred in unchanged modules/tests; no test failures were reported.
Only validation-note wording and result records changed after the numerical
source and native tests were finalized.

## Remaining limits

Larger k, the finite-df gap 100<df<100000, k<2, and the existing affine
restrictions remain unsupported. Finite k and df through 100 remain the
current targets. The independent log target does not promise relative
accuracy after exponentiating very large negative log values or last-bit
monotonicity across numerical switches. Other compilers were not tested.
