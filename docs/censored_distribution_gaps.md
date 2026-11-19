# Censored distribution coverage gaps

Inventory updated: 2026-11-18.

For complete native Fortran coverage, the eventual target is censored
likelihoods for all named univariate continuous families, with explicit
interval conventions for angular distributions. SciPy's
[`CensoredData`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
API targets univariate continuous fitting.

The current SciFort fitting registry in [src/fit.f90](../src/fit.f90)
contains 109 continuous families. One hundred seven already have family-specific censored
likelihoods in [src/censored.f90](../src/censored.f90): `normal`,
`exponential`, `weibull`, `lognormal`, `gamma`, `beta`, `logistic`,
`gumbel_r`, `gumbel_l`, `laplace`, `cauchy`, `halfcauchy`, `halflogistic`,
`rayleigh`, `weibull_max`, `erlang`, `gibrat`, `halfnorm`, `hypsecant`, `laplace_asymmetric`, `loglaplace`, `pareto`, `lomax`, `fisk`, `burr`, `burr12`, `uniform`, `powerlaw`, `arcsine`, `genlogistic`, `gompertz`, `invweibull`, `moyal`, `levy`, `levy_l`, `bradford`, `triang`, `trapezoid`, `dgamma`, `dweibull`, `loggamma`, `invgamma`, `betaprime`, `maxwell`, `powernorm`, `powerlognorm`, `skewcauchy`, `exponpow`, `exponweib`, `mielke`, `fatiguelife`, `johnsonsu`, `genexpon`, `cosine`, `anglit`, `semicircular`, `truncexpon`, `johnsonsb`, `kappa3`, `loguniform`, `truncpareto`, `genhalflogistic`, `foldnorm`, `foldcauchy`, `genpareto`, `gennorm`, `halfgennorm`, `rdist`, `chi2`, `chi`, `nakagami`, `t`, `f`, `jf_skew_t`, `invgauss`, `wald`, `recipinvgauss`, `gengamma`, `truncnorm`, `truncweibull_min`, `alpha`, `genextreme`, `crystalball`, `argus`, `pearson3`, `kappa4`, `gausshyper`, `irwinhall`, `vonmises_line`, `exponnorm`, `skewnorm`, `dpareto_lognorm`, `rice`, `ncx2`, `ncf`, `kstwobign`, `wrapcauchy`, `rel_breitwigner`, `geninvgauss`, `norminvgauss`, `genhyperbolic`, `nct`, `landau`, `kstwo`, `tukeylambda`, `ksone`, and `vonmises` (fixed-branch SciFort extension). That leaves 2 candidates below.

Names follow SciFort's API. The count covers named distribution families in
the current fitting registry, rather than histogram-backed or custom models.
Existing support does not establish full SciPy parity: lognormal, Gibrat, gamma, and Erlang
censored fits require fixed location, beta fits require fixed location and
scale, and gamma/Erlang/beta fits exclude exact support endpoints. See
[censored likelihoods and fitting](censored_fitting.md) for their contracts.

Additional numerical coverage remains within implemented families:
KStwo censoring requires integer 1<=n<=32, normal scale<=huge/2,
abs(loc)/scale<=100, and finite noncollapsed affine support. Larger sample
sizes remain unsupported by that censored method.
noncentral-t censoring requires 0.5<=df<=100 and abs(nc)<=10;
infinite df and wider shape ranges remain unsupported by that method.
Landau censoring uses an independent saddle-contour left tail to avoid
loss of relative accuracy in the existing scalar left-tail switch.
These likelihoods target scaled log accuracy, rather than relative accuracy
of probabilities whose logarithms have very large magnitude.
Tukey lambda currently requires abs(lam)<=5 and either lam=0 or
abs(lam)>=1e-6, with documented positive-shape affine support limits;
KSone requires integer 1<=n<=100 and documented affine support limits;
von Mises censoring is a SciFort extension on a fixed affine [-pi,pi]
branch, with the same numerical limits as vonmises_line. Crossing arcs
remain unsupported and are separate from named-family coverage.
GIG censoring currently requires abs(p)<=5 and 0.5<=b<=10;
normal-inverse Gaussian/generalized hyperbolic require 0.5<=a<=10,
abs(b)<=0.8*a, and (for GH) abs(p)<=5. The GH boundary abs(b)=a
and wider shape ranges remain unsupported by these censored methods.
Hyperbolic overflowed standardized coordinates currently use infinite limits;
retaining representable log tails in those cases remains.
wrapped-Cauchy censoring currently requires 0<c<=0.99, normal scale,
abs(loc)/scale<=10, and a finite noncollapsed affine [0,2*pi] support;
relativistic Breit-Wigner censoring requires 0.1<=rho<=100;
Rice censoring currently requires 0<=b<=10; noncentral chi-squared/F require
each degree-of-freedom parameter in [0.1,2000] and 0<=nc<=100.
Mixtures return NaN if their remainder
bound is not met within 10,001 terms; this particularly limits very remote
noncentral radial tails. Noncentral F requires normal representable
scale*dfd/dfn. Extending these limits remains;
EMG censoring requires 0.01<=k<=100; skew-normal requires abs(a)<=10;
double Pareto-lognormal requires abs(u)<=100, 0.01<=s<=10, and 0.01<=a,b<=100;
overflowed EMG/skew-normal standard coordinates retain infinite-coordinate behavior;
Gauss hypergeometric censoring requires 0.1<=a,b<=20, abs(c)<=10, -0.5<=z<=10;
Irwin-Hall requires integer 1<=n<=32; line von Mises requires 0<=kappa<=100;
these three require normal scale, abs(loc)/scale<=1e4, and representable
affine supports with scale*extent<=huge/2;
ARGUS censoring currently requires 1e-6<=chi<=100;
Pearson III uses the existing normal approximation for abs(skew)<1.6e-5,
otherwise requires 0.002<=abs(skew)<=2000;
kappa4 requires abs(h),abs(k)<=10, with each zero or at least 1e-6 in
magnitude, and k*log(h)<=log(10) when h>0;
alpha censoring currently requires 1e-6<=a<=100;
Crystal Ball requires 1e-6<=beta<=100 and 1<m<=1e6;
GEV requires abs(c)<=100 and either c=0 or abs(c)>=1e-6;
generalized gamma requires tiny<=a<=1e6 and 1e-6<=abs(c)<=1e6;
truncated normal has finite standardized bounds within [-100,100];
truncated Weibull requires 1e-6<=c<=1e6 and finite ordered truncation bounds;
truncation requires normal nonzero scaled/physical endpoints and
abs(loc)/(scale*(b-a))<=1e4; very narrow windows additionally require
zero loc and a power-of-two scale;
Student t/F/Jones-Faddy censoring requires normal positive beta parameters
with sum at most 1e6; Student t also accepts positive-infinite df;
chi/chi-squared censoring requires a positive representable half-shape;
chi/chi-squared currently require df<=2e6 and Nakagami nu<=1e6;
generalized-normal censoring requires a representable reciprocal power;
R censoring currently requires shape at most 1e6, a representable positive
half-shape, and a finite noncollapsed affine support. These restrictions
are separate from the remaining-family count below. See the API guide.

## Remaining continuous families

The groups organize related work; they are not a claim that every member
has the same numerical difficulty or can use the same fitting algorithm.

| Work group | Count | Remaining distributions |
| --- | ---: | --- |
| Numerically demanding or specialist families | 2 | `levy_stable`, `studentized_range` |
| **Total** | **2** | |

For angular distributions, distinguish intervals in a fixed coordinate range
from arcs crossing the angular boundary. `vonmises_line` is already covered with
ordinary bounded intervals: SciPy explicitly defines it on the real line, while
`vonmises` has a circular CDF convention. Its new censored method explicitly
uses the same fixed branch as vonmises_line and is labeled a SciFort extension. See
[SciPy's von Mises documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises.html)
and [wrapped Cauchy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wrapcauchy.html).
Wrapped Cauchy now has ordinary censored likelihoods on its fixed
[loc,loc+2*pi*scale] coordinate interval. Crossing-boundary arcs are not
encoded as reversed intervals; a separate arc representation remains
outside the implemented method.

## Discrete distributions: separate SciFort extension

Discrete censored likelihoods are mathematically useful, but would be a
**SciFort extension** requiring discrete endpoint conventions. They are
separate from the continuous `CensoredData` compatibility work above.
The current registry contains these 21 discrete fitting families:

`bernoulli`, `poisson`, `geometric`, `binomial`, `negative_binomial`,
`randint`, `planck`, `dlaplace`, `logser`, `betabinom`, `hypergeom`,
`nhypergeom`, `boltzmann`, `betanbinom`, `yulesimon`, `zipf`, `zipfian`,
`skellam`, `nchypergeom_fisher`, `nchypergeom_wallenius`, and `poisson_binom`.

Multivariate censoring would need a separate observation representation and
joint probability calculations; it is outside this univariate inventory.

## Suggested implementation order

1. Extend noncentral radial mixture convergence ranges and the documented
   generalized-gamma and truncation parameter ranges; Rice now has censoring support.
2. Extend the bounded-family numerical ranges and retain explicit support
   and interval conventions; all three formerly remaining bounded names
   now have censored likelihoods and fits.
3. Extend Tukey-lambda and KSone numerical ranges; their censored methods
   now include log-probability inversion and finite-sample endpoint formulas.
4. Continue through the other groups, prioritizing reusable stable kernels
   and numerical correctness over raw distribution counts.

Several additions can share validated kernels, so the remaining work need
not mean separate implementations for every method. Likelihood support and
fitting restrictions should be tracked separately as families are added.
Native Fortran coverage takes priority; C ABI and wrapper expansion remain
deferred. See the broader [SciPy gaps inventory](scipy_gaps.md).

## References

- [SciPy distribution inventory](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [SciPy CensoredData](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
- [SciFort censored likelihoods and fitting](censored_fitting.md)
- [SciFort roadmap](../ROADMAP.md)
