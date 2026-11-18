# Censored distribution coverage gaps

Inventory updated: 2026-11-17.

For complete native Fortran coverage, the eventual target is censored
likelihoods for all named univariate continuous families, with explicit
interval conventions for angular distributions. SciPy's
[`CensoredData`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
API targets univariate continuous fitting.

The current SciFort fitting registry in [src/fit.f90](../src/fit.f90)
contains 109 continuous families. Eighty-six already have family-specific censored
likelihoods in [src/censored.f90](../src/censored.f90): `normal`,
`exponential`, `weibull`, `lognormal`, `gamma`, `beta`, `logistic`,
`gumbel_r`, `gumbel_l`, `laplace`, `cauchy`, `halfcauchy`, `halflogistic`,
`rayleigh`, `weibull_max`, `erlang`, `gibrat`, `halfnorm`, `hypsecant`, `laplace_asymmetric`, `loglaplace`, `pareto`, `lomax`, `fisk`, `burr`, `burr12`, `uniform`, `powerlaw`, `arcsine`, `genlogistic`, `gompertz`, `invweibull`, `moyal`, `levy`, `levy_l`, `bradford`, `triang`, `trapezoid`, `dgamma`, `dweibull`, `loggamma`, `invgamma`, `betaprime`, `maxwell`, `powernorm`, `powerlognorm`, `skewcauchy`, `exponpow`, `exponweib`, `mielke`, `fatiguelife`, `johnsonsu`, `genexpon`, `cosine`, `anglit`, `semicircular`, `truncexpon`, `johnsonsb`, `kappa3`, `loguniform`, `truncpareto`, `genhalflogistic`, `foldnorm`, `foldcauchy`, `genpareto`, `gennorm`, `halfgennorm`, `rdist`, `chi2`, `chi`, `nakagami`, `t`, `f`, `jf_skew_t`, `invgauss`, `wald`, `recipinvgauss`, `gengamma`, `truncnorm`, `truncweibull_min`, `alpha`, `genextreme`, `crystalball`, `argus`, `pearson3`, and `kappa4`. That leaves 23 candidates below.

Names follow SciFort's API. The count covers named distribution families in
the current fitting registry, rather than histogram-backed or custom models.
Existing support does not establish full SciPy parity: lognormal, Gibrat, gamma, and Erlang
censored fits require fixed location, beta fits require fixed location and
scale, and gamma/Erlang/beta fits exclude exact support endpoints. See
[censored likelihoods and fitting](censored_fitting.md) for their contracts.

Additional numerical coverage remains within implemented families:
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
| Lifetime and tail models | 1 | `rice` |
| Bounded and truncated distributions | 3 | `gausshyper`, `irwinhall`, `vonmises_line` |
| Additional shape and skew families | 4 | `exponnorm`, `skewnorm`, `tukeylambda`, `dpareto_lognorm` |
| Numerically demanding or specialist families | 13 | `ncx2`, `ncf`, `nct`, `geninvgauss`, `norminvgauss`, `genhyperbolic`, `landau`, `levy_stable`, `rel_breitwigner`, `studentized_range`, `ksone`, `kstwo`, `kstwobign` |
| Angular distributions requiring explicit interval conventions | 2 | `vonmises`, `wrapcauchy` |
| **Total** | **23** | |

For angular distributions, distinguish intervals in a fixed coordinate range
from arcs crossing the angular boundary. `vonmises_line` belongs in the
ordinary bounded group: SciPy explicitly defines it on the real line, while
`vonmises` has a circular CDF convention. See
[SciPy's von Mises documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises.html)
and [wrapped Cauchy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wrapcauchy.html).

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

1. Add Rice with stable noncentral radial tails and extend the documented
   generalized-gamma and truncation parameter ranges.
2. Extend the remaining bounded and truncated models with explicit support
   and interval conventions.
3. Continue the remaining shape and skew models through reusable stable tails.
4. Continue through the other groups, prioritizing reusable stable kernels
   and numerical correctness over raw distribution counts.

Several additions can share validated kernels, so the remaining work need
not mean 23 independent numerical implementations. Likelihood support and
fitting restrictions should be tracked separately as families are added.
Native Fortran coverage takes priority; C ABI and wrapper expansion remain
deferred. See the broader [SciPy gaps inventory](scipy_gaps.md).

## References

- [SciPy distribution inventory](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [SciPy CensoredData](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
- [SciFort censored likelihoods and fitting](censored_fitting.md)
- [SciFort roadmap](../ROADMAP.md)
