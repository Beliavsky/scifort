# Censored distribution coverage gaps

Inventory date: 2026-11-12.

For complete native Fortran coverage, the eventual target is censored
likelihoods for all named univariate continuous families, with explicit
interval conventions for angular distributions. SciPy's
[`CensoredData`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
API targets univariate continuous fitting.

The current SciFort fitting registry in [src/fit.f90](../src/fit.f90)
contains 109 continuous families. Nine already have family-specific censored
likelihoods in [src/censored.f90](../src/censored.f90): `normal`,
`exponential`, `weibull`, `lognormal`, `gamma`, `beta`, `logistic`,
`gumbel_r`, and `gumbel_l`. That leaves 100 candidates below.

Names follow SciFort's API. The count covers named distribution families in
the current fitting registry, rather than histogram-backed or custom models.
Existing support does not establish full SciPy parity: lognormal and gamma
censored fits require fixed location, beta fits require fixed location and
scale, and gamma/beta fits exclude exact support endpoints. See
[censored likelihoods and fitting](censored_fitting.md) for their contracts.

## Remaining continuous families

The groups organize related work; they are not a claim that every member
has the same numerical difficulty or can use the same fitting algorithm.

| Work group | Count | Remaining distributions |
| --- | ---: | --- |
| Basic location/scale and skew models | 7 | `laplace`, `cauchy`, `t`, `hypsecant`, `moyal`, `laplace_asymmetric`, `skewcauchy` |
| Related families and transformations | 20 | `chi2`, `chi`, `erlang`, `rayleigh`, `maxwell`, `halfnorm`, `halfcauchy`, `halflogistic`, `gibrat`, `wald`, `levy`, `levy_l`, `weibull_max`, `invweibull`, `loglaplace`, `invgamma`, `betaprime`, `f`, `foldnorm`, `foldcauchy` |
| Lifetime and tail models | 17 | `pareto`, `lomax`, `genpareto`, `fisk`, `burr`, `burr12`, `mielke`, `fatiguelife`, `invgauss`, `recipinvgauss`, `nakagami`, `rice`, `gompertz`, `genexpon`, `exponpow`, `exponweib`, `gengamma` |
| Bounded and truncated distributions | 20 | `uniform`, `powerlaw`, `triang`, `trapezoid`, `arcsine`, `bradford`, `truncexpon`, `loguniform`, `truncnorm`, `truncpareto`, `truncweibull_min`, `rdist`, `cosine`, `semicircular`, `anglit`, `argus`, `johnsonsb`, `gausshyper`, `irwinhall`, `vonmises_line` |
| Additional shape and skew families | 21 | `alpha`, `genlogistic`, `gennorm`, `halfgennorm`, `genhalflogistic`, `loggamma`, `powernorm`, `powerlognorm`, `dgamma`, `dweibull`, `exponnorm`, `johnsonsu`, `genextreme`, `kappa3`, `kappa4`, `crystalball`, `jf_skew_t`, `pearson3`, `skewnorm`, `tukeylambda`, `dpareto_lognorm` |
| Numerically demanding or specialist families | 13 | `ncx2`, `ncf`, `nct`, `geninvgauss`, `norminvgauss`, `genhyperbolic`, `landau`, `levy_stable`, `rel_breitwigner`, `studentized_range`, `ksone`, `kstwo`, `kstwobign` |
| Angular distributions requiring explicit interval conventions | 2 | `vonmises`, `wrapcauchy` |
| **Total** | **100** | |

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

1. Add Laplace censored likelihoods and bounded fitting.
2. Extend related families through validated transformations and special cases.
3. Add Pareto, Lomax, and generalized Pareto censoring support.
4. Continue through the other groups, prioritizing reusable stable kernels
   and numerical correctness over raw distribution counts.

Several additions can share validated kernels, so the remaining work need
not mean 100 independent numerical implementations. Likelihood support and
fitting restrictions should be tracked separately as families are added.
Native Fortran coverage takes priority; C ABI and wrapper expansion remain
deferred. See the broader [SciPy gaps inventory](scipy_gaps.md).

## References

- [SciPy distribution inventory](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [SciPy CensoredData](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
- [SciFort censored likelihoods and fitting](censored_fitting.md)
- [SciFort roadmap](../ROADMAP.md)
