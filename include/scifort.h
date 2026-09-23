/* SPDX-License-Identifier: MIT */
/* Copyright (c) 2026 SciFort contributors */

#ifndef SCIFORT_H
#define SCIFORT_H

#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

enum {
    SCIFORT_STATUS_OK = 0,
    SCIFORT_STATUS_INVALID_ARGUMENT = 1
};

int scifort_version_major(void);
int scifort_version_minor(void);
int scifort_version_patch(void);

/* Normal distribution. */
double scifort_normal_pdf_f64(double x, double loc, double scale);
double scifort_normal_cdf_f64(double x, double loc, double scale);
double scifort_normal_ppf_f64(double p, double loc, double scale);

/* Uniform distribution. */
double scifort_uniform_pdf_f64(double x, double loc, double scale);
double scifort_uniform_cdf_f64(double x, double loc, double scale);
double scifort_uniform_ppf_f64(double p, double loc, double scale);

/* Exponential distribution. */
double scifort_exponential_pdf_f64(double x, double loc, double scale);
double scifort_exponential_cdf_f64(double x, double loc, double scale);
double scifort_exponential_ppf_f64(double p, double loc, double scale);

/* Laplace distribution. */
double scifort_laplace_pdf_f64(double x, double loc, double scale);
double scifort_laplace_cdf_f64(double x, double loc, double scale);
double scifort_laplace_ppf_f64(double p, double loc, double scale);

/* Logistic distribution. */
double scifort_logistic_pdf_f64(double x, double loc, double scale);
double scifort_logistic_cdf_f64(double x, double loc, double scale);
double scifort_logistic_ppf_f64(double p, double loc, double scale);

/* Cauchy distribution. */
double scifort_cauchy_pdf_f64(double x, double loc, double scale);
double scifort_cauchy_cdf_f64(double x, double loc, double scale);
double scifort_cauchy_ppf_f64(double p, double loc, double scale);

/* Rayleigh distribution. */
double scifort_rayleigh_pdf_f64(double x, double loc, double scale);
double scifort_rayleigh_cdf_f64(double x, double loc, double scale);
double scifort_rayleigh_ppf_f64(double p, double loc, double scale);

/* Gamma distribution. */
double scifort_gamma_pdf_f64(double x, double a, double loc, double scale);
double scifort_gamma_cdf_f64(double x, double a, double loc, double scale);
double scifort_gamma_ppf_f64(double p, double a, double loc, double scale);

/* Chi2 distribution. */
double scifort_chi2_pdf_f64(double x, double df, double loc, double scale);
double scifort_chi2_cdf_f64(double x, double df, double loc, double scale);
double scifort_chi2_ppf_f64(double p, double df, double loc, double scale);

/* T distribution. */
double scifort_t_pdf_f64(double x, double df, double loc, double scale);
double scifort_t_cdf_f64(double x, double df, double loc, double scale);
double scifort_t_ppf_f64(double p, double df, double loc, double scale);

/* Lognormal distribution. */
double scifort_lognormal_pdf_f64(double x, double s, double loc, double scale);
double scifort_lognormal_cdf_f64(double x, double s, double loc, double scale);
double scifort_lognormal_ppf_f64(double p, double s, double loc, double scale);

/* Weibull distribution. */
double scifort_weibull_pdf_f64(double x, double c, double loc, double scale);
double scifort_weibull_cdf_f64(double x, double c, double loc, double scale);
double scifort_weibull_ppf_f64(double p, double c, double loc, double scale);

/* Pareto distribution. */
double scifort_pareto_pdf_f64(double x, double b, double loc, double scale);
double scifort_pareto_cdf_f64(double x, double b, double loc, double scale);
double scifort_pareto_ppf_f64(double p, double b, double loc, double scale);

/* Beta distribution. */
double scifort_beta_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_beta_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_beta_ppf_f64(double p, double a, double b, double loc, double scale);

/* F distribution. */
double scifort_f_pdf_f64(double x, double dfn, double dfd, double loc, double scale);
double scifort_f_cdf_f64(double x, double dfn, double dfd, double loc, double scale);
double scifort_f_ppf_f64(double p, double dfn, double dfd, double loc, double scale);

/* Right-skewed Gumbel distribution. */
double scifort_gumbel_r_pdf_f64(double x, double loc, double scale);
double scifort_gumbel_r_cdf_f64(double x, double loc, double scale);
double scifort_gumbel_r_ppf_f64(double p, double loc, double scale);

/* Left-skewed Gumbel distribution. */
double scifort_gumbel_l_pdf_f64(double x, double loc, double scale);
double scifort_gumbel_l_cdf_f64(double x, double loc, double scale);
double scifort_gumbel_l_ppf_f64(double p, double loc, double scale);

/* Power-function distribution. */
double scifort_powerlaw_pdf_f64(double x, double a, double loc, double scale);
double scifort_powerlaw_cdf_f64(double x, double a, double loc, double scale);
double scifort_powerlaw_ppf_f64(double p, double a, double loc, double scale);

/* Triangular distribution. */
double scifort_triang_pdf_f64(double x, double c, double loc, double scale);
double scifort_triang_cdf_f64(double x, double c, double loc, double scale);
double scifort_triang_ppf_f64(double p, double c, double loc, double scale);

/* Generalized Pareto distribution. */
double scifort_genpareto_pdf_f64(double x, double c, double loc, double scale);
double scifort_genpareto_cdf_f64(double x, double c, double loc, double scale);
double scifort_genpareto_ppf_f64(double p, double c, double loc, double scale);

/* Arcsine distribution. */
double scifort_arcsine_pdf_f64(double x, double loc, double scale);
double scifort_arcsine_cdf_f64(double x, double loc, double scale);
double scifort_arcsine_ppf_f64(double p, double loc, double scale);

/* Half-normal distribution. */
double scifort_halfnorm_pdf_f64(double x, double loc, double scale);
double scifort_halfnorm_cdf_f64(double x, double loc, double scale);
double scifort_halfnorm_ppf_f64(double p, double loc, double scale);
double scifort_halfcauchy_pdf_f64(double x, double loc, double scale);
double scifort_halfcauchy_cdf_f64(double x, double loc, double scale);
double scifort_halfcauchy_ppf_f64(double p, double loc, double scale);
double scifort_lomax_pdf_f64(double x, double c, double loc, double scale);
double scifort_lomax_cdf_f64(double x, double c, double loc, double scale);
double scifort_lomax_ppf_f64(double p, double c, double loc, double scale);

/* Chi distribution. */
double scifort_chi_pdf_f64(double x, double df, double loc, double scale);
double scifort_chi_cdf_f64(double x, double df, double loc, double scale);
double scifort_chi_ppf_f64(double p, double df, double loc, double scale);

/* Maxwell distribution. */
double scifort_maxwell_pdf_f64(double x, double loc, double scale);
double scifort_maxwell_cdf_f64(double x, double loc, double scale);
double scifort_maxwell_ppf_f64(double p, double loc, double scale);

/* Cosine distribution. */
double scifort_cosine_pdf_f64(double x, double loc, double scale);
double scifort_cosine_cdf_f64(double x, double loc, double scale);
double scifort_cosine_ppf_f64(double p, double loc, double scale);

/* Semicircular distribution. */
double scifort_semicircular_pdf_f64(double x, double loc, double scale);
double scifort_semicircular_cdf_f64(double x, double loc, double scale);
double scifort_semicircular_ppf_f64(double p, double loc, double scale);

/* Anglit distribution. */
double scifort_anglit_pdf_f64(double x, double loc, double scale);
double scifort_anglit_cdf_f64(double x, double loc, double scale);
double scifort_anglit_ppf_f64(double p, double loc, double scale);

/* Moyal distribution. */
double scifort_moyal_pdf_f64(double x, double loc, double scale);
double scifort_moyal_cdf_f64(double x, double loc, double scale);
double scifort_moyal_ppf_f64(double p, double loc, double scale);

/* Hyperbolic-secant distribution. */
double scifort_hypsecant_pdf_f64(double x, double loc, double scale);
double scifort_hypsecant_cdf_f64(double x, double loc, double scale);
double scifort_hypsecant_ppf_f64(double p, double loc, double scale);

/* Half-logistic distribution. */
double scifort_halflogistic_pdf_f64(double x, double loc, double scale);
double scifort_halflogistic_cdf_f64(double x, double loc, double scale);
double scifort_halflogistic_ppf_f64(double p, double loc, double scale);

/* Inverse-gamma distribution. */
double scifort_invgamma_pdf_f64(double x, double a, double loc, double scale);
double scifort_invgamma_cdf_f64(double x, double a, double loc, double scale);
double scifort_invgamma_ppf_f64(double p, double a, double loc, double scale);

/* Inverse-Gaussian distribution. */
double scifort_invgauss_pdf_f64(double x, double mu_shape, double loc, double scale);
double scifort_invgauss_cdf_f64(double x, double mu_shape, double loc, double scale);
double scifort_invgauss_ppf_f64(double p, double mu_shape, double loc, double scale);

/* Levy distribution. */
double scifort_levy_pdf_f64(double x, double loc, double scale);
double scifort_levy_cdf_f64(double x, double loc, double scale);
double scifort_levy_ppf_f64(double p, double loc, double scale);

/* Log-Laplace distribution. */
double scifort_loglaplace_pdf_f64(double x, double c, double loc, double scale);
double scifort_loglaplace_cdf_f64(double x, double c, double loc, double scale);
double scifort_loglaplace_ppf_f64(double p, double c, double loc, double scale);
double scifort_bradford_pdf_f64(double x, double c, double loc, double scale);
double scifort_bradford_cdf_f64(double x, double c, double loc, double scale);
double scifort_bradford_ppf_f64(double p, double c, double loc, double scale);
double scifort_truncexpon_pdf_f64(double x, double b, double loc, double scale);
double scifort_truncexpon_cdf_f64(double x, double b, double loc, double scale);
double scifort_truncexpon_ppf_f64(double p, double b, double loc, double scale);
double scifort_fisk_pdf_f64(double x, double c, double loc, double scale);
double scifort_fisk_cdf_f64(double x, double c, double loc, double scale);
double scifort_fisk_ppf_f64(double p, double c, double loc, double scale);
double scifort_dweibull_pdf_f64(double x, double c, double loc, double scale);
double scifort_dweibull_cdf_f64(double x, double c, double loc, double scale);
double scifort_dweibull_ppf_f64(double p, double c, double loc, double scale);


/* Alpha distribution. */
double scifort_alpha_pdf_f64(double x, double a, double loc, double scale);
double scifort_alpha_cdf_f64(double x, double a, double loc, double scale);
double scifort_alpha_ppf_f64(double p, double a, double loc, double scale);

/* Birnbaum-Saunders fatigue-life distribution. */
double scifort_fatiguelife_pdf_f64(double x, double c, double loc, double scale);
double scifort_fatiguelife_cdf_f64(double x, double c, double loc, double scale);
double scifort_fatiguelife_ppf_f64(double p, double c, double loc, double scale);

/* Generalized logistic distribution. */
double scifort_genlogistic_pdf_f64(double x, double c, double loc, double scale);
double scifort_genlogistic_cdf_f64(double x, double c, double loc, double scale);
double scifort_genlogistic_ppf_f64(double p, double c, double loc, double scale);

/* Generalized normal distribution. */
double scifort_gennorm_pdf_f64(double x, double beta, double loc, double scale);
double scifort_gennorm_cdf_f64(double x, double beta, double loc, double scale);
double scifort_gennorm_ppf_f64(double p, double beta, double loc, double scale);
double scifort_nakagami_pdf_f64(double x, double nu, double loc, double scale);
double scifort_nakagami_cdf_f64(double x, double nu, double loc, double scale);
double scifort_nakagami_ppf_f64(double p, double nu, double loc, double scale);
double scifort_powernorm_pdf_f64(double x, double c, double loc, double scale);
double scifort_powernorm_cdf_f64(double x, double c, double loc, double scale);
double scifort_powernorm_ppf_f64(double p, double c, double loc, double scale);
double scifort_loggamma_pdf_f64(double x, double c, double loc, double scale);
double scifort_loggamma_cdf_f64(double x, double c, double loc, double scale);
double scifort_loggamma_ppf_f64(double p, double c, double loc, double scale);
double scifort_wald_pdf_f64(double x, double loc, double scale);
double scifort_wald_cdf_f64(double x, double loc, double scale);
double scifort_wald_ppf_f64(double p, double loc, double scale);
double scifort_gompertz_pdf_f64(double x, double c, double loc, double scale);
double scifort_gompertz_cdf_f64(double x, double c, double loc, double scale);
double scifort_gompertz_ppf_f64(double p, double c, double loc, double scale);
double scifort_invweibull_pdf_f64(double x, double c, double loc, double scale);
double scifort_invweibull_cdf_f64(double x, double c, double loc, double scale);
double scifort_invweibull_ppf_f64(double p, double c, double loc, double scale);
double scifort_betaprime_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_betaprime_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_betaprime_ppf_f64(double p, double a, double b, double loc, double scale);
double scifort_burr12_pdf_f64(double x, double c, double d, double loc, double scale);
double scifort_burr12_cdf_f64(double x, double c, double d, double loc, double scale);
double scifort_burr12_ppf_f64(double p, double c, double d, double loc, double scale);

/* Generalized half-logistic distribution. */
double scifort_genhalflogistic_pdf_f64(double x, double c, double loc, double scale);
double scifort_genhalflogistic_cdf_f64(double x, double c, double loc, double scale);
double scifort_genhalflogistic_ppf_f64(double p, double c, double loc, double scale);

/* Exponential-power distribution. */
double scifort_exponpow_pdf_f64(double x, double b, double loc, double scale);
double scifort_exponpow_cdf_f64(double x, double b, double loc, double scale);
double scifort_exponpow_ppf_f64(double p, double b, double loc, double scale);

/* Exponentiated Weibull distribution. */
double scifort_exponweib_pdf_f64(double x, double a, double c, double loc, double scale);
double scifort_exponweib_cdf_f64(double x, double a, double c, double loc, double scale);
double scifort_exponweib_ppf_f64(double p, double a, double c, double loc, double scale);

/* Power log-normal distribution. */
double scifort_powerlognorm_pdf_f64(double x, double c, double s, double loc, double scale);
double scifort_powerlognorm_cdf_f64(double x, double c, double s, double loc, double scale);
double scifort_powerlognorm_ppf_f64(double p, double c, double s, double loc, double scale);

/* Left-Levy distribution. */
double scifort_levy_l_pdf_f64(double x, double loc, double scale);
double scifort_levy_l_cdf_f64(double x, double loc, double scale);
double scifort_levy_l_ppf_f64(double p, double loc, double scale);

/* Weibull maximum distribution. */
double scifort_weibull_max_pdf_f64(double x, double c, double loc, double scale);
double scifort_weibull_max_cdf_f64(double x, double c, double loc, double scale);
double scifort_weibull_max_ppf_f64(double p, double c, double loc, double scale);

/* R distribution. */
double scifort_rdist_pdf_f64(double x, double c, double loc, double scale);
double scifort_rdist_cdf_f64(double x, double c, double loc, double scale);
double scifort_rdist_ppf_f64(double p, double c, double loc, double scale);

/* Skew-Cauchy distribution. */
double scifort_skewcauchy_pdf_f64(double x, double a, double loc, double scale);
double scifort_skewcauchy_cdf_f64(double x, double a, double loc, double scale);
double scifort_skewcauchy_ppf_f64(double p, double a, double loc, double scale);

/* Double-gamma distribution. */
double scifort_dgamma_pdf_f64(double x, double a, double loc, double scale);
double scifort_dgamma_cdf_f64(double x, double a, double loc, double scale);
double scifort_dgamma_ppf_f64(double p, double a, double loc, double scale);

/* Asymmetric Laplace distribution. */
double scifort_laplace_asymmetric_pdf_f64(double x, double kappa, double loc, double scale);
double scifort_laplace_asymmetric_cdf_f64(double x, double kappa, double loc, double scale);
double scifort_laplace_asymmetric_ppf_f64(double p, double kappa, double loc, double scale);

/* Truncated normal distribution with standardized finite bounds a < b. */
double scifort_truncnorm_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_truncnorm_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_truncnorm_ppf_f64(double p, double a, double b, double loc, double scale);

/* Log-uniform / reciprocal distribution. */
double scifort_loguniform_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_loguniform_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_loguniform_ppf_f64(double p, double a, double b, double loc, double scale);


/* Folded normal distribution. */
double scifort_foldnorm_pdf_f64(double x, double c, double loc, double scale);
double scifort_foldnorm_cdf_f64(double x, double c, double loc, double scale);
double scifort_foldnorm_ppf_f64(double p, double c, double loc, double scale);

/* Folded Cauchy distribution. */
double scifort_foldcauchy_pdf_f64(double x, double c, double loc, double scale);
double scifort_foldcauchy_cdf_f64(double x, double c, double loc, double scale);
double scifort_foldcauchy_ppf_f64(double p, double c, double loc, double scale);

/* Reciprocal inverse-Gaussian distribution. */
double scifort_recipinvgauss_pdf_f64(double x, double mu, double loc, double scale);
double scifort_recipinvgauss_cdf_f64(double x, double mu, double loc, double scale);
double scifort_recipinvgauss_ppf_f64(double p, double mu, double loc, double scale);

/* Truncated Pareto distribution. */
double scifort_truncpareto_pdf_f64(double x, double b, double c, double loc, double scale);
double scifort_truncpareto_cdf_f64(double x, double b, double c, double loc, double scale);
double scifort_truncpareto_ppf_f64(double p, double b, double c, double loc, double scale);

/* Bernoulli distribution. */
double scifort_bernoulli_pmf_f64(double k, double p, double loc);
double scifort_bernoulli_cdf_f64(double k, double p, double loc);
double scifort_bernoulli_ppf_f64(double probability, double p, double loc);

/* Poisson distribution. */
double scifort_poisson_pmf_f64(double k, double mu, double loc);
double scifort_poisson_cdf_f64(double k, double mu, double loc);
double scifort_poisson_ppf_f64(double p, double mu, double loc);

/* Geometric distribution. */
double scifort_geometric_pmf_f64(double k, double p, double loc);
double scifort_geometric_cdf_f64(double k, double p, double loc);
double scifort_geometric_ppf_f64(double probability, double p, double loc);

/* Binomial distribution. */
double scifort_binomial_pmf_f64(double k, double n, double p, double loc);
double scifort_binomial_cdf_f64(double k, double n, double p, double loc);
double scifort_binomial_ppf_f64(double probability, double n, double p, double loc);

/* Negative Binomial distribution. */
double scifort_negative_binomial_pmf_f64(double k, double n, double p, double loc);
double scifort_negative_binomial_cdf_f64(double k, double n, double p, double loc);
double scifort_negative_binomial_ppf_f64(double probability, double n, double p, double loc);

/* Bulk normal routines; status is SCIFORT_STATUS_OK or SCIFORT_STATUS_INVALID_ARGUMENT. */
void scifort_normal_pdf_vec_f64(
    size_t n,
    const double *x,
    double loc,
    double scale,
    double *y,
    int *status
);

void scifort_normal_cdf_vec_f64(
    size_t n,
    const double *x,
    double loc,
    double scale,
    double *y,
    int *status
);

void scifort_normal_ppf_vec_f64(
    size_t n,
    const double *p,
    double loc,
    double scale,
    double *x,
    int *status
);


/* Exponentially modified normal distribution. */
double scifort_exponnorm_pdf_f64(double x, double k, double loc, double scale);
double scifort_exponnorm_cdf_f64(double x, double k, double loc, double scale);
double scifort_exponnorm_ppf_f64(double p, double k, double loc, double scale);

/* Johnson SB distribution. */
double scifort_johnsonsb_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_johnsonsb_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_johnsonsb_ppf_f64(double p, double a, double b, double loc, double scale);

/* Johnson SU distribution. */
double scifort_johnsonsu_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_johnsonsu_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_johnsonsu_ppf_f64(double p, double a, double b, double loc, double scale);

/* Trapezoidal distribution. */
double scifort_trapezoid_pdf_f64(double x, double c, double d, double loc, double scale);
double scifort_trapezoid_cdf_f64(double x, double c, double d, double loc, double scale);
double scifort_trapezoid_ppf_f64(double p, double c, double d, double loc, double scale);


/* Burr Type III distribution. */
double scifort_burr_pdf_f64(double x, double c, double d, double loc, double scale);
double scifort_burr_cdf_f64(double x, double c, double d, double loc, double scale);
double scifort_burr_ppf_f64(double p, double c, double d, double loc, double scale);

/* Mielke beta-kappa / Dagum distribution. */
double scifort_mielke_pdf_f64(double x, double k, double s, double loc, double scale);
double scifort_mielke_cdf_f64(double x, double k, double s, double loc, double scale);
double scifort_mielke_ppf_f64(double p, double k, double s, double loc, double scale);

/* Gibrat distribution. */
double scifort_gibrat_pdf_f64(double x, double loc, double scale);
double scifort_gibrat_cdf_f64(double x, double loc, double scale);
double scifort_gibrat_ppf_f64(double p, double loc, double scale);

/* Wrapped Cauchy distribution. */
double scifort_wrapcauchy_pdf_f64(double x, double c, double loc, double scale);
double scifort_wrapcauchy_cdf_f64(double x, double c, double loc, double scale);
double scifort_wrapcauchy_ppf_f64(double p, double c, double loc, double scale);

/* Generalized extreme-value distribution (SciPy sign convention). */
double scifort_genextreme_pdf_f64(double x, double c, double loc, double scale);
double scifort_genextreme_cdf_f64(double x, double c, double loc, double scale);
double scifort_genextreme_ppf_f64(double p, double c, double loc, double scale);

/* Kappa-3 distribution. */
double scifort_kappa3_pdf_f64(double x, double a, double loc, double scale);
double scifort_kappa3_cdf_f64(double x, double a, double loc, double scale);
double scifort_kappa3_ppf_f64(double p, double a, double loc, double scale);

/* Kappa-4 distribution. */
double scifort_kappa4_pdf_f64(double x, double h, double k, double loc, double scale);
double scifort_kappa4_cdf_f64(double x, double h, double k, double loc, double scale);
double scifort_kappa4_ppf_f64(double p, double h, double k, double loc, double scale);

/* Doubly truncated Weibull-minimum distribution. */
double scifort_truncweibull_min_pdf_f64(double x, double c, double a, double b, double loc, double scale);
double scifort_truncweibull_min_cdf_f64(double x, double c, double a, double b, double loc, double scale);
double scifort_truncweibull_min_ppf_f64(double p, double c, double a, double b, double loc, double scale);


/* Generalized gamma distribution. */
double scifort_gengamma_pdf_f64(double x, double a, double c, double loc, double scale);
double scifort_gengamma_cdf_f64(double x, double a, double c, double loc, double scale);
double scifort_gengamma_ppf_f64(double p, double a, double c, double loc, double scale);

/* Half generalized normal distribution. */
double scifort_halfgennorm_pdf_f64(double x, double beta, double loc, double scale);
double scifort_halfgennorm_cdf_f64(double x, double beta, double loc, double scale);
double scifort_halfgennorm_ppf_f64(double p, double beta, double loc, double scale);

/* ARGUS distribution. */
double scifort_argus_pdf_f64(double x, double chi, double loc, double scale);
double scifort_argus_cdf_f64(double x, double chi, double loc, double scale);
double scifort_argus_ppf_f64(double p, double chi, double loc, double scale);

/* Erlang distribution. */
double scifort_erlang_pdf_f64(double x, double a, double loc, double scale);
double scifort_erlang_cdf_f64(double x, double a, double loc, double scale);
double scifort_erlang_ppf_f64(double p, double a, double loc, double scale);


/* Crystal Ball distribution. */
double scifort_crystalball_pdf_f64(double x, double beta, double m, double loc, double scale);
double scifort_crystalball_cdf_f64(double x, double beta, double m, double loc, double scale);
double scifort_crystalball_ppf_f64(double p, double beta, double m, double loc, double scale);

/* Jones-Faddy skew-t distribution. */
double scifort_jf_skew_t_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_jf_skew_t_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_jf_skew_t_ppf_f64(double p, double a, double b, double loc, double scale);

/* Pearson type III distribution. */
double scifort_pearson3_pdf_f64(double x, double skew, double loc, double scale);
double scifort_pearson3_cdf_f64(double x, double skew, double loc, double scale);
double scifort_pearson3_ppf_f64(double p, double skew, double loc, double scale);

/* Relativistic Breit-Wigner distribution. */
double scifort_rel_breitwigner_pdf_f64(double x, double rho, double loc, double scale);
double scifort_rel_breitwigner_cdf_f64(double x, double rho, double loc, double scale);
double scifort_rel_breitwigner_ppf_f64(double p, double rho, double loc, double scale);

/* Generalized exponential distribution. */
double scifort_genexpon_pdf_f64(double x, double a, double b, double c, double loc, double scale);
double scifort_genexpon_cdf_f64(double x, double a, double b, double c, double loc, double scale);
double scifort_genexpon_ppf_f64(double p, double a, double b, double c, double loc, double scale);

/* Skew-normal distribution. */
double scifort_skewnorm_pdf_f64(double x, double a, double loc, double scale);
double scifort_skewnorm_cdf_f64(double x, double a, double loc, double scale);
double scifort_skewnorm_ppf_f64(double p, double a, double loc, double scale);

/* Tukey lambda distribution. */
double scifort_tukeylambda_pdf_f64(double x, double lam, double loc, double scale);
double scifort_tukeylambda_cdf_f64(double x, double lam, double loc, double scale);
double scifort_tukeylambda_ppf_f64(double p, double lam, double loc, double scale);

/* Rice distribution. */
double scifort_rice_pdf_f64(double x, double b, double loc, double scale);
double scifort_rice_cdf_f64(double x, double b, double loc, double scale);
double scifort_rice_ppf_f64(double p, double b, double loc, double scale);

/* Double-Pareto lognormal distribution. */
double scifort_dpareto_lognorm_pdf_f64(double x, double u, double s, double a, double b, double loc, double scale);
double scifort_dpareto_lognorm_cdf_f64(double x, double u, double s, double a, double b, double loc, double scale);
double scifort_dpareto_lognorm_ppf_f64(double p, double u, double s, double a, double b, double loc, double scale);

/* Circular von Mises distribution. */
double scifort_vonmises_pdf_f64(double x, double kappa, double loc, double scale);
double scifort_vonmises_cdf_f64(double x, double kappa, double loc, double scale);
double scifort_vonmises_ppf_f64(double p, double kappa, double loc, double scale);

/* Von Mises distribution on [-pi, pi]. */
double scifort_vonmises_line_pdf_f64(double x, double kappa, double loc, double scale);
double scifort_vonmises_line_cdf_f64(double x, double kappa, double loc, double scale);
double scifort_vonmises_line_ppf_f64(double p, double kappa, double loc, double scale);

/* Asymptotic two-sided Kolmogorov distribution. */
double scifort_kstwobign_pdf_f64(double x, double loc, double scale);
double scifort_kstwobign_cdf_f64(double x, double loc, double scale);
double scifort_kstwobign_ppf_f64(double p, double loc, double scale);

/* Irwin-Hall distribution. */
double scifort_irwinhall_pdf_f64(double x, double n, double loc, double scale);
double scifort_irwinhall_cdf_f64(double x, double n, double loc, double scale);
double scifort_irwinhall_ppf_f64(double p, double n, double loc, double scale);

/* One-sided finite-sample Kolmogorov-Smirnov distribution. */
double scifort_ksone_pdf_f64(double x, double n, double loc, double scale);
double scifort_ksone_cdf_f64(double x, double n, double loc, double scale);
double scifort_ksone_ppf_f64(double p, double n, double loc, double scale);
double scifort_kstwo_pdf_f64(double x, double n, double loc, double scale);
double scifort_kstwo_cdf_f64(double x, double n, double loc, double scale);
double scifort_kstwo_ppf_f64(double p, double n, double loc, double scale);

/* Noncentral chi-square distribution. */
double scifort_ncx2_pdf_f64(double x, double df, double nc, double loc, double scale);
double scifort_ncx2_cdf_f64(double x, double df, double nc, double loc, double scale);
double scifort_ncx2_ppf_f64(double p, double df, double nc, double loc, double scale);

/* Noncentral F distribution. */
double scifort_ncf_pdf_f64(double x, double dfn, double dfd, double nc, double loc, double scale);
double scifort_ncf_cdf_f64(double x, double dfn, double dfd, double nc, double loc, double scale);
double scifort_ncf_ppf_f64(double p, double dfn, double dfd, double nc, double loc, double scale);

/* Discrete uniform distribution on low, ..., high - 1. */
double scifort_randint_pmf_f64(double k, double low, double high, double loc);
double scifort_randint_cdf_f64(double k, double low, double high, double loc);
double scifort_randint_ppf_f64(double p, double low, double high, double loc);

/* Planck discrete exponential distribution. */
double scifort_planck_pmf_f64(double k, double lambda, double loc);
double scifort_planck_cdf_f64(double k, double lambda, double loc);
double scifort_planck_ppf_f64(double p, double lambda, double loc);

/* Discrete Laplace distribution. */
double scifort_dlaplace_pmf_f64(double k, double a, double loc);
double scifort_dlaplace_cdf_f64(double k, double a, double loc);
double scifort_dlaplace_ppf_f64(double p, double a, double loc);

/* Logarithmic-series distribution. */
double scifort_logser_pmf_f64(double k, double p, double loc);
double scifort_logser_cdf_f64(double k, double p, double loc);
double scifort_logser_ppf_f64(double q, double p, double loc);

/* Beta-binomial distribution. */
double scifort_betabinom_pmf_f64(double k, double n, double a, double b, double loc);
double scifort_betabinom_cdf_f64(double k, double n, double a, double b, double loc);
double scifort_betabinom_ppf_f64(double q, double n, double a, double b, double loc);

/* Hypergeometric distribution. */
double scifort_hypergeom_pmf_f64(double k, double M, double n, double N, double loc);
double scifort_hypergeom_cdf_f64(double k, double M, double n, double N, double loc);
double scifort_hypergeom_ppf_f64(double q, double M, double n, double N, double loc);

/* Negative hypergeometric distribution. */
double scifort_nhypergeom_pmf_f64(double k, double M, double n, double r, double loc);
double scifort_nhypergeom_cdf_f64(double k, double M, double n, double r, double loc);
double scifort_nhypergeom_ppf_f64(double q, double M, double n, double r, double loc);

/* Boltzmann truncated discrete exponential distribution. */
double scifort_boltzmann_pmf_f64(double k, double lambda, double N, double loc);
double scifort_boltzmann_cdf_f64(double k, double lambda, double N, double loc);
double scifort_boltzmann_ppf_f64(double q, double lambda, double N, double loc);

double scifort_betanbinom_pmf_f64(double k, double n, double a, double b, double loc);
double scifort_betanbinom_cdf_f64(double k, double n, double a, double b, double loc);
double scifort_betanbinom_ppf_f64(double q, double n, double a, double b, double loc);
double scifort_yulesimon_pmf_f64(double k, double alpha, double loc);
double scifort_yulesimon_cdf_f64(double k, double alpha, double loc);
double scifort_yulesimon_ppf_f64(double q, double alpha, double loc);
double scifort_zipf_pmf_f64(double k, double a, double loc);
double scifort_zipf_cdf_f64(double k, double a, double loc);
double scifort_zipf_ppf_f64(double q, double a, double loc);
double scifort_zipfian_pmf_f64(double k, double a, double n, double loc);
double scifort_zipfian_cdf_f64(double k, double a, double n, double loc);
double scifort_zipfian_ppf_f64(double q, double a, double n, double loc);

/* Generalized inverse Gaussian distribution. */
double scifort_geninvgauss_pdf_f64(double x, double p, double b, double loc, double scale);
double scifort_geninvgauss_cdf_f64(double x, double p, double b, double loc, double scale);
double scifort_geninvgauss_ppf_f64(double q, double p, double b, double loc, double scale);

/* Normal inverse Gaussian distribution. */
double scifort_norminvgauss_pdf_f64(double x, double a, double b, double loc, double scale);
double scifort_norminvgauss_cdf_f64(double x, double a, double b, double loc, double scale);
double scifort_norminvgauss_ppf_f64(double q, double a, double b, double loc, double scale);

/* Skellam distribution. */
double scifort_skellam_pmf_f64(double k, double mu1, double mu2, double loc);
double scifort_skellam_cdf_f64(double k, double mu1, double mu2, double loc);
double scifort_skellam_ppf_f64(double q, double mu1, double mu2, double loc);

/* Generalized hyperbolic distribution. */
double scifort_genhyperbolic_pdf_f64(double x, double p, double a, double b, double loc, double scale);
double scifort_genhyperbolic_cdf_f64(double x, double p, double a, double b, double loc, double scale);
double scifort_genhyperbolic_ppf_f64(double q, double p, double a, double b, double loc, double scale);

/* Fisher noncentral hypergeometric distribution. */
double scifort_nchypergeom_fisher_pmf_f64(double k, double M, double n, double N, double odds, double loc);
double scifort_nchypergeom_fisher_cdf_f64(double k, double M, double n, double N, double odds, double loc);
double scifort_nchypergeom_fisher_ppf_f64(double q, double M, double n, double N, double odds, double loc);

/* Noncentral Student t distribution. */
double scifort_nct_pdf_f64(double x, double df, double nc, double loc, double scale);
double scifort_nct_cdf_f64(double x, double df, double nc, double loc, double scale);
double scifort_nct_ppf_f64(double q, double df, double nc, double loc, double scale);

/* Gauss hypergeometric distribution. */
double scifort_gausshyper_pdf_f64(double x, double a, double b, double c, double z, double loc, double scale);
double scifort_gausshyper_cdf_f64(double x, double a, double b, double c, double z, double loc, double scale);
double scifort_gausshyper_ppf_f64(double q, double a, double b, double c, double z, double loc, double scale);
double scifort_landau_pdf_f64(double x, double loc, double scale);
double scifort_landau_cdf_f64(double x, double loc, double scale);
double scifort_landau_ppf_f64(double q, double loc, double scale);


/* Levy-stable distribution (SciPy default S1 parameterization). */
double scifort_levy_stable_pdf_f64(double x, double alpha, double beta, double loc, double scale);
double scifort_levy_stable_cdf_f64(double x, double alpha, double beta, double loc, double scale);
double scifort_levy_stable_ppf_f64(double q, double alpha, double beta, double loc, double scale);

/* Studentized-range distribution. */
double scifort_studentized_range_pdf_f64(double x, double k, double df, double loc, double scale);
double scifort_studentized_range_cdf_f64(double x, double k, double df, double loc, double scale);
double scifort_studentized_range_ppf_f64(double q, double k, double df, double loc, double scale);

/* Wallenius noncentral hypergeometric distribution. */
double scifort_nchypergeom_wallenius_pmf_f64(double k, double M, double n, double N, double odds, double loc);
double scifort_nchypergeom_wallenius_cdf_f64(double k, double M, double n, double N, double odds, double loc);
double scifort_nchypergeom_wallenius_ppf_f64(double q, double M, double n, double N, double odds, double loc);

/* Poisson-binomial distribution. */
double scifort_poisson_binom_pmf_f64(double k, const double *p, size_t n, double loc);
double scifort_poisson_binom_cdf_f64(double k, const double *p, size_t n, double loc);
double scifort_poisson_binom_ppf_f64(double q, const double *p, size_t n, double loc);

#ifdef __cplusplus
}
#endif

#endif
