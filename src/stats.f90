! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_stats
    use scifort_beta, only : beta_cdf, beta_isf, beta_logcdf, beta_logpdf, &
        beta_logsf, beta_pdf, beta_ppf, beta_sf
    use scifort_binomial, only : binomial_cdf, binomial_isf, binomial_logcdf, binomial_logpmf, &
        binomial_logsf, binomial_pmf, binomial_ppf, binomial_sf
    use scifort_cauchy, only : cauchy_cdf, cauchy_isf, cauchy_logcdf, &
        cauchy_logpdf, cauchy_logsf, cauchy_pdf, cauchy_ppf, cauchy_sf
    use scifort_chi2, only : chi2_cdf, chi2_isf, chi2_logcdf, chi2_logpdf, &
        chi2_logsf, chi2_pdf, chi2_ppf, chi2_sf
    use scifort_exponential, only : exponential_cdf, exponential_isf, &
        exponential_logcdf, exponential_logpdf, exponential_logsf, &
        exponential_pdf, exponential_ppf, exponential_sf
    use scifort_f_distribution, only : f_cdf, f_isf, f_logcdf, f_logpdf, &
        f_logsf, f_pdf, f_ppf, f_sf
    use scifort_gamma, only : gamma_cdf, gamma_isf, gamma_logcdf, gamma_logpdf, &
        gamma_logsf, gamma_pdf, gamma_ppf, gamma_sf
    use scifort_kinds, only : dp
    use scifort_laplace, only : laplace_cdf, laplace_isf, laplace_logcdf, &
        laplace_logpdf, laplace_logsf, laplace_pdf, laplace_ppf, laplace_sf
    use scifort_logistic, only : logistic_cdf, logistic_isf, logistic_logcdf, &
        logistic_logpdf, logistic_logsf, logistic_pdf, logistic_ppf, logistic_sf
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf
    use scifort_poisson, only : poisson_cdf, poisson_isf, poisson_logcdf, poisson_logpmf, &
        poisson_logsf, poisson_pmf, poisson_ppf, poisson_sf
    use scifort_student_t, only : t_cdf, t_isf, t_logcdf, t_logpdf, &
        t_logsf, t_pdf, t_ppf, t_sf
    use scifort_uniform, only : uniform_cdf, uniform_isf, uniform_logcdf, &
        uniform_logpdf, uniform_logsf, uniform_pdf, uniform_ppf, uniform_sf
    implicit none
    private

    public :: beta_cdf
    public :: beta_isf
    public :: beta_logcdf
    public :: beta_logpdf
    public :: beta_logsf
    public :: beta_pdf
    public :: beta_ppf
    public :: beta_sf
    public :: binomial_cdf
    public :: binomial_isf
    public :: binomial_logcdf
    public :: binomial_logpmf
    public :: binomial_logsf
    public :: binomial_pmf
    public :: binomial_ppf
    public :: binomial_sf
    public :: cauchy_cdf
    public :: cauchy_isf
    public :: cauchy_logcdf
    public :: cauchy_logpdf
    public :: cauchy_logsf
    public :: cauchy_pdf
    public :: cauchy_ppf
    public :: cauchy_sf
    public :: chi2_cdf
    public :: chi2_isf
    public :: chi2_logcdf
    public :: chi2_logpdf
    public :: chi2_logsf
    public :: chi2_pdf
    public :: chi2_ppf
    public :: chi2_sf
    public :: dp
    public :: exponential_cdf
    public :: exponential_isf
    public :: exponential_logcdf
    public :: exponential_logpdf
    public :: exponential_logsf
    public :: exponential_pdf
    public :: exponential_ppf
    public :: exponential_sf
    public :: f_cdf
    public :: f_isf
    public :: f_logcdf
    public :: f_logpdf
    public :: f_logsf
    public :: f_pdf
    public :: f_ppf
    public :: f_sf
    public :: gamma_cdf
    public :: gamma_isf
    public :: gamma_logcdf
    public :: gamma_logpdf
    public :: gamma_logsf
    public :: gamma_pdf
    public :: gamma_ppf
    public :: gamma_sf
    public :: laplace_cdf
    public :: laplace_isf
    public :: laplace_logcdf
    public :: laplace_logpdf
    public :: laplace_logsf
    public :: laplace_pdf
    public :: laplace_ppf
    public :: laplace_sf
    public :: logistic_cdf
    public :: logistic_isf
    public :: logistic_logcdf
    public :: logistic_logpdf
    public :: logistic_logsf
    public :: logistic_pdf
    public :: logistic_ppf
    public :: logistic_sf
    public :: normal_cdf
    public :: normal_isf
    public :: normal_logcdf
    public :: normal_logpdf
    public :: normal_logsf
    public :: normal_pdf
    public :: normal_ppf
    public :: normal_sf
    public :: poisson_cdf
    public :: poisson_isf
    public :: poisson_logcdf
    public :: poisson_logpmf
    public :: poisson_logsf
    public :: poisson_pmf
    public :: poisson_ppf
    public :: poisson_sf
    public :: t_cdf
    public :: t_isf
    public :: t_logcdf
    public :: t_logpdf
    public :: t_logsf
    public :: t_pdf
    public :: t_ppf
    public :: t_sf
    public :: uniform_cdf
    public :: uniform_isf
    public :: uniform_logcdf
    public :: uniform_logpdf
    public :: uniform_logsf
    public :: uniform_pdf
    public :: uniform_ppf
    public :: uniform_sf
end module scifort_stats
