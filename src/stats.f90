! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_stats
    use scifort_cauchy, only : cauchy_cdf, cauchy_isf, cauchy_logcdf, &
        cauchy_logpdf, cauchy_logsf, cauchy_pdf, cauchy_ppf, cauchy_sf
    use scifort_exponential, only : exponential_cdf, exponential_isf, &
        exponential_logcdf, exponential_logpdf, exponential_logsf, &
        exponential_pdf, exponential_ppf, exponential_sf
    use scifort_kinds, only : dp
    use scifort_laplace, only : laplace_cdf, laplace_isf, laplace_logcdf, &
        laplace_logpdf, laplace_logsf, laplace_pdf, laplace_ppf, laplace_sf
    use scifort_logistic, only : logistic_cdf, logistic_isf, logistic_logcdf, &
        logistic_logpdf, logistic_logsf, logistic_pdf, logistic_ppf, logistic_sf
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf
    use scifort_uniform, only : uniform_cdf, uniform_isf, uniform_logcdf, &
        uniform_logpdf, uniform_logsf, uniform_pdf, uniform_ppf, uniform_sf
    implicit none
    private

    public :: cauchy_cdf
    public :: cauchy_isf
    public :: cauchy_logcdf
    public :: cauchy_logpdf
    public :: cauchy_logsf
    public :: cauchy_pdf
    public :: cauchy_ppf
    public :: cauchy_sf
    public :: dp
    public :: exponential_cdf
    public :: exponential_isf
    public :: exponential_logcdf
    public :: exponential_logpdf
    public :: exponential_logsf
    public :: exponential_pdf
    public :: exponential_ppf
    public :: exponential_sf
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
    public :: uniform_cdf
    public :: uniform_isf
    public :: uniform_logcdf
    public :: uniform_logpdf
    public :: uniform_logsf
    public :: uniform_pdf
    public :: uniform_ppf
    public :: uniform_sf
end module scifort_stats
