! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Maxwell distribution with location loc and scale > 0. It is the chi
! distribution with three degrees of freedom and matches scipy.stats.maxwell.

module scifort_maxwell
    use scifort_chi, only : chi_cdf, chi_isf, chi_logcdf, chi_logpdf, &
        chi_logsf, chi_pdf, chi_ppf, chi_sf
    use scifort_kinds, only : dp
    implicit none
    private

    public :: maxwell_cdf
    public :: maxwell_isf
    public :: maxwell_logcdf
    public :: maxwell_logpdf
    public :: maxwell_logsf
    public :: maxwell_pdf
    public :: maxwell_ppf
    public :: maxwell_sf

contains

    pure elemental function maxwell_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_pdf(x, 3.0_dp, loc, scale)
    end function maxwell_pdf

    pure elemental function maxwell_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_logpdf(x, 3.0_dp, loc, scale)
    end function maxwell_logpdf

    pure elemental function maxwell_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_cdf(x, 3.0_dp, loc, scale)
    end function maxwell_cdf

    pure elemental function maxwell_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_sf(x, 3.0_dp, loc, scale)
    end function maxwell_sf

    pure elemental function maxwell_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_logcdf(x, 3.0_dp, loc, scale)
    end function maxwell_logcdf

    pure elemental function maxwell_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        y = chi_logsf(x, 3.0_dp, loc, scale)
    end function maxwell_logsf

    pure elemental function maxwell_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        x = chi_ppf(p, 3.0_dp, loc, scale)
    end function maxwell_ppf

    pure elemental function maxwell_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        x = chi_isf(p, 3.0_dp, loc, scale)
    end function maxwell_isf

end module scifort_maxwell
