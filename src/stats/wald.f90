! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Wald distribution, the scipy.stats.invgauss special case with shape mu=1.
module scifort_wald
    use scifort_invgauss, only : invgauss_cdf, invgauss_isf, invgauss_logcdf, &
        invgauss_logpdf, invgauss_logsf, invgauss_pdf, invgauss_ppf, invgauss_sf
    use scifort_kinds, only : dp
    implicit none
    private
    public :: wald_pdf, wald_logpdf, wald_cdf, wald_sf
    public :: wald_logcdf, wald_logsf, wald_ppf, wald_isf
contains
    pure elemental function wald_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_pdf(x, 1.0_dp, loc, scale)
    end function wald_pdf
    pure elemental function wald_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_logpdf(x, 1.0_dp, loc, scale)
    end function wald_logpdf
    pure elemental function wald_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_cdf(x, 1.0_dp, loc, scale)
    end function wald_cdf
    pure elemental function wald_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_sf(x, 1.0_dp, loc, scale)
    end function wald_sf
    pure elemental function wald_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_logcdf(x, 1.0_dp, loc, scale)
    end function wald_logcdf
    pure elemental function wald_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = invgauss_logsf(x, 1.0_dp, loc, scale)
    end function wald_logsf
    pure elemental function wald_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = invgauss_ppf(p, 1.0_dp, loc, scale)
    end function wald_ppf
    pure elemental function wald_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = invgauss_isf(p, 1.0_dp, loc, scale)
    end function wald_isf
end module scifort_wald
