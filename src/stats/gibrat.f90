! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Gibrat distribution matching scipy.stats.gibrat. This is the lognormal
! distribution with fixed log-shape s=1.
module scifort_gibrat
    use scifort_kinds, only : dp
    use scifort_lognormal, only : lognormal_cdf, lognormal_isf, lognormal_logcdf, &
        lognormal_logpdf, lognormal_logsf, lognormal_pdf, lognormal_ppf, lognormal_sf
    implicit none
    private

    public :: gibrat_cdf, gibrat_isf, gibrat_logcdf, gibrat_logpdf
    public :: gibrat_logsf, gibrat_pdf, gibrat_ppf, gibrat_sf

contains

    pure elemental function gibrat_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_pdf(x, 1.0_dp, loc, scale)
    end function gibrat_pdf

    pure elemental function gibrat_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_logpdf(x, 1.0_dp, loc, scale)
    end function gibrat_logpdf

    pure elemental function gibrat_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_cdf(x, 1.0_dp, loc, scale)
    end function gibrat_cdf

    pure elemental function gibrat_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_sf(x, 1.0_dp, loc, scale)
    end function gibrat_sf

    pure elemental function gibrat_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_logcdf(x, 1.0_dp, loc, scale)
    end function gibrat_logcdf

    pure elemental function gibrat_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = lognormal_logsf(x, 1.0_dp, loc, scale)
    end function gibrat_logsf

    pure elemental function gibrat_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = lognormal_ppf(p, 1.0_dp, loc, scale)
    end function gibrat_ppf

    pure elemental function gibrat_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = lognormal_isf(p, 1.0_dp, loc, scale)
    end function gibrat_isf

end module scifort_gibrat
