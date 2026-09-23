! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Erlang distribution matching scipy.stats.erlang computational behavior.
! SciPy warns for noninteger shape values but evaluates them as gamma values;
! this numerical layer likewise accepts every finite a > 0.
module scifort_erlang
    use scifort_gamma, only : gamma_cdf, gamma_isf, gamma_logcdf, gamma_logpdf, &
        gamma_logsf, gamma_pdf, gamma_ppf, gamma_sf
    use scifort_kinds, only : dp
    implicit none
    private

    public :: erlang_cdf, erlang_isf, erlang_logcdf, erlang_logpdf
    public :: erlang_logsf, erlang_pdf, erlang_ppf, erlang_sf

contains

    pure elemental function erlang_pdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_pdf(x, a, loc, scale)
    end function erlang_pdf

    pure elemental function erlang_logpdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_logpdf(x, a, loc, scale)
    end function erlang_logpdf

    pure elemental function erlang_cdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_cdf(x, a, loc, scale)
    end function erlang_cdf

    pure elemental function erlang_sf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_sf(x, a, loc, scale)
    end function erlang_sf

    pure elemental function erlang_logcdf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_logcdf(x, a, loc, scale)
    end function erlang_logcdf

    pure elemental function erlang_logsf(x, a, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        y = gamma_logsf(x, a, loc, scale)
    end function erlang_logsf

    pure elemental function erlang_ppf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = gamma_ppf(p, a, loc, scale)
    end function erlang_ppf

    pure elemental function erlang_isf(p, a, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        x = gamma_isf(p, a, loc, scale)
    end function erlang_isf

end module scifort_erlang
