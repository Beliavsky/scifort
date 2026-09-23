! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized Pareto distribution with finite shape c. In standardized
! coordinates z = (x - loc) / scale, f(z) = (1 + c*z)**(-1 - 1/c),
! with the c = 0 limit equal to the exponential distribution.

module scifort_genpareto
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: genpareto_cdf
    public :: genpareto_isf
    public :: genpareto_logcdf
    public :: genpareto_logpdf
    public :: genpareto_logsf
    public :: genpareto_pdf
    public :: genpareto_ppf
    public :: genpareto_sf

contains

    pure elemental function genpareto_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: exponent
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 0.0_dp) then
            y = negative_infinity(x)
            return
        end if
        if (c == 0.0_dp) then
            y = -z - log(sigma)
            return
        end if

        t = 1.0_dp + c * z
        if (t < 0.0_dp) then
            y = negative_infinity(x)
        else if (t == 0.0_dp) then
            exponent = -1.0_dp - 1.0_dp / c
            if (exponent > 0.0_dp) then
                y = negative_infinity(x)
            else if (exponent < 0.0_dp) then
                y = positive_infinity(x)
            else
                y = -log(sigma)
            end if
        else
            y = -(1.0_dp + 1.0_dp / c) * log1p_safe(c * z) - log(sigma)
        end if
    end function genpareto_logpdf

    pure elemental function genpareto_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logp

        logp = genpareto_logpdf(x, c, loc, scale)
        if (ieee_is_nan(logp)) then
            y = logp
        else if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else if (logp == positive_infinity(logp)) then
            y = positive_infinity(logp)
        else if (logp > log(huge(1.0_dp))) then
            y = positive_infinity(logp)
        else
            y = exp(logp)
        end if
    end function genpareto_pdf

    pure elemental function genpareto_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 0.0_dp) then
            y = 0.0_dp
        else if (c == 0.0_dp) then
            y = -z
        else
            t = 1.0_dp + c * z
            if (t <= 0.0_dp) then
                y = negative_infinity(x)
            else
                y = -log1p_safe(c * z) / c
            end if
        end if
    end function genpareto_logsf

    pure elemental function genpareto_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = genpareto_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == negative_infinity(logsf)) then
            y = 0.0_dp
        else
            y = exp(logsf)
        end if
    end function genpareto_sf

    pure elemental function genpareto_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = genpareto_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == 0.0_dp) then
            y = 0.0_dp
        else if (logsf == negative_infinity(logsf)) then
            y = 1.0_dp
        else
            y = -expm1_safe(logsf)
        end if
    end function genpareto_cdf

    pure elemental function genpareto_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: logsf

        logsf = genpareto_logsf(x, c, loc, scale)
        if (ieee_is_nan(logsf)) then
            y = logsf
        else if (logsf == 0.0_dp) then
            y = negative_infinity(x)
        else if (logsf == negative_infinity(logsf)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(logsf))
        end if
    end function genpareto_logcdf

    pure elemental function genpareto_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            if (c < 0.0_dp) then
                x = mu - sigma / c
            else
                x = positive_infinity(p)
            end if
        else if (c == 0.0_dp) then
            x = mu - sigma * log1p_safe(-p)
        else
            q = -c * log1p_safe(-p)
            if (c > 0.0_dp .and. q > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = expm1_safe(q) / c
                x = mu + sigma * z
            end if
        end if
    end function genpareto_ppf

    pure elemental function genpareto_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            if (c < 0.0_dp) then
                x = mu - sigma / c
            else
                x = positive_infinity(p)
            end if
        else if (c == 0.0_dp) then
            x = mu - sigma * log(p)
        else
            q = -c * log(p)
            if (c > 0.0_dp .and. q > log(huge(1.0_dp))) then
                x = positive_infinity(p)
            else
                z = expm1_safe(q) / c
                x = mu + sigma * z
            end if
        end if
    end function genpareto_isf

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_genpareto
