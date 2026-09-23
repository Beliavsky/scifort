! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Truncated normal distribution with finite standardized limits a < b.
! The normal law is truncated on [a,b] before applying loc and scale,
! matching scipy.stats.truncnorm for finite a and b.

module scifort_truncnorm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_logpdf, &
        normal_logsf, normal_pdf, normal_sf
    implicit none
    private

    public :: truncnorm_cdf
    public :: truncnorm_isf
    public :: truncnorm_logcdf
    public :: truncnorm_logpdf
    public :: truncnorm_logsf
    public :: truncnorm_pdf
    public :: truncnorm_ppf
    public :: truncnorm_sf

contains

    pure elemental function truncnorm_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y

        y = exp(truncnorm_logpdf(x, a, b, loc, scale))
    end function truncnorm_pdf

    pure elemental function truncnorm_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z < a .or. z > b .or. .not. ieee_is_finite(z)) then
            y = negative_infinity(x)
        else
            y = normal_logpdf(z) - log_gauss_mass(a, b) - log(sigma)
        end if
    end function truncnorm_logpdf

    pure elemental function truncnorm_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y
        y = exp(truncnorm_logcdf(x, a, b, loc, scale))
    end function truncnorm_cdf

    pure elemental function truncnorm_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y
        y = exp(truncnorm_logsf(x, a, b, loc, scale))
    end function truncnorm_sf

    pure elemental function truncnorm_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = negative_infinity(x)
            else if (z >= b) then
                y = 0.0_dp
            else
                y = log_gauss_mass(a, z) - log_gauss_mass(a, b)
            end if
        end if
    end function truncnorm_logcdf

    pure elemental function truncnorm_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= a) then
                y = 0.0_dp
            else if (z >= b) then
                y = negative_infinity(x)
            else
                y = log_gauss_mass(z, b) - log_gauss_mass(a, b)
            end if
        end if
    end function truncnorm_logsf

    pure elemental function truncnorm_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: x

        integer :: iteration
        real(dp) :: hi
        real(dp) :: logp
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            x = quiet_nan(p)
            return
        end if
        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        end if
        if (p <= 0.0_dp) then
            x = mu + sigma * a
            return
        else if (p >= 1.0_dp) then
            x = mu + sigma * b
            return
        end if
        logp = log(p)
        lo = a
        hi = b
        do iteration = 1, 90
            mid = lo + 0.5_dp * (hi - lo)
            if (truncnorm_logcdf(mid, a, b) < logp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = mu + sigma * (lo + 0.5_dp * (hi - lo))
    end function truncnorm_ppf

    pure elemental function truncnorm_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: x

        integer :: iteration
        real(dp) :: hi
        real(dp) :: logp
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shapes(a, b))) then
            x = quiet_nan(p)
            return
        end if
        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
            return
        end if
        if (p >= 1.0_dp) then
            x = mu + sigma * a
            return
        else if (p <= 0.0_dp) then
            x = mu + sigma * b
            return
        end if
        logp = log(p)
        lo = a
        hi = b
        do iteration = 1, 90
            mid = lo + 0.5_dp * (hi - lo)
            if (truncnorm_logsf(mid, a, b) > logp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = mu + sigma * (lo + 0.5_dp * (hi - lo))
    end function truncnorm_isf

    pure elemental function log_gauss_mass(left, right) result(value)
        real(dp), intent(in) :: left !! finite left endpoint
        real(dp), intent(in) :: right !! finite right endpoint greater than left
        real(dp) :: value

        real(dp) :: log_large
        real(dp) :: log_small

        if (right <= 0.0_dp) then
            log_large = normal_logcdf(right)
            log_small = normal_logcdf(left)
        else if (left >= 0.0_dp) then
            log_large = normal_logsf(left)
            log_small = normal_logsf(right)
        else
            value = log(normal_cdf(right) - normal_cdf(left))
            return
        end if
        value = log_large + log1p_safe(-exp(log_small - log_large))
    end function log_gauss_mass

    pure elemental logical function valid_shapes(a, b) result(valid)
        real(dp), intent(in) :: a !! lower truncation point to validate
        real(dp), intent(in) :: b !! upper truncation point to validate
        valid = ieee_is_finite(a) .and. ieee_is_finite(b) .and. a < b
    end function valid_shapes

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location argument
        real(dp), intent(in), optional :: scale !! optional scale argument
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_truncnorm
