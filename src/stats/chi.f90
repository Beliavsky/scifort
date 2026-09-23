! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Chi distribution with df > 0 degrees of freedom, location loc, and scale > 0.
! With z = (x - loc) / scale >= 0 and a = df / 2,
! f(x) = z**(df - 1) exp(-z**2 / 2) /
!        (2**(df / 2 - 1) Gamma(df / 2) scale).
! This matches scipy.stats.chi.

module scifort_chi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_sqrt_two
    use scifort_gamma, only : gamma_logpdf
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    implicit none
    private

    public :: chi_cdf
    public :: chi_isf
    public :: chi_logcdf
    public :: chi_logpdf
    public :: chi_logsf
    public :: chi_pdf
    public :: chi_ppf
    public :: chi_sf

contains

    pure elemental function chi_pdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: log_density

        log_density = chi_logpdf(x, df, loc, scale)
        if (ieee_is_nan(log_density)) then
            y = log_density
        else if (log_density == negative_infinity(log_density)) then
            y = 0.0_dp
        else if (log_density == positive_infinity(log_density)) then
            y = positive_infinity(log_density)
        else
            y = exp(log_density)
        end if
    end function chi_pdf

    pure elemental function chi_logpdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp) then
            y = negative_infinity(z)
        else if (z == 0.0_dp) then
            if (df < 1.0_dp) then
                y = positive_infinity(z)
            else if (df == 1.0_dp) then
                y = 0.5_dp * log(2.0_dp / acos(-1.0_dp)) - log(sigma)
            else
                y = negative_infinity(z)
            end if
        else if (z > max_standard_z()) then
            y = negative_infinity(z)
        else
            a = 0.5_dp * df
            t = (0.5_dp * z) * z
            y = log(z) + gamma_logpdf(t, a) - log(sigma)
        end if
    end function chi_logpdf

    pure elemental function chi_cdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z > max_standard_z()) then
            y = 1.0_dp
        else
            a = 0.5_dp * df
            t = (0.5_dp * z) * z
            y = gammainc(a, t)
        end if
    end function chi_cdf

    pure elemental function chi_sf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 1.0_dp
        else if (z > max_standard_z()) then
            y = 0.0_dp
        else
            a = 0.5_dp * df
            t = (0.5_dp * z) * z
            y = gammaincc(a, t)
        end if
    end function chi_sf

    pure elemental function chi_logcdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = negative_infinity(z)
        else if (z > max_standard_z()) then
            y = 0.0_dp
        else
            a = 0.5_dp * df
            t = (0.5_dp * z) * z
            y = log_gammainc(a, t)
        end if
    end function chi_logcdf

    pure elemental function chi_logsf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z > max_standard_z()) then
            y = negative_infinity(z)
        else
            a = 0.5_dp * df
            t = (0.5_dp * z) * z
            y = log_gammaincc(a, t)
        end if
    end function chi_logsf

    pure elemental function chi_ppf(p, df, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            t = gammaincinv(0.5_dp * df, p)
            z = scifort_sqrt_two * sqrt(t)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function chi_ppf

    pure elemental function chi_isf(p, df, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            t = gammainccinv(0.5_dp * df, p)
            z = scifort_sqrt_two * sqrt(t)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function chi_isf

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! finite positive scale parameter
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity when needed
        real(dp) :: x

        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive

    pure elemental logical function valid_df(df) result(valid)
        real(dp), intent(in) :: df !! degrees of freedom to check

        valid = ieee_is_finite(df) .and. df > 0.0_dp
    end function valid_df

    pure elemental function max_standard_z() result(zmax)
        real(dp) :: zmax

        zmax = scifort_sqrt_two * sqrt(huge(1.0_dp))
    end function max_standard_z

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

end module scifort_chi
