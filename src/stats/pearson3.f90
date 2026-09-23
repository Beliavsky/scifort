! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Pearson type III distribution matching scipy.stats.pearson3.
module scifort_pearson3
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_gamma, only : gamma_cdf, gamma_isf, gamma_logcdf, gamma_logpdf, &
        gamma_logsf, gamma_ppf, gamma_sf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_ppf, normal_sf
    implicit none
    private

    real(dp), parameter :: normal_transition = 1.6e-5_dp

    public :: pearson3_cdf, pearson3_isf, pearson3_logcdf, pearson3_logpdf
    public :: pearson3_logsf, pearson3_pdf, pearson3_ppf, pearson3_sf

contains

    pure elemental function pearson3_logpdf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (abs(skew) < normal_transition) then
            y = normal_logpdf(x, mu, sigma)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            beta = 2.0_dp / skew
            alpha = beta * beta
            t = alpha + beta * z
            y = log(abs(beta)) + gamma_logpdf(t, alpha) - log(sigma)
        end if
    end function pearson3_logpdf

    pure elemental function pearson3_pdf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = pearson3_logpdf(x, skew, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function pearson3_pdf

    pure elemental function pearson3_cdf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (abs(skew) < normal_transition) then
            y = normal_cdf(x, mu, sigma)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = 1.0_dp
        else
            z = (x - mu) / sigma
            beta = 2.0_dp / skew
            alpha = beta * beta
            t = alpha + beta * z
            if (beta > 0.0_dp) then
                y = gamma_cdf(t, alpha)
            else
                y = gamma_sf(t, alpha)
            end if
        end if
    end function pearson3_cdf

    pure elemental function pearson3_sf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (abs(skew) < normal_transition) then
            y = normal_sf(x, mu, sigma)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 1.0_dp
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            beta = 2.0_dp / skew
            alpha = beta * beta
            t = alpha + beta * z
            if (beta > 0.0_dp) then
                y = gamma_sf(t, alpha)
            else
                y = gamma_cdf(t, alpha)
            end if
        end if
    end function pearson3_sf

    pure elemental function pearson3_logcdf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (abs(skew) < normal_transition) then
            y = normal_logcdf(x, mu, sigma)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            beta = 2.0_dp / skew
            alpha = beta * beta
            t = alpha + beta * z
            if (beta > 0.0_dp) then
                y = gamma_logcdf(t, alpha)
            else
                y = gamma_logsf(t, alpha)
            end if
        end if
    end function pearson3_logcdf

    pure elemental function pearson3_logsf(x, skew, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (abs(skew) < normal_transition) then
            y = normal_logsf(x, mu, sigma)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            beta = 2.0_dp / skew
            alpha = beta * beta
            t = alpha + beta * z
            if (beta > 0.0_dp) then
                y = gamma_logsf(t, alpha)
            else
                y = gamma_logcdf(t, alpha)
            end if
        end if
    end function pearson3_logsf

    pure elemental function pearson3_ppf(p, skew, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (abs(skew) < normal_transition) then
            x = normal_ppf(p, mu, sigma)
        else if (p == 0.0_dp) then
            if (skew > 0.0_dp) then
                x = mu - 2.0_dp * sigma / skew
            else
                x = negative_infinity(p)
            end if
        else if (p == 1.0_dp) then
            if (skew < 0.0_dp) then
                x = mu - 2.0_dp * sigma / skew
            else
                x = positive_infinity(p)
            end if
        else
            beta = 2.0_dp / skew
            alpha = beta * beta
            if (beta > 0.0_dp) then
                t = gamma_ppf(p, alpha)
            else
                t = gamma_isf(p, alpha)
            end if
            z = t / beta - beta
            x = mu + sigma * z
        end if
    end function pearson3_ppf

    pure elemental function pearson3_isf(p, skew, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, alpha, beta, mu, sigma, t, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(skew) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (abs(skew) < normal_transition) then
            x = normal_isf(p, mu, sigma)
        else if (p == 0.0_dp) then
            if (skew < 0.0_dp) then
                x = mu - 2.0_dp * sigma / skew
            else
                x = positive_infinity(p)
            end if
        else if (p == 1.0_dp) then
            if (skew > 0.0_dp) then
                x = mu - 2.0_dp * sigma / skew
            else
                x = negative_infinity(p)
            end if
        else
            beta = 2.0_dp / skew
            alpha = beta * beta
            if (beta > 0.0_dp) then
                t = gamma_isf(p, alpha)
            else
                t = gamma_ppf(p, alpha)
            end if
            z = t / beta - beta
            x = mu + sigma * z
        end if
    end function pearson3_isf

    pure elemental logical function valid_shape(skew) result(valid)
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        valid = ieee_is_finite(skew)
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location, default 0
        real(dp), intent(out) :: sigma !! resolved scale, default 1
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_pearson3
