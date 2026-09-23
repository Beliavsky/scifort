! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Crystal Ball distribution matching scipy.stats.crystalball parameterization.
module scifort_crystalball
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logsf
    implicit none
    private

    public :: crystalball_cdf, crystalball_isf, crystalball_logcdf
    public :: crystalball_logpdf, crystalball_logsf, crystalball_pdf
    public :: crystalball_ppf, crystalball_sf

contains

    pure elemental function crystalball_logpdf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: d, logmnorm, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta, m) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            logmnorm = log_normalizer_denominator(beta, m)
            if (z > -beta) then
                y = -logmnorm - 0.5_dp * z * z - log(sigma)
            else
                d = m / beta - beta - z
                y = -logmnorm + m * (log(m) - log(beta)) - &
                    0.5_dp * beta * beta - m * log(d) - log(sigma)
            end if
        end if
    end function crystalball_logpdf

    pure elemental function crystalball_pdf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = crystalball_logpdf(x, beta, m, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function crystalball_pdf

    pure elemental function crystalball_logcdf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: d, logmnorm, logq, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta, m) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            logmnorm = log_normalizer_denominator(beta, m)
            if (z <= -beta) then
                d = m / beta - beta - z
                y = -logmnorm + m * (log(m) - log(beta)) - &
                    0.5_dp * beta * beta - log(m - 1.0_dp) + &
                    (1.0_dp - m) * log(d)
            else
                logq = scifort_log_sqrt_two_pi + normal_logsf(z) - logmnorm
                y = log1mexp_local(logq)
            end if
        end if
    end function crystalball_logcdf

    pure elemental function crystalball_logsf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logmnorm, logp, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta, m) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = 0.0_dp
        else if (x == positive_infinity(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            logmnorm = log_normalizer_denominator(beta, m)
            if (z > -beta) then
                y = scifort_log_sqrt_two_pi + normal_logsf(z) - logmnorm
            else
                logp = crystalball_logcdf(x, beta, m, mu, sigma)
                y = log1mexp_local(logp)
            end if
        end if
    end function crystalball_logsf

    pure elemental function crystalball_cdf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = crystalball_logcdf(x, beta, m, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function crystalball_cdf

    pure elemental function crystalball_sf(x, beta, m, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = crystalball_logsf(x, beta, m, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function crystalball_sf

    pure elemental function crystalball_ppf(p, beta, m, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: logmnorm, logpbeta, logd, mu, qn, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta, m) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            logmnorm = log_normalizer_denominator(beta, m)
            logpbeta = log_tail_integral(beta, m) - logmnorm
            if (log(p) < logpbeta) then
                logd = (log(p) + logmnorm - m * (log(m) - log(beta)) + &
                    0.5_dp * beta * beta + log(m - 1.0_dp)) / (1.0_dp - m)
                z = m / beta - beta - exp(logd)
            else
                qn = exp(log1p_safe(-p) + logmnorm - scifort_log_sqrt_two_pi)
                qn = min(1.0_dp, max(0.0_dp, qn))
                z = normal_isf(qn)
            end if
            x = mu + sigma * z
        end if
    end function crystalball_ppf

    pure elemental function crystalball_isf(p, beta, m, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: logmnorm, logqbeta, logd, mu, qn, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(beta, m) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else
            logmnorm = log_normalizer_denominator(beta, m)
            logqbeta = scifort_log_sqrt_two_pi + normal_logcdf(beta) - logmnorm
            if (log(p) <= logqbeta) then
                qn = exp(log(p) + logmnorm - scifort_log_sqrt_two_pi)
                qn = min(1.0_dp, max(0.0_dp, qn))
                z = normal_isf(qn)
            else
                logd = (log1p_safe(-p) + logmnorm - m * (log(m) - log(beta)) + &
                    0.5_dp * beta * beta + log(m - 1.0_dp)) / (1.0_dp - m)
                z = m / beta - beta - exp(logd)
            end if
            x = mu + sigma * z
        end if
    end function crystalball_isf

    pure elemental function log_normalizer_denominator(beta, m) result(y)
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp) :: y, a, b, hi

        a = log_tail_integral(beta, m)
        b = scifort_log_sqrt_two_pi + normal_logcdf(beta)
        hi = max(a, b)
        y = hi + log(exp(a - hi) + exp(b - hi))
    end function log_normalizer_denominator

    pure elemental function log_tail_integral(beta, m) result(y)
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp) :: y

        y = log(m) - log(beta) - log(m - 1.0_dp) - 0.5_dp * beta * beta
    end function log_tail_integral

    pure elemental function log1mexp_local(logp) result(y)
        real(dp), intent(in) :: logp !! logarithm of a probability in (-inf,0]
        real(dp) :: y

        if (logp == negative_infinity(logp)) then
            y = 0.0_dp
        else if (logp >= 0.0_dp) then
            y = negative_infinity(logp)
        else if (logp < -0.693147180559945309417232121458176568_dp) then
            y = log1p_safe(-exp(logp))
        else
            y = log(-expm1_safe(logp))
        end if
    end function log1mexp_local

    pure elemental logical function valid_shape(beta, m) result(valid)
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        valid = beta > 0.0_dp .and. m > 1.0_dp .and. &
            ieee_is_finite(beta) .and. ieee_is_finite(m)
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

end module scifort_crystalball
