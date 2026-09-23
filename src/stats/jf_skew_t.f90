! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Jones-Faddy skew-t distribution matching scipy.stats.jf_skew_t.
module scifort_jf_skew_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_incomplete_beta, only : betainccinv, betaincinv, incomplete_beta_xy
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: jf_skew_t_cdf, jf_skew_t_isf, jf_skew_t_logcdf
    public :: jf_skew_t_logpdf, jf_skew_t_logsf, jf_skew_t_pdf
    public :: jf_skew_t_ppf, jf_skew_t_sf

contains

    pure elemental function jf_skew_t_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logt, logu, mu, s, sigma, t, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            s = a + b
            call beta_coordinates(z, s, t, u, logt, logu)
            y = (a + 0.5_dp) * logt + (b + 0.5_dp) * logu + &
                2.0_dp * scifort_log_two - log_beta(a, b) - &
                0.5_dp * log(s) - log(sigma)
        end if
    end function jf_skew_t_logpdf

    pure elemental function jf_skew_t_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = jf_skew_t_logpdf(x, a, b, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function jf_skew_t_pdf

    pure elemental subroutine probabilities(x, a, b, loc, scale, p, q, logp, logq)
        real(dp), intent(in) :: x !! point at which probabilities are evaluated
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp), intent(out) :: p !! lower-tail probability
        real(dp), intent(out) :: q !! upper-tail probability
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: logt, logu, mu, s, sigma, t, u, nanv, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
                ieee_is_nan(x)) then
            nanv = quiet_nan(x)
            p = nanv
            q = nanv
            logp = nanv
            logq = nanv
        else if (x == negative_infinity(x)) then
            p = 0.0_dp
            q = 1.0_dp
            logp = negative_infinity(x)
            logq = 0.0_dp
        else if (x == positive_infinity(x)) then
            p = 1.0_dp
            q = 0.0_dp
            logp = 0.0_dp
            logq = negative_infinity(x)
        else
            s = a + b
            z = (x - mu) / sigma
            call beta_coordinates(z, s, t, u, logt, logu)
            call incomplete_beta_xy(a, b, t, u, p, q, logp, logq)
        end if
    end subroutine probabilities

    pure elemental function jf_skew_t_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, q, lp, lq
        call probabilities(x, a, b, loc, scale, y, q, lp, lq)
    end function jf_skew_t_cdf

    pure elemental function jf_skew_t_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, lp, lq
        call probabilities(x, a, b, loc, scale, p, y, lp, lq)
    end function jf_skew_t_sf

    pure elemental function jf_skew_t_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q, lq
        call probabilities(x, a, b, loc, scale, p, q, y, lq)
    end function jf_skew_t_logcdf

    pure elemental function jf_skew_t_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p, q, lp
        call probabilities(x, a, b, loc, scale, p, q, lp, y)
    end function jf_skew_t_logsf

    pure elemental function jf_skew_t_ppf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, s, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = negative_infinity(p)
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            s = a + b
            q = betaincinv(a, b, p)
            z = (2.0_dp * q - 1.0_dp) * sqrt(s) / &
                (2.0_dp * sqrt(q) * sqrt(1.0_dp - q))
            x = mu + sigma * z
        end if
    end function jf_skew_t_ppf

    pure elemental function jf_skew_t_isf(p, a, b, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, q, s, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shape(a, b) .or. .not. valid_loc_scale(mu, sigma) .or. &
                .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = negative_infinity(p)
        else
            s = a + b
            q = betainccinv(a, b, p)
            z = (2.0_dp * q - 1.0_dp) * sqrt(s) / &
                (2.0_dp * sqrt(q) * sqrt(1.0_dp - q))
            x = mu + sigma * z
        end if
    end function jf_skew_t_isf

    pure elemental subroutine beta_coordinates(z, s, t, u, logt, logu)
        real(dp), intent(in) :: z !! standardized variate
        real(dp), intent(in) :: s !! positive sum a+b
        real(dp), intent(out) :: t !! lower incomplete-beta coordinate
        real(dp), intent(out) :: u !! accurately formed complement 1-t
        real(dp), intent(out) :: logt !! logarithm of t
        real(dp), intent(out) :: logu !! logarithm of u
        real(dp) :: az, h, ratio, root_s

        root_s = sqrt(s)
        az = abs(z)
        if (az == 0.0_dp) then
            t = 0.5_dp
            u = 0.5_dp
            logt = -scifort_log_two
            logu = -scifort_log_two
            return
        end if

        if (az > root_s) then
            ratio = root_s / az
            h = az * sqrt(1.0_dp + ratio * ratio)
        else
            ratio = az / root_s
            h = root_s * sqrt(1.0_dp + ratio * ratio)
        end if

        if (z > 0.0_dp) then
            logu = -scifort_log_two + log(s) - 2.0_dp * log(h) - &
                log1p_safe(z / h)
            u = exp(logu)
            t = 1.0_dp - u
            logt = log1p_safe(-u)
        else
            logt = -scifort_log_two + log(s) - 2.0_dp * log(h) - &
                log1p_safe(-z / h)
            t = exp(logt)
            u = 1.0_dp - t
            logu = log1p_safe(-t)
        end if
    end subroutine beta_coordinates

    pure elemental logical function valid_shape(a, b) result(valid)
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        valid = a > 0.0_dp .and. b > 0.0_dp .and. &
            ieee_is_finite(a) .and. ieee_is_finite(b)
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

end module scifort_jf_skew_t
