! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Beta distribution with shapes a, b > 0, location loc, and scale > 0 on
! [loc, loc + scale]: f(x) = z**(a - 1) (1 - z)**(b - 1) / (B(a, b) scale),
! z = (x - loc) / scale. Matches scipy.stats.beta.

module scifort_beta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : beta_inverse_xy, incomplete_beta_xy, &
        log_beta_kernel
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta, stirling_threshold
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: beta_cdf
    public :: beta_isf
    public :: beta_logcdf
    public :: beta_logpdf
    public :: beta_logsf
    public :: beta_pdf
    public :: beta_ppf
    public :: beta_sf

contains

    pure elemental function beta_pdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > 1.0_dp) then
            y = 0.0_dp
        else if (z == 0.0_dp) then
            y = endpoint_density(a, b) / sigma
        else if (z == 1.0_dp) then
            y = endpoint_density(b, a) / sigma
        else
            y = exp(standard_logpdf(a, b, z, 1.0_dp - z)) / sigma
        end if
    end function beta_pdf

    pure elemental function beta_logpdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > 1.0_dp) then
            y = negative_infinity(z)
        else if (z == 0.0_dp) then
            y = log_nonnegative(endpoint_density(a, b)) - log(sigma)
        else if (z == 1.0_dp) then
            y = log_nonnegative(endpoint_density(b, a)) - log(sigma)
        else
            y = standard_logpdf(a, b, z, 1.0_dp - z) - log(sigma)
        end if
    end function beta_logpdf

    pure elemental function beta_cdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q
        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = 1.0_dp
        else
            call incomplete_beta_xy(a, b, z, 1.0_dp - z, p, q, logp, logq)
            y = p
        end if
    end function beta_cdf

    pure elemental function beta_sf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q
        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 1.0_dp
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            call incomplete_beta_xy(a, b, z, 1.0_dp - z, p, q, logp, logq)
            y = q
        end if
    end function beta_sf

    pure elemental function beta_logcdf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q
        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = negative_infinity(z)
        else if (z >= 1.0_dp) then
            y = 0.0_dp
        else
            call incomplete_beta_xy(a, b, z, 1.0_dp - z, p, q, logp, logq)
            y = logp
        end if
    end function beta_logcdf

    pure elemental function beta_logsf(x, a, b, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q
        real(dp) :: z
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z <= 0.0_dp) then
            y = 0.0_dp
        else if (z >= 1.0_dp) then
            y = negative_infinity(z)
        else
            call incomplete_beta_xy(a, b, z, 1.0_dp - z, p, q, logp, logq)
            y = logq
        end if
    end function beta_logsf

    pure elemental function beta_ppf(p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: u
        real(dp) :: v
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (p <= 0.0_dp) then
            y = mu
        else if (p >= 1.0_dp) then
            y = mu + sigma
        else
            call beta_inverse_xy(a, b, p, 1.0_dp - p, u, v)
            y = mu + sigma * u
        end if
    end function beta_ppf

    pure elemental function beta_isf(p, a, b, loc, scale) result(y)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! width of the support, > 0 (default 1)
        real(dp) :: y

        real(dp) :: u
        real(dp) :: v
        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(a) .and. valid_shape(b))) then
            y = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (p <= 0.0_dp) then
            y = mu + sigma
        else if (p >= 1.0_dp) then
            y = mu
        else
            call beta_inverse_xy(a, b, 1.0_dp - p, p, u, v)
            y = mu + sigma * u
        end if
    end function beta_isf

    ! log of the standard density for 0 < z < 1 with w = 1 - z. When both
    ! shapes are large, the density is formed from the stable kernel
    ! z**a w**b / B(a, b) of scifort_incomplete_beta.
    pure elemental function standard_logpdf(a, b, z, w) result(y)
        real(dp), intent(in) :: a !! first shape, finite and > 0
        real(dp), intent(in) :: b !! second shape, finite and > 0
        real(dp), intent(in) :: z !! standardized value, 0 < z < 1
        real(dp), intent(in) :: w !! 1 - z
        real(dp) :: y

        if (min(a, b) >= stirling_threshold) then
            y = log_beta_kernel(a, b, z, w) - log(z) - log(w)
        else
            y = (a - 1.0_dp) * log(z) + (b - 1.0_dp) * log1p_safe(-z) - log_beta(a, b)
        end if
    end function standard_logpdf

    ! Standard density at the endpoint where the first shape applies:
    ! infinite for s < 1, t for s = 1 (since 1 / B(1, t) = t), zero above.
    pure elemental function endpoint_density(s, t) result(y)
        real(dp), intent(in) :: s !! shape attached to the endpoint
        real(dp), intent(in) :: t !! other shape
        real(dp) :: y

        if (s < 1.0_dp) then
            y = positive_infinity(s)
        else if (s == 1.0_dp) then
            y = t
        else
            y = 0.0_dp
        end if
    end function endpoint_density

    pure elemental function log_nonnegative(p) result(y)
        real(dp), intent(in) :: p !! value, >= 0, to take the log of
        real(dp) :: y

        if (p > 0.0_dp) then
            y = log(p)
        else
            y = negative_infinity(p)
        end if
    end function log_nonnegative

    pure elemental logical function valid_shape(a) result(valid)
        real(dp), intent(in) :: a !! shape to check

        valid = ieee_is_finite(a) .and. a > 0.0_dp
    end function valid_shape

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

end module scifort_beta
