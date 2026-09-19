! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Student t distribution with df > 0 degrees of freedom, location loc, and
! scale > 0. Matches scipy.stats.t, including df = +infinity, which gives the
! normal distribution. For t > 0 the upper tail is
! (1/2) I_x(df/2, 1/2) with x = df / (df + t**2); x and 1 - x are formed
! separately so that neither tail loses accuracy.

module scifort_student_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_incomplete_beta, only : beta_inverse_xy, incomplete_beta_xy
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, normal_logpdf, &
        normal_logsf, normal_pdf, normal_ppf, normal_sf
    implicit none
    private

    public :: t_cdf
    public :: t_isf
    public :: t_logcdf
    public :: t_logpdf
    public :: t_logsf
    public :: t_pdf
    public :: t_ppf
    public :: t_sf

    ! Below this value of x = df / (df + t**2) the leading term of the tail
    ! is exact in binary64 and is used in logarithmic form so that x does not
    ! underflow.
    real(dp), parameter :: tiny_ratio = 1.0e-300_dp

contains

    pure elemental function t_pdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_pdf(x, mu, sigma)
        else
            y = exp(standard_logpdf((x - mu) / sigma, df)) / sigma
        end if
    end function t_pdf

    pure elemental function t_logpdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_logpdf(x, mu, sigma)
        else
            y = standard_logpdf((x - mu) / sigma, df) - log(sigma)
        end if
    end function t_logpdf

    pure elemental function t_cdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: far
        real(dp) :: log_far
        real(dp) :: log_near
        real(dp) :: mu
        real(dp) :: near
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_cdf(x, mu, sigma)
        else
            z = (x - mu) / sigma
            call standard_tails(z, df, far, near, log_far, log_near)
            if (z < 0.0_dp) then
                y = far
            else
                y = near
            end if
        end if
    end function t_cdf

    pure elemental function t_sf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: far
        real(dp) :: log_far
        real(dp) :: log_near
        real(dp) :: mu
        real(dp) :: near
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_sf(x, mu, sigma)
        else
            z = (x - mu) / sigma
            call standard_tails(z, df, far, near, log_far, log_near)
            if (z > 0.0_dp) then
                y = far
            else
                y = near
            end if
        end if
    end function t_sf

    pure elemental function t_logcdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: far
        real(dp) :: log_far
        real(dp) :: log_near
        real(dp) :: mu
        real(dp) :: near
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_logcdf(x, mu, sigma)
        else
            z = (x - mu) / sigma
            call standard_tails(z, df, far, near, log_far, log_near)
            if (z < 0.0_dp) then
                y = log_far
            else
                y = log_near
            end if
        end if
    end function t_logcdf

    pure elemental function t_logsf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: far
        real(dp) :: log_far
        real(dp) :: log_near
        real(dp) :: mu
        real(dp) :: near
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_logsf(x, mu, sigma)
        else
            z = (x - mu) / sigma
            call standard_tails(z, df, far, near, log_far, log_near)
            if (z > 0.0_dp) then
                y = log_far
            else
                y = log_near
            end if
        end if
    end function t_logsf

    pure elemental function t_ppf(p, df, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (.not. ieee_is_finite(df)) then
            y = normal_ppf(p, mu, sigma)
        else if (p <= 0.5_dp) then
            y = mu - sigma * standard_upper_quantile(p, df)
        else
            y = mu + sigma * standard_upper_quantile(1.0_dp - p, df)
        end if
    end function t_ppf

    pure elemental function t_isf(p, df, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else if (.not. ieee_is_finite(df)) then
            y = normal_isf(p, mu, sigma)
        else if (p <= 0.5_dp) then
            y = mu + sigma * standard_upper_quantile(p, df)
        else
            y = mu - sigma * standard_upper_quantile(1.0_dp - p, df)
        end if
    end function t_isf

    ! log density: -(df + 1)/2 log(1 + z**2/df) - log(df)/2 - log(B(df/2, 1/2)).
    pure elemental function standard_logpdf(z, df) result(y)
        real(dp), intent(in) :: z
        real(dp), intent(in) :: df
        real(dp) :: y

        real(dp) :: log_term
        real(dp) :: r

        if (ieee_is_nan(z)) then
            y = z
            return
        else if (.not. ieee_is_finite(z)) then
            y = negative_infinity(z)
            return
        end if

        r = abs(z) / sqrt(df)
        if (r <= 1.0_dp) then
            log_term = log1p_safe(r * r)
        else
            log_term = 2.0_dp * log(r) + log1p_safe(1.0_dp / r / r)
        end if
        y = -0.5_dp * (df + 1.0_dp) * log_term - 0.5_dp * log(df) - &
            log_beta(0.5_dp * df, 0.5_dp)
    end function standard_logpdf

    ! Tail probabilities for the standardized value z. far is the tail beyond
    ! |z| (the lower tail for z < 0, the upper tail for z > 0) and near is its
    ! complement, with their logarithms.
    pure elemental subroutine standard_tails(z, df, far, near, log_far, log_near)
        real(dp), intent(in) :: z
        real(dp), intent(in) :: df
        real(dp), intent(out) :: far
        real(dp), intent(out) :: near
        real(dp), intent(out) :: log_far
        real(dp), intent(out) :: log_near

        real(dp) :: a
        real(dp) :: log_constant
        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: p
        real(dp) :: q
        real(dp) :: r
        real(dp) :: r2
        real(dp) :: x
        real(dp) :: y

        if (ieee_is_nan(z)) then
            far = z
            near = z
            log_far = z
            log_near = z
            return
        else if (z == 0.0_dp) then
            far = 0.5_dp
            near = 0.5_dp
            log_far = -scifort_log_two
            log_near = -scifort_log_two
            return
        else if (.not. ieee_is_finite(z)) then
            far = 0.0_dp
            near = 1.0_dp
            log_far = negative_infinity(z)
            log_near = 0.0_dp
            return
        end if

        a = 0.5_dp * df
        if (abs(z) <= sqrt(df)) then
            r = abs(z) / sqrt(df)
            r2 = r * r
            x = 1.0_dp / (1.0_dp + r2)
            y = r2 / (1.0_dp + r2)
        else
            r = sqrt(df) / abs(z)
            r2 = r * r
            x = r2 / (1.0_dp + r2)
            y = 1.0_dp / (1.0_dp + r2)
        end if

        if (x < tiny_ratio) then
            ! (1/2) I_x(a, 1/2) = r**df / (2 a B(a, 1/2)) (1 + O(x)) with
            ! x = r**2. The power is formed directly when it is representable so
            ! that the result does not inherit the rounding of a large logarithm.
            log_constant = -scifort_log_two - log(a) - log_beta(a, 0.5_dp)
            log_far = df * log(r) + log_constant
            if (log_far > log(tiny(1.0_dp)) .and. df * log(r) > log(tiny(1.0_dp))) then
                far = r**df * exp(log_constant)
            else
                far = exp(log_far)
            end if
            near = 1.0_dp - far
            log_near = -far
        else if (y <= 0.0_dp) then
            ! |z| is so small relative to sqrt(df) that r**2 underflows.
            far = 0.5_dp
            near = 0.5_dp
            log_far = -scifort_log_two
            log_near = -scifort_log_two
        else
            call incomplete_beta_xy(a, 0.5_dp, x, y, p, q, logp, logq)
            far = 0.5_dp * p
            log_far = -scifort_log_two + logp
            near = 0.5_dp + 0.5_dp * q
            log_near = log1p_safe(-far)
        end if
    end subroutine standard_tails

    ! The value t >= 0 with upper tail probability p, 0 <= p <= 1/2, for the
    ! standard distribution: solve I_x(df/2, 1/2) = 2p for x and y = 1 - x,
    ! then t = sqrt(df y / x).
    pure elemental function standard_upper_quantile(p, df) result(t)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: df
        real(dp) :: t

        real(dp) :: a
        real(dp) :: log_x
        real(dp) :: x
        real(dp) :: y

        if (p <= 0.0_dp) then
            t = positive_infinity(p)
            return
        else if (p >= 0.5_dp) then
            t = 0.0_dp
            return
        end if

        a = 0.5_dp * df
        call beta_inverse_xy(a, 0.5_dp, 2.0_dp * p, 1.0_dp - 2.0_dp * p, x, y)
        if (x >= tiny_ratio) then
            t = sqrt(df) * sqrt(y) / sqrt(x)
        else
            ! Invert the leading term 2p = x**a / (a B(a, 1/2)):
            ! t = sqrt(df / x) = sqrt(df) (2p)**(-1/df) (a B(a, 1/2))**(-1/df).
            ! The power of p is formed directly; overflow to infinity is correct.
            log_x = (log(2.0_dp * p) + log(a) + log_beta(a, 0.5_dp)) / a
            if (-0.5_dp * log_x < log(huge(1.0_dp)) - 0.5_dp * log(df)) then
                t = sqrt(df) * (2.0_dp * p)**(-1.0_dp / df) * &
                    exp(-(log(a) + log_beta(a, 0.5_dp)) / df)
            else
                t = positive_infinity(p)
            end if
        end if
    end function standard_upper_quantile

    pure elemental logical function valid_df(df) result(valid)
        real(dp), intent(in) :: df

        valid = df > 0.0_dp
    end function valid_df

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp), intent(out) :: mu
        real(dp), intent(out) :: sigma

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_student_t
