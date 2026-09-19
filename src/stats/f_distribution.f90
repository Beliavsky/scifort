! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! F distribution with dfn, dfd > 0 degrees of freedom, location loc, and
! scale > 0. Matches scipy.stats.f. With a = dfn/2, b = dfd/2, and
! standardized z >= 0, the lower tail is I_x(a, b) with
! x = dfn z / (dfn z + dfd) and 1 - x = dfd / (dfn z + dfd); both are formed
! directly so that neither tail loses accuracy.

module scifort_f_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_beta, only : beta_inverse_xy, incomplete_beta_xy, &
        log_beta_kernel
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta, stirling_threshold
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: f_cdf
    public :: f_isf
    public :: f_logcdf
    public :: f_logpdf
    public :: f_logsf
    public :: f_pdf
    public :: f_ppf
    public :: f_sf

    ! Below this value a coordinate of the incomplete beta function would
    ! lose precision to underflow; the leading term is used instead.
    real(dp), parameter :: tiny_coordinate = 1.0e-300_dp

contains

    pure elemental function f_pdf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > huge(z)) then
            y = 0.0_dp
        else if (z == 0.0_dp) then
            if (dfn < 2.0_dp) then
                y = positive_infinity(z)
            else if (dfn == 2.0_dp) then
                y = 1.0_dp / sigma
            else
                y = 0.0_dp
            end if
        else
            y = exp(standard_logpdf(z, dfn, dfd)) / sigma
        end if
    end function f_pdf

    pure elemental function f_logpdf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (ieee_is_nan(z)) then
            y = z
        else if (z < 0.0_dp .or. z > huge(z)) then
            y = negative_infinity(z)
        else if (z == 0.0_dp) then
            if (dfn < 2.0_dp) then
                y = positive_infinity(z)
            else if (dfn == 2.0_dp) then
                y = -log(sigma)
            else
                y = negative_infinity(z)
            end if
        else
            y = standard_logpdf(z, dfn, dfd) - log(sigma)
        end if
    end function f_logpdf

    pure elemental function f_cdf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: mu
        real(dp) :: q
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
        else
            call standard_tails((x - mu) / sigma, dfn, dfd, y, q, logp, logq)
        end if
    end function f_cdf

    pure elemental function f_sf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: logq
        real(dp) :: mu
        real(dp) :: p
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
        else
            call standard_tails((x - mu) / sigma, dfn, dfd, p, y, logp, logq)
        end if
    end function f_sf

    pure elemental function f_logcdf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logq
        real(dp) :: mu
        real(dp) :: p
        real(dp) :: q
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
        else
            call standard_tails((x - mu) / sigma, dfn, dfd, p, q, y, logq)
        end if
    end function f_logcdf

    pure elemental function f_logsf(x, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: logp
        real(dp) :: mu
        real(dp) :: p
        real(dp) :: q
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(x)
        else
            call standard_tails((x - mu) / sigma, dfn, dfd, p, q, logp, y)
        end if
    end function f_logsf

    pure elemental function f_ppf(p, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else
            y = mu + sigma * standard_quantile(p, 1.0_dp - p, dfn, dfd)
        end if
    end function f_ppf

    pure elemental function f_isf(p, dfn, dfd, loc, scale) result(y)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, finite and > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(dfn) .and. valid_df(dfd))) then
            y = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            y = quiet_nan(p)
        else
            y = mu + sigma * standard_quantile(1.0_dp - p, p, dfn, dfd)
        end if
    end function f_isf

    ! log density for finite z > 0. For moderate degrees of freedom:
    ! a log(w) - (a + b) log(1 + w) - log(z) - log(B(a, b)), w = dfn z / dfd.
    ! When both a and b are large the stable beta kernel is used:
    ! log(x**a (1 - x)**b / B(a, b)) - log(z).
    pure elemental function standard_logpdf(z, dfn, dfd) result(y)
        real(dp), intent(in) :: z !! standardized value, finite and > 0
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp) :: y

        real(dp) :: a
        real(dp) :: b
        real(dp) :: log_term
        real(dp) :: ratio
        real(dp) :: w
        real(dp) :: x
        real(dp) :: x_complement

        a = 0.5_dp * dfn
        b = 0.5_dp * dfd
        if (min(a, b) >= stirling_threshold) then
            call beta_coordinates(z, dfn, dfd, x, x_complement)
            y = log_beta_kernel(a, b, x, x_complement) - log(z)
        else
            ! With w = dfn z / dfd the density is
            ! w**a (1 + w)**(-a - b) / (z B(a, b)). For w > 1 the terms
            ! a log(w) - (a + b) log(1 + w) are rewritten as
            ! -a log1p(1/w) - b log(1 + w) so that no large terms cancel.
            ratio = dfn / dfd
            w = ratio * z
            if (w <= 1.0_dp) then
                y = a * (log(ratio) + log(z)) - (a + b) * log1p_safe(w)
            else if (ieee_is_finite(w)) then
                log_term = log(w) + log1p_safe(1.0_dp / w)
                y = -a * log1p_safe(1.0_dp / w) - b * log_term
            else
                y = -b * (log(ratio) + log(z))
            end if
            y = y - log(z) - log_beta(a, b)
        end if
    end function standard_logpdf

    ! x = dfn z / (dfn z + dfd) and its complement dfd / (dfn z + dfd),
    ! each formed without cancellation or overflow.
    pure elemental subroutine beta_coordinates(z, dfn, dfd, x, x_complement)
        real(dp), intent(in) :: z !! standardized value, > 0
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp), intent(out) :: x !! dfn z / (dfn z + dfd)
        real(dp), intent(out) :: x_complement !! dfd / (dfn z + dfd)

        real(dp) :: v
        real(dp) :: w

        if (z <= dfd / dfn) then
            w = (dfn / dfd) * z
            x = w / (1.0_dp + w)
            x_complement = 1.0_dp / (1.0_dp + w)
        else
            v = (dfd / dfn) / z
            x = 1.0_dp / (1.0_dp + v)
            x_complement = v / (1.0_dp + v)
        end if
    end subroutine beta_coordinates

    ! Lower and upper tails and their logarithms at standardized z.
    pure elemental subroutine standard_tails(z, dfn, dfd, p, q, logp, logq)
        real(dp), intent(in) :: z !! standardized value
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp), intent(out) :: p !! lower-tail probability
        real(dp), intent(out) :: q !! upper-tail probability
        real(dp), intent(out) :: logp !! log(p)
        real(dp), intent(out) :: logq !! log(q)

        real(dp) :: a
        real(dp) :: b
        real(dp) :: x
        real(dp) :: x_complement

        if (ieee_is_nan(z)) then
            p = z
            q = z
            logp = z
            logq = z
            return
        else if (z <= 0.0_dp) then
            p = 0.0_dp
            q = 1.0_dp
            logp = negative_infinity(z)
            logq = 0.0_dp
            return
        else if (.not. ieee_is_finite(z)) then
            p = 1.0_dp
            q = 0.0_dp
            logp = 0.0_dp
            logq = negative_infinity(z)
            return
        end if

        a = 0.5_dp * dfn
        b = 0.5_dp * dfd
        call beta_coordinates(z, dfn, dfd, x, x_complement)
        if (x < tiny_coordinate) then
            ! I_x(a, b) = x**a / (a B(a, b)) (1 + O(x)), with x = dfn z / dfd.
            call leading_term(a, b, log(dfn / dfd) + log(z), logp, p)
            q = 1.0_dp - p
            logq = -p
        else if (x_complement < tiny_coordinate) then
            call leading_term(b, a, log(dfd / dfn) - log(z), logq, q)
            p = 1.0_dp - q
            logp = -q
        else
            call incomplete_beta_xy(a, b, x, x_complement, p, q, logp, logq)
        end if
    end subroutine standard_tails

    ! The leading term c**s / (s B(s, t)) of I_c(s, t) for tiny c, given
    ! log(c), and its logarithm. The power is formed directly when it is
    ! representable so that the value does not inherit the rounding of a
    ! large logarithm.
    pure elemental subroutine leading_term(s, t, log_c, log_value, value)
        real(dp), intent(in) :: s !! shape attached to the tiny coordinate c
        real(dp), intent(in) :: t !! other shape
        real(dp), intent(in) :: log_c !! log(c)
        real(dp), intent(out) :: log_value !! log(value)
        real(dp), intent(out) :: value !! leading term c**s / (s B(s, t))

        real(dp) :: log_constant

        log_constant = -log(s) - log_beta(s, t)
        log_value = s * log_c + log_constant
        if (s * log_c > log(tiny(1.0_dp)) .and. log_value > log(tiny(1.0_dp))) then
            value = exp(log_c)**s * exp(log_constant)
        else
            value = exp(log_value)
        end if
    end subroutine leading_term

    ! Standardized quantile for lower tail p and upper tail q = 1 - p:
    ! z = (dfd / dfn) x / (1 - x) with I_x(dfn/2, dfd/2) = p.
    pure elemental function standard_quantile(p, q, dfn, dfd) result(z)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: q !! upper-tail probability, 1 - p
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp) :: z

        real(dp) :: x
        real(dp) :: x_complement

        if (p <= 0.0_dp) then
            z = 0.0_dp
        else if (q <= 0.0_dp) then
            z = positive_infinity(p)
        else
            call beta_inverse_xy(0.5_dp * dfn, 0.5_dp * dfd, p, q, x, x_complement)
            if (x_complement <= 0.0_dp) then
                z = positive_infinity(p)
            else
                z = (dfd / dfn) * (x / x_complement)
            end if
        end if
    end function standard_quantile

    pure elemental logical function valid_df(df) result(valid)
        real(dp), intent(in) :: df !! degrees of freedom to check

        valid = ieee_is_finite(df) .and. df > 0.0_dp
    end function valid_df

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

end module scifort_f_distribution
