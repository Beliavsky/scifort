! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Noncentral F distribution matching scipy.stats.ncf.
module scifort_ncf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_digamma, only : digamma_positive
    use scifort_incomplete_beta, only : betainc
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_beta, log_gamma_one_plus
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: ncf_cdf, ncf_isf, ncf_logcdf, ncf_logpdf
    public :: ncf_logsf, ncf_pdf, ncf_ppf, ncf_sf
    public :: ncf_logpdf_derivatives

contains

    pure elemental function ncf_logpdf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, logf, dd1, dd2, dnc, dx

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(dfn, dfd, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z < 0.0_dp .or. (.not. ieee_is_finite(z) .and. z > 0.0_dp)) then
            y = negative_infinity(z)
        else if (z == 0.0_dp) then
            if (dfn < 2.0_dp) then
                y = positive_infinity(z)
            else if (dfn == 2.0_dp) then
                y = -0.5_dp * nc - log(sigma)
            else
                y = negative_infinity(z)
            end if
        else
            call ncf_logpdf_derivatives(z, dfn, dfd, nc, logf, dd1, dd2, dnc, dx)
            y = logf - log(sigma)
        end if
    end function ncf_logpdf

    pure elemental function ncf_pdf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = ncf_logpdf(x, dfn, dfd, nc, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function ncf_pdf

    pure elemental function ncf_cdf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(dfn, dfd, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 0.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 1.0_dp
            else
                y = ncf_tail(z, dfn, dfd, nc, .false.)
            end if
        end if
    end function ncf_cdf

    pure elemental function ncf_sf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(dfn, dfd, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 0.0_dp) then
                y = 1.0_dp
            else if (.not. ieee_is_finite(z)) then
                y = 0.0_dp
            else
                y = ncf_tail(z, dfn, dfd, nc, .true.)
            end if
        end if
    end function ncf_sf

    pure elemental function ncf_logcdf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = ncf_cdf(x, dfn, dfd, nc, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function ncf_logcdf

    pure elemental function ncf_logsf(x, dfn, dfd, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = ncf_sf(x, dfn, dfd, nc, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function ncf_logsf

    pure elemental function ncf_ppf(p, dfn, dfd, nc, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(dfn, dfd, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = mu
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else if (p > 0.5_dp) then
            x = ncf_isf(1.0_dp - p, dfn, dfd, nc, mu, sigma)
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do while (ncf_tail(hi, dfn, dfd, nc, .false.) < p)
                hi = 2.0_dp * hi
                if (.not. ieee_is_finite(hi)) exit
            end do
            do iter = 1, 180
                mid = 0.5_dp * (lo + hi)
                if (ncf_tail(mid, dfn, dfd, nc, .false.) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ncf_ppf

    pure elemental function ncf_isf(p, dfn, dfd, nc, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(dfn, dfd, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p >= 1.0_dp) then
            x = mu
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = 1.0_dp
            do while (ncf_tail(hi, dfn, dfd, nc, .true.) > p)
                hi = 2.0_dp * hi
                if (.not. ieee_is_finite(hi)) exit
            end do
            do iter = 1, 180
                mid = 0.5_dp * (lo + hi)
                if (ncf_tail(mid, dfn, dfd, nc, .true.) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ncf_isf

    pure subroutine ncf_logpdf_derivatives(x, dfn, dfd, nc, logf, ddfn, ddfd, dnc, dx)
        real(dp), intent(in) :: x !! positive standardized observation
        real(dp), intent(in) :: dfn !! numerator degrees of freedom
        real(dp), intent(in) :: dfd !! denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: ddfn !! derivative with respect to dfn
        real(dp), intent(out) :: ddfd !! derivative with respect to dfd
        real(dp), intent(out) :: dnc !! derivative with respect to nc
        real(dp), intent(out) :: dx !! derivative with respect to x
        real(dp) :: m, logw, log_component, logtol, maxlog, s0, sd1, sd2, snc, sx
        integer :: k0, k, margin

        if (.not. valid_shapes(dfn, dfd, nc) .or. x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
            logf = quiet_nan(x)
            ddfn = logf
            ddfd = logf
            dnc = logf
            dx = logf
            return
        end if
        if (nc == 0.0_dp) then
            call component_values(0, 0.0_dp, x, dfn, dfd, logf, ddfn, ddfd, dx)
            dnc = 0.5_dp * (((dfn + dfd) / dfn) * f_coordinate(x, dfn, dfd) - 1.0_dp)
            return
        end if

        m = 0.5_dp * nc
        k0 = int(floor(m))
        logw = -m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp))
        maxlog = -huge(1.0_dp)
        s0 = 0.0_dp
        sd1 = 0.0_dp
        sd2 = 0.0_dp
        snc = 0.0_dp
        sx = 0.0_dp
        call add_component(k0, logw, x, dfn, dfd, nc, maxlog, s0, sd1, sd2, snc, sx, log_component)
        margin = int(10.0_dp * sqrt(max(1.0_dp, m)) + 20.0_dp)
        logtol = log(4.0_dp * epsilon(1.0_dp))

        do k = k0 - 1, 0, -1
            logw = logw + log(real(k + 1, dp)) - log(m)
            call add_component(k, logw, x, dfn, dfd, nc, maxlog, s0, sd1, sd2, snc, sx, log_component)
            if (k < k0 - margin .and. log_component < maxlog + logtol) exit
        end do

        logw = -m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp))
        do k = k0 + 1, k0 + 100000
            logw = logw + log(m) - log(real(k, dp))
            call add_component(k, logw, x, dfn, dfd, nc, maxlog, s0, sd1, sd2, snc, sx, log_component)
            if (k > k0 + margin .and. log_component < maxlog + logtol) exit
        end do

        if (s0 <= 0.0_dp) then
            logf = negative_infinity(x)
            ddfn = quiet_nan(x)
            ddfd = ddfn
            dnc = ddfn
            dx = ddfn
        else
            logf = maxlog + log(s0)
            ddfn = sd1 / s0
            ddfd = sd2 / s0
            dnc = snc / s0
            dx = sx / s0
        end if
    end subroutine ncf_logpdf_derivatives

    pure subroutine add_component(k, logw, x, dfn, dfd, nc, maxlog, s0, sd1, sd2, snc, sx, log_component)
        integer, intent(in) :: k !! Poisson mixture index
        real(dp), intent(in) :: logw !! log Poisson weight
        real(dp), intent(in) :: x !! positive standardized observation
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! positive noncentrality parameter
        real(dp), intent(inout) :: maxlog !! running logarithmic scale
        real(dp), intent(inout) :: s0 !! scaled density accumulator
        real(dp), intent(inout) :: sd1 !! scaled numerator-df derivative accumulator
        real(dp), intent(inout) :: sd2 !! scaled denominator-df derivative accumulator
        real(dp), intent(inout) :: snc !! scaled noncentrality-derivative accumulator
        real(dp), intent(inout) :: sx !! scaled x-derivative accumulator
        real(dp), intent(out) :: log_component !! logarithm of this mixture-density term
        real(dp) :: logt, hd1, hd2, hx, fac, rescale

        call component_values(k, logw, x, dfn, dfd, logt, hd1, hd2, hx)
        log_component = logt
        if (logt > maxlog) then
            if (maxlog > -huge(1.0_dp) / 2.0_dp) then
                rescale = exp(maxlog - logt)
                s0 = s0 * rescale
                sd1 = sd1 * rescale
                sd2 = sd2 * rescale
                snc = snc * rescale
                sx = sx * rescale
            end if
            maxlog = logt
        end if
        fac = exp(logt - maxlog)
        s0 = s0 + fac
        sd1 = sd1 + fac * hd1
        sd2 = sd2 + fac * hd2
        snc = snc + fac * (real(k, dp) / nc - 0.5_dp)
        sx = sx + fac * hx
    end subroutine add_component

    pure subroutine component_values(k, logw, x, dfn, dfd, logt, hd1, hd2, hx)
        integer, intent(in) :: k !! mixture index
        real(dp), intent(in) :: logw !! logarithmic mixture weight
        real(dp), intent(in) :: x !! positive standardized observation
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(out) :: logt !! logarithm of weighted density component
        real(dp), intent(out) :: hd1 !! component log-derivative with respect to dfn
        real(dp), intent(out) :: hd2 !! component log-derivative with respect to dfd
        real(dp), intent(out) :: hx !! component log-derivative with respect to x
        real(dp) :: a, b, y, omy, d
        real(dp) :: dlogy_d1, dlogomy_d1, dlogj_d1
        real(dp) :: dlogy_d2, dlogomy_d2, dlogj_d2

        a = 0.5_dp * dfn + real(k, dp)
        b = 0.5_dp * dfd
        d = dfd + dfn * x
        y = dfn * x / d
        omy = dfd / d
        logt = logw - log_beta(a, b) + (a - 1.0_dp) * log(y) + &
            (b - 1.0_dp) * log(omy) + log(dfn) + log(dfd) - 2.0_dp * log(d)

        dlogy_d1 = omy / dfn
        dlogomy_d1 = -y / dfn
        dlogj_d1 = (1.0_dp - 2.0_dp * y) / dfn
        hd1 = 0.5_dp * (-digamma_positive(a) + digamma_positive(a + b) + log(y)) + &
            (a - 1.0_dp) * dlogy_d1 + (b - 1.0_dp) * dlogomy_d1 + dlogj_d1

        dlogy_d2 = -omy / dfd
        dlogomy_d2 = y / dfd
        dlogj_d2 = (2.0_dp * y - 1.0_dp) / dfd
        hd2 = 0.5_dp * (-digamma_positive(b) + digamma_positive(a + b) + log(omy)) + &
            (a - 1.0_dp) * dlogy_d2 + (b - 1.0_dp) * dlogomy_d2 + dlogj_d2

        hx = ((a - 1.0_dp) * omy - (b - 1.0_dp) * y - 2.0_dp * y) / x
    end subroutine component_values

    pure elemental function f_coordinate(x, dfn, dfd) result(y)
        real(dp), intent(in) :: x !! positive standardized F variate
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp) :: y
        y = dfn * x / (dfd + dfn * x)
    end function f_coordinate

    pure function ncf_tail(x, dfn, dfd, nc, upper) result(prob)
        real(dp), intent(in) :: x !! positive standardized F variate
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality
        logical, intent(in) :: upper !! true for SF, false for CDF
        real(dp) :: prob
        real(dp) :: y, omy, m, weight, term, total, tol
        integer :: k0, k, margin

        y = f_coordinate(x, dfn, dfd)
        omy = dfd / (dfd + dfn * x)
        if (nc == 0.0_dp) then
            if (upper) then
                prob = betainc(0.5_dp * dfd, 0.5_dp * dfn, omy)
            else
                prob = betainc(0.5_dp * dfn, 0.5_dp * dfd, y)
            end if
            return
        end if
        m = 0.5_dp * nc
        k0 = int(floor(m))
        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        if (upper) then
            term = weight * betainc(0.5_dp * dfd, 0.5_dp * dfn + real(k0, dp), omy)
        else
            term = weight * betainc(0.5_dp * dfn + real(k0, dp), 0.5_dp * dfd, y)
        end if
        total = term
        tol = 8.0_dp * epsilon(1.0_dp)
        margin = int(10.0_dp * sqrt(max(1.0_dp, m)) + 20.0_dp)

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 - 1, 0, -1
            weight = weight * real(k + 1, dp) / m
            if (upper) then
                term = weight * betainc(0.5_dp * dfd, 0.5_dp * dfn + real(k, dp), omy)
            else
                term = weight * betainc(0.5_dp * dfn + real(k, dp), 0.5_dp * dfd, y)
            end if
            total = total + term
            if (k < k0 - margin .and. weight < tol * max(total, tiny(1.0_dp))) exit
        end do

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 + 1, k0 + 100000
            weight = weight * m / real(k, dp)
            if (upper) then
                term = weight * betainc(0.5_dp * dfd, 0.5_dp * dfn + real(k, dp), omy)
            else
                term = weight * betainc(0.5_dp * dfn + real(k, dp), 0.5_dp * dfd, y)
            end if
            total = total + term
            if (k > k0 + margin .and. weight < tol * max(total, tiny(1.0_dp))) exit
        end do
        prob = max(0.0_dp, min(1.0_dp, total))
    end function ncf_tail

    pure elemental function valid_shapes(dfn, dfd, nc) result(ok)
        real(dp), intent(in) :: dfn !! candidate numerator degrees of freedom
        real(dp), intent(in) :: dfd !! candidate denominator degrees of freedom
        real(dp), intent(in) :: nc !! candidate noncentrality parameter
        logical :: ok
        ok = ieee_is_finite(dfn) .and. dfn > 0.0_dp .and. &
            ieee_is_finite(dfd) .and. dfd > 0.0_dp .and. &
            ieee_is_finite(nc) .and. nc >= 0.0_dp
    end function valid_shapes

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location, absent for zero
        real(dp), intent(in), optional :: scale !! requested scale, absent for one
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_ncf
