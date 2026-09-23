! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Noncentral chi-square distribution matching scipy.stats.ncx2.
module scifort_ncx2
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_incomplete_gamma, only : gammainc, gammaincc
    use scifort_kinds, only : dp
    use scifort_log_gamma, only : log_gamma_one_plus
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    use scifort_digamma, only : digamma_positive
    implicit none
    private

    public :: ncx2_cdf, ncx2_isf, ncx2_logcdf, ncx2_logpdf
    public :: ncx2_logsf, ncx2_pdf, ncx2_ppf, ncx2_sf
    public :: ncx2_logpdf_derivatives

contains

    pure elemental function ncx2_logpdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: mu, sigma, z, ddf, dnc, dx, logf

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if
        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        z = (x - mu) / sigma
        if (z < 0.0_dp .or. .not. ieee_is_finite(z)) then
            if (z > 0.0_dp) then
                y = negative_infinity(z)
            else
                y = negative_infinity(z)
            end if
        else if (z == 0.0_dp) then
            if (df < 2.0_dp) then
                y = positive_infinity(z)
            else if (df == 2.0_dp) then
                y = -0.5_dp * nc - log(2.0_dp) - log(sigma)
            else
                y = negative_infinity(z)
            end if
        else
            call ncx2_logpdf_derivatives(z, df, nc, logf, ddf, dnc, dx)
            y = logf - log(sigma)
        end if
    end function ncx2_logpdf

    pure elemental function ncx2_pdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = ncx2_logpdf(x, df, nc, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function ncx2_pdf

    pure elemental function ncx2_cdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
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
                y = ncx2_tail(z, df, nc, .false.)
            end if
        end if
    end function ncx2_cdf

    pure elemental function ncx2_sf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
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
                y = ncx2_tail(z, df, nc, .true.)
            end if
        end if
    end function ncx2_sf

    pure elemental function ncx2_logcdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = ncx2_cdf(x, df, nc, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function ncx2_logcdf

    pure elemental function ncx2_logsf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, p
        p = ncx2_sf(x, df, nc, loc, scale)
        if (ieee_is_nan(p)) then
            y = p
        else if (p <= 0.0_dp) then
            y = negative_infinity(p)
        else
            y = log(p)
        end if
    end function ncx2_logsf

    pure elemental function ncx2_ppf(p, df, nc, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(df, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p <= 0.0_dp) then
            x = mu
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
        else if (p > 0.5_dp) then
            x = ncx2_isf(1.0_dp - p, df, nc, mu, sigma)
        else
            lo = 0.0_dp
            hi = max(1.0_dp, df + nc + 10.0_dp * sqrt(2.0_dp * (df + 2.0_dp * nc)) + 10.0_dp)
            do while (ncx2_tail(hi, df, nc, .false.) < p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 220
                mid = 0.5_dp * (lo + hi)
                if (ncx2_tail(mid, df, nc, .false.) < p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ncx2_ppf

    pure elemental function ncx2_isf(p, df, nc, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, sigma, lo, hi, mid
        integer :: iter
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(df, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p >= 1.0_dp) then
            x = mu
        else if (p <= 0.0_dp) then
            x = positive_infinity(p)
        else
            lo = 0.0_dp
            hi = max(1.0_dp, df + nc + 10.0_dp * sqrt(2.0_dp * (df + 2.0_dp * nc)) + 10.0_dp)
            do while (ncx2_tail(hi, df, nc, .true.) > p)
                hi = 2.0_dp * hi
            end do
            do iter = 1, 220
                mid = 0.5_dp * (lo + hi)
                if (ncx2_tail(mid, df, nc, .true.) > p) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = mu + sigma * 0.5_dp * (lo + hi)
        end if
    end function ncx2_isf

    pure subroutine ncx2_logpdf_derivatives(x, df, nc, logf, ddf, dnc, dx)
        real(dp), intent(in) :: x !! positive standardized observation
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: ddf !! derivative of log density with respect to df
        real(dp), intent(out) :: dnc !! derivative of log density with respect to nc
        real(dp), intent(out) :: dx !! derivative of log density with respect to x
        real(dp) :: m, logw, log_component, logtol
        integer :: k0, k, margin
        real(dp) :: maxlog, s0, sdf, snc, sx

        if (.not. valid_shapes(df, nc) .or. x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
            logf = quiet_nan(x)
            ddf = logf
            dnc = logf
            dx = logf
            return
        end if
        if (nc == 0.0_dp) then
            logf = central_logpdf(x, 0.5_dp * df)
            ddf = 0.5_dp * (log(0.5_dp * x) - digamma_positive(0.5_dp * df))
            dnc = 0.5_dp * (x / df - 1.0_dp)
            dx = (0.5_dp * df - 1.0_dp) / x - 0.5_dp
            return
        end if

        m = 0.5_dp * nc
        k0 = int(floor(m))
        logw = -m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp))
        maxlog = -huge(1.0_dp)
        s0 = 0.0_dp
        sdf = 0.0_dp
        snc = 0.0_dp
        sx = 0.0_dp
        call add_pdf_component(k0, logw, x, df, nc, maxlog, s0, sdf, snc, sx, log_component)

        margin = int(10.0_dp * sqrt(max(1.0_dp, m)) + 20.0_dp)
        logtol = log(4.0_dp * epsilon(1.0_dp))
        logw = -m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp))
        do k = k0 - 1, 0, -1
            logw = logw + log(real(k + 1, dp)) - log(m)
            call add_pdf_component(k, logw, x, df, nc, maxlog, s0, sdf, snc, sx, log_component)
            if (k < k0 - margin .and. log_component < maxlog + logtol) exit
        end do

        logw = -m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp))
        do k = k0 + 1, k0 + 100000
            logw = logw + log(m) - log(real(k, dp))
            call add_pdf_component(k, logw, x, df, nc, maxlog, s0, sdf, snc, sx, log_component)
            if (k > k0 + margin .and. log_component < maxlog + logtol) exit
        end do

        if (s0 <= 0.0_dp) then
            logf = negative_infinity(x)
            ddf = quiet_nan(x)
            dnc = ddf
            dx = ddf
        else
            logf = maxlog + log(s0)
            ddf = sdf / s0
            dnc = snc / s0
            dx = sx / s0
        end if
    end subroutine ncx2_logpdf_derivatives

    pure subroutine add_pdf_component(k, logw, x, df, nc, maxlog, s0, sdf, snc, sx, log_component)
        integer, intent(in) :: k !! Poisson mixture index
        real(dp), intent(in) :: logw !! logarithm of Poisson weight
        real(dp), intent(in) :: x !! positive standardized observation
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in) :: nc !! positive noncentrality
        real(dp), intent(inout) :: maxlog !! running logarithmic scale
        real(dp), intent(inout) :: s0 !! scaled density accumulator
        real(dp), intent(inout) :: sdf !! scaled df-derivative accumulator
        real(dp), intent(inout) :: snc !! scaled noncentrality-derivative accumulator
        real(dp), intent(inout) :: sx !! scaled x-derivative accumulator
        real(dp), intent(out) :: log_component !! logarithm of this mixture-density term
        real(dp) :: a, logt, fac, hdf, hnc, hx, rescale

        a = 0.5_dp * df + real(k, dp)
        logt = logw + central_logpdf(x, a)
        log_component = logt
        if (logt > maxlog) then
            if (maxlog > -huge(1.0_dp) / 2.0_dp) then
                rescale = exp(maxlog - logt)
                s0 = s0 * rescale
                sdf = sdf * rescale
                snc = snc * rescale
                sx = sx * rescale
            end if
            maxlog = logt
        end if
        fac = exp(logt - maxlog)
        hdf = 0.5_dp * (log(0.5_dp * x) - digamma_positive(a))
        hnc = real(k, dp) / nc - 0.5_dp
        hx = (a - 1.0_dp) / x - 0.5_dp
        s0 = s0 + fac
        sdf = sdf + fac * hdf
        snc = snc + fac * hnc
        sx = sx + fac * hx
    end subroutine add_pdf_component

    pure elemental function central_logpdf(x, a) result(y)
        real(dp), intent(in) :: x !! positive chi-square variate
        real(dp), intent(in) :: a !! half degrees of freedom, positive
        real(dp) :: y
        y = (a - 1.0_dp) * log(0.5_dp * x) - 0.5_dp * x - log_gamma(a) - log(2.0_dp)
    end function central_logpdf

    pure function ncx2_tail(x, df, nc, upper) result(prob)
        real(dp), intent(in) :: x !! positive standardized variate
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality
        logical, intent(in) :: upper !! true for SF, false for CDF
        real(dp) :: prob
        real(dp) :: m, t, weight, term, total, tol
        integer :: k0, k, margin

        t = 0.5_dp * x
        if (nc == 0.0_dp) then
            if (upper) then
                prob = gammaincc(0.5_dp * df, t)
            else
                prob = gammainc(0.5_dp * df, t)
            end if
            return
        end if
        m = 0.5_dp * nc
        k0 = int(floor(m))
        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        if (upper) then
            term = weight * gammaincc(0.5_dp * df + real(k0, dp), t)
        else
            term = weight * gammainc(0.5_dp * df + real(k0, dp), t)
        end if
        total = term
        tol = 8.0_dp * epsilon(1.0_dp)
        margin = int(10.0_dp * sqrt(max(1.0_dp, m)) + 20.0_dp)

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 - 1, 0, -1
            weight = weight * real(k + 1, dp) / m
            if (upper) then
                term = weight * gammaincc(0.5_dp * df + real(k, dp), t)
            else
                term = weight * gammainc(0.5_dp * df + real(k, dp), t)
            end if
            total = total + term
            if (k < k0 - margin .and. weight < tol * max(total, tiny(1.0_dp))) exit
        end do

        weight = exp(-m + real(k0, dp) * log(m) - log_gamma_one_plus(real(k0, dp)))
        do k = k0 + 1, k0 + 100000
            weight = weight * m / real(k, dp)
            if (upper) then
                term = weight * gammaincc(0.5_dp * df + real(k, dp), t)
            else
                term = weight * gammainc(0.5_dp * df + real(k, dp), t)
            end if
            total = total + term
            if (k > k0 + margin .and. weight < tol * max(total, tiny(1.0_dp))) exit
        end do
        prob = max(0.0_dp, min(1.0_dp, total))
    end function ncx2_tail

    pure elemental function valid_shapes(df, nc) result(ok)
        real(dp), intent(in) :: df !! candidate degrees of freedom
        real(dp), intent(in) :: nc !! candidate noncentrality
        logical :: ok
        ok = ieee_is_finite(df) .and. df > 0.0_dp .and. &
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

end module scifort_ncx2
