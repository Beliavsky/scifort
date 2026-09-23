! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Noncentral Student t distribution matching scipy.stats.nct.
!
! If U ~ Gamma(df/2, 1), then X = (Z + nc) / sqrt(2 U / df),
! with Z standard normal.  CDFs and densities are evaluated as expectations
! over U.  A double-exponential change of variable is applied to the gamma
! probability coordinate, so the quadrature is stable for df both below and
! above one and does not integrate the singular gamma density directly.
module scifort_nct
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_pi
    use scifort_digamma, only : digamma_positive
    use scifort_incomplete_gamma, only : gammainccinv, gammaincinv
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf
    use scifort_special_elementary, only : log_ndtr
    implicit none
    private

    public :: nct_cdf, nct_isf, nct_logcdf, nct_logpdf
    public :: nct_logsf, nct_pdf, nct_ppf, nct_sf
    public :: nct_logpdf_derivatives

    integer, parameter :: quadrature_max_half_steps = 256
    real(dp), parameter :: quadrature_h_coarse = 0.0625_dp
    real(dp), parameter :: quadrature_h_medium = 0.03125_dp
    real(dp), parameter :: quadrature_h_fine = 0.015625_dp

contains

    pure elemental function nct_logpdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_logpdf(x, mu + sigma * nc, sigma)
        else
            z = (x - mu) / sigma
            y = standard_logpdf(z, df, nc) - log(sigma)
        end if
    end function nct_logpdf

    pure elemental function nct_pdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly

        ly = nct_logpdf(x, df, nc, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function nct_pdf

    pure elemental function nct_logcdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = negative_infinity(x)
            else
                y = 0.0_dp
            end if
        else if (.not. ieee_is_finite(df)) then
            y = normal_logcdf(x, mu + sigma * nc, sigma)
        else
            z = (x - mu) / sigma
            y = standard_logtail(z, df, nc, .false.)
        end if
    end function nct_logcdf

    pure elemental function nct_logsf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(x)
            end if
        else if (.not. ieee_is_finite(df)) then
            y = normal_logsf(x, mu + sigma * nc, sigma)
        else
            z = (x - mu) / sigma
            y = standard_logtail(z, df, nc, .true.)
        end if
    end function nct_logsf

    pure elemental function nct_cdf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of lower-tail probability P(X <= x)
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf, mu, sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_cdf(x, mu + sigma * nc, sigma)
        else
            lcdf = nct_logcdf(x, df, nc, mu, sigma)
            if (ieee_is_nan(lcdf)) then
                y = lcdf
            else if (lcdf < -log(2.0_dp)) then
                y = exp(lcdf)
            else
                lsf = nct_logsf(x, df, nc, mu, sigma)
                y = -expm1_safe(lsf)
            end if
        end if
    end function nct_cdf

    pure elemental function nct_sf(x, df, nc, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of upper-tail probability P(X > x)
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lcdf, lsf, mu, sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(df)) then
            y = normal_sf(x, mu + sigma * nc, sigma)
        else
            lsf = nct_logsf(x, df, nc, mu, sigma)
            if (ieee_is_nan(lsf)) then
                y = lsf
            else if (lsf < -log(2.0_dp)) then
                y = exp(lsf)
            else
                lcdf = nct_logcdf(x, df, nc, mu, sigma)
                y = -expm1_safe(lcdf)
            end if
        end if
    end function nct_sf

    pure elemental function nct_ppf(probability, df, nc, loc, scale) result(x)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability <= 0.0_dp) then
            x = negative_infinity(probability)
        else if (probability >= 1.0_dp) then
            x = positive_infinity(probability)
        else if (.not. ieee_is_finite(df)) then
            x = normal_ppf(probability, mu + sigma * nc, sigma)
        else
            z = standard_quantile(probability, df, nc, .false.)
            x = mu + sigma * z
        end if
    end function nct_ppf

    pure elemental function nct_isf(probability, df, nc, loc, scale) result(x)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: df !! positive degrees of freedom; +infinity gives a shifted normal
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_parameters(df, nc) .or. .not. valid_loc_scale(mu, sigma) .or. &
            .not. valid_probability(probability)) then
            x = quiet_nan(probability)
        else if (probability >= 1.0_dp) then
            x = negative_infinity(probability)
        else if (probability <= 0.0_dp) then
            x = positive_infinity(probability)
        else if (.not. ieee_is_finite(df)) then
            x = normal_isf(probability, mu + sigma * nc, sigma)
        else
            z = standard_quantile(probability, df, nc, .true.)
            x = mu + sigma * z
        end if
    end function nct_isf

    pure subroutine nct_logpdf_derivatives(z, df, nc, logf, ddf, dnc, dz)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(out) :: logf !! standardized log density
        real(dp), intent(out) :: ddf !! derivative with respect to df
        real(dp), intent(out) :: dnc !! derivative with respect to nc
        real(dp), intent(out) :: dz !! derivative with respect to standardized observation
        real(dp) :: a, arg, factor_df, li, logq, logqc, logjac, maxlog
        real(dp) :: s, sum_df, sum_nc, sum_w, sum_z, t, u, v
        real(dp) :: terms(-quadrature_max_half_steps:quadrature_max_half_steps)
        real(dp) :: fdf(-quadrature_max_half_steps:quadrature_max_half_steps)
        real(dp) :: fnc(-quadrature_max_half_steps:quadrature_max_half_steps)
        real(dp) :: fz(-quadrature_max_half_steps:quadrature_max_half_steps)
        real(dp) :: h
        integer :: j, half_steps

        if (.not. valid_parameters(df, nc) .or. .not. ieee_is_finite(z)) then
            logf = quiet_nan(z); ddf = logf; dnc = logf; dz = logf
            return
        end if
        if (.not. ieee_is_finite(df)) then
            logf = normal_logpdf(z, nc, 1.0_dp)
            ddf = 0.0_dp; dnc = z - nc; dz = nc - z
            return
        end if

        a = 0.5_dp * df
        call quadrature_grid(df, h, half_steps)
        maxlog = negative_infinity(z)
        do j = -half_steps, half_steps
            t = h * real(j, dp)
            call gamma_probability_node(t, a, u, logq, logqc, logjac)
            if (u <= 0.0_dp .or. .not. ieee_is_finite(u)) then
                terms(j) = negative_infinity(z)
                fdf(j) = 0.0_dp; fnc(j) = 0.0_dp; fz(j) = 0.0_dp
            else
                s = sqrt(2.0_dp * u / df)
                arg = z * s - nc
                li = logjac + log(s) - 0.5_dp * arg * arg - scifort_log_sqrt_two_pi
                if (abs(j) == half_steps) li = li - log(2.0_dp)
                terms(j) = li
                factor_df = 0.5_dp * (log(u) - digamma_positive(a) - 0.5_dp / a + &
                    arg * z * s / (2.0_dp * a))
                fdf(j) = factor_df
                fnc(j) = arg
                fz(j) = -arg * s
                maxlog = max(maxlog, li)
            end if
        end do
        if (.not. ieee_is_finite(maxlog)) then
            logf = negative_infinity(z); ddf = 0.0_dp; dnc = 0.0_dp; dz = 0.0_dp
            return
        end if
        sum_w = 0.0_dp; sum_df = 0.0_dp; sum_nc = 0.0_dp; sum_z = 0.0_dp
        do j = -half_steps, half_steps
            if (ieee_is_finite(terms(j))) then
                v = exp(terms(j) - maxlog)
                sum_w = sum_w + v
                sum_df = sum_df + v * fdf(j)
                sum_nc = sum_nc + v * fnc(j)
                sum_z = sum_z + v * fz(j)
            end if
        end do
        logf = maxlog + log(sum_w) + log(h)
        ddf = sum_df / sum_w
        dnc = sum_nc / sum_w
        dz = sum_z / sum_w
    end subroutine nct_logpdf_derivatives

    pure function standard_logpdf(z, df, nc) result(logf)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp) :: logf, ddf, dnc, dz
        call nct_logpdf_derivatives(z, df, nc, logf, ddf, dnc, dz)
    end function standard_logpdf

    pure function standard_logtail(z, df, nc, upper) result(logp)
        real(dp), intent(in) :: z !! standardized observation
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        logical, intent(in) :: upper !! true for upper tail and false for lower tail
        real(dp) :: logp, a, arg, h, li, logjac, logq, logqc, t, u
        integer :: j, half_steps

        a = 0.5_dp * df
        call quadrature_grid(df, h, half_steps)
        logp = negative_infinity(z)
        do j = -half_steps, half_steps
            t = h * real(j, dp)
            call gamma_probability_node(t, a, u, logq, logqc, logjac)
            if (u <= 0.0_dp) then
                arg = -nc
            else if (.not. ieee_is_finite(u)) then
                if (z > 0.0_dp) then
                    arg = positive_infinity(z)
                else if (z < 0.0_dp) then
                    arg = negative_infinity(z)
                else
                    arg = -nc
                end if
            else
                arg = z * sqrt(2.0_dp * u / df) - nc
            end if
            if (upper) arg = -arg
            li = logjac + log_ndtr(arg)
            if (abs(j) == half_steps) li = li - log(2.0_dp)
            logp = logaddexp_pair(logp, li)
        end do
        logp = min(0.0_dp, logp + log(h))
    end function standard_logtail

    pure function standard_quantile(probability, df, nc, upper) result(z)
        real(dp), intent(in) :: probability !! requested tail probability in (0,1)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        logical, intent(in) :: upper !! true when probability is an upper tail
        real(dp) :: z, lo, hi, mid, logtarget, value
        integer :: iter

        logtarget = log(probability)
        lo = min(-1.0_dp, nc - 1.0_dp)
        hi = max(1.0_dp, nc + 1.0_dp)
        if (.not. upper) then
            do iter = 1, 80
                if (standard_logtail(lo, df, nc, .false.) <= logtarget) exit
                lo = 2.0_dp * lo - 1.0_dp
            end do
            do iter = 1, 80
                if (standard_logtail(hi, df, nc, .false.) >= logtarget) exit
                hi = 2.0_dp * hi + 1.0_dp
            end do
            do iter = 1, 70
                mid = 0.5_dp * (lo + hi)
                value = standard_logtail(mid, df, nc, .false.)
                if (value < logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        else
            do iter = 1, 80
                if (standard_logtail(lo, df, nc, .true.) >= logtarget) exit
                lo = 2.0_dp * lo - 1.0_dp
            end do
            do iter = 1, 80
                if (standard_logtail(hi, df, nc, .true.) <= logtarget) exit
                hi = 2.0_dp * hi + 1.0_dp
            end do
            do iter = 1, 70
                mid = 0.5_dp * (lo + hi)
                value = standard_logtail(mid, df, nc, .true.)
                if (value > logtarget) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
        end if
        z = 0.5_dp * (lo + hi)
    end function standard_quantile

    pure subroutine quadrature_grid(df, h, half_steps)
        real(dp), intent(in) :: df !! positive finite degrees of freedom
        real(dp), intent(out) :: h !! quadrature spacing on the double-exponential coordinate
        integer, intent(out) :: half_steps !! number of positive and negative quadrature steps

        if (df < 0.3_dp) then
            h = quadrature_h_fine
            half_steps = 256
        else if (df < 1.0_dp) then
            h = quadrature_h_medium
            half_steps = 128
        else
            h = quadrature_h_coarse
            half_steps = 64
        end if
    end subroutine quadrature_grid

    pure subroutine gamma_probability_node(t, a, u, logq, logqc, logjac)
        real(dp), intent(in) :: t !! double-exponential integration coordinate
        real(dp), intent(in) :: a !! gamma shape df/2
        real(dp), intent(out) :: u !! gamma quantile at the transformed probability
        real(dp), intent(out) :: logq !! logarithm of lower gamma probability coordinate
        real(dp), intent(out) :: logqc !! logarithm of upper gamma probability coordinate
        real(dp), intent(out) :: logjac !! logarithm of dq/dt
        real(dp) :: e, ld, q, qc, v

        v = scifort_pi * sinh(t)
        if (v >= 0.0_dp) then
            e = exp(-v)
            ld = log1p_safe(e)
            logq = -ld
            logqc = -v - ld
            qc = exp(logqc)
            u = gammainccinv(a, qc)
        else
            e = exp(v)
            ld = log1p_safe(e)
            logq = v - ld
            logqc = -ld
            q = exp(logq)
            u = gammaincinv(a, q)
        end if
        logjac = log(scifort_pi) + log(cosh(t)) + logq + logqc
    end subroutine gamma_probability_node

    pure elemental function logaddexp_pair(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, m
        if (.not. ieee_is_finite(a)) then
            y = b
        else if (.not. ieee_is_finite(b)) then
            y = a
        else
            m = max(a, b)
            y = m + log1p_safe(exp(-abs(a - b)))
        end if
    end function logaddexp_pair

    pure elemental logical function valid_parameters(df, nc)
        real(dp), intent(in) :: df !! degrees of freedom to validate
        real(dp), intent(in) :: nc !! noncentrality to validate
        valid_parameters = df > 0.0_dp .and. .not. ieee_is_nan(df) .and. ieee_is_finite(nc)
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability to validate
        valid_probability = ieee_is_finite(p) .and. p >= 0.0_dp .and. p <= 1.0_dp
    end function valid_probability

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! optional location parameter
        real(dp), intent(in), optional :: scale !! optional positive scale parameter
        real(dp), intent(out) :: mu !! resolved location parameter
        real(dp), intent(out) :: sigma !! resolved positive scale parameter
        mu = 0.0_dp; sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_nct
