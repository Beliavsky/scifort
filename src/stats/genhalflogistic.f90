! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized half-logistic distribution matching scipy.stats.genhalflogistic.
! Standard support is 0 <= z <= 1/c with c > 0.
module scifort_genhalflogistic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private
    public :: genhalflogistic_pdf, genhalflogistic_logpdf, genhalflogistic_cdf
    public :: genhalflogistic_sf, genhalflogistic_logcdf, genhalflogistic_logsf
    public :: genhalflogistic_ppf, genhalflogistic_isf
contains
    pure elemental function genhalflogistic_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lt, lu, mu, sigma, t, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu .or. .not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z > 1.0_dp / c) then
                y = negative_infinity(x)
            else if (z == 1.0_dp / c) then
                if (c < 1.0_dp) then
                    y = negative_infinity(x)
                else if (c == 1.0_dp) then
                    y = scifort_log_two - log(sigma)
                else
                    y = positive_infinity(x)
                end if
            else
                t = 1.0_dp - c * z
                lt = log(t)
                lu = lt / c
                y = scifort_log_two + (1.0_dp / c - 1.0_dp) * lt - &
                    2.0_dp * log1p_safe(exp(lu)) - log(sigma)
            end if
        end if
    end function genhalflogistic_logpdf

    pure elemental function genhalflogistic_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = genhalflogistic_logpdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function genhalflogistic_pdf

    pure elemental function genhalflogistic_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: lt, lu, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x < mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z >= 1.0_dp / c) then
                y = negative_infinity(x)
            else
                lt = log1p_safe(-c * z)
                lu = lt / c
                y = scifort_log_two + lu - log1p_safe(exp(lu))
            end if
        end if
    end function genhalflogistic_logsf

    pure elemental function genhalflogistic_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = genhalflogistic_logsf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function genhalflogistic_sf

    pure elemental function genhalflogistic_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = genhalflogistic_logsf(x, c, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == negative_infinity(ls)) then
            y = 1.0_dp
        else
            y = -expm1_safe(ls)
        end if
    end function genhalflogistic_cdf

    pure elemental function genhalflogistic_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ls
        ls = genhalflogistic_logsf(x, c, loc, scale)
        if (ieee_is_nan(ls)) then
            y = ls
        else if (ls == 0.0_dp) then
            y = negative_infinity(x)
        else if (ls == negative_infinity(ls)) then
            y = 0.0_dp
        else
            y = log(-expm1_safe(ls))
        end if
    end function genhalflogistic_logcdf

    pure elemental function genhalflogistic_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: lratio, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            z = 1.0_dp / c
            x = affine(mu, sigma, z, p)
        else
            lratio = log1p_safe(-p) - log1p_safe(p)
            z = -expm1_safe(c * lratio) / c
            x = affine(mu, sigma, z, p)
        end if
    end function genhalflogistic_ppf

    pure elemental function genhalflogistic_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: lratio, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            z = 1.0_dp / c
            x = affine(mu, sigma, z, p)
        else
            lratio = log(p) - log(2.0_dp - p)
            z = -expm1_safe(c * lratio) / c
            x = affine(mu, sigma, z, p)
        end if
    end function genhalflogistic_isf

    pure elemental function valid_shape(c) result(ok)
        real(dp), intent(in) :: c !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location
        real(dp), intent(in), optional :: scale !! requested scale
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized quantile
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine
end module scifort_genhalflogistic
