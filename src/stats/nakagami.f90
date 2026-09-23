! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Nakagami distribution matching scipy.stats.nakagami.
! In standardized coordinates z=(x-loc)/scale >= 0,
! f(z)=2*nu**nu/Gamma(nu) * z**(2*nu-1) * exp(-nu*z**2), nu>0.
module scifort_nakagami
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_incomplete_gamma, only : gammaincinv, gammainccinv, log_gammainc, log_gammaincc
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private
    public :: nakagami_pdf, nakagami_logpdf, nakagami_cdf, nakagami_sf
    public :: nakagami_logcdf, nakagami_logsf, nakagami_ppf, nakagami_isf
contains
    pure elemental function nakagami_logpdf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, a
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(nu))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z < 0.0_dp) then
                y = negative_infinity(x)
            else if (z == 0.0_dp) then
                a = 2.0_dp * nu - 1.0_dp
                if (a > 0.0_dp) then
                    y = negative_infinity(x)
                else if (a < 0.0_dp) then
                    y = positive_infinity(x)
                else
                    y = scifort_log_two + nu * log(nu) - log_gamma(nu) - log(sigma)
                end if
            else
                y = scifort_log_two + nu * log(nu) - log_gamma(nu) + &
                    (2.0_dp * nu - 1.0_dp) * log(z) - nu * z * z - log(sigma)
            end if
        end if
    end function nakagami_logpdf

    pure elemental function nakagami_pdf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = nakagami_logpdf(x, nu, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function nakagami_pdf

    pure elemental function nakagami_logcdf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, t
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(nu))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = negative_infinity(x)
        else if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (z > sqrt(huge(1.0_dp) / nu)) then
                y = 0.0_dp
            else
                t = nu * z * z
                y = log_gammainc(nu, t)
            end if
        end if
    end function nakagami_logcdf

    pure elemental function nakagami_logsf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z, t
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(nu))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x <= mu) then
            y = 0.0_dp
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (z > sqrt(huge(1.0_dp) / nu)) then
                y = negative_infinity(x)
            else
                t = nu * z * z
                y = log_gammaincc(nu, t)
            end if
        end if
    end function nakagami_logsf

    pure elemental function nakagami_cdf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = nakagami_logcdf(x, nu, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function nakagami_cdf

    pure elemental function nakagami_sf(x, nu, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = nakagami_logsf(x, nu, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function nakagami_sf

    pure elemental function nakagami_ppf(p, nu, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(nu)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            z = sqrt(gammaincinv(nu, p) / nu)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function nakagami_ppf

    pure elemental function nakagami_isf(p, nu, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(nu)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else
            z = sqrt(gammainccinv(nu, p) / nu)
            x = affine_positive(mu, sigma, z, p)
        end if
    end function nakagami_isf

    pure elemental function valid_shape(nu) result(ok)
        real(dp), intent(in) :: nu !! candidate shape parameter
        logical :: ok
        ok = ieee_is_finite(nu) .and. nu > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        if (present(loc)) mu = loc
        sigma = 1.0_dp
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure elemental function affine_positive(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! nonnegative standardized quantile
        real(dp), intent(in) :: seed !! value used to form NaN if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            x = positive_infinity(seed)
        else if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = positive_infinity(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_positive
end module scifort_nakagami
