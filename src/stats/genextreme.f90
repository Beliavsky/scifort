! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Generalized extreme-value distribution matching scipy.stats.genextreme.
! SciPy uses the sign convention F(z)=exp(-(1-c*z)**(1/c)).
module scifort_genextreme
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: genextreme_cdf, genextreme_isf, genextreme_logcdf, genextreme_logpdf
    public :: genextreme_logsf, genextreme_pdf, genextreme_ppf, genextreme_sf

contains

    pure elemental function genextreme_logpdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logt, logu, mu, sigma, t, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x)) then
            y = negative_infinity(x)
        else
            z = (x - mu) / sigma
            if (c == 0.0_dp) then
                logt = -z
                if (logt > log(huge(1.0_dp))) then
                    y = negative_infinity(x)
                else
                    t = exp(logt)
                    y = -t + logt - log(sigma)
                end if
            else
                u = 1.0_dp - c * z
                if (u < 0.0_dp) then
                    y = negative_infinity(x)
                else if (u == 0.0_dp) then
                    if (c < 0.0_dp .or. c < 1.0_dp) then
                        y = negative_infinity(x)
                    else if (c == 1.0_dp) then
                        y = -log(sigma)
                    else
                        y = positive_infinity(x)
                    end if
                else
                    logu = log1p_safe(-c * z)
                    logt = logu / c
                    if (logt > log(huge(1.0_dp))) then
                        y = negative_infinity(x)
                    else
                        if (logt < log(tiny(1.0_dp))) then
                            t = 0.0_dp
                        else
                            t = exp(logt)
                        end if
                        y = -t + logt - logu - log(sigma)
                    end if
                end if
            end if
        end if
    end function genextreme_logpdf

    pure elemental function genextreme_pdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = genextreme_logpdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else if (ly == positive_infinity(ly)) then
            y = positive_infinity(ly)
        else
            y = exp(ly)
        end if
    end function genextreme_pdf

    pure elemental function genextreme_logcdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y
        real(dp) :: logt, mu, sigma, u, z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c))) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (x == negative_infinity(x)) then
            y = negative_infinity(x)
        else if (x == positive_infinity(x)) then
            y = 0.0_dp
        else
            z = (x - mu) / sigma
            if (c == 0.0_dp) then
                logt = -z
            else
                u = 1.0_dp - c * z
                if (u <= 0.0_dp) then
                    if (c > 0.0_dp) then
                        y = 0.0_dp
                    else
                        y = negative_infinity(x)
                    end if
                    return
                end if
                logt = log1p_safe(-c * z) / c
            end if
            if (logt > log(huge(1.0_dp))) then
                y = negative_infinity(x)
            else if (logt < log(tiny(1.0_dp))) then
                y = 0.0_dp
            else
                y = -exp(logt)
            end if
        end if
    end function genextreme_logcdf

    pure elemental function genextreme_cdf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, ly
        ly = genextreme_logcdf(x, c, loc, scale)
        if (ieee_is_nan(ly)) then
            y = ly
        else if (ly == negative_infinity(ly)) then
            y = 0.0_dp
        else
            y = exp(ly)
        end if
    end function genextreme_cdf

    pure elemental function genextreme_sf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = genextreme_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 1.0_dp
        else if (lc == 0.0_dp) then
            y = 0.0_dp
        else
            y = -expm1_safe(lc)
        end if
    end function genextreme_sf

    pure elemental function genextreme_logsf(x, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, lc
        lc = genextreme_logcdf(x, c, loc, scale)
        if (ieee_is_nan(lc)) then
            y = lc
        else if (lc == negative_infinity(lc)) then
            y = 0.0_dp
        else if (lc == 0.0_dp) then
            y = negative_infinity(x)
        else
            y = log(-expm1_safe(lc))
        end if
    end function genextreme_logsf

    pure elemental function genextreme_ppf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0,1]
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, g, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            call lower_endpoint(c, mu, sigma, p, x)
        else if (p == 1.0_dp) then
            call upper_endpoint(c, mu, sigma, p, x)
        else
            g = -log(-log(p))
            z = transform_gumbel(g, c, p)
            x = affine_value(mu, sigma, z, p)
        end if
    end function genextreme_ppf

    pure elemental function genextreme_isf(p, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, g, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. ieee_is_finite(c)) .or. &
            .not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            call lower_endpoint(c, mu, sigma, p, x)
        else if (p == 0.0_dp) then
            call upper_endpoint(c, mu, sigma, p, x)
        else
            g = -log(-log1p_safe(-p))
            z = transform_gumbel(g, c, p)
            x = affine_value(mu, sigma, z, p)
        end if
    end function genextreme_isf

    pure elemental function transform_gumbel(g, c, seed) result(z)
        real(dp), intent(in) :: g !! standard Gumbel quantile
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in) :: seed !! probability used to construct infinities
        real(dp) :: z, t
        if (c == 0.0_dp) then
            z = g
        else
            t = -c * g
            if (t > log(huge(1.0_dp))) then
                if (c > 0.0_dp) then
                    z = negative_infinity(seed)
                else
                    z = positive_infinity(seed)
                end if
            else
                z = -expm1_safe(t) / c
            end if
        end if
    end function transform_gumbel

    pure elemental subroutine lower_endpoint(c, mu, sigma, seed, x)
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: seed !! probability used to construct infinity
        real(dp), intent(out) :: x !! lower support endpoint
        if (c < 0.0_dp) then
            x = affine_value(mu, sigma, 1.0_dp / c, seed)
        else
            x = negative_infinity(seed)
        end if
    end subroutine lower_endpoint

    pure elemental subroutine upper_endpoint(c, mu, sigma, seed, x)
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: seed !! probability used to construct infinity
        real(dp), intent(out) :: x !! upper support endpoint
        if (c > 0.0_dp) then
            x = affine_value(mu, sigma, 1.0_dp / c, seed)
        else
            x = positive_infinity(seed)
        end if
    end subroutine upper_endpoint

    pure elemental function affine_value(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location
        real(dp), intent(in) :: sigma !! positive scale
        real(dp), intent(in) :: z !! standardized value
        real(dp), intent(in) :: seed !! value used to construct infinity if needed
        real(dp) :: x
        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else if (abs(z) > (huge(1.0_dp) - abs(mu)) / sigma) then
            if (z < 0.0_dp) then
                x = negative_infinity(seed)
            else
                x = positive_infinity(seed)
            end if
        else
            x = mu + sigma * z
        end if
    end function affine_value

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_genextreme
