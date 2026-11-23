! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright 1984, 1987, 1988, 1989 by Stephen L. Moshier
! Copyright 1984, 1995, 2000 by Stephen L. Moshier
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
! All rights reserved.
!
! Inverse-normal rational approximations and Horner evaluation adapted from
! SciPy v1.15.3 xsf/cephes/ndtri.h, polevl.h and const.h (Cephes).
! See CODE_PROVENANCE.md and LICENSES/scipy-cephes-ndtri-BSD-3-Clause.txt.

module scifort_normal
    use scifort_constants, only : scifort_inv_sqrt_two_pi, &
        scifort_log_sqrt_two_pi, scifort_sqrt_two
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, &
        quiet_nan, valid_loc_scale
    implicit none
    private

    public :: normal_cdf
    public :: normal_isf
    public :: normal_logcdf
    public :: normal_logpdf
    public :: normal_logsf
    public :: normal_pdf
    public :: normal_ppf
    public :: normal_sf

contains

    pure elemental function normal_pdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = scifort_inv_sqrt_two_pi * exp(-0.5_dp * z * z) / sigma
    end function normal_pdf

    pure elemental function normal_logpdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = -0.5_dp * z * z - log(sigma) - scifort_log_sqrt_two_pi
    end function normal_logpdf

    pure elemental function normal_cdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = standard_normal_cdf(z)
    end function normal_cdf

    pure elemental function normal_sf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = standard_normal_sf(z)
    end function normal_sf

    pure elemental function normal_logcdf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = standard_normal_logcdf(z)
    end function normal_logcdf

    pure elemental function normal_logsf(x, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        y = standard_normal_logcdf(-z)
    end function normal_logsf

    pure elemental function normal_ppf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu + sigma * standard_normal_ppf(p)
        end if
    end function normal_ppf

    pure elemental function normal_isf(p, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
            return
        end if

        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else
            x = mu - sigma * standard_normal_ppf(p)
        end if
    end function normal_isf

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

    pure elemental function standard_normal_cdf(z) result(p)
        real(dp), intent(in) :: z !! standardized value
        real(dp) :: p

        p = 0.5_dp * erfc(-z / scifort_sqrt_two)
    end function standard_normal_cdf

    pure elemental function standard_normal_sf(z) result(p)
        real(dp), intent(in) :: z !! standardized value
        real(dp) :: p

        p = 0.5_dp * erfc(z / scifort_sqrt_two)
    end function standard_normal_sf

    pure elemental function standard_normal_logcdf(z) result(y)
        real(dp), intent(in) :: z !! standardized value
        real(dp) :: y

        real(dp) :: inv_z2
        real(dp) :: series
        real(dp) :: tail

        if (z < -20.0_dp) then
            if (-z > sqrt(2.0_dp) * sqrt(huge(1.0_dp))) then
                y = negative_infinity(z)
                return
            else if (-z > sqrt(huge(1.0_dp))) then
                inv_z2 = 0.0_dp
            else
                inv_z2 = 1.0_dp / (z * z)
            end if
            series = 1.0_dp - inv_z2 * (1.0_dp - inv_z2 * &
                (3.0_dp - inv_z2 * (15.0_dp - inv_z2 * &
                (105.0_dp - inv_z2 * (945.0_dp - 10395.0_dp * inv_z2)))))
            y = (-0.5_dp * z) * z - log(-z) - scifort_log_sqrt_two_pi + log(series)
        else if (z <= 0.0_dp) then
            y = log(standard_normal_cdf(z))
        else
            tail = standard_normal_sf(z)
            y = log1p_safe(-tail)
        end if
    end function standard_normal_logcdf

    pure elemental function standard_normal_ppf(p) result(x)
        implicit none
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp) :: x
        ! Cephes ndtri coefficients, in descending powers; see provenance.
        real(dp), parameter :: p0(5) = [ &
            -5.99633501014107895267e1_dp, &
            9.80010754185999661536e1_dp, &
            -5.66762857469070293439e1_dp, &
            1.39312609387279679503e1_dp, &
            -1.23916583867381258016e0_dp ]
        real(dp), parameter :: q0(8) = [ &
            1.95448858338141759834e0_dp, &
            4.67627912898881538453e0_dp, &
            8.63602421390890590575e1_dp, &
            -2.25462687854119370527e2_dp, &
            2.00260212380060660359e2_dp, &
            -8.20372256168333339912e1_dp, &
            1.59056225126211695515e1_dp, &
            -1.18331621121330003142e0_dp ]
        real(dp), parameter :: p1(9) = [ &
            4.05544892305962419923e0_dp, &
            3.15251094599893866154e1_dp, &
            5.71628192246421288162e1_dp, &
            4.40805073893200834700e1_dp, &
            1.46849561928858024014e1_dp, &
            2.18663306850790267539e0_dp, &
            -1.40256079171354495875e-1_dp, &
            -3.50424626827848203418e-2_dp, &
            -8.57456785154685413611e-4_dp ]
        real(dp), parameter :: q1(8) = [ &
            1.57799883256466749731e1_dp, &
            4.53907635128879210584e1_dp, &
            4.13172038254672030440e1_dp, &
            1.50425385692907503408e1_dp, &
            2.50464946208309415979e0_dp, &
            -1.42182922854787788574e-1_dp, &
            -3.80806407691578277194e-2_dp, &
            -9.33259480895457427372e-4_dp ]
        real(dp), parameter :: p2(9) = [ &
            3.23774891776946035970e0_dp, &
            6.91522889068984211695e0_dp, &
            3.93881025292474443415e0_dp, &
            1.33303460815807542389e0_dp, &
            2.01485389549179081538e-1_dp, &
            1.23716634817820021358e-2_dp, &
            3.01581553508235416007e-4_dp, &
            2.65806974686737550832e-6_dp, &
            6.23974539184983293730e-9_dp ]
        real(dp), parameter :: q2(8) = [ &
            6.02427039364742014255e0_dp, &
            3.67983563856160859403e0_dp, &
            1.37702099489081330271e0_dp, &
            2.16236993594496635890e-1_dp, &
            1.34204006088543189037e-2_dp, &
            3.28014464682127739104e-4_dp, &
            2.89247864745380683936e-6_dp, &
            6.79019408009981274425e-9_dp ]
        real(dp), parameter :: tail_cutoff = 0.13533528323661269189_dp
        real(dp), parameter :: sqrt_two_pi = 2.50662827463100050242_dp
        real(dp) :: y, z, y2, leading, correction
        logical :: upper_tail

        if (p <= 0.0_dp) then
            x = negative_infinity(p)
            return
        else if (p >= 1.0_dp) then
            x = positive_infinity(p)
            return
        end if

        y = p
        upper_tail = .false.
        if (y > 1.0_dp-tail_cutoff) then
            y = 1.0_dp-y
            upper_tail = .true.
        end if
        if (y > tail_cutoff) then
            y = y-0.5_dp
            y2 = y*y
            x = (y+y*(y2*normal_polynomial(y2,p0)/normal_monic_polynomial(y2,q0))) &
                *sqrt_two_pi
            return
        end if

        z = sqrt(-2.0_dp*log(y))
        leading = z-log(z)/z
        y = 1.0_dp/z
        if (z < 8.0_dp) then
            correction = y*normal_polynomial(y,p1)/normal_monic_polynomial(y,q1)
        else
            correction = y*normal_polynomial(y,p2)/normal_monic_polynomial(y,q2)
        end if
        x = leading-correction
        if (.not. upper_tail) x = -x
    end function standard_normal_ppf

    pure function normal_polynomial(x, coefficients) result(y)
        implicit none
        real(dp), intent(in) :: x !! polynomial argument
        real(dp), intent(in) :: coefficients(:) !! coefficients in descending powers, nonempty
        real(dp) :: y
        integer :: i
        y = coefficients(1)
        do i=2,size(coefficients)
            y = y*x+coefficients(i)
        end do
    end function normal_polynomial

    pure function normal_monic_polynomial(x, coefficients) result(y)
        implicit none
        real(dp), intent(in) :: x !! polynomial argument
        real(dp), intent(in) :: coefficients(:) !! descending coefficients excluding leading 1
        real(dp) :: y
        integer :: i
        y = x+coefficients(1)
        do i=2,size(coefficients)
            y = y*x+coefficients(i)
        end do
    end function normal_monic_polynomial

end module scifort_normal
