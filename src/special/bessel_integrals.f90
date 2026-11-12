! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Integrals involving order-zero cylindrical and modified Bessel functions.
! The implementation is independent and follows the defining integrals and
! convergent series in NIST DLMF 10.22.39-10.22.40 and 10.43.5, with direct
! quadrature and asymptotic tail integration away from the origin.
module scifort_bessel_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_classic, only : j0, j1, k0, y0, y1
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431_dp
    real(dp), parameter :: gl16_nodes(8) = [ &
        0.09501250983763744051_dp, &
        0.2816035507792589154_dp, &
        0.4580167776572273697_dp, &
        0.6178762444026437706_dp, &
        0.7554044083550029987_dp, &
        0.8656312023878317552_dp, &
        0.9445750230732326003_dp, &
        0.9894009349916499385_dp]
    real(dp), parameter :: gl16_weights(8) = [ &
        0.1894506104550686409_dp, &
        0.1826034150449236393_dp, &
        0.1691565193950026469_dp, &
        0.1495959888165767082_dp, &
        0.1246289712555340712_dp, &
        0.09515851168249260528_dp, &
        0.06225352393864745648_dp, &
        0.02715245941175417621_dp]
    integer, parameter :: max_series_terms = 2048
    real(dp), parameter :: oscillatory_asymptotic_cutoff = 40.0_dp

    public :: it2i0k0
    public :: it2j0y0
    public :: iti0k0
    public :: itj0y0

contains

    pure elemental subroutine itj0y0(x, ij0, iy0)
        real(dp), intent(in) :: x !! real upper integration limit
        real(dp), intent(out) :: ij0 !! integral of J0(t) from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0(t) from zero to x; NaN for x < 0

        real(dp) :: ax

        if (ieee_is_nan(x)) then
            ij0 = x
            iy0 = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                ij0 = 1.0_dp
                iy0 = 0.0_dp
            else
                ij0 = -1.0_dp
                iy0 = quiet_nan(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            ij0 = x
            iy0 = x
            return
        end if

        ax = abs(x)
        call j0_y0_integrals_positive(ax, ij0, iy0)
        if (x < 0.0_dp) then
            ij0 = -ij0
            iy0 = quiet_nan(x)
        end if
    end subroutine itj0y0

    pure elemental subroutine it2j0y0(x, ij0, iy0)
        real(dp), intent(in) :: x !! real argument; second integral requires x > 0
        real(dp), intent(out) :: ij0 !! integral of (1-J0(t))/t from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0(t)/t from x to infinity

        real(dp) :: ax

        if (ieee_is_nan(x)) then
            ij0 = x
            iy0 = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            ij0 = positive_infinity(x)
            if (x > 0.0_dp) then
                iy0 = 0.0_dp
            else
                iy0 = quiet_nan(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            ij0 = 0.0_dp
            iy0 = -positive_infinity(x)
            return
        end if

        ax = abs(x)
        call j0_y0_related_integrals_positive(ax, ij0, iy0)
        if (x < 0.0_dp) iy0 = quiet_nan(x)
    end subroutine it2j0y0

    pure elemental subroutine iti0k0(x, ii0, ik0)
        real(dp), intent(in) :: x !! real upper integration limit
        real(dp), intent(out) :: ii0 !! integral of I0(t) from zero to x
        real(dp), intent(out) :: ik0 !! integral of K0(t) from zero to x; NaN for x < 0

        real(dp) :: ax
        real(dp) :: dummy

        if (ieee_is_nan(x)) then
            ii0 = x
            ik0 = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                ii0 = positive_infinity(x)
                ik0 = 0.5_dp * scifort_pi
            else
                ii0 = -positive_infinity(x)
                ik0 = quiet_nan(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            ii0 = x
            ik0 = x
            return
        end if

        ax = abs(x)
        call i0_related_series(ax, ii0, dummy)
        if (x < 0.0_dp) ii0 = -ii0

        if (x < 0.0_dp) then
            ik0 = quiet_nan(x)
        else if (ax <= 6.0_dp) then
            ik0 = k0_integral_series(ax)
        else if (ax >= oscillatory_asymptotic_cutoff) then
            ik0 = 0.5_dp * scifort_pi
        else
            call k0_tail_integrals(ax, ik0, dummy)
            ik0 = 0.5_dp * scifort_pi - ik0
        end if
    end subroutine iti0k0

    pure elemental subroutine it2i0k0(x, ii0, ik0)
        real(dp), intent(in) :: x !! real argument; second integral requires x > 0
        real(dp), intent(out) :: ii0 !! integral of (I0(t)-1)/t from zero to x
        real(dp), intent(out) :: ik0 !! integral of K0(t)/t from x to infinity

        real(dp) :: ax
        real(dp) :: dummy

        if (ieee_is_nan(x)) then
            ii0 = x
            ik0 = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            ii0 = positive_infinity(x)
            if (x > 0.0_dp) then
                ik0 = 0.0_dp
            else
                ik0 = quiet_nan(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            ii0 = 0.0_dp
            ik0 = positive_infinity(x)
            return
        end if

        ax = abs(x)
        call i0_related_series(ax, dummy, ii0)

        if (x < 0.0_dp) then
            ik0 = quiet_nan(x)
        else if (ax <= 8.0_dp) then
            ik0 = k0_over_t_tail_series(ax)
        else
            call k0_tail_integrals(ax, dummy, ik0)
        end if
    end subroutine it2i0k0

    pure subroutine j0_y0_integrals_positive(x, ij0, iy0)
        real(dp), intent(in) :: x !! positive finite upper integration limit
        real(dp), intent(out) :: ij0 !! integral of J0 from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0 from zero to x

        real(dp) :: tail_j
        real(dp) :: tail_y

        if (x <= 1.0_dp) then
            call j0_y0_integral_series(x, ij0, iy0)
        else if (x < oscillatory_asymptotic_cutoff) then
            call j0_y0_integral_series(1.0_dp, ij0, iy0)
            call add_j0_y0_quadrature(1.0_dp, x, ij0, iy0)
        else
            call bessel_zero_order_tail(x, tail_j, tail_y)
            ij0 = 1.0_dp - tail_j
            iy0 = -tail_y
        end if
    end subroutine j0_y0_integrals_positive

    pure subroutine j0_y0_integral_series(x, ij0, iy0)
        real(dp), intent(in) :: x !! positive finite argument no larger than one
        real(dp), intent(out) :: ij0 !! integral of J0 from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0 from zero to x

        real(dp) :: add_j
        real(dp) :: add_y
        real(dp) :: harmonic
        real(dp) :: log_term
        real(dp) :: term
        real(dp) :: z
        integer :: k

        z = 0.25_dp * x * x
        log_term = log(0.5_dp * x) + euler_gamma
        term = 1.0_dp
        harmonic = 0.0_dp
        ij0 = x
        iy0 = (2.0_dp / scifort_pi) * x * (log_term - 1.0_dp)

        do k = 1, max_series_terms
            term = -term * z / real(k * k, dp)
            harmonic = harmonic + 1.0_dp / real(k, dp)
            add_j = x * term / real(2 * k + 1, dp)
            add_y = (2.0_dp / scifort_pi) * x * term / real(2 * k + 1, dp) * &
                (log_term - 1.0_dp / real(2 * k + 1, dp) - harmonic)
            ij0 = ij0 + add_j
            iy0 = iy0 + add_y
            if (abs(add_j) <= epsilon(1.0_dp) * max(1.0_dp, abs(ij0)) .and. &
                    abs(add_y) <= epsilon(1.0_dp) * max(1.0_dp, abs(iy0))) exit
        end do
    end subroutine j0_y0_integral_series

    pure subroutine add_j0_y0_quadrature(a, b, ij0, iy0)
        real(dp), intent(in) :: a !! positive lower integration limit
        real(dp), intent(in) :: b !! upper integration limit, >= a
        real(dp), intent(inout) :: ij0 !! accumulated J0 integral
        real(dp), intent(inout) :: iy0 !! accumulated Y0 integral

        real(dp) :: half
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: panel_width
        real(dp) :: right
        real(dp) :: t1
        real(dp) :: t2
        integer :: i
        integer :: p
        integer :: panels

        if (b <= a) return
        panels = max(1, ceiling(b - a))
        panel_width = (b - a) / real(panels, dp)
        do p = 0, panels - 1
            left = a + real(p, dp) * panel_width
            right = left + panel_width
            midpoint = 0.5_dp * (left + right)
            half = 0.5_dp * panel_width
            do i = 1, size(gl16_nodes)
                t1 = midpoint - half * gl16_nodes(i)
                t2 = midpoint + half * gl16_nodes(i)
                ij0 = ij0 + half * gl16_weights(i) * (j0(t1) + j0(t2))
                iy0 = iy0 + half * gl16_weights(i) * (y0(t1) + y0(t2))
            end do
        end do
    end subroutine add_j0_y0_quadrature

    pure subroutine j0_y0_related_integrals_positive(x, ij0, iy0)
        real(dp), intent(in) :: x !! positive finite argument
        real(dp), intent(out) :: ij0 !! integral of (1-J0)/t from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0/t from x to infinity

        real(dp) :: anchor_j
        real(dp) :: anchor_y
        real(dp) :: tail_j
        real(dp) :: tail_y

        if (x <= 10.0_dp) then
            call j0_y0_related_series(x, ij0, iy0)
        else if (x < oscillatory_asymptotic_cutoff) then
            call j0_y0_related_series(10.0_dp, anchor_j, anchor_y)
            ij0 = anchor_j
            iy0 = anchor_y
            call add_j0_y0_related_quadrature(10.0_dp, x, ij0, iy0)
        else
            call bessel_zero_order_over_t_tail(x, tail_j, tail_y)
            ij0 = euler_gamma + log(0.5_dp * x) + tail_j
            iy0 = tail_y
        end if
    end subroutine j0_y0_related_integrals_positive

    pure subroutine j0_y0_related_series(x, ij0, iy0)
        real(dp), intent(in) :: x !! positive finite argument no larger than ten
        real(dp), intent(out) :: ij0 !! integral of (1-J0)/t from zero to x
        real(dp), intent(out) :: iy0 !! integral of Y0/t from x to infinity

        real(dp) :: add_j
        real(dp) :: add_y
        real(dp) :: harmonic
        real(dp) :: log_half
        real(dp) :: q
        real(dp) :: sign_value
        real(dp) :: z
        integer :: k

        z = 0.25_dp * x * x
        log_half = log(0.5_dp * x)
        ij0 = 0.0_dp
        iy0 = -(log_half + euler_gamma) ** 2 / scifort_pi + scifort_pi / 6.0_dp
        q = 1.0_dp
        harmonic = 0.0_dp
        sign_value = -1.0_dp

        do k = 1, max_series_terms
            q = q * z / real(k * k, dp)
            harmonic = harmonic + 1.0_dp / real(k, dp)
            sign_value = -sign_value
            add_j = sign_value * q / real(2 * k, dp)
            add_y = -(2.0_dp / scifort_pi) * sign_value * &
                (harmonic - euler_gamma + 1.0_dp / real(2 * k, dp) - log_half) * &
                q / real(2 * k, dp)
            ij0 = ij0 + add_j
            iy0 = iy0 + add_y
            if (k > 8 .and. abs(add_j) <= epsilon(1.0_dp) * max(1.0_dp, abs(ij0)) .and. &
                    abs(add_y) <= epsilon(1.0_dp) * max(1.0_dp, abs(iy0))) exit
        end do
    end subroutine j0_y0_related_series

    pure subroutine add_j0_y0_related_quadrature(a, b, ij0, iy0)
        real(dp), intent(in) :: a !! positive lower integration limit
        real(dp), intent(in) :: b !! upper integration limit, >= a
        real(dp), intent(inout) :: ij0 !! accumulated (1-J0)/t integral
        real(dp), intent(inout) :: iy0 !! accumulated Y0/t tail integral

        real(dp) :: half
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: panel_width
        real(dp) :: right
        real(dp) :: t1
        real(dp) :: t2
        integer :: i
        integer :: p
        integer :: panels

        if (b <= a) return
        panels = max(1, ceiling(b - a))
        panel_width = (b - a) / real(panels, dp)
        do p = 0, panels - 1
            left = a + real(p, dp) * panel_width
            right = left + panel_width
            midpoint = 0.5_dp * (left + right)
            half = 0.5_dp * panel_width
            do i = 1, size(gl16_nodes)
                t1 = midpoint - half * gl16_nodes(i)
                t2 = midpoint + half * gl16_nodes(i)
                ij0 = ij0 + half * gl16_weights(i) * &
                    ((1.0_dp - j0(t1)) / t1 + (1.0_dp - j0(t2)) / t2)
                iy0 = iy0 - half * gl16_weights(i) * (y0(t1) / t1 + y0(t2) / t2)
            end do
        end do
    end subroutine add_j0_y0_related_quadrature

    pure subroutine bessel_zero_order_tail(x, tail_j, tail_y)
        real(dp), intent(in) :: x !! large positive finite argument
        real(dp), intent(out) :: tail_j !! integral of J0 from x to infinity
        real(dp), intent(out) :: tail_y !! integral of Y0 from x to infinity

        real(dp) :: acoef
        real(dp) :: bcoef
        real(dp) :: coefficient
        real(dp) :: invx
        real(dp) :: invx_power
        integer :: k
        integer :: n

        invx = 1.0_dp / x
        coefficient = -1.0_dp
        bcoef = coefficient
        acoef = invx
        invx_power = 1.0_dp
        do k = 1, 10
            n = 2 * k
            coefficient = -real((n - 1) * (n - 1), dp) * coefficient
            invx_power = invx_power * invx * invx
            bcoef = bcoef + coefficient * invx_power
            acoef = acoef - real(n + 1, dp) * coefficient * invx_power * invx
        end do
        tail_j = acoef * j0(x) + bcoef * j1(x)
        tail_y = acoef * y0(x) + bcoef * y1(x)
    end subroutine bessel_zero_order_tail

    pure subroutine bessel_zero_order_over_t_tail(x, tail_j, tail_y)
        real(dp), intent(in) :: x !! large positive finite argument
        real(dp), intent(out) :: tail_j !! integral of J0(t)/t from x to infinity
        real(dp), intent(out) :: tail_y !! integral of Y0(t)/t from x to infinity

        real(dp) :: acoef
        real(dp) :: bcoef
        real(dp) :: coefficient
        real(dp) :: invx
        real(dp) :: invx_power
        integer :: k
        integer :: n

        invx = 1.0_dp / x
        coefficient = -1.0_dp
        invx_power = invx
        bcoef = coefficient * invx_power
        acoef = 2.0_dp * invx * invx
        do k = 1, 10
            n = 2 * k + 1
            coefficient = -real((n - 1) * (n - 1), dp) * coefficient
            invx_power = invx_power * invx * invx
            bcoef = bcoef + coefficient * invx_power
            acoef = acoef - real(n + 1, dp) * coefficient * invx_power * invx
        end do
        tail_j = acoef * j0(x) + bcoef * j1(x)
        tail_y = acoef * y0(x) + bcoef * y1(x)
    end subroutine bessel_zero_order_over_t_tail

    pure subroutine i0_related_series(x, ii0, ii0_over_t)
        real(dp), intent(in) :: x !! nonnegative finite argument
        real(dp), intent(out) :: ii0 !! integral of I0 from zero to x
        real(dp), intent(out) :: ii0_over_t !! integral of (I0-1)/t from zero to x

        real(dp) :: add_i
        real(dp) :: add_over_t
        real(dp) :: q
        real(dp) :: z
        integer :: k

        if (x == 0.0_dp) then
            ii0 = 0.0_dp
            ii0_over_t = 0.0_dp
            return
        end if

        z = 0.25_dp * x * x
        q = 1.0_dp
        ii0 = x
        ii0_over_t = 0.0_dp
        do k = 1, max_series_terms
            q = q * z / real(k * k, dp)
            if (.not. ieee_is_finite(q)) then
                ii0 = positive_infinity(x)
                ii0_over_t = positive_infinity(x)
                return
            end if
            add_i = x * q / real(2 * k + 1, dp)
            add_over_t = q / real(2 * k, dp)
            ii0 = ii0 + add_i
            ii0_over_t = ii0_over_t + add_over_t
            if (.not. ieee_is_finite(ii0) .or. .not. ieee_is_finite(ii0_over_t)) then
                ii0 = positive_infinity(x)
                ii0_over_t = positive_infinity(x)
                return
            end if
            if (real(k, dp) > 0.5_dp * x .and. &
                    add_i <= epsilon(1.0_dp) * max(1.0_dp, ii0) .and. &
                    add_over_t <= epsilon(1.0_dp) * max(1.0_dp, ii0_over_t)) exit
        end do
    end subroutine i0_related_series

    pure function k0_integral_series(x) result(value)
        real(dp), intent(in) :: x !! positive finite argument no larger than six
        real(dp) :: value

        real(dp) :: add
        real(dp) :: harmonic
        real(dp) :: log_half
        real(dp) :: q
        real(dp) :: sum_value
        real(dp) :: z
        integer :: k

        log_half = log(0.5_dp * x)
        z = 0.25_dp * x * x
        sum_value = -(log_half + euler_gamma - 1.0_dp)
        q = 1.0_dp
        harmonic = 0.0_dp
        do k = 1, max_series_terms
            q = q * z / real(k * k, dp)
            harmonic = harmonic + 1.0_dp / real(k, dp)
            add = q / real(2 * k + 1, dp) * &
                (harmonic - log_half - euler_gamma + 1.0_dp / real(2 * k + 1, dp))
            sum_value = sum_value + add
            if (k > 8 .and. abs(add) <= epsilon(1.0_dp) * max(1.0_dp, abs(sum_value))) exit
        end do
        value = x * sum_value
    end function k0_integral_series

    pure function k0_over_t_tail_series(x) result(value)
        real(dp), intent(in) :: x !! positive finite argument no larger than eight
        real(dp) :: value

        real(dp) :: add
        real(dp) :: harmonic
        real(dp) :: log_half
        real(dp) :: q
        real(dp) :: z
        integer :: k

        log_half = log(0.5_dp * x)
        z = 0.25_dp * x * x
        value = 0.5_dp * (log_half + euler_gamma) ** 2 + scifort_pi ** 2 / 24.0_dp
        q = 1.0_dp
        harmonic = 0.0_dp
        do k = 1, max_series_terms
            q = q * z / real(k * k, dp)
            harmonic = harmonic + 1.0_dp / real(k, dp)
            add = -(harmonic - euler_gamma + 1.0_dp / real(2 * k, dp) - log_half) * &
                q / real(2 * k, dp)
            value = value + add
            if (k > 8 .and. abs(add) <= epsilon(1.0_dp) * max(1.0_dp, abs(value))) exit
        end do
    end function k0_over_t_tail_series

    pure subroutine k0_tail_integrals(x, tail_k0, tail_k0_over_t)
        real(dp), intent(in) :: x !! positive finite lower integration limit
        real(dp), intent(out) :: tail_k0 !! integral of K0 from x to infinity
        real(dp), intent(out) :: tail_k0_over_t !! integral of K0(t)/t from x to infinity

        real(dp) :: half
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: right
        real(dp) :: t1
        real(dp) :: t2
        real(dp) :: v1
        real(dp) :: v2
        integer :: i
        integer :: p

        tail_k0 = 0.0_dp
        tail_k0_over_t = 0.0_dp
        do p = 0, 39
            left = x + real(p, dp)
            right = left + 1.0_dp
            midpoint = 0.5_dp * (left + right)
            half = 0.5_dp
            do i = 1, size(gl16_nodes)
                t1 = midpoint - half * gl16_nodes(i)
                t2 = midpoint + half * gl16_nodes(i)
                v1 = k0(t1)
                v2 = k0(t2)
                tail_k0 = tail_k0 + half * gl16_weights(i) * (v1 + v2)
                tail_k0_over_t = tail_k0_over_t + half * gl16_weights(i) * &
                    (v1 / t1 + v2 / t2)
            end do
        end do
    end subroutine k0_tail_integrals

end module scifort_bessel_integrals
