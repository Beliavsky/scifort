! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Scaled modified-Bessel-I helpers for nonnegative real order and argument.
! The power series is DLMF 10.25.2; the large-argument expansion is
! DLMF 10.40.1. These routines are intended primarily for directional
! statistics, where the order is integer or half-integer.

module scifort_bessel_i
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    implicit none
    private

    public :: besseli_ratio
    public :: log_besseli_scaled

contains

    pure elemental function log_besseli_scaled(nu, x) result(y)
        real(dp), intent(in) :: nu !! modified-Bessel order, finite and >= 0
        real(dp), intent(in) :: x !! nonnegative finite argument
        real(dp) :: y

        integer, parameter :: max_asym_terms = 40
        integer, parameter :: max_series_terms = 200000
        integer :: k
        real(dp) :: log_sum_scaled
        real(dp) :: log_term
        real(dp) :: log_x2
        real(dp) :: max_log
        real(dp) :: mu
        real(dp) :: odd
        real(dp) :: rel
        real(dp) :: series_sum
        real(dp) :: term

        if (ieee_is_nan(nu) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(nu) .or. &
                .not. ieee_is_finite(x) .or. nu < 0.0_dp .or. x < 0.0_dp) then
            y = quiet_nan(x)
            return
        end if
        if (x == 0.0_dp) then
            if (nu == 0.0_dp) then
                y = 0.0_dp
            else
                y = negative_infinity(x)
            end if
            return
        end if

        ! For x large compared with the order, the exponentially scaled
        ! asymptotic series is both fast and cancellation-free.
        if (x > max(50.0_dp, 2.0_dp * (nu + 1.0_dp) * (nu + 1.0_dp))) then
            mu = 4.0_dp * nu * nu
            term = 1.0_dp
            series_sum = 1.0_dp
            do k = 1, max_asym_terms
                odd = real(2 * k - 1, dp)
                term = term * (-(mu - odd * odd)) / (real(k, dp) * 8.0_dp * x)
                if (abs(term) > abs(series_sum) .and. k > 4) exit
                series_sum = series_sum + term
                if (abs(term) <= 4.0_dp * epsilon(1.0_dp) * abs(series_sum)) exit
            end do
            if (series_sum > 0.0_dp .and. ieee_is_finite(series_sum)) then
                y = -0.5_dp * log(2.0_dp * scifort_pi * x) + log(series_sum)
                return
            end if
        end if

        ! I_nu(x) = (x/2)^nu/Gamma(nu+1) * sum_k
        ! (x^2/4)^k/(k! Gamma(nu+k+1)/Gamma(nu+1)).  Accumulate the
        ! positive sum by log-sum-exp so the series remains safe even when
        ! its unscaled terms would overflow.
        log_term = 0.0_dp
        max_log = 0.0_dp
        log_sum_scaled = 0.0_dp
        series_sum = 1.0_dp
        log_x2 = 2.0_dp * log(x) - log(4.0_dp)
        do k = 1, max_series_terms
            log_term = log_term + log_x2 - log(real(k, dp)) - log(nu + real(k, dp))
            if (log_term > max_log) then
                series_sum = series_sum * exp(max_log - log_term) + 1.0_dp
                max_log = log_term
            else
                series_sum = series_sum + exp(log_term - max_log)
            end if
            log_sum_scaled = max_log + log(series_sum)
            rel = exp(log_term - log_sum_scaled)
            if (k > max(10, int(0.5_dp * x)) .and. &
                    rel <= 4.0_dp * epsilon(1.0_dp)) exit
        end do

        y = nu * log(0.5_dp * x) - log_gamma(nu + 1.0_dp) - x + log_sum_scaled
    end function log_besseli_scaled

    pure elemental function besseli_ratio(nu, x) result(r)
        real(dp), intent(in) :: nu !! denominator order in I_(nu+1)(x)/I_nu(x), >= 0
        real(dp), intent(in) :: x !! nonnegative finite argument
        real(dp) :: r

        real(dp) :: a
        real(dp) :: b

        if (ieee_is_nan(nu) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(nu) .or. &
                .not. ieee_is_finite(x) .or. nu < 0.0_dp .or. x < 0.0_dp) then
            r = quiet_nan(x)
        else if (x == 0.0_dp) then
            r = 0.0_dp
        else if (x < sqrt(tiny(1.0_dp))) then
            r = x / (2.0_dp * (nu + 1.0_dp))
        else
            a = log_besseli_scaled(nu + 1.0_dp, x)
            b = log_besseli_scaled(nu, x)
            r = exp(a - b)
        end if
    end function besseli_ratio

end module scifort_bessel_i
