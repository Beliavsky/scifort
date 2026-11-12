! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Wright's generalized Bessel function for nonnegative real parameters.
! The implementation is original SciFort code based on the defining series
! in NIST DLMF 10.46.1, evaluated around its largest term in the log domain.

module scifort_wright_bessel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_constants, only : scifort_pi
    use scifort_gamma_combinatorics, only : polygamma
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_special_elementary, only : digamma, gammaln
    implicit none
    private

    public :: log_wright_bessel
    public :: wright_bessel

    real(dp), parameter :: log_huge_dp = log(huge(1.0_dp))
    real(dp), parameter :: log_relative_cutoff = log(epsilon(1.0_dp)) - 4.0_dp
    real(dp), parameter :: direct_mode_limit = 1.0e10_dp

contains

    pure elemental function wright_bessel(a, b, x) result(value)
        real(dp), intent(in) :: a !! first Wright parameter, >= 0
        real(dp), intent(in) :: b !! second Wright parameter, >= 0
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: value

        real(dp) :: log_value

        log_value = log_wright_bessel(a, b, x)
        if (ieee_is_nan(log_value)) then
            value = quiet_nan(x)
        else if (.not. ieee_is_finite(log_value)) then
            if (log_value > 0.0_dp) then
                value = positive_infinity(x)
            else
                value = 0.0_dp
            end if
        else if (log_value > log_huge_dp) then
            value = positive_infinity(x)
        else
            value = exp(log_value)
        end if
    end function wright_bessel

    pure elemental function log_wright_bessel(a, b, x) result(value)
        real(dp), intent(in) :: a !! first Wright parameter, >= 0
        real(dp), intent(in) :: b !! second Wright parameter, >= 0
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: value

        real(dp) :: mode

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(x)) then
            value = quiet_nan(x)
            return
        end if
        if (a < 0.0_dp .or. b < 0.0_dp .or. x < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            value = positive_infinity(x)
            return
        end if
        if (x <= 0.0_dp) then
            if (b <= 0.0_dp) then
                value = negative_infinity(x)
            else
                value = -gammaln(b)
            end if
            return
        end if
        if (a <= 0.0_dp) then
            if (b <= 0.0_dp) then
                value = negative_infinity(x)
            else
                value = x - gammaln(b)
            end if
            return
        end if

        mode = wright_mode(a, b, x)
        if (.not. ieee_is_finite(mode)) then
            value = quiet_nan(x)
        else if (mode > direct_mode_limit) then
            value = wright_log_saddle(a, b, x, mode)
        else
            value = wright_log_sum(a, b, x, mode)
        end if
    end function log_wright_bessel

    pure function wright_mode(a, b, x) result(mode)
        real(dp), intent(in) :: a !! strictly positive first parameter
        real(dp), intent(in) :: b !! nonnegative second parameter
        real(dp), intent(in) :: x !! strictly positive finite argument
        real(dp) :: mode

        integer :: iteration
        real(dp) :: high
        real(dp) :: low
        real(dp) :: midpoint
        real(dp) :: log_x

        log_x = log(x)
        if (b <= 0.0_dp) then
            low = 1.0_dp
        else
            low = 0.0_dp
        end if

        if (wright_score(low, a, b, log_x) <= 0.0_dp) then
            mode = low
            return
        end if

        high = max(1.0_dp, 2.0_dp * low)
        do iteration = 1, 2048
            if (wright_score(high, a, b, log_x) <= 0.0_dp) exit
            low = high
            if (high >= 0.25_dp * huge(1.0_dp)) then
                mode = quiet_nan(x)
                return
            end if
            high = 2.0_dp * high
        end do
        if (wright_score(high, a, b, log_x) > 0.0_dp) then
            mode = quiet_nan(x)
            return
        end if

        do iteration = 1, 80
            midpoint = 0.5_dp * (low + high)
            if (wright_score(midpoint, a, b, log_x) > 0.0_dp) then
                low = midpoint
            else
                high = midpoint
            end if
        end do
        mode = 0.5_dp * (low + high)
    end function wright_mode

    pure function wright_score(t, a, b, log_x) result(score)
        real(dp), intent(in) :: t !! continuous summation index
        real(dp), intent(in) :: a !! strictly positive first parameter
        real(dp), intent(in) :: b !! nonnegative second parameter
        real(dp), intent(in) :: log_x !! natural logarithm of the argument
        real(dp) :: score

        score = log_x - digamma(t + 1.0_dp) - a * digamma(a * t + b)
    end function wright_score

    pure function wright_log_term(t, a, b, log_x) result(log_term)
        real(dp), intent(in) :: t !! nonnegative summation index
        real(dp), intent(in) :: a !! strictly positive first parameter
        real(dp), intent(in) :: b !! nonnegative second parameter
        real(dp), intent(in) :: log_x !! natural logarithm of the argument
        real(dp) :: log_term

        real(dp) :: gamma_argument

        gamma_argument = a * t + b
        if (gamma_argument <= 0.0_dp) then
            log_term = negative_infinity(t)
        else
            log_term = t * log_x - gammaln(t + 1.0_dp) - gammaln(gamma_argument)
        end if
    end function wright_log_term

    pure function wright_log_sum(a, b, x, mode) result(log_sum)
        real(dp), intent(in) :: a !! strictly positive first parameter
        real(dp), intent(in) :: b !! nonnegative second parameter
        real(dp), intent(in) :: x !! strictly positive finite argument
        real(dp), intent(in) :: mode !! continuous location of the largest term
        real(dp) :: log_sum

        integer :: offset
        integer(int64) :: k
        integer(int64) :: k_center
        integer(int64) :: k_min
        integer(int64) :: k_peak
        integer(int64) :: trial
        real(dp) :: best_log
        real(dp) :: compensation
        real(dp) :: log_term
        real(dp) :: log_x
        real(dp) :: scaled_sum
        real(dp) :: scaled_term
        real(dp) :: sum_increment
        real(dp) :: sum_next

        if (b <= 0.0_dp) then
            k_min = 1_int64
        else
            k_min = 0_int64
        end if

        k_center = max(k_min, int(mode, int64))
        k_peak = k_center
        log_x = log(x)
        best_log = negative_infinity(x)
        do offset = -3, 4
            if (offset < 0 .and. k_center < k_min - int(offset, int64)) cycle
            trial = k_center + int(offset, int64)
            if (trial < k_min) cycle
            log_term = wright_log_term(real(trial, dp), a, b, log_x)
            if (log_term > best_log) then
                best_log = log_term
                k_peak = trial
            end if
        end do

        if (.not. ieee_is_finite(best_log)) then
            log_sum = best_log
            return
        end if

        scaled_sum = 1.0_dp
        compensation = 0.0_dp

        k = k_peak - 1_int64
        do while (k >= k_min)
            log_term = wright_log_term(real(k, dp), a, b, log_x) - best_log
            if (log_term < log_relative_cutoff) exit
            scaled_term = exp(log_term)
            sum_increment = scaled_term - compensation
            sum_next = scaled_sum + sum_increment
            compensation = (sum_next - scaled_sum) - sum_increment
            scaled_sum = sum_next
            if (k == k_min) exit
            k = k - 1_int64
        end do

        k = k_peak + 1_int64
        do
            log_term = wright_log_term(real(k, dp), a, b, log_x) - best_log
            if (log_term < log_relative_cutoff) exit
            scaled_term = exp(log_term)
            sum_increment = scaled_term - compensation
            sum_next = scaled_sum + sum_increment
            compensation = (sum_next - scaled_sum) - sum_increment
            scaled_sum = sum_next
            if (k == huge(k)) then
                log_sum = wright_log_saddle(a, b, x, mode)
                return
            end if
            k = k + 1_int64
        end do

        log_sum = best_log + log(scaled_sum)
    end function wright_log_sum

    pure function wright_log_saddle(a, b, x, mode) result(log_sum)
        real(dp), intent(in) :: a !! strictly positive first parameter
        real(dp), intent(in) :: b !! nonnegative second parameter
        real(dp), intent(in) :: x !! strictly positive finite argument
        real(dp), intent(in) :: mode !! continuous saddle point
        real(dp) :: log_sum

        real(dp) :: correction
        real(dp) :: curvature
        real(dp) :: derivative3
        real(dp) :: derivative4
        real(dp) :: gamma_argument
        real(dp) :: log_peak

        gamma_argument = a * mode + b
        log_peak = wright_log_term(mode, a, b, log(x))
        curvature = polygamma(1, mode + 1.0_dp) + a * a * polygamma(1, gamma_argument)
        derivative3 = -polygamma(2, mode + 1.0_dp) - a**3 * polygamma(2, gamma_argument)
        derivative4 = -polygamma(3, mode + 1.0_dp) - a**4 * polygamma(3, gamma_argument)

        if (curvature <= 0.0_dp .or. .not. ieee_is_finite(curvature)) then
            log_sum = quiet_nan(x)
            return
        end if

        correction = derivative4 / (8.0_dp * curvature**2) + &
            5.0_dp * derivative3**2 / (24.0_dp * curvature**3)
        if (correction <= -1.0_dp) then
            log_sum = quiet_nan(x)
            return
        end if

        log_sum = log_peak + 0.5_dp * log(2.0_dp * scifort_pi / curvature) + log(1.0_dp + correction)
    end function wright_log_saddle

end module scifort_wright_bessel
