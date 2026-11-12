! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real-axis Bessel derivatives from the standard shifted-order identities.
! See DLMF 10.6.7 for J/Y and DLMF 10.29.5 for modified I/K.

module scifort_bessel_derivatives
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_jy, only : jv, yv
    use scifort_bessel_order, only : iv, kv
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter :: family_j = 1
    integer, parameter :: family_y = 2
    integer, parameter :: family_i = 3
    integer, parameter :: family_k = 4
    integer, parameter :: recurrence_weight_limit = 1000

    public :: ivp
    public :: jvp
    public :: kvp
    public :: yvp

contains

    pure elemental function jvp(v, x, n) result(y)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by jv for this order
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: y
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            y = quiet_nan(x)
        else if (derivative_order == 0) then
            y = jv(v, x)
        else
            y = bessel_derivative_sum(v, x, derivative_order, family_j, -1)
        end if
    end function jvp

    pure elemental function yvp(v, x, n) result(y)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! nonnegative real argument accepted by yv
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: y
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            y = quiet_nan(x)
        else if (derivative_order == 0) then
            y = yv(v, x)
        else
            y = bessel_derivative_sum(v, x, derivative_order, family_y, -1)
        end if
    end function yvp

    pure elemental function ivp(v, x, n) result(y)
        real(dp), intent(in) :: v !! real modified-Bessel-I order
        real(dp), intent(in) :: x !! real argument accepted by iv for this order
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: y
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            y = quiet_nan(x)
        else if (derivative_order == 0) then
            y = iv(v, x)
        else
            y = bessel_derivative_sum(v, x, derivative_order, family_i, 1)
        end if
    end function ivp

    pure elemental function kvp(v, x, n) result(y)
        real(dp), intent(in) :: v !! real modified-Bessel-K order
        real(dp), intent(in) :: x !! nonnegative real argument accepted by kv
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: y
        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n
        if (derivative_order < 0) then
            y = quiet_nan(x)
        else if (derivative_order == 0) then
            y = kv(v, x)
        else
            y = bessel_derivative_sum(v, x, derivative_order, family_k, 1)
            if (mod(derivative_order, 2) /= 0) y = -y
        end if
    end function kvp

    pure function bessel_derivative_sum(v, x, n, family, phase) result(y)
        real(dp), intent(in) :: v !! real Bessel order before differentiation
        real(dp), intent(in) :: x !! real Bessel argument
        integer, intent(in) :: n !! positive derivative order
        integer, intent(in) :: family !! internal Bessel-family selector
        integer, intent(in) :: phase !! shifted-order coefficient phase, -1 or +1
        real(dp) :: y

        logical :: nan_seen
        logical :: negative_infinity_seen
        logical :: positive_infinity_seen
        integer :: k
        real(dp) :: base_value
        real(dp) :: log_weight
        real(dp) :: max_term
        real(dp) :: order
        real(dp) :: scaled_sum
        real(dp) :: term
        real(dp) :: weight

        nan_seen = .false.
        negative_infinity_seen = .false.
        positive_infinity_seen = .false.
        max_term = 0.0_dp
        scaled_sum = 0.0_dp

        if (n <= recurrence_weight_limit) then
            weight = scale(1.0_dp, -n)
            do k = 0, n
                if (k > 0) then
                    weight = real(phase, dp) * weight * real(n - k + 1, dp) / real(k, dp)
                end if
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = bessel_family_value(family, order, x)
                call accumulate_weighted_value(weight, base_value, max_term, scaled_sum, &
                    positive_infinity_seen, negative_infinity_seen, nan_seen)
            end do
        else
            do k = 0, n
                log_weight = log_gamma(real(n, dp) + 1.0_dp) - &
                    log_gamma(real(k, dp) + 1.0_dp) - &
                    log_gamma(real(n - k, dp) + 1.0_dp) - real(n, dp) * log(2.0_dp)
                weight = exp(log_weight)
                if (phase < 0 .and. mod(k, 2) /= 0) weight = -weight
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = bessel_family_value(family, order, x)
                call accumulate_weighted_value(weight, base_value, max_term, scaled_sum, &
                    positive_infinity_seen, negative_infinity_seen, nan_seen)
            end do
        end if

        if (nan_seen .or. (positive_infinity_seen .and. negative_infinity_seen)) then
            y = quiet_nan(x)
        else if (positive_infinity_seen) then
            y = positive_infinity(x)
        else if (negative_infinity_seen) then
            y = negative_infinity(x)
        else if (max_term == 0.0_dp .or. scaled_sum == 0.0_dp) then
            y = 0.0_dp
        else if (abs(scaled_sum) > 1.0_dp) then
            if (max_term > huge(y) / abs(scaled_sum)) then
                if (scaled_sum > 0.0_dp) then
                    y = positive_infinity(x)
                else
                    y = negative_infinity(x)
                end if
            else
                y = max_term * scaled_sum
            end if
        else
            y = max_term * scaled_sum
        end if
    end function bessel_derivative_sum

    pure function bessel_family_value(family, order, x) result(y)
        integer, intent(in) :: family !! internal Bessel-family selector
        real(dp), intent(in) :: order !! shifted real Bessel order
        real(dp), intent(in) :: x !! real Bessel argument
        real(dp) :: y

        select case (family)
        case (family_j)
            y = jv(order, x)
        case (family_y)
            y = yv(order, x)
        case (family_i)
            y = iv(order, x)
        case (family_k)
            y = kv(order, x)
        case default
            y = quiet_nan(x)
        end select
    end function bessel_family_value

    pure subroutine accumulate_weighted_value(weight, value, max_term, scaled_sum, &
            positive_infinity_seen, negative_infinity_seen, nan_seen)
        real(dp), intent(in) :: weight !! normalized signed binomial coefficient
        real(dp), intent(in) :: value !! shifted-order Bessel value
        real(dp), intent(inout) :: max_term !! largest absolute finite term accumulated
        real(dp), intent(inout) :: scaled_sum !! finite sum divided by max_term
        logical, intent(inout) :: positive_infinity_seen !! whether a +infinity term occurred
        logical, intent(inout) :: negative_infinity_seen !! whether a -infinity term occurred
        logical, intent(inout) :: nan_seen !! whether a NaN term occurred

        real(dp) :: abs_term
        real(dp) :: term

        if (ieee_is_nan(value)) then
            nan_seen = .true.
            return
        end if
        if (.not. ieee_is_finite(value)) then
            if ((value > 0.0_dp) .eqv. (weight > 0.0_dp)) then
                positive_infinity_seen = .true.
            else
                negative_infinity_seen = .true.
            end if
            return
        end if

        term = weight * value
        if (term == 0.0_dp) return
        abs_term = abs(term)
        if (max_term == 0.0_dp) then
            max_term = abs_term
            scaled_sum = sign(1.0_dp, term)
        else if (abs_term > max_term) then
            scaled_sum = scaled_sum * (max_term / abs_term) + sign(1.0_dp, term)
            max_term = abs_term
        else
            scaled_sum = scaled_sum + term / max_term
        end if
    end subroutine accumulate_weighted_value

end module scifort_bessel_derivatives
