! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Parabolic-cylinder D_v(x), V_v(x), W(a,x), derivatives, and SciPy-style sequences.
!
! The implementation is independent and follows NIST DLMF Chapter 12.
! Positive D_v arguments use the DLMF 12.7.14 relation to Kummer U. The
! independent V_v solution uses exact reciprocal-gamma origin data, ODE
! continuation, and the DLMF large-positive-argument expansion; negative
! arguments use the DLMF connection formula. Integer D_n values use the
! Hermite recurrence, which also supplies the complex pbdn_seq sequence.
! W(a,x) uses the DLMF origin values and Taylor continuation of its ODE.

module scifort_parabolic_cylinder
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_gamma_family, only : gammasgn, loggamma, rgamma
    use scifort_hypergeometric, only : hyperu
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: parabolic_cylinder_status_invalid_argument
    public :: parabolic_cylinder_status_success
    public :: pbdn_seq
    public :: pbdv
    public :: pbdv_seq
    public :: pbwa
    public :: pbvv
    public :: pbvv_seq

    real(dp), parameter :: log_two = 0.693147180559945309417232121458176568_dp
    real(dp), parameter :: log_sqrt_pi = 0.572364942924700087071713675676529356_dp
    integer, parameter :: parabolic_cylinder_status_success = 0
    integer, parameter :: parabolic_cylinder_status_invalid_argument = 1
    integer, parameter :: parabolic_cylinder_sequence_limit = 100000
    real(dp), parameter :: pbdv_step_size = 0.25_dp
    integer, parameter :: pbdv_max_steps = 10000
    integer, parameter :: pbdv_max_taylor_terms = 120
    integer, parameter :: pbdv_integer_recurrence_limit = 10000
    real(dp), parameter :: pbwa_legacy_domain_limit = 5.0_dp

contains

    pure elemental subroutine pbwa(a, x, w, wprime)
        real(dp), intent(in) :: a !! real parameter, with abs(a) <= 5 for SciPy compatibility
        real(dp), intent(in) :: x !! real argument, with abs(x) <= 5 for SciPy compatibility
        real(dp), intent(out) :: w !! parabolic-cylinder W(a,x)
        real(dp), intent(out) :: wprime !! derivative of W(a,x) with respect to x

        real(dp) :: h
        real(dp) :: remaining
        real(dp) :: x0
        integer :: step

        if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
            w = quiet_nan(x)
            wprime = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(x)) then
            w = quiet_nan(x)
            wprime = quiet_nan(x)
            return
        end if
        if (abs(a) > pbwa_legacy_domain_limit .or. abs(x) > pbwa_legacy_domain_limit) then
            w = quiet_nan(x)
            wprime = quiet_nan(x)
            return
        end if

        call pbwa_origin(a, w, wprime)
        if (x == 0.0_dp) return

        x0 = 0.0_dp
        do step = 1, pbdv_max_steps
            remaining = x - x0
            if (abs(remaining) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))) return
            h = sign(min(abs(remaining), pbdv_step_size), remaining)
            call pbwa_taylor_step(a, x0, h, w, wprime)
            x0 = x0 + h
            if (.not. ieee_is_finite(w) .or. .not. ieee_is_finite(wprime)) return
        end do

        w = quiet_nan(x)
        wprime = quiet_nan(x)
    end subroutine pbwa

    pure subroutine pbwa_origin(a, w, wprime)
        real(dp), intent(in) :: a !! finite real W-function parameter
        real(dp), intent(out) :: w !! W(a,0)
        real(dp), intent(out) :: wprime !! derivative W'(a,0)

        complex(dp) :: lg_denominator
        complex(dp) :: lg_numerator
        real(dp) :: log_ratio

        lg_numerator = loggamma(cmplx(0.25_dp, 0.5_dp * a, kind=dp))
        lg_denominator = loggamma(cmplx(0.75_dp, 0.5_dp * a, kind=dp))
        log_ratio = real(lg_numerator - lg_denominator, dp)

        w = exp(-0.75_dp * log_two + 0.5_dp * log_ratio)
        wprime = -exp(-0.25_dp * log_two - 0.5_dp * log_ratio)
    end subroutine pbwa_origin

    pure subroutine pbwa_taylor_step(a, x0, h, w, wprime)
        real(dp), intent(in) :: a !! finite real W-function parameter
        real(dp), intent(in) :: x0 !! current expansion center
        real(dp), intent(in) :: h !! nonzero real continuation step
        real(dp), intent(inout) :: w !! input/output W-function value
        real(dp), intent(inout) :: wprime !! input/output derivative value

        real(dp) :: coefficients(0:pbdv_max_taylor_terms - 1)
        real(dp) :: c_next
        real(dp) :: derivative_new
        real(dp) :: derivative_term
        real(dp) :: h_power
        real(dp) :: q0
        real(dp) :: q1
        real(dp) :: value_new
        real(dp) :: value_term
        integer :: n

        coefficients = 0.0_dp
        coefficients(0) = w
        coefficients(1) = wprime
        value_new = w + wprime * h
        derivative_new = wprime
        h_power = h
        q0 = 0.25_dp * x0 * x0 - a
        q1 = 0.5_dp * x0

        coefficients(2) = -0.5_dp * q0 * coefficients(0)
        h_power = h_power * h
        value_new = value_new + coefficients(2) * h_power
        derivative_new = derivative_new + 2.0_dp * coefficients(2) * (h_power / h)

        coefficients(3) = (-q0 * coefficients(1) - q1 * coefficients(0)) / 6.0_dp
        h_power = h_power * h
        value_new = value_new + coefficients(3) * h_power
        derivative_new = derivative_new + 3.0_dp * coefficients(3) * (h_power / h)

        do n = 2, pbdv_max_taylor_terms - 3
            c_next = -q0 * coefficients(n) - q1 * coefficients(n - 1) - &
                0.25_dp * coefficients(n - 2)
            c_next = c_next / real((n + 2) * (n + 1), dp)
            coefficients(n + 2) = c_next

            h_power = h_power * h
            value_term = c_next * h_power
            derivative_term = real(n + 2, dp) * c_next * (h_power / h)
            value_new = value_new + value_term
            derivative_new = derivative_new + derivative_term

            if (n >= 10 .and. abs(value_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value_new)) .and. &
                    abs(derivative_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative_new))) exit
        end do

        w = value_new
        wprime = derivative_new
    end subroutine pbwa_taylor_step

    pure elemental subroutine pbdv(v, x, d, dprime)
        real(dp), intent(in) :: v !! real order
        real(dp), intent(in) :: x !! real argument
        real(dp), intent(out) :: d !! parabolic-cylinder D_v(x)
        real(dp), intent(out) :: dprime !! derivative of D_v(x) with respect to x

        real(dp) :: sign_value

        if (ieee_is_nan(v) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(v)) then
            d = quiet_nan(x)
            dprime = quiet_nan(x)
            return
        end if

        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                d = 0.0_dp
                dprime = -0.0_dp
            else if (is_nonnegative_integer(v)) then
                if (modulo(int(v), 2) == 0) then
                    d = 0.0_dp
                    dprime = 0.0_dp
                else
                    d = -0.0_dp
                    dprime = -0.0_dp
                end if
            else
                sign_value = gammasgn(-v)
                d = sign(positive_infinity(0.0_dp), sign_value)
                dprime = sign(positive_infinity(0.0_dp), -sign_value)
            end if
            return
        end if

        if (x == 0.0_dp) then
            call pbdv_origin(v, d, dprime)
        else if (x > 0.0_dp) then
            d = pbdv_positive_value(v, x)
            dprime = 0.5_dp * x * d - pbdv_positive_value(v + 1.0_dp, x)
        else if (is_small_nonnegative_integer(v)) then
            call pbdv_integer_order(int(v), x, d, dprime)
        else
            call pbdv_negative_continuation(v, x, d, dprime)
        end if
    end subroutine pbdv

    pure elemental subroutine pbvv(v, x, vv, vprime)
        real(dp), intent(in) :: v !! real order in SciPy's V_v(x) notation
        real(dp), intent(in) :: x !! real argument
        real(dp), intent(out) :: vv !! parabolic-cylinder V_v(x)
        real(dp), intent(out) :: vprime !! derivative of V_v(x) with respect to x

        real(dp) :: coefficient
        real(dp) :: cosine
        real(dp) :: d
        real(dp) :: dprime
        real(dp) :: positive_value
        real(dp) :: positive_derivative
        real(dp) :: sine
        real(dp) :: t

        if (ieee_is_nan(v) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(v)) then
            vv = quiet_nan(x)
            vprime = quiet_nan(x)
            return
        end if

        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                vv = positive_infinity(x)
                vprime = positive_infinity(x)
            else
                cosine = cospi_reduced(v)
                if (cosine /= 0.0_dp) then
                    vv = sign(positive_infinity(x), -cosine)
                    vprime = sign(positive_infinity(x), cosine)
                else
                    coefficient = -sinpi_reduced(v) * rgamma(v + 1.0_dp)
                    vv = sign(0.0_dp, coefficient)
                    vprime = sign(0.0_dp, -coefficient)
                end if
            end if
            return
        end if

        if (x == 0.0_dp) then
            call pbvv_origin(v, vv, vprime)
            return
        end if

        if (x > 0.0_dp) then
            call pbvv_positive(v, x, vv, vprime)
            return
        end if

        t = -x
        call pbvv_positive(v, t, positive_value, positive_derivative)
        call pbdv(v, t, d, dprime)
        sine = sinpi_reduced(v)
        cosine = cospi_reduced(v)
        coefficient = sine * rgamma(v + 1.0_dp)
        vv = -coefficient * d - cosine * positive_value
        vprime = coefficient * dprime + cosine * positive_derivative
    end subroutine pbvv

    pure subroutine pbdv_seq(v, x, dv, dpdx, status)
        real(dp), intent(in) :: v !! highest requested order following SciPy sequence semantics
        real(dp), intent(in) :: x !! real evaluation argument
        real(dp), allocatable, intent(out) :: dv(:) !! D values in increasing sequence order
        real(dp), allocatable, intent(out) :: dpdx(:) !! corresponding derivatives with respect to x
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        real(dp) :: start_order
        real(dp) :: value
        real(dp) :: derivative
        integer :: k
        integer :: n
        integer :: n1

        if (.not. valid_sequence_order(v)) then
            allocate(dv(0), dpdx(0))
            if (present(status)) status = parabolic_cylinder_status_invalid_argument
            return
        end if

        n = int(v)
        n1 = max(n, 1)
        if (n1 > parabolic_cylinder_sequence_limit) then
            allocate(dv(0), dpdx(0))
            if (present(status)) status = parabolic_cylinder_status_invalid_argument
            return
        end if
        start_order = modulo(real(n1, dp) + v - real(n, dp), 1.0_dp)
        allocate(dv(n1 + 1), dpdx(n1 + 1))
        do k = 0, n1
            call pbdv(start_order + real(k, dp), x, value, derivative)
            dv(k + 1) = value
            dpdx(k + 1) = derivative
        end do
        if (present(status)) status = parabolic_cylinder_status_success
    end subroutine pbdv_seq

    pure subroutine pbvv_seq(v, x, vv, vp, status)
        real(dp), intent(in) :: v !! highest requested order following SciPy sequence semantics
        real(dp), intent(in) :: x !! real evaluation argument
        real(dp), allocatable, intent(out) :: vv(:) !! V values in increasing sequence order
        real(dp), allocatable, intent(out) :: vp(:) !! corresponding derivatives with respect to x
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        real(dp) :: start_order
        real(dp) :: value
        real(dp) :: derivative
        integer :: k
        integer :: n
        integer :: n1

        if (.not. valid_sequence_order(v)) then
            allocate(vv(0), vp(0))
            if (present(status)) status = parabolic_cylinder_status_invalid_argument
            return
        end if

        n = int(v)
        n1 = max(n, 1)
        if (n1 > parabolic_cylinder_sequence_limit) then
            allocate(vv(0), vp(0))
            if (present(status)) status = parabolic_cylinder_status_invalid_argument
            return
        end if
        start_order = modulo(real(n1, dp) + v - real(n, dp), 1.0_dp)
        allocate(vv(n1 + 1), vp(n1 + 1))
        do k = 0, n1
            call pbvv(start_order + real(k, dp), x, value, derivative)
            vv(k + 1) = value
            vp(k + 1) = derivative
        end do
        if (present(status)) status = parabolic_cylinder_status_success
    end subroutine pbvv_seq

    pure subroutine pbdn_seq(n, z, dn, dpdz, status)
        integer, intent(in) :: n !! highest nonnegative integer order; zero still returns D_0 and D_1
        complex(dp), intent(in) :: z !! complex evaluation argument
        complex(dp), allocatable, intent(out) :: dn(:) !! D_k(z), starting at k=0
        complex(dp), allocatable, intent(out) :: dpdz(:) !! corresponding derivatives with respect to z
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        complex(dp), allocatable :: values(:)
        integer :: k
        integer :: n1

        if (n < 0 .or. n > parabolic_cylinder_sequence_limit) then
            allocate(dn(0), dpdz(0))
            if (present(status)) status = parabolic_cylinder_status_invalid_argument
            return
        end if

        n1 = max(n, 1)
        allocate(values(0:n1 + 1), dn(n1 + 1), dpdz(n1 + 1))
        values(0) = exp(-0.25_dp * z * z)
        values(1) = z * values(0)
        do k = 1, n1
            values(k + 1) = z * values(k) - real(k, dp) * values(k - 1)
        end do
        do k = 0, n1
            dn(k + 1) = values(k)
            dpdz(k + 1) = 0.5_dp * z * values(k) - values(k + 1)
        end do
        if (present(status)) status = parabolic_cylinder_status_success
    end subroutine pbdn_seq

    pure subroutine pbvv_positive(v, x, vv, vprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp), intent(out) :: vv !! V_v(x)
        real(dp), intent(out) :: vprime !! derivative V_v'(x)

        real(dp) :: h
        real(dp) :: remaining
        real(dp) :: threshold
        real(dp) :: x0
        integer :: step

        threshold = max(12.0_dp, 2.0_dp * sqrt(abs(v) + 1.0_dp) + 8.0_dp)
        if (x >= threshold) then
            call pbvv_positive_asymptotic(v, x, vv, vprime)
            return
        end if

        call pbvv_origin(v, vv, vprime)
        x0 = 0.0_dp
        do step = 1, pbdv_max_steps
            remaining = x - x0
            if (abs(remaining) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, x)) exit
            h = min(remaining, pbdv_step_size)
            call parabolic_taylor_step(v, x0, h, vv, vprime)
            x0 = x0 + h
            if (.not. ieee_is_finite(vv) .or. .not. ieee_is_finite(vprime)) exit
        end do
    end subroutine pbvv_positive

    pure subroutine pbvv_origin(v, vv, vprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(out) :: vv !! V_v(0)
        real(dp), intent(out) :: vprime !! derivative V_v'(0)

        real(dp) :: r1
        real(dp) :: r2
        real(dp) :: log_abs

        r1 = rgamma(1.0_dp + 0.5_dp * v)
        r2 = rgamma(-0.5_dp * v)
        if (r1 == 0.0_dp .or. r2 == 0.0_dp) then
            vv = 0.0_dp
        else
            log_abs = log(scifort_pi) - 0.5_dp * v * log_two + &
                2.0_dp * log(abs(r1)) + log(abs(r2))
            vv = signed_exp(log_abs, r2)
        end if

        r1 = rgamma(0.5_dp * (v + 1.0_dp))
        r2 = rgamma(0.5_dp * (1.0_dp - v))
        if (r1 == 0.0_dp .or. r2 == 0.0_dp) then
            vprime = 0.0_dp
        else
            log_abs = log(scifort_pi) + 0.5_dp * (1.0_dp - v) * log_two + &
                2.0_dp * log(abs(r1)) + log(abs(r2))
            vprime = signed_exp(log_abs, r2)
        end if
    end subroutine pbvv_origin

    pure subroutine pbvv_positive_asymptotic(v, x, vv, vprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! sufficiently large positive finite argument
        real(dp), intent(out) :: vv !! asymptotic V_v(x)
        real(dp), intent(out) :: vprime !! asymptotic derivative V_v'(x)

        integer, parameter :: max_terms = 320
        real(dp) :: bracket
        real(dp) :: log_prefactor
        real(dp) :: next_term
        real(dp) :: ratio
        real(dp) :: sum_value
        real(dp) :: term
        real(dp) :: weighted_sum
        integer :: s

        term = 1.0_dp
        sum_value = 1.0_dp
        weighted_sum = 0.0_dp
        do s = 0, max_terms - 2
            ratio = ((v + 1.0_dp + 2.0_dp * real(s, dp)) / x) * &
                ((v + 2.0_dp + 2.0_dp * real(s, dp)) / x) / &
                (2.0_dp * real(s + 1, dp))
            next_term = term * ratio
            if (s >= 2 .and. abs(next_term) > abs(term)) exit
            term = next_term
            sum_value = sum_value + term
            weighted_sum = weighted_sum + real(s + 1, dp) * term
            if (abs(term) <= 8.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(sum_value))) exit
        end do

        log_prefactor = 0.5_dp * log(2.0_dp / scifort_pi) + 0.25_dp * x * x - &
            (v + 1.0_dp) * log(x)
        if (sum_value == 0.0_dp) then
            vv = 0.0_dp
        else
            vv = signed_exp(log_prefactor + log(abs(sum_value)), sum_value)
        end if
        bracket = (0.5_dp * x - (v + 1.0_dp) / x) * sum_value - &
            2.0_dp * weighted_sum / x
        if (bracket == 0.0_dp) then
            vprime = 0.0_dp
        else
            vprime = signed_exp(log_prefactor + log(abs(bracket)), bracket)
        end if
    end subroutine pbvv_positive_asymptotic

    pure function pbdv_positive_value(v, x) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite argument
        real(dp) :: value

        real(dp) :: log_scale
        real(dp) :: scale
        real(dp) :: t
        real(dp) :: u_value

        t = 0.5_dp * x * x
        u_value = hyperu(-0.5_dp * v, 0.5_dp, t)
        if (u_value == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        log_scale = 0.5_dp * v * log_two - 0.25_dp * x * x
        if (.not. ieee_is_finite(u_value)) then
            if (log_scale < log(tiny(1.0_dp))) then
                value = quiet_nan(x)
            else
                value = sign(positive_infinity(x), u_value)
            end if
        else if (log_scale > log(huge(1.0_dp))) then
            value = sign(positive_infinity(x), u_value)
        else if (log_scale < log(tiny(1.0_dp))) then
            if (log(abs(u_value)) + log_scale < log(tiny(1.0_dp))) then
                value = sign(0.0_dp, u_value)
            else
                value = sign(1.0_dp, u_value) * exp(log(abs(u_value)) + log_scale)
            end if
        else
            scale = exp(log_scale)
            value = scale * u_value
        end if
    end function pbdv_positive_value

    pure subroutine pbdv_origin(v, d, dprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(out) :: d !! D_v(0)
        real(dp), intent(out) :: dprime !! derivative D_v'(0)

        real(dp) :: denominator
        real(dp) :: log_abs

        denominator = 0.5_dp * (1.0_dp - v)
        if (is_gamma_pole(denominator)) then
            d = 0.0_dp
        else
            log_abs = 0.5_dp * v * log_two + log_sqrt_pi - gammaln(denominator)
            d = signed_exp(log_abs, gammasgn(denominator))
        end if

        denominator = -0.5_dp * v
        if (is_gamma_pole(denominator)) then
            dprime = 0.0_dp
        else
            log_abs = 0.5_dp * (v + 1.0_dp) * log_two + log_sqrt_pi - &
                gammaln(denominator)
            dprime = signed_exp(log_abs, -gammasgn(denominator))
        end if
    end subroutine pbdv_origin

    pure subroutine pbdv_negative_continuation(v, x, d, dprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! negative finite target argument
        real(dp), intent(out) :: d !! D_v(x)
        real(dp), intent(out) :: dprime !! derivative D_v'(x)

        real(dp) :: h
        real(dp) :: remaining
        real(dp) :: x0
        integer :: step

        call pbdv_origin(v, d, dprime)
        x0 = 0.0_dp

        do step = 1, pbdv_max_steps
            remaining = x - x0
            if (abs(remaining) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))) exit
            h = -min(abs(remaining), pbdv_step_size)
            call parabolic_taylor_step(v, x0, h, d, dprime)
            x0 = x0 + h
            if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(dprime)) exit
        end do
    end subroutine pbdv_negative_continuation

    pure subroutine parabolic_taylor_step(v, x0, h, d, dprime)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x0 !! current expansion center
        real(dp), intent(in) :: h !! nonzero real continuation step
        real(dp), intent(inout) :: d !! input/output D_v value
        real(dp), intent(inout) :: dprime !! input/output derivative value

        real(dp) :: coefficients(0:pbdv_max_taylor_terms - 1)
        real(dp) :: c_next
        real(dp) :: derivative_new
        real(dp) :: derivative_term
        real(dp) :: h_power
        real(dp) :: q0
        real(dp) :: q1
        real(dp) :: value_new
        real(dp) :: value_term
        integer :: n

        coefficients = 0.0_dp
        coefficients(0) = d
        coefficients(1) = dprime
        value_new = d + dprime * h
        derivative_new = dprime
        h_power = h
        q0 = v + 0.5_dp - 0.25_dp * x0 * x0
        q1 = -0.5_dp * x0

        do n = 0, pbdv_max_taylor_terms - 3
            c_next = -q0 * coefficients(n)
            if (n >= 1) c_next = c_next - q1 * coefficients(n - 1)
            if (n >= 2) c_next = c_next + 0.25_dp * coefficients(n - 2)
            c_next = c_next / real((n + 2) * (n + 1), dp)
            coefficients(n + 2) = c_next

            h_power = h_power * h
            value_term = c_next * h_power
            derivative_term = real(n + 2, dp) * c_next * (h_power / h)
            value_new = value_new + value_term
            derivative_new = derivative_new + derivative_term

            if (n >= 10 .and. abs(value_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value_new)) .and. &
                    abs(derivative_term) <= &
                    16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(derivative_new))) exit
        end do

        d = value_new
        dprime = derivative_new
    end subroutine parabolic_taylor_step

    pure subroutine pbdv_integer_order(n, x, d, dprime)
        integer, intent(in) :: n !! nonnegative integer order
        real(dp), intent(in) :: x !! finite real argument
        real(dp), intent(out) :: d !! D_n(x)
        real(dp), intent(out) :: dprime !! derivative D_n'(x)

        real(dp) :: current
        real(dp) :: d_next
        real(dp) :: previous
        real(dp) :: start
        integer :: k

        start = exp(-0.25_dp * x * x)
        if (n == 0) then
            d = start
            dprime = -0.5_dp * x * start
            return
        end if

        previous = start
        current = x * start
        if (n > 1) then
            do k = 1, n - 1
                d_next = x * current - real(k, dp) * previous
                previous = current
                current = d_next
            end do
        end if
        d = current
        d_next = x * current - real(n, dp) * previous
        dprime = 0.5_dp * x * d - d_next
    end subroutine pbdv_integer_order

    pure elemental function is_gamma_pole(x) result(is_pole)
        real(dp), intent(in) :: x !! Gamma argument tested for a nonpositive-integer pole
        logical :: is_pole

        is_pole = ieee_is_finite(x) .and. x <= 0.0_dp .and. x == anint(x)
    end function is_gamma_pole

    pure elemental function is_nonnegative_integer(x) result(is_integer)
        real(dp), intent(in) :: x !! value tested for a nonnegative integer
        logical :: is_integer

        is_integer = ieee_is_finite(x) .and. x >= 0.0_dp .and. x == anint(x) .and. &
            x <= real(huge(1), dp)
    end function is_nonnegative_integer

    pure elemental function is_small_nonnegative_integer(x) result(is_integer)
        real(dp), intent(in) :: x !! value tested for recurrence-sized nonnegative integer
        logical :: is_integer

        is_integer = is_nonnegative_integer(x) .and. &
            x <= real(pbdv_integer_recurrence_limit, dp)
    end function is_small_nonnegative_integer

    pure elemental function valid_sequence_order(v) result(valid)
        real(dp), intent(in) :: v !! proposed real sequence order
        logical :: valid

        valid = ieee_is_finite(v) .and. abs(v) <= real(parabolic_cylinder_sequence_limit, dp)
    end function valid_sequence_order

    pure elemental function sinpi_reduced(x) result(value)
        real(dp), intent(in) :: x !! finite real argument of sin(pi*x)
        real(dp) :: value

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        if (r == 0.0_dp .or. abs(r) == 1.0_dp) then
            value = sign(0.0_dp, r)
        else if (abs(r) == 0.5_dp) then
            value = sign(1.0_dp, r)
        else
            value = sin(scifort_pi * r)
        end if
    end function sinpi_reduced

    pure elemental function cospi_reduced(x) result(value)
        real(dp), intent(in) :: x !! finite real argument of cos(pi*x)
        real(dp) :: value

        real(dp) :: r

        r = modulo(x, 2.0_dp)
        if (r > 1.0_dp) r = r - 2.0_dp
        if (r == 0.0_dp) then
            value = 1.0_dp
        else if (abs(r) == 1.0_dp) then
            value = -1.0_dp
        else if (abs(r) == 0.5_dp) then
            value = 0.0_dp
        else
            value = cos(scifort_pi * r)
        end if
    end function cospi_reduced

    pure function signed_exp(log_abs, sign_value) result(value)
        real(dp), intent(in) :: log_abs !! logarithm of the absolute result
        real(dp), intent(in) :: sign_value !! sign carrier
        real(dp) :: value

        if (log_abs > log(huge(1.0_dp))) then
            value = sign(positive_infinity(log_abs), sign_value)
        else if (log_abs < log(tiny(1.0_dp))) then
            value = sign(0.0_dp, sign_value)
        else
            value = sign(1.0_dp, sign_value) * exp(log_abs)
        end if
    end function signed_exp

end module scifort_parabolic_cylinder
