! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real/complex Hankel functions and derivatives for real order.
! Real-axis values reuse SciFort's cylindrical Bessel J/Y kernels.  Complex
! values use DLMF 10.27.8 to connect the recessive half-plane solution to a
! scaled K function, angular ODE continuation of the scaled Bessel equations,
! and DLMF 10.11.3-10.11.5 across the negative-real branch cut.

module scifort_hankel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_derivatives, only : jvp, yvp
    use scifort_bessel_jy, only : jv, yv
    use scifort_bessel_order, only : kve
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)
    integer, parameter :: arc_equation_hankel = 1
    integer, parameter :: arc_equation_kve = 2
    integer, parameter :: arc_max_steps = 100000
    integer, parameter :: derivative_weight_limit = 1000
    real(dp), parameter :: arc_ode_tolerance = 1.0e-15_dp

    interface hankel1
        module procedure hankel1_real
        module procedure hankel1_complex
    end interface hankel1

    interface hankel1e
        module procedure hankel1e_real
        module procedure hankel1e_complex
    end interface hankel1e

    interface hankel2
        module procedure hankel2_real
        module procedure hankel2_complex
    end interface hankel2

    interface hankel2e
        module procedure hankel2e_real
        module procedure hankel2e_complex
    end interface hankel2e

    interface h1vp
        module procedure h1vp_real
        module procedure h1vp_complex
    end interface h1vp

    interface h2vp
        module procedure h2vp_real
        module procedure h2vp_complex
    end interface h2vp

    public :: h1vp
    public :: h2vp
    public :: hankel1
    public :: hankel1e
    public :: hankel2
    public :: hankel2e

contains

    pure elemental function hankel1_real(v, x) result(value)
        real(dp), intent(in) :: v !! real order of the Hankel function
        real(dp), intent(in) :: x !! finite real argument; zero is singular
        complex(dp) :: value

        complex(dp) :: phase_minus

        if (invalid_real_inputs(v, x)) then
            value = complex_nan()
        else if (is_zero(x)) then
            value = complex_nan()
        else if (x > 0.0_dp) then
            value = positive_hankel(v, x, 1)
        else
            phase_minus = phase_pi(v, -1)
            value = -phase_minus * positive_hankel(v, -x, 2)
        end if
    end function hankel1_real

    pure elemental function hankel2_real(v, x) result(value)
        real(dp), intent(in) :: v !! real order of the Hankel function
        real(dp), intent(in) :: x !! finite real argument; zero is singular
        complex(dp) :: value

        complex(dp) :: h1_positive
        complex(dp) :: h2_positive
        complex(dp) :: phase_plus
        real(dp) :: cosine_pi_v

        if (invalid_real_inputs(v, x)) then
            value = complex_nan()
        else if (is_zero(x)) then
            if (is_zero(v)) then
                value = cmplx(quiet_nan(x), positive_infinity(x), kind=dp)
            else
                value = complex_nan()
            end if
        else if (x > 0.0_dp) then
            value = positive_hankel(v, x, 2)
        else
            h1_positive = positive_hankel(v, -x, 1)
            h2_positive = positive_hankel(v, -x, 2)
            phase_plus = phase_pi(v, 1)
            cosine_pi_v = cospi_value(v)
            value = phase_plus * h1_positive + 2.0_dp * cosine_pi_v * h2_positive
        end if
    end function hankel2_real

    pure elemental function hankel1e_real(v, x) result(value)
        real(dp), intent(in) :: v !! real order of exp(-i*x)*H1_v(x)
        real(dp), intent(in) :: x !! finite nonzero real argument
        complex(dp) :: value

        if (invalid_real_inputs(v, x) .or. is_zero(x)) then
            value = complex_nan()
        else
            value = hankel1_real(v, x) * unit_phase(-x)
        end if
    end function hankel1e_real

    pure elemental function hankel2e_real(v, x) result(value)
        real(dp), intent(in) :: v !! real order of exp(+i*x)*H2_v(x)
        real(dp), intent(in) :: x !! finite nonzero real argument
        complex(dp) :: value

        if (invalid_real_inputs(v, x) .or. is_zero(x)) then
            value = complex_nan()
        else
            value = hankel2_real(v, x) * unit_phase(x)
        end if
    end function hankel2e_real

    pure elemental function h1vp_real(v, x, n) result(value)
        real(dp), intent(in) :: v !! real Hankel-function order before differentiation
        real(dp), intent(in) :: x !! finite real argument; zero is singular
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        complex(dp) :: value

        complex(dp) :: derivative_positive
        complex(dp) :: phase_minus
        integer :: derivative_order
        real(dp) :: parity

        derivative_order = 1
        if (present(n)) derivative_order = n

        if (derivative_order < 0 .or. invalid_real_inputs(v, x)) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = hankel1_real(v, x)
        else if (is_zero(x)) then
            value = complex_nan()
        else if (x > 0.0_dp) then
            value = positive_hankel_derivative(v, x, derivative_order, 1)
        else
            derivative_positive = positive_hankel_derivative(v, -x, derivative_order, 2)
            phase_minus = phase_pi(v, -1)
            parity = 1.0_dp
            if (mod(derivative_order, 2) == 0) parity = -1.0_dp
            value = parity * phase_minus * derivative_positive
        end if
    end function h1vp_real

    pure elemental function h2vp_real(v, x, n) result(value)
        real(dp), intent(in) :: v !! real Hankel-function order before differentiation
        real(dp), intent(in) :: x !! finite real argument; zero is singular
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        complex(dp) :: value

        complex(dp) :: h1_derivative
        complex(dp) :: h2_derivative
        complex(dp) :: phase_plus
        integer :: derivative_order
        real(dp) :: cosine_pi_v
        real(dp) :: parity

        derivative_order = 1
        if (present(n)) derivative_order = n

        if (derivative_order < 0 .or. invalid_real_inputs(v, x)) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = hankel2_real(v, x)
        else if (is_zero(x)) then
            value = complex_nan()
        else if (x > 0.0_dp) then
            value = positive_hankel_derivative(v, x, derivative_order, 2)
        else
            h1_derivative = positive_hankel_derivative(v, -x, derivative_order, 1)
            h2_derivative = positive_hankel_derivative(v, -x, derivative_order, 2)
            phase_plus = phase_pi(v, 1)
            cosine_pi_v = cospi_value(v)
            parity = 1.0_dp
            if (mod(derivative_order, 2) /= 0) parity = -1.0_dp
            value = parity * (phase_plus * h1_derivative + &
                2.0_dp * cosine_pi_v * h2_derivative)
        end if
    end function h2vp_real

    pure elemental function hankel1_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of the Hankel function
        complex(dp), intent(in) :: z !! finite complex argument; zero is singular
        complex(dp) :: value

        if (invalid_complex_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = hankel1_real(v, real(z, dp))
        else
            value = hankel1e_complex(v, z) * exp(imaginary_unit * z)
        end if
    end function hankel1_complex

    pure elemental function hankel2_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of the Hankel function
        complex(dp), intent(in) :: z !! finite complex argument; zero is singular
        complex(dp) :: value

        if (invalid_complex_inputs(v, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            if (is_zero(v)) then
                value = cmplx(quiet_nan(0.0_dp), positive_infinity(0.0_dp), kind=dp)
            else
                value = complex_nan()
            end if
        else if (aimag(z) == 0.0_dp) then
            value = hankel2_real(v, real(z, dp))
        else
            value = hankel2e_complex(v, z) * exp(-imaginary_unit * z)
        end if
    end function hankel2_complex

    pure elemental function hankel1e_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(-i*z)*H1_v(z)
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        if (invalid_complex_inputs(v, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = hankel1e_real(v, real(z, dp))
        else
            value = complex_hankel_scaled(v, z, 1)
        end if
    end function hankel1e_complex

    pure elemental function hankel2e_complex(v, z) result(value)
        real(dp), intent(in) :: v !! real order of exp(+i*z)*H2_v(z)
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        if (invalid_complex_inputs(v, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = hankel2e_real(v, real(z, dp))
        else
            value = complex_hankel_scaled(v, z, 2)
        end if
    end function hankel2e_complex

    pure elemental function h1vp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real Hankel-function order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument; zero is singular
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n

        if (derivative_order < 0 .or. invalid_complex_inputs(v, z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = h1vp_real(v, real(z, dp), derivative_order)
        else if (complex_is_zero(z)) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = hankel1_complex(v, z)
        else
            value = complex_hankel_derivative(v, z, derivative_order, 1)
        end if
    end function h1vp_complex

    pure elemental function h2vp_complex(v, z, n) result(value)
        real(dp), intent(in) :: v !! real Hankel-function order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument; zero is singular
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        integer :: derivative_order

        derivative_order = 1
        if (present(n)) derivative_order = n

        if (derivative_order < 0 .or. invalid_complex_inputs(v, z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = h2vp_real(v, real(z, dp), derivative_order)
        else if (complex_is_zero(z)) then
            value = complex_nan()
        else if (derivative_order == 0) then
            value = hankel2_complex(v, z)
        else
            value = complex_hankel_derivative(v, z, derivative_order, 2)
        end if
    end function h2vp_complex

    pure function complex_hankel_scaled(v, z, kind_index) result(value)
        real(dp), intent(in) :: v !! finite real Hankel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        complex(dp) :: h1_scaled
        complex(dp) :: h2_scaled
        complex(dp) :: phase_minus
        complex(dp) :: phase_plus
        complex(dp) :: w
        real(dp) :: cosine_pi_v

        if (real(z, dp) >= 0.0_dp) then
            value = right_half_hankel_scaled(v, z, kind_index)
            return
        end if

        w = -z
        h1_scaled = right_half_hankel_scaled(v, w, 1)
        h2_scaled = right_half_hankel_scaled(v, w, 2)
        if (.not. complex_is_finite(h1_scaled) .or. .not. complex_is_finite(h2_scaled)) then
            value = complex_nan()
            return
        end if

        phase_minus = phase_pi(v, -1)
        phase_plus = phase_pi(v, 1)
        cosine_pi_v = cospi_value(v)

        if (aimag(z) > 0.0_dp) then
            if (kind_index == 1) then
                value = -phase_minus * h2_scaled
            else
                value = phase_plus * h1_scaled + 2.0_dp * cosine_pi_v * &
                    exp(-2.0_dp * imaginary_unit * w) * h2_scaled
            end if
        else
            if (kind_index == 1) then
                value = 2.0_dp * cosine_pi_v * exp(2.0_dp * imaginary_unit * w) * &
                    h1_scaled + phase_minus * h2_scaled
            else
                value = -phase_plus * h1_scaled
            end if
        end if
    end function complex_hankel_scaled

    pure function right_half_hankel_scaled(v, z, kind_index) result(value)
        real(dp), intent(in) :: v !! finite real Hankel order
        complex(dp), intent(in) :: z !! nonzero argument with nonnegative real part
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        complex(dp) :: k_scaled
        complex(dp) :: w
        real(dp), parameter :: two_over_pi = 2.0_dp / acos(-1.0_dp)

        if (aimag(z) == 0.0_dp) then
            if (kind_index == 1) then
                value = hankel1e_real(v, real(z, dp))
            else
                value = hankel2e_real(v, real(z, dp))
            end if
        else if (kind_index == 1 .and. aimag(z) > 0.0_dp) then
            w = -imaginary_unit * z
            k_scaled = complex_kve_right(v, w)
            value = -two_over_pi * imaginary_unit * phase_pi(0.5_dp * v, -1) * &
                k_scaled
        else if (kind_index == 2 .and. aimag(z) < 0.0_dp) then
            w = imaginary_unit * z
            k_scaled = complex_kve_right(v, w)
            value = two_over_pi * imaginary_unit * phase_pi(0.5_dp * v, 1) * &
                k_scaled
        else
            value = continued_hankel_scaled(v, z, kind_index)
        end if
    end function right_half_hankel_scaled

    pure function complex_kve_right(v, z) result(value)
        real(dp), intent(in) :: v !! finite real modified-Bessel order
        complex(dp), intent(in) :: z !! nonzero argument with nonnegative real part
        complex(dp) :: value

        complex(dp) :: derivative
        logical :: success
        real(dp) :: radius
        real(dp) :: theta

        radius = abs(z)
        if (radius == 0.0_dp) then
            value = complex_nan()
            return
        end if
        if (aimag(z) == 0.0_dp) then
            value = cmplx(kve(v, real(z, dp)), 0.0_dp, kind=dp)
            return
        end if

        theta = atan(aimag(z), real(z, dp))
        value = cmplx(kve(v, radius), 0.0_dp, kind=dp)
        derivative = cmplx(kve(v, radius) - &
            0.5_dp * (kve(v - 1.0_dp, radius) + kve(v + 1.0_dp, radius)), &
            0.0_dp, kind=dp)
        if (.not. complex_is_finite(value) .or. .not. complex_is_finite(derivative)) then
            value = complex_nan()
            return
        end if

        call integrate_scaled_arc(v, radius, theta, arc_equation_kve, 0, &
            value, derivative, success)
        if (.not. success) value = complex_nan()
    end function complex_kve_right

    pure function continued_hankel_scaled(v, z, kind_index) result(value)
        real(dp), intent(in) :: v !! finite real Hankel order
        complex(dp), intent(in) :: z !! nonreal right-half-plane argument
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        complex(dp) :: derivative
        complex(dp) :: lower_scaled
        complex(dp) :: upper_scaled
        logical :: success
        integer :: sign_index
        real(dp) :: radius
        real(dp) :: theta

        radius = abs(z)
        theta = atan(aimag(z), real(z, dp))
        if (kind_index == 1) then
            sign_index = 1
            value = hankel1e_real(v, radius)
            lower_scaled = hankel1e_real(v - 1.0_dp, radius)
            upper_scaled = hankel1e_real(v + 1.0_dp, radius)
        else
            sign_index = -1
            value = hankel2e_real(v, radius)
            lower_scaled = hankel2e_real(v - 1.0_dp, radius)
            upper_scaled = hankel2e_real(v + 1.0_dp, radius)
        end if
        derivative = 0.5_dp * (lower_scaled - upper_scaled) - &
            real(sign_index, dp) * imaginary_unit * value
        if (.not. complex_is_finite(value) .or. .not. complex_is_finite(derivative)) then
            value = complex_nan()
            return
        end if

        call integrate_scaled_arc(v, radius, theta, arc_equation_hankel, &
            sign_index, value, derivative, success)
        if (.not. success) value = complex_nan()
    end function continued_hankel_scaled

    pure subroutine integrate_scaled_arc(v, radius, theta_end, equation, sign_index, &
            value, derivative, success)
        ! Dormand-Prince embedded RK5(4) coefficients from
        ! J. Comput. Appl. Math. 6 (1980), 19-26,
        ! doi:10.1016/0771-050X(80)90013-3.
        real(dp), intent(in) :: v !! finite real Bessel order
        real(dp), intent(in) :: radius !! positive constant radius of continuation arc
        real(dp), intent(in) :: theta_end !! principal target phase in radians
        integer, intent(in) :: equation !! internal scaled-equation selector
        integer, intent(in) :: sign_index !! +1 for H1e, -1 for H2e, unused for Kve
        complex(dp), intent(inout) :: value !! scaled function at start on entry and target on return
        complex(dp), intent(inout) :: derivative !! derivative in z at start and target
        logical, intent(out) :: success !! true when adaptive continuation reaches theta_end

        complex(dp) :: k11, k12, k21, k22, k31, k32, k41, k42
        complex(dp) :: k51, k52, k61, k62, k71, k72
        complex(dp) :: t1, t2
        complex(dp) :: y1
        complex(dp) :: y1_fourth
        complex(dp) :: y1_fifth
        complex(dp) :: y2
        complex(dp) :: y2_fourth
        complex(dp) :: y2_fifth
        integer :: step_count
        real(dp) :: direction
        real(dp) :: error_norm
        real(dp) :: factor
        real(dp) :: h
        real(dp) :: remaining
        real(dp) :: scale1
        real(dp) :: scale2
        real(dp) :: theta

        if (theta_end == 0.0_dp) then
            success = .true.
            return
        end if

        direction = sign(1.0_dp, theta_end)
        theta = 0.0_dp
        y1 = value
        y2 = derivative
        h = direction * min(abs(theta_end), 0.05_dp)
        success = .false.

        do step_count = 1, arc_max_steps
            remaining = theta_end - theta
            if (direction * remaining <= 0.0_dp) exit
            if (abs(h) > abs(remaining)) h = remaining

            call scaled_arc_rhs(v, radius, theta, equation, sign_index, y1, y2, &
                k11, k12)
            t1 = y1 + h * (1.0_dp / 5.0_dp) * k11
            t2 = y2 + h * (1.0_dp / 5.0_dp) * k12
            call scaled_arc_rhs(v, radius, theta + h / 5.0_dp, equation, &
                sign_index, t1, t2, k21, k22)

            t1 = y1 + h * (3.0_dp * k11 / 40.0_dp + 9.0_dp * k21 / 40.0_dp)
            t2 = y2 + h * (3.0_dp * k12 / 40.0_dp + 9.0_dp * k22 / 40.0_dp)
            call scaled_arc_rhs(v, radius, theta + 3.0_dp * h / 10.0_dp, equation, &
                sign_index, t1, t2, k31, k32)

            t1 = y1 + h * (44.0_dp * k11 / 45.0_dp - 56.0_dp * k21 / 15.0_dp + &
                32.0_dp * k31 / 9.0_dp)
            t2 = y2 + h * (44.0_dp * k12 / 45.0_dp - 56.0_dp * k22 / 15.0_dp + &
                32.0_dp * k32 / 9.0_dp)
            call scaled_arc_rhs(v, radius, theta + 4.0_dp * h / 5.0_dp, equation, &
                sign_index, t1, t2, k41, k42)

            t1 = y1 + h * (19372.0_dp * k11 / 6561.0_dp - &
                25360.0_dp * k21 / 2187.0_dp + 64448.0_dp * k31 / 6561.0_dp - &
                212.0_dp * k41 / 729.0_dp)
            t2 = y2 + h * (19372.0_dp * k12 / 6561.0_dp - &
                25360.0_dp * k22 / 2187.0_dp + 64448.0_dp * k32 / 6561.0_dp - &
                212.0_dp * k42 / 729.0_dp)
            call scaled_arc_rhs(v, radius, theta + 8.0_dp * h / 9.0_dp, equation, &
                sign_index, t1, t2, k51, k52)

            t1 = y1 + h * (9017.0_dp * k11 / 3168.0_dp - &
                355.0_dp * k21 / 33.0_dp + 46732.0_dp * k31 / 5247.0_dp + &
                49.0_dp * k41 / 176.0_dp - 5103.0_dp * k51 / 18656.0_dp)
            t2 = y2 + h * (9017.0_dp * k12 / 3168.0_dp - &
                355.0_dp * k22 / 33.0_dp + 46732.0_dp * k32 / 5247.0_dp + &
                49.0_dp * k42 / 176.0_dp - 5103.0_dp * k52 / 18656.0_dp)
            call scaled_arc_rhs(v, radius, theta + h, equation, sign_index, &
                t1, t2, k61, k62)

            y1_fifth = y1 + h * (35.0_dp * k11 / 384.0_dp + &
                500.0_dp * k31 / 1113.0_dp + 125.0_dp * k41 / 192.0_dp - &
                2187.0_dp * k51 / 6784.0_dp + 11.0_dp * k61 / 84.0_dp)
            y2_fifth = y2 + h * (35.0_dp * k12 / 384.0_dp + &
                500.0_dp * k32 / 1113.0_dp + 125.0_dp * k42 / 192.0_dp - &
                2187.0_dp * k52 / 6784.0_dp + 11.0_dp * k62 / 84.0_dp)
            call scaled_arc_rhs(v, radius, theta + h, equation, sign_index, &
                y1_fifth, y2_fifth, k71, k72)

            y1_fourth = y1 + h * (5179.0_dp * k11 / 57600.0_dp + &
                7571.0_dp * k31 / 16695.0_dp + 393.0_dp * k41 / 640.0_dp - &
                92097.0_dp * k51 / 339200.0_dp + 187.0_dp * k61 / 2100.0_dp + &
                k71 / 40.0_dp)
            y2_fourth = y2 + h * (5179.0_dp * k12 / 57600.0_dp + &
                7571.0_dp * k32 / 16695.0_dp + 393.0_dp * k42 / 640.0_dp - &
                92097.0_dp * k52 / 339200.0_dp + 187.0_dp * k62 / 2100.0_dp + &
                k72 / 40.0_dp)

            if (.not. complex_is_finite(y1_fifth) .or. &
                    .not. complex_is_finite(y2_fifth) .or. &
                    .not. complex_is_finite(y1_fourth) .or. &
                    .not. complex_is_finite(y2_fourth)) return

            scale1 = max(1.0_dp, abs(y1), abs(y1_fifth))
            scale2 = max(1.0_dp, abs(y2), abs(y2_fifth))
            error_norm = max(abs(y1_fifth - y1_fourth) / scale1, &
                abs(y2_fifth - y2_fourth) / scale2) / arc_ode_tolerance
            if (.not. ieee_is_finite(error_norm)) return

            if (error_norm <= 1.0_dp) then
                theta = theta + h
                y1 = y1_fifth
                y2 = y2_fifth
                if (theta == theta_end) exit
            end if

            if (error_norm == 0.0_dp) then
                factor = 5.0_dp
            else
                factor = 0.9_dp * error_norm ** (-0.2_dp)
                factor = min(5.0_dp, max(0.2_dp, factor))
            end if
            h = h * factor
            if (abs(h) < spacing(theta)) return
        end do

        if (direction * (theta_end - theta) > 8.0_dp * spacing(theta_end)) return
        value = y1
        derivative = y2
        success = .true.
    end subroutine integrate_scaled_arc

    pure subroutine scaled_arc_rhs(v, radius, theta, equation, sign_index, &
            value, derivative, rhs_value, rhs_derivative)
        real(dp), intent(in) :: v !! finite real Bessel order
        real(dp), intent(in) :: radius !! positive radius of the continuation arc
        real(dp), intent(in) :: theta !! current polar angle in radians
        integer, intent(in) :: equation !! internal scaled-equation selector
        integer, intent(in) :: sign_index !! +1 for H1e, -1 for H2e, unused for Kve
        complex(dp), intent(in) :: value !! scaled function value at current point
        complex(dp), intent(in) :: derivative !! derivative in z at current point
        complex(dp), intent(out) :: rhs_value !! derivative of value with respect to theta
        complex(dp), intent(out) :: rhs_derivative !! derivative of z-derivative with respect to theta

        complex(dp) :: z
        real(dp) :: order_squared

        z = radius * cmplx(cos(theta), sin(theta), kind=dp)
        order_squared = v * v
        rhs_value = imaginary_unit * z * derivative
        select case (equation)
        case (arc_equation_hankel)
            rhs_derivative = (2.0_dp * real(sign_index, dp) * z - imaginary_unit) * &
                derivative + (real(sign_index, dp) + imaginary_unit * order_squared / z) * &
                value
        case (arc_equation_kve)
            rhs_derivative = (2.0_dp * imaginary_unit * z - imaginary_unit) * &
                derivative + imaginary_unit * (1.0_dp + order_squared / z) * value
        case default
            rhs_value = complex_nan()
            rhs_derivative = complex_nan()
        end select
    end subroutine scaled_arc_rhs

    pure function complex_hankel_derivative(v, z, n, kind_index) result(value)
        real(dp), intent(in) :: v !! real Hankel order before differentiation
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        integer, intent(in) :: n !! positive derivative order
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        complex(dp) :: base_value
        complex(dp) :: scaled_sum
        complex(dp) :: term
        integer :: k
        integer :: sign_index
        real(dp) :: abs_term
        real(dp) :: log_weight
        real(dp) :: max_term
        real(dp) :: order
        real(dp) :: weight

        max_term = 0.0_dp
        scaled_sum = (0.0_dp, 0.0_dp)

        if (n <= derivative_weight_limit) then
            weight = scale(1.0_dp, -n)
            do k = 0, n
                if (k > 0) then
                    weight = -weight * real(n - k + 1, dp) / real(k, dp)
                end if
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = complex_hankel_scaled(order, z, kind_index)
                if (.not. complex_is_finite(base_value)) then
                    value = complex_nan()
                    return
                end if
                term = weight * base_value
                abs_term = abs(term)
                if (abs_term == 0.0_dp) cycle
                if (max_term == 0.0_dp) then
                    max_term = abs_term
                    scaled_sum = term / abs_term
                else if (abs_term > max_term) then
                    scaled_sum = scaled_sum * (max_term / abs_term) + term / abs_term
                    max_term = abs_term
                else
                    scaled_sum = scaled_sum + term / max_term
                end if
            end do
        else
            do k = 0, n
                log_weight = log_gamma(real(n, dp) + 1.0_dp) - &
                    log_gamma(real(k, dp) + 1.0_dp) - &
                    log_gamma(real(n - k, dp) + 1.0_dp) - real(n, dp) * log(2.0_dp)
                weight = exp(log_weight)
                if (mod(k, 2) /= 0) weight = -weight
                if (weight == 0.0_dp) cycle
                order = v - real(n, dp) + 2.0_dp * real(k, dp)
                base_value = complex_hankel_scaled(order, z, kind_index)
                if (.not. complex_is_finite(base_value)) then
                    value = complex_nan()
                    return
                end if
                term = weight * base_value
                abs_term = abs(term)
                if (abs_term == 0.0_dp) cycle
                if (max_term == 0.0_dp) then
                    max_term = abs_term
                    scaled_sum = term / abs_term
                else if (abs_term > max_term) then
                    scaled_sum = scaled_sum * (max_term / abs_term) + term / abs_term
                    max_term = abs_term
                else
                    scaled_sum = scaled_sum + term / max_term
                end if
            end do
        end if

        if (max_term == 0.0_dp) then
            value = (0.0_dp, 0.0_dp)
            return
        end if
        scaled_sum = max_term * scaled_sum
        if (kind_index == 1) then
            sign_index = 1
        else
            sign_index = -1
        end if
        value = scaled_sum * exp(real(sign_index, dp) * imaginary_unit * z)
    end function complex_hankel_derivative

    pure elemental function positive_hankel(v, x, kind_index) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite real argument
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        real(dp) :: j_value
        real(dp) :: y_value

        j_value = jv(v, x)
        y_value = yv(v, x)
        if (kind_index == 1) then
            value = cmplx(j_value, y_value, kind=dp)
        else
            value = cmplx(j_value, -y_value, kind=dp)
        end if
    end function positive_hankel

    pure elemental function positive_hankel_derivative(v, x, n, kind_index) result(value)
        real(dp), intent(in) :: v !! finite real order
        real(dp), intent(in) :: x !! positive finite real argument
        integer, intent(in) :: n !! positive derivative order
        integer, intent(in) :: kind_index !! Hankel kind, 1 or 2
        complex(dp) :: value

        real(dp) :: j_value
        real(dp) :: y_value

        j_value = jvp(v, x, n)
        y_value = yvp(v, x, n)
        if (kind_index == 1) then
            value = cmplx(j_value, y_value, kind=dp)
        else
            value = cmplx(j_value, -y_value, kind=dp)
        end if
    end function positive_hankel_derivative

    pure elemental logical function complex_is_zero(z) result(zero)
        complex(dp), intent(in) :: z !! finite complex value to test for exact zero

        zero = is_zero(real(z, dp)) .and. is_zero(aimag(z))
    end function complex_is_zero

    pure elemental logical function complex_is_finite(z) result(finite)
        complex(dp), intent(in) :: z !! complex value to classify

        finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure elemental logical function invalid_complex_inputs(v, z) result(invalid)
        real(dp), intent(in) :: v !! Hankel order to validate
        complex(dp), intent(in) :: z !! complex argument to validate

        invalid = ieee_is_nan(v) .or. ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
        if (invalid) return
        invalid = .not. ieee_is_finite(v) .or. .not. complex_is_finite(z)
    end function invalid_complex_inputs

    pure elemental logical function is_zero(x) result(zero)
        real(dp), intent(in) :: x !! finite real value to test for either signed zero

        zero = .not. (x < 0.0_dp .or. x > 0.0_dp)
    end function is_zero

    pure elemental logical function invalid_real_inputs(v, x) result(invalid)
        real(dp), intent(in) :: v !! Hankel order to validate
        real(dp), intent(in) :: x !! real argument to validate

        invalid = ieee_is_nan(v) .or. ieee_is_nan(x)
        if (invalid) return
        invalid = .not. ieee_is_finite(v) .or. .not. ieee_is_finite(x)
    end function invalid_real_inputs

    pure elemental function phase_pi(v, sign_imaginary) result(value)
        real(dp), intent(in) :: v !! finite real multiplier of pi
        integer, intent(in) :: sign_imaginary !! +1 for exp(+i*pi*v), -1 for exp(-i*pi*v)
        complex(dp) :: value

        value = cmplx(cospi_value(v), &
            real(sign_imaginary, dp) * sinpi_value(v), kind=dp)
    end function phase_pi

    pure elemental function unit_phase(theta) result(value)
        real(dp), intent(in) :: theta !! finite real phase in radians
        complex(dp) :: value

        value = cmplx(cos(theta), sin(theta), kind=dp)
    end function unit_phase

    pure elemental function sinpi_value(v) result(value)
        real(dp), intent(in) :: v !! finite real argument of sin(pi*v)
        real(dp) :: value

        real(dp) :: r
        real(dp), parameter :: pi = acos(-1.0_dp)

        r = modulo(v, 2.0_dp)
        if (r <= 0.5_dp) then
            value = sin(pi * r)
        else if (r <= 1.5_dp) then
            value = sin(pi * (1.0_dp - r))
        else
            value = sin(pi * (r - 2.0_dp))
        end if
    end function sinpi_value

    pure elemental function cospi_value(v) result(value)
        real(dp), intent(in) :: v !! finite real argument of cos(pi*v)
        real(dp) :: value

        real(dp) :: r
        real(dp), parameter :: pi = acos(-1.0_dp)

        r = modulo(v, 2.0_dp)
        if (r <= 1.0_dp) then
            value = cos(pi * r)
        else
            value = cos(pi * (2.0_dp - r))
        end if
    end function cospi_value

    pure elemental function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

end module scifort_hankel
