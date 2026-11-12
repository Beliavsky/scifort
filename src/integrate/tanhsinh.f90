! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent double-exponential tanh-sinh quadrature implementation based on
! H. Takahasi and M. Mori, Publ. RIMS Kyoto Univ. 9 (1974), 721-741,
! DOI 10.2977/PRIMS/1195192451. No SciPy implementation code is copied here.

module scifort_integrate_tanhsinh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_integrate_quadrature, only : integrate_status_invalid_argument, &
        integrate_status_numerical_failure, integrate_status_success
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: tanhsinh_status_maxlevel = 4

    public :: tanhsinh

    abstract interface
        function tanhsinh_integrand(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! finite integration coordinate
            real(dp) :: value
        end function tanhsinh_integrand
    end interface

contains

    subroutine tanhsinh(f, a, b, result, abserr, maxlevel, minlevel, atol, rtol, &
            nfev, status)
        procedure(tanhsinh_integrand) :: f !! scalar real integrand callback
        real(dp), intent(in) :: a !! lower integration limit; may be infinite
        real(dp), intent(in) :: b !! upper integration limit; may be infinite
        real(dp), intent(out) :: result !! double-exponential quadrature estimate
        real(dp), intent(out) :: abserr !! absolute difference between the last two refinement levels
        integer, intent(in), optional :: maxlevel !! maximum dyadic refinement level; default 10
        integer, intent(in), optional :: minlevel !! earliest convergence test level; default 2
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance; default zero
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance; default epsilon**0.75
        integer, intent(out), optional :: nfev !! total integrand evaluations across all levels
        integer, intent(out), optional :: status !! zero, max-level, invalid-argument, or numerical-failure status

        integer :: evaluations
        integer :: level
        integer :: selected_maxlevel
        integer :: selected_minlevel
        integer :: interval_kind
        integer :: level_status
        real(dp) :: absolute_tolerance
        real(dp) :: current
        real(dp) :: error_estimate
        real(dp) :: lower
        real(dp) :: orientation
        real(dp) :: previous
        real(dp) :: relative_tolerance
        real(dp) :: tail_error
        real(dp) :: upper

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        evaluations = 0
        if (present(nfev)) nfev = evaluations
        if (present(status)) status = integrate_status_invalid_argument

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
        selected_maxlevel = 10
        selected_minlevel = 2
        absolute_tolerance = 0.0_dp
        relative_tolerance = epsilon(1.0_dp)**0.75_dp
        if (present(maxlevel)) selected_maxlevel = maxlevel
        if (present(minlevel)) selected_minlevel = minlevel
        if (present(atol)) absolute_tolerance = atol
        if (present(rtol)) relative_tolerance = rtol
        if (selected_maxlevel < 0 .or. selected_minlevel < 0 .or. &
                selected_minlevel > selected_maxlevel) return
        if (.not. ieee_is_finite(absolute_tolerance) .or. &
                .not. ieee_is_finite(relative_tolerance)) return
        if (absolute_tolerance < 0.0_dp .or. relative_tolerance < 0.0_dp) return

        if (a == b) then
            result = 0.0_dp
            abserr = 0.0_dp
            if (present(status)) status = integrate_status_success
            return
        end if

        call normalize_interval(a, b, lower, upper, orientation, interval_kind, level_status)
        if (level_status /= integrate_status_success) return

        previous = quiet_nan(0.0_dp)
        error_estimate = huge(1.0_dp)
        do level = 0, selected_maxlevel
            call evaluate_level(f, lower, upper, interval_kind, orientation, level, &
                current, tail_error, evaluations, level_status)
            if (level_status /= integrate_status_success) then
                result = quiet_nan(0.0_dp)
                abserr = quiet_nan(0.0_dp)
                if (present(nfev)) nfev = evaluations
                if (present(status)) status = level_status
                return
            end if
            if (level > 0) then
                error_estimate = max(abs(current - previous), tail_error)
                if (level >= selected_minlevel) then
                    if (error_estimate <= max(absolute_tolerance, &
                            relative_tolerance * abs(current))) then
                        result = current
                        abserr = error_estimate
                        if (present(nfev)) nfev = evaluations
                        if (present(status)) status = integrate_status_success
                        return
                    end if
                end if
            end if
            previous = current
        end do

        result = current
        abserr = error_estimate
        if (present(nfev)) nfev = evaluations
        if (present(status)) status = tanhsinh_status_maxlevel
    end subroutine tanhsinh

    subroutine normalize_interval(a, b, lower, upper, orientation, interval_kind, status)
        real(dp), intent(in) :: a !! requested lower integration limit
        real(dp), intent(in) :: b !! requested upper integration limit
        real(dp), intent(out) :: lower !! normalized lower integration limit
        real(dp), intent(out) :: upper !! normalized upper integration limit
        real(dp), intent(out) :: orientation !! plus or minus one for normalized orientation
        integer, intent(out) :: interval_kind !! finite, right-infinite, left-infinite, or whole-line code
        integer, intent(out) :: status !! integrate_status_* result

        integer, parameter :: finite_interval = 1
        integer, parameter :: right_infinite_interval = 2
        integer, parameter :: left_infinite_interval = 3
        integer, parameter :: whole_line_interval = 4

        status = integrate_status_success
        orientation = 1.0_dp
        lower = a
        upper = b
        interval_kind = 0

        if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
            interval_kind = finite_interval
            if (b < a) then
                lower = b
                upper = a
                orientation = -1.0_dp
            end if
        elseif (ieee_is_finite(a) .and. .not. ieee_is_finite(b) .and. b > 0.0_dp) then
            interval_kind = right_infinite_interval
        elseif (.not. ieee_is_finite(a) .and. a < 0.0_dp .and. ieee_is_finite(b)) then
            interval_kind = left_infinite_interval
        elseif (.not. ieee_is_finite(a) .and. a < 0.0_dp .and. &
                .not. ieee_is_finite(b) .and. b > 0.0_dp) then
            interval_kind = whole_line_interval
        elseif (.not. ieee_is_finite(a) .and. a > 0.0_dp .and. ieee_is_finite(b)) then
            interval_kind = right_infinite_interval
            lower = b
            upper = a
            orientation = -1.0_dp
        elseif (ieee_is_finite(a) .and. .not. ieee_is_finite(b) .and. b < 0.0_dp) then
            interval_kind = left_infinite_interval
            lower = b
            upper = a
            orientation = -1.0_dp
        elseif (.not. ieee_is_finite(a) .and. a > 0.0_dp .and. &
                .not. ieee_is_finite(b) .and. b < 0.0_dp) then
            interval_kind = whole_line_interval
            lower = b
            upper = a
            orientation = -1.0_dp
        else
            status = integrate_status_invalid_argument
        end if
    end subroutine normalize_interval

    subroutine evaluate_level(f, lower, upper, interval_kind, orientation, level, &
            value, tail_error, evaluations, status)
        procedure(tanhsinh_integrand) :: f !! scalar integrand callback
        real(dp), intent(in) :: lower !! normalized lower limit
        real(dp), intent(in) :: upper !! normalized upper limit
        integer, intent(in) :: interval_kind !! normalized interval-type code
        real(dp), intent(in) :: orientation !! plus or minus one for original orientation
        integer, intent(in) :: level !! dyadic refinement level
        real(dp), intent(out) :: value !! quadrature estimate at this refinement level
        real(dp), intent(out) :: tail_error !! transformed-integrand magnitude at the truncation boundary
        integer, intent(inout) :: evaluations !! cumulative callback evaluation count
        integer, intent(out) :: status !! integrate_status_* result

        real(dp), parameter :: pi_over_two = 1.57079632679489661923132169163975144_dp
        real(dp), parameter :: t_limit = 6.0_dp
        integer, parameter :: finite_interval = 1
        integer, parameter :: right_infinite_interval = 2
        integer, parameter :: left_infinite_interval = 3
        integer, parameter :: whole_line_interval = 4

        integer :: k
        integer :: kmax
        real(dp) :: ch
        real(dp) :: exp_term
        real(dp) :: fvalue
        real(dp) :: h
        real(dp) :: half_width
        real(dp) :: jacobian
        real(dp) :: q
        real(dp) :: sech_squared
        real(dp) :: summation
        real(dp) :: t
        real(dp) :: term
        real(dp) :: x
        real(dp) :: y

        status = integrate_status_success
        value = 0.0_dp
        tail_error = 0.0_dp
        h = scale(1.0_dp, -level)
        kmax = ceiling(t_limit / h)
        summation = 0.0_dp
        half_width = 0.0_dp
        if (interval_kind == finite_interval) then
            half_width = 0.5_dp * upper - 0.5_dp * lower
        end if

        do k = -kmax, kmax
            t = real(k, dp) * h
            y = pi_over_two * sinh(t)
            ch = cosh(t)
            select case (interval_kind)
            case (finite_interval)
                if (y >= 0.0_dp) then
                    q = exp(-2.0_dp * y)
                    x = upper - half_width * (2.0_dp * q / (1.0_dp + q))
                else
                    q = exp(2.0_dp * y)
                    x = lower + half_width * (2.0_dp * q / (1.0_dp + q))
                end if
                sech_squared = 4.0_dp * q / (1.0_dp + q)**2
                jacobian = half_width * pi_over_two * ch * sech_squared
                if (x == lower .or. x == upper) cycle
            case (right_infinite_interval)
                exp_term = exp(2.0_dp * y)
                x = lower + exp_term
                jacobian = 2.0_dp * pi_over_two * ch * exp_term
            case (left_infinite_interval)
                exp_term = exp(-2.0_dp * y)
                x = upper - exp_term
                jacobian = 2.0_dp * pi_over_two * ch * exp_term
            case (whole_line_interval)
                x = 0.5_dp * sinh(2.0_dp * y)
                jacobian = pi_over_two * ch * cosh(2.0_dp * y)
            case default
                status = integrate_status_invalid_argument
                return
            end select

            if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(jacobian)) cycle
            if (jacobian == 0.0_dp) cycle
            fvalue = f(x)
            evaluations = evaluations + 1
            if (.not. ieee_is_finite(fvalue)) then
                status = integrate_status_numerical_failure
                return
            end if
            if (abs(jacobian) > 1.0_dp) then
                if (abs(fvalue) > huge(1.0_dp) / abs(jacobian)) then
                    status = integrate_status_numerical_failure
                    return
                end if
            end if
            term = fvalue * jacobian
            if (abs(k) == kmax) tail_error = max(tail_error, abs(term))
            if (term > 0.0_dp .and. summation > huge(1.0_dp) - term) then
                status = integrate_status_numerical_failure
                return
            end if
            if (term < 0.0_dp .and. summation < -huge(1.0_dp) - term) then
                status = integrate_status_numerical_failure
                return
            end if
            summation = summation + term
            if (.not. ieee_is_finite(summation)) then
                status = integrate_status_numerical_failure
                return
            end if
        end do
        value = orientation * h * summation
        if (.not. ieee_is_finite(value)) status = integrate_status_numerical_failure
    end subroutine evaluate_level

end module scifort_integrate_tanhsinh
