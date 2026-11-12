! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_integrate_nested
    use scifort_integrate, only : dblquad, dp, nquad, quad_status_invalid_argument, &
        quad_status_success, tplquad
    implicit none

    real(dp) :: abserr
    real(dp) :: result
    integer :: status

    call dblquad(double_polynomial, 0.0_dp, 2.0_dp, zero_bound, one_bound, &
        result, abserr, status=status)
    call check_true('dblquad box status', status == quad_status_success)
    call check_close('dblquad box', result, 2.0_dp / 3.0_dp, 2.0e-11_dp)

    call dblquad(double_variable, 0.0_dp, 1.0_dp, lower_y, upper_y, &
        result, abserr, status=status)
    call check_true('dblquad variable status', status == quad_status_success)
    call check_close('dblquad variable', result, 1.0_dp / 3.0_dp, 2.0e-11_dp)

    call tplquad(triple_polynomial, 0.0_dp, 1.0_dp, zero_bound, one_bound, &
        zero_surface, one_surface, result, abserr, status=status)
    call check_true('tplquad box status', status == quad_status_success)
    call check_close('tplquad box', result, 0.125_dp, 3.0e-11_dp)

    call tplquad(triple_polynomial, 0.0_dp, 1.0_dp, zero_bound, upper_middle, &
        zero_surface, upper_inner, result, abserr, status=status)
    call check_true('tplquad variable status', status == quad_status_success)
    call check_close('tplquad variable', result, 0.05416666666666668_dp, 1.0e-9_dp)

    call nquad(nd_product, 2, unit_ranges, result, abserr, status=status)
    call check_true('nquad unit square status', status == quad_status_success)
    call check_close('nquad unit square', result, 0.25_dp, 3.0e-11_dp)

    call nquad(nd_product, 3, variable_ranges, result, abserr, status=status)
    call check_true('nquad variable status', status == quad_status_success)
    call check_close('nquad variable', result, 0.05416666666666668_dp, 1.0e-9_dp)

    call nquad(nd_product, 0, unit_ranges, result, abserr, status=status)
    call check_true('nquad invalid dimension', status == quad_status_invalid_argument)

    print *, 'integrate nested tests passed'

contains

    function double_polynomial(y, x) result(value)
        real(dp), intent(in) :: y !! inner integration coordinate
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = x * y * y
    end function double_polynomial

    function double_variable(y, x) result(value)
        real(dp), intent(in) :: y !! inner integration coordinate
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = x * y
    end function double_variable

    function triple_polynomial(z, y, x) result(value)
        real(dp), intent(in) :: z !! innermost integration coordinate
        real(dp), intent(in) :: y !! middle integration coordinate
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = x * y * z
    end function triple_polynomial

    function zero_bound(x) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate, unused
        real(dp) :: value

        value = 0.0_dp + 0.0_dp * x
    end function zero_bound

    function one_bound(x) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate, unused
        real(dp) :: value

        value = 1.0_dp + 0.0_dp * x
    end function one_bound

    function lower_y(x) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = x
    end function lower_y

    function upper_y(x) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = 2.0_dp - x
    end function upper_y

    function upper_middle(x) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp) :: value

        value = 1.0_dp - 2.0_dp * x
    end function upper_middle

    function zero_surface(x, y) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate, unused
        real(dp), intent(in) :: y !! middle integration coordinate, unused
        real(dp) :: value

        value = 0.0_dp + 0.0_dp * (x + y)
    end function zero_surface

    function one_surface(x, y) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate, unused
        real(dp), intent(in) :: y !! middle integration coordinate, unused
        real(dp) :: value

        value = 1.0_dp + 0.0_dp * (x + y)
    end function one_surface

    function upper_inner(x, y) result(value)
        real(dp), intent(in) :: x !! outer integration coordinate
        real(dp), intent(in) :: y !! middle integration coordinate
        real(dp) :: value

        value = 1.0_dp - x - 2.0_dp * y
    end function upper_inner

    function nd_product(x) result(value)
        real(dp), intent(in) :: x(:) !! integration coordinates, innermost first
        real(dp) :: value

        value = product(x)
    end function nd_product

    subroutine unit_ranges(axis, x, lower, upper)
        integer, intent(in) :: axis !! integration axis, unused
        real(dp), intent(in) :: x(:) !! integration coordinate workspace, unused
        real(dp), intent(out) :: lower !! lower axis bound
        real(dp), intent(out) :: upper !! upper axis bound

        lower = 0.0_dp + 0.0_dp * real(axis, dp) + 0.0_dp * sum(x)
        upper = 1.0_dp
    end subroutine unit_ranges

    subroutine variable_ranges(axis, x, lower, upper)
        integer, intent(in) :: axis !! integration axis, 1=z, 2=y, 3=x
        real(dp), intent(in) :: x(:) !! coordinate workspace with outer coordinates filled
        real(dp), intent(out) :: lower !! lower axis bound
        real(dp), intent(out) :: upper !! upper axis bound

        lower = 0.0_dp
        select case (axis)
        case (3)
            upper = 1.0_dp
        case (2)
            upper = 1.0_dp - 2.0_dp * x(3)
        case (1)
            upper = 1.0_dp - x(3) - 2.0_dp * x(2)
        case default
            upper = 0.0_dp
        end select
    end subroutine variable_ranges

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition that must be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_integrate_nested
