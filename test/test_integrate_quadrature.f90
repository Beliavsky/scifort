! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_integrate_quadrature
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
    use scifort_integrate, only : cumulative_simpson, cumulative_trapezoid, dp, fixed_quad, &
        integrate_status_invalid_argument, newton_cotes, quad, quad_status_invalid_argument, &
        quad_status_success, romb, simpson, trapezoid
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: weights(:)
    real(dp) :: abserr
    real(dp) :: error_coeff
    real(dp) :: result
    real(dp) :: xneg_inf
    real(dp) :: xpos_inf
    integer :: i
    integer :: neval
    integer :: status

    call check_close('trapezoid irregular', &
        trapezoid([0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp], &
        [0.0_dp, 0.5_dp, 2.0_dp, 3.0_dp]), 10.5_dp, 2.0e-15_dp)
    call check_close('trapezoid decreasing x', &
        trapezoid([1.0_dp, 4.0_dp, 9.0_dp], [3.0_dp, 2.0_dp, 0.0_dp]), &
        -15.5_dp, 2.0e-15_dp)

    values = cumulative_trapezoid([0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp], &
        [0.0_dp, 0.5_dp, 2.0_dp, 3.0_dp], initial=0.0_dp, status=status)
    call check_true('cumulative trapezoid status', status == 0)
    call check_array('cumulative trapezoid irregular', values, &
        [0.0_dp, 0.25_dp, 4.0_dp, 10.5_dp], 2.0e-15_dp)
    values = cumulative_trapezoid([1.0_dp, 2.0_dp], initial=1.0_dp, status=status)
    call check_true('cumulative trapezoid nonzero initial status', &
        status == integrate_status_invalid_argument .and. size(values) == 0)

    call check_close('simpson irregular even sample count', &
        simpson([0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp], &
        [0.0_dp, 0.5_dp, 2.0_dp, 3.0_dp]), 10.3_dp, 3.0e-15_dp)
    call check_close('simpson irregular odd sample count', &
        simpson([1.0_dp, -2.0_dp, 0.5_dp, 4.0_dp, 3.0_dp], &
        [0.0_dp, 0.3_dp, 1.1_dp, 2.0_dp, 4.0_dp]), &
        9.40488425925926_dp, 5.0e-15_dp)
    call check_close('simpson repeated x compatibility', &
        simpson([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
        [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]), &
        4.333333333333333_dp, 3.0e-15_dp)

    values = cumulative_simpson([0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp], &
        [0.0_dp, 0.5_dp, 2.0_dp, 3.0_dp], initial=2.5_dp, status=status)
    call check_true('cumulative simpson status', status == 0)
    call check_array('cumulative simpson irregular', values, &
        [2.5_dp, 2.75_dp, 6.5_dp, 12.8_dp], 4.0e-15_dp)
    values = cumulative_simpson([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
        dx=0.0_dp, initial=3.0_dp, status=status)
    call check_true('cumulative simpson zero dx status', status == 0)
    call check_array('cumulative simpson zero dx', values, &
        [3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp], 0.0_dp)
    values = cumulative_simpson([1.0_dp, 2.0_dp, 4.0_dp], &
        [0.0_dp, 1.0_dp, 1.0_dp], status=status)
    call check_true('cumulative simpson requires increasing x', &
        status == integrate_status_invalid_argument .and. size(values) == 0)

    call check_close('romb 5 samples', &
        romb(exp([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]), dx=0.25_dp), &
        1.7182826879247572_dp, 3.0e-15_dp)
    call check_close('romb 9 samples', &
        romb(exp([0.0_dp, 0.125_dp, 0.25_dp, 0.375_dp, 0.5_dp, 0.625_dp, &
        0.75_dp, 0.875_dp, 1.0_dp]), dx=0.125_dp), &
        1.7182818287945303_dp, 4.0e-15_dp)

    result = fixed_quad(polynomial8, 0.0_dp, 1.0_dp, n=5, status=status)
    call check_true('fixed quad status', status == 0)
    call check_close('fixed quad degree 8', result, 1.0_dp / 9.0_dp, 4.0e-15_dp)

    call newton_cotes(4, weights, error_coeff, status)
    call check_true('newton cotes order status', status == 0)
    call check_array('newton cotes order 4 weights', weights, &
        [0.3111111111111111_dp, 1.4222222222222223_dp, 0.5333333333333333_dp, &
        1.4222222222222223_dp, 0.3111111111111111_dp], 3.0e-15_dp)
    call check_close('newton cotes order 4 error', error_coeff, &
        -0.008465608465608466_dp, 2.0e-17_dp)
    call newton_cotes([0.0_dp, 0.4_dp, 1.5_dp, 3.0_dp], weights, error_coeff, status=status)
    call check_true('newton cotes irregular status', status == 0)
    call check_array('newton cotes irregular weights', weights, &
        [0.5_dp, 0.0_dp, 2.0_dp, 0.5_dp], 2.0e-15_dp)
    call check_close('newton cotes irregular error', error_coeff, &
        -0.0843750000000001_dp, 3.0e-16_dp)

    call newton_cotes(20, weights, error_coeff, status)
    call check_true('newton cotes order 20 status', status == 0)
    call check_close('newton cotes order 20 weight sum', sum(weights), 20.0_dp, 1.0e-10_dp)
    call check_close('newton cotes order 20 normalized moment', &
        sum(weights * [(2.0_dp * real(i, dp) / 20.0_dp - 1.0_dp, i = 0, 20)] ** 20), &
        20.0_dp / 21.0_dp, 2.0e-11_dp)
    call check_true('newton cotes order 20 symmetry', &
        maxval(abs(weights - weights(size(weights):1:-1))) < 1.0e-8_dp)

    call quad(sine, 0.0_dp, pi, result, abserr, neval=neval, status=status)
    call check_true('quad finite status', status == quad_status_success .and. neval > 0)
    call check_close('quad finite sine', result, 2.0_dp, 2.0e-10_dp)

    call quad(sine, pi, 0.0_dp, result, abserr, status=status)
    call check_true('quad reversed finite status', status == quad_status_success)
    call check_close('quad reversed finite', result, -2.0_dp, 2.0e-10_dp)

    xpos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    xneg_inf = ieee_value(0.0_dp, ieee_negative_inf)
    call quad(exponential_tail, 0.0_dp, xpos_inf, result, abserr, status=status)
    call check_true('quad semi-infinite status', status == quad_status_success)
    call check_close('quad semi-infinite exponential', result, 1.0_dp, 2.0e-8_dp)

    call quad(gaussian, xneg_inf, xpos_inf, result, abserr, status=status)
    call check_true('quad whole-line status', status == quad_status_success)
    call check_close('quad whole-line gaussian', result, sqrt(pi), 2.0e-8_dp)

    call quad(gaussian, xpos_inf, xneg_inf, result, abserr, status=status)
    call check_true('quad reversed whole-line status', status == quad_status_success)
    call check_close('quad reversed whole-line', result, -sqrt(pi), 2.0e-8_dp)

    call quad(sine, 0.0_dp, 1.0_dp, result, abserr, limit=0, status=status)
    call check_true('quad invalid limit status', status == quad_status_invalid_argument)

    print *, 'integrate quadrature tests passed'

contains

    function polynomial8(x) result(value)
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        value = x ** 8
    end function polynomial8

    function sine(x) result(value)
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        value = sin(x)
    end function sine

    function exponential_tail(x) result(value)
        real(dp), intent(in) :: x !! nonnegative evaluation point
        real(dp) :: value

        value = exp(-x)
    end function exponential_tail

    function gaussian(x) result(value)
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        value = exp(-x * x)
    end function gaussian

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', name, ' size ', size(actual), size(expected)
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_array

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
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_integrate_quadrature
