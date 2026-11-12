! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_integrate_quadpack_modes
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
    use scifort_integrate, only : dp, quad, quad_status_invalid_argument, quad_status_success
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: abserr
    real(dp) :: result
    real(dp) :: xneg_inf
    real(dp) :: xpos_inf
    integer :: neval
    integer :: status

    xpos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    xneg_inf = ieee_value(0.0_dp, ieee_negative_inf)

    call quad(cusp, 0.0_dp, 1.0_dp, result, abserr, points=[0.3_dp], &
        neval=neval, status=status)
    call check_true('quad points status', status == quad_status_success .and. neval > 0)
    call check_close('quad points cusp', result, 0.29_dp, 2.0e-8_dp)

    call quad(cusp, 1.0_dp, 0.0_dp, result, abserr, &
        points=[0.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], status=status)
    call check_true('quad reversed points status', status == quad_status_success)
    call check_close('quad reversed points cusp', result, -0.29_dp, 2.0e-8_dp)

    call quad(cusp, 0.0_dp, 1.0_dp, result, abserr, points=[1.2_dp], status=status)
    call check_true('quad points outside interval', status == quad_status_invalid_argument)

    call quad(exponential_decay, 0.0_dp, 2.0_dp, result, abserr, &
        weight='cos', wvar=3.0_dp, status=status)
    call check_true('quad finite cos status', status == quad_status_success)
    call check_close('quad finite cos', result, 0.07566107555324406_dp, 3.0e-8_dp)

    call quad(exponential_decay, 0.0_dp, 2.0_dp, result, abserr, &
        weight='sin', wvar=3.0_dp, status=status)
    call check_true('quad finite sin status', status == quad_status_success)
    call check_close('quad finite sin', result, 0.2647980022491830_dp, 3.0e-8_dp)

    call quad(exponential_decay, 0.0_dp, xpos_inf, result, abserr, &
        weight='cos', wvar=2.0_dp, status=status)
    call check_true('quad infinite cos status', status == quad_status_success)
    call check_close('quad infinite cos', result, 0.2_dp, 2.0e-7_dp)

    call quad(exponential_decay, 0.0_dp, xpos_inf, result, abserr, &
        weight='sin', wvar=2.0_dp, status=status)
    call check_true('quad infinite sin status', status == quad_status_success)
    call check_close('quad infinite sin', result, 0.4_dp, 2.0e-7_dp)

    call quad(exponential_growth, xneg_inf, 0.0_dp, result, abserr, &
        weight='cos', wvar=2.0_dp, status=status)
    call check_true('quad left-infinite cos status', status == quad_status_success)
    call check_close('quad left-infinite cos', result, 0.2_dp, 2.0e-7_dp)

    call quad(exponential_growth, xneg_inf, 0.0_dp, result, abserr, &
        weight='sin', wvar=2.0_dp, status=status)
    call check_true('quad left-infinite sin status', status == quad_status_success)
    call check_close('quad left-infinite sin', result, -0.4_dp, 2.0e-7_dp)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='alg', &
        wvar=0.5_dp, wvar2=0.25_dp, status=status)
    call check_true('quad alg status', status == quad_status_success)
    call check_close('quad alg', result, 0.49943953415088005_dp, 5.0e-6_dp)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='alg-loga', &
        wvar=0.5_dp, wvar2=0.25_dp, status=status)
    call check_true('quad alg-loga status', status == quad_status_success)
    call check_close('quad alg-loga', result, -0.39076701082437015_dp, 8.0e-6_dp)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='alg-logb', &
        wvar=0.5_dp, wvar2=0.25_dp, status=status)
    call check_true('quad alg-logb status', status == quad_status_success)
    call check_close('quad alg-logb', result, -0.5225908331798218_dp, 8.0e-6_dp)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='alg-log', &
        wvar=0.5_dp, wvar2=0.25_dp, status=status)
    call check_true('quad alg-log status', status == quad_status_success)
    call check_close('quad alg-log', result, 0.19034045795155038_dp, 1.2e-5_dp)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='cauchy', &
        wvar=0.25_dp, status=status)
    call check_true('quad cauchy status', status == quad_status_success)
    call check_close('quad cauchy', result, log(3.0_dp), 3.0e-8_dp)

    call quad(one, 1.0_dp, 0.0_dp, result, abserr, weight='cauchy', &
        wvar=0.25_dp, status=status)
    call check_true('quad reversed cauchy status', status == quad_status_success)
    call check_close('quad reversed cauchy', result, -log(3.0_dp), 3.0e-8_dp)

    call quad(one, xneg_inf, xpos_inf, result, abserr, weight='cos', &
        wvar=2.0_dp, status=status)
    call check_true('quad whole-line weighted invalid', status == quad_status_invalid_argument)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, points=[0.5_dp], &
        weight='cos', wvar=2.0_dp, status=status)
    call check_true('quad points and weight invalid', status == quad_status_invalid_argument)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='cauchy', &
        wvar=0.0_dp, status=status)
    call check_true('quad cauchy endpoint invalid', status == quad_status_invalid_argument)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='alg', &
        wvar=-1.0_dp, wvar2=0.0_dp, status=status)
    call check_true('quad alg exponent invalid', status == quad_status_invalid_argument)

    call quad(one, 0.0_dp, 1.0_dp, result, abserr, weight='unknown', &
        wvar=1.0_dp, status=status)
    call check_true('quad weight name invalid', status == quad_status_invalid_argument)

    call quad(exponential_decay, 0.0_dp, xpos_inf, result, abserr, &
        weight='cos', wvar=2.0_dp, limlst=2, status=status)
    call check_true('quad limlst invalid', status == quad_status_invalid_argument)

    call quad(exponential_decay, 0.0_dp, pi, result, abserr, &
        weight='cos', wvar=2.0_dp, maxp1=0, status=status)
    call check_true('quad maxp1 invalid', status == quad_status_invalid_argument)

    print *, 'integrate quadpack mode tests passed'

contains

    function cusp(x) result(value)
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        value = abs(x - 0.3_dp)
    end function cusp

    function exponential_decay(x) result(value)
        real(dp), intent(in) :: x !! nonnegative evaluation point
        real(dp) :: value

        value = exp(-x)
    end function exponential_decay

    function exponential_growth(x) result(value)
        real(dp), intent(in) :: x !! nonpositive evaluation point
        real(dp) :: value

        value = exp(x)
    end function exponential_growth

    function one(x) result(value)
        real(dp), intent(in) :: x !! evaluation point, unused
        real(dp) :: value

        value = 1.0_dp + 0.0_dp * x
    end function one

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

end program test_integrate_quadpack_modes
