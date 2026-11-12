! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_integrate_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, &
        ieee_positive_inf, ieee_value
    use scifort_integrate, only : dp, integrate_status_success, qmc_quad, quad_status_success, &
        quad_vec, tanhsinh, tanhsinh_status_maxlevel
    use scifort_random, only : rng_seed, rng_state
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: abserr
    real(dp) :: qmc_error
    real(dp) :: qmc_result
    real(dp) :: result
    real(dp) :: vector_result(3)
    real(dp) :: xneg_inf
    real(dp) :: xpos_inf
    type(rng_state) :: state
    integer :: neval
    integer :: status

    xpos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    xneg_inf = ieee_value(0.0_dp, ieee_negative_inf)

    call tanhsinh(square, 0.0_dp, 1.0_dp, result, abserr, nfev=neval, status=status)
    call check_true('tanhsinh polynomial status', status == integrate_status_success .and. neval > 0)
    call check_close('tanhsinh polynomial', result, 1.0_dp / 3.0_dp, 3.0e-14_dp)

    call tanhsinh(inverse_sqrt, 0.0_dp, 1.0_dp, result, abserr, status=status)
    call check_true('tanhsinh endpoint singular status', status == integrate_status_success)
    call check_close('tanhsinh endpoint singular', result, 2.0_dp, 8.0e-13_dp)

    call tanhsinh(exponential_decay, 0.0_dp, xpos_inf, result, abserr, status=status)
    call check_true('tanhsinh semi-infinite status', status == integrate_status_success)
    call check_close('tanhsinh semi-infinite', result, 1.0_dp, 8.0e-13_dp)

    call tanhsinh(gaussian, xneg_inf, xpos_inf, result, abserr, status=status)
    call check_true('tanhsinh whole-line status', status == integrate_status_success)
    call check_close('tanhsinh whole-line', result, sqrt(pi), 8.0e-13_dp)

    call tanhsinh(gaussian, xpos_inf, xneg_inf, result, abserr, status=status)
    call check_true('tanhsinh reversed whole-line status', status == integrate_status_success)
    call check_close('tanhsinh reversed whole-line', result, -sqrt(pi), 8.0e-13_dp)


    call tanhsinh(slow_tail, 1.0_dp, xpos_inf, result, abserr, status=status)
    call check_true('tanhsinh slow-tail maxlevel', status == tanhsinh_status_maxlevel)

    call quad_vec(vector_function, 0.0_dp, 1.0_dp, vector_result, abserr, &
        norm='2', neval=neval, status=status)
    call check_true('quad vec status', status == quad_status_success .and. neval > 0)
    call check_close('quad vec sine', vector_result(1), 1.0_dp - cos(1.0_dp), 2.0e-10_dp)
    call check_close('quad vec square', vector_result(2), 1.0_dp / 3.0_dp, 2.0e-10_dp)
    call check_close('quad vec exponential', vector_result(3), 1.0_dp - exp(-1.0_dp), 2.0e-10_dp)

    call rng_seed(state, 20260928)
    call qmc_quad(product_integrand, [0.0_dp, 0.0_dp], [1.0_dp, 2.0_dp], &
        qmc_result, qmc_error, state, n_estimates=8, n_points=2048, neval=neval, status=status)
    call check_true('qmc quad status', status == integrate_status_success .and. neval == 16384)
    call check_close('qmc quad product', qmc_result, 1.0_dp, 2.0e-3_dp)
    call check_true('qmc quad standard error', qmc_error >= 0.0_dp .and. qmc_error < 2.0e-3_dp)

    call rng_seed(state, 17)
    call qmc_quad(one_nd, [1.0_dp, 0.0_dp], [0.0_dp, 2.0_dp], &
        qmc_result, qmc_error, state, n_estimates=4, n_points=64, status=status)
    call check_true('qmc reversed bounds status', status == integrate_status_success)
    call check_close('qmc reversed bounds', qmc_result, -2.0_dp, 2.0e-15_dp)
    call check_close('qmc constant standard error', qmc_error, 0.0_dp, 0.0_dp)

    call rng_seed(state, 17)
    call qmc_quad(one_nd, [0.0_dp], [1.0_dp], qmc_result, qmc_error, state, &
        n_estimates=1, n_points=16, status=status)
    call check_true('qmc one estimate standard error', &
        status == integrate_status_success .and. ieee_is_nan(qmc_error))

    print *, 'advanced integrate tests passed'

contains

    function square(x) result(value)
        real(dp), intent(in) :: x !! evaluation point
        real(dp) :: value

        value = x * x
    end function square

    function inverse_sqrt(x) result(value)
        real(dp), intent(in) :: x !! strictly positive evaluation point near an endpoint
        real(dp) :: value

        value = 1.0_dp / sqrt(x)
    end function inverse_sqrt

    function exponential_decay(x) result(value)
        real(dp), intent(in) :: x !! nonnegative evaluation point
        real(dp) :: value

        value = exp(-x)
    end function exponential_decay


    function slow_tail(x) result(value)
        real(dp), intent(in) :: x !! point on the slowly decaying semi-infinite tail
        real(dp) :: value

        value = x**(-1.01_dp)
    end function slow_tail

    function gaussian(x) result(value)
        real(dp), intent(in) :: x !! real evaluation point
        real(dp) :: value

        if (abs(x) > 30.0_dp) then
            value = 0.0_dp
        else
            value = exp(-x * x)
        end if
    end function gaussian

    subroutine vector_function(x, values)
        real(dp), intent(in) :: x !! scalar evaluation point
        real(dp), intent(out) :: values(:) !! three output components

        values(1) = sin(x)
        values(2) = x * x
        values(3) = exp(-x)
    end subroutine vector_function

    function product_integrand(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional integration point
        real(dp) :: value

        value = x(1) * x(2)
    end function product_integrand

    function one_nd(x) result(value)
        real(dp), intent(in) :: x(:) !! multidimensional integration point, unused
        real(dp) :: value

        value = 1.0_dp + 0.0_dp * sum(x)
    end function one_nd

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

end program test_integrate_advanced
