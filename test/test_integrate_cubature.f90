! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_integrate_cubature
    use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
    use scifort_integrate, only : cubature, cubature_status_max_subdivisions, dp, &
        integrate_status_success
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: abserr
    real(dp) :: points(2, 1)
    real(dp) :: result
    real(dp) :: vector_error(3)
    real(dp) :: vector_result(3)
    real(dp) :: xneg_inf
    real(dp) :: xpos_inf
    integer :: neval
    integer :: status
    integer :: subdivisions

    xpos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    xneg_inf = ieee_value(0.0_dp, ieee_negative_inf)

    call cubature(product_2d, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], result, abserr, &
        rule='genz-malik', neval=neval, subdivisions=subdivisions, status=status)
    call check_true('genz-malik product status', status == integrate_status_success)
    call check_close('genz-malik product', result, 0.25_dp, 2.0e-15_dp)
    call check_true('genz-malik product evaluations', neval == 17 .and. subdivisions == 0)

    call cubature(product_2d, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], result, abserr, &
        rule='gk21', neval=neval, subdivisions=subdivisions, status=status)
    call check_true('gk21 product status', status == integrate_status_success)
    call check_close('gk21 product', result, 0.25_dp, 2.0e-15_dp)
    call check_true('gk21 product evaluations', neval == 541 .and. subdivisions == 0)

    call cubature(vector_2d, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], &
        vector_result, vector_error, rule='genz-malik', status=status)
    call check_true('vector cubature status', status == integrate_status_success)
    call check_close('vector constant', vector_result(1), 1.0_dp, 2.0e-15_dp)
    call check_close('vector linear', vector_result(2), 1.0_dp, 2.0e-15_dp)
    call check_close('vector product', vector_result(3), 0.25_dp, 2.0e-15_dp)

    call cubature(product_2d, [1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp], result, abserr, &
        rule='genz-malik', status=status)
    call check_true('reversed bound status', status == integrate_status_success)
    call check_close('reversed bound sign', result, -0.25_dp, 2.0e-15_dp)

    points(:, 1) = [0.25_dp, 0.75_dp]
    call cubature(cusp_2d, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], result, abserr, &
        rule='genz-malik', rtol=1.0e-10_dp, points=points, &
        subdivisions=subdivisions, status=status)
    call check_true('break point status', status == integrate_status_success)
    call check_close('break point cusp', result, 0.625_dp, 2.0e-15_dp)
    call check_true('break point avoids refinement', subdivisions == 0)

    call cubature(gaussian_1d, [xneg_inf], [xpos_inf], result, abserr, &
        rule='gk15', rtol=1.0e-10_dp, status=status)
    call check_true('whole-line status', status == integrate_status_success)
    call check_close('whole-line gaussian', result, sqrt(pi), 5.0e-12_dp)

    call cubature(exponential_1d, [xneg_inf], [0.0_dp], result, abserr, &
        rule='gk21', rtol=1.0e-10_dp, status=status)
    call check_true('left-infinite status', status == integrate_status_success)
    call check_close('left-infinite exponential', result, 1.0_dp, 2.0e-13_dp)

    call cubature(cosine_sum_3d, [0.0_dp, 0.0_dp, 0.0_dp], &
        [1.0_dp, 1.0_dp, 1.0_dp], result, abserr, rule='genz-malik', &
        rtol=1.0e-15_dp, max_subdivisions=0, subdivisions=subdivisions, status=status)
    call check_true('max subdivisions status', &
        status == cubature_status_max_subdivisions .and. subdivisions == 0)
    call check_close('max subdivisions estimate', result, 3.0_dp * sin(1.0_dp), 4.0e-7_dp)

    print *, 'cubature integration tests passed'

contains

    function product_2d(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional integration point
        real(dp) :: value

        value = x(1) * x(2)
    end function product_2d

    subroutine vector_2d(x, values)
        real(dp), intent(in) :: x(:) !! two-dimensional integration point
        real(dp), intent(out) :: values(:) !! constant, linear, and product outputs

        values(1) = 1.0_dp
        values(2) = x(1) + x(2)
        values(3) = x(1) * x(2)
    end subroutine vector_2d

    function cusp_2d(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional integration point
        real(dp) :: value

        value = abs(x(1) - 0.25_dp) + abs(x(2) - 0.75_dp)
    end function cusp_2d

    function gaussian_1d(x) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional integration point
        real(dp) :: value

        value = exp(-x(1) * x(1))
    end function gaussian_1d

    function exponential_1d(x) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional integration point
        real(dp) :: value

        value = exp(x(1))
    end function exponential_1d

    function cosine_sum_3d(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional integration point
        real(dp) :: value

        value = sum(cos(x))
    end function cosine_sum_3d

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

end program test_integrate_cubature
