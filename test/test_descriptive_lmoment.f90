! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_lmoment
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : lmoment, lmoment_all, lmoment_keepdims
    implicit none

    real(dp), parameter :: tol = 8.0e-13_dp
    real(dp) :: x(6)
    real(dp) :: matrix(3, 4)
    real(dp) :: qnan
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: matrix_values(:, :)
    real(dp), allocatable :: scalar_axis(:)
    real(dp), allocatable :: keep_default(:, :, :)
    real(dp), allocatable :: keep_orders(:, :, :)
    real(dp), allocatable :: keep_scalar(:, :)
    real(dp) :: rank15(1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: values15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: keep15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    integer :: i
    integer :: shape15(15)

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]

    values = lmoment(x)
    call check_array('lmoment default', values, &
        [5.5_dp, 3.166666666666666_dp, 0.5368421052631575_dp, &
            0.2631578947368491_dp])
    call check_close('lmoment scalar order1', lmoment(x, 1), 5.5_dp)
    call check_close('lmoment scalar order2', lmoment(x, 2), 3.166666666666666_dp)
    call check_close('lmoment scalar order5', lmoment(x, 5), -0.05263157894739012_dp)
    call check_close('lmoment scalar order6', lmoment(x, 6), 0.26315789473684237_dp)
    call check_nan('lmoment scalar order7 insufficient', lmoment(x, 7))

    values = lmoment(x, [1, 3, 5])
    call check_array('lmoment vector orders', values, &
        [5.5_dp, 0.5368421052631575_dp, -0.05263157894739012_dp])
    values = lmoment(x, standardize=.false.)
    call check_array('lmoment raw', values, &
        [5.5_dp, 3.166666666666666_dp, 1.6999999999999984_dp, &
            0.8333333333333552_dp])
    values = lmoment([3.0_dp, 1.0_dp, 2.0_dp], sorted=.true.)
    call check_array_nan('lmoment sorted true', values, &
        [2.0_dp, -0.33333333333333326_dp, -3.0_dp, qnan])

    matrix = reshape([(real(i, dp), i = 1, 12)], shape(matrix))
    matrix_values = lmoment(matrix, axis=1)
    if (any(shape(matrix_values) /= [4, 4])) then
        print *, 'FAIL lmoment matrix axis1 shape', shape(matrix_values)
        error stop 1
    end if
    call check_array_nan('lmoment matrix axis1', reshape(matrix_values, [size(matrix_values)]), &
        [2.0_dp, 0.6666666666666665_dp, 0.0_dp, qnan, &
         5.0_dp, 0.6666666666666670_dp, 0.0_dp, qnan, &
         8.0_dp, 0.6666666666666661_dp, 0.0_dp, qnan, &
         11.0_dp, 0.6666666666666661_dp, 0.0_dp, qnan])

    scalar_axis = lmoment(matrix, 2, axis=2)
    call check_array('lmoment matrix scalar axis2', scalar_axis, [2.5_dp, 2.5_dp, 2.5_dp])
    matrix_values = lmoment(matrix, [1, 3], axis=2)
    if (any(shape(matrix_values) /= [2, 3])) then
        print *, 'FAIL lmoment matrix vector axis2 shape', shape(matrix_values)
        error stop 1
    end if
    call check_array('lmoment matrix vector axis2', reshape(matrix_values, [size(matrix_values)]), &
        [5.5_dp, 0.0_dp, 6.5_dp, 0.0_dp, 7.5_dp, 0.0_dp])

    keep_default = lmoment_keepdims(matrix, axis=1)
    if (any(shape(keep_default) /= [1, 4, 4])) then
        print *, 'FAIL lmoment default keepdims shape', shape(keep_default)
        error stop 1
    end if
    keep_orders = lmoment_keepdims(matrix, [1, 3], axis=2)
    if (any(shape(keep_orders) /= [2, 3, 1])) then
        print *, 'FAIL lmoment vector keepdims shape', shape(keep_orders)
        error stop 1
    end if
    keep_scalar = lmoment_keepdims(matrix, 2, axis=2)
    if (any(shape(keep_scalar) /= [3, 1])) then
        print *, 'FAIL lmoment scalar keepdims shape', shape(keep_scalar)
        error stop 1
    end if
    call check_array('lmoment scalar keepdims values', reshape(keep_scalar, [size(keep_scalar)]), &
        [2.5_dp, 2.5_dp, 2.5_dp])

    values = lmoment_all(matrix)
    call check_array('lmoment all default', values, &
        [6.5_dp, 2.166666666666666_dp, 0.0_dp, 8.198570028001158e-16_dp])
    call check_close('lmoment all scalar', lmoment_all(matrix, 2), 2.166666666666666_dp)
    values = lmoment_all(matrix, [1, 3])
    call check_array('lmoment all vector', values, [6.5_dp, 0.0_dp])

    values = lmoment([1.0_dp, qnan, 3.0_dp], nan_policy='omit')
    call check_array_nan('lmoment omit', values, [2.0_dp, 1.0_dp, qnan, qnan])
    values = lmoment([1.0_dp, qnan, 3.0_dp], nan_policy='propagate')
    if (.not. all(ieee_is_nan(values))) then
        print *, 'FAIL lmoment propagate'
        error stop 1
    end if
    values = lmoment([1.0_dp, qnan, 3.0_dp], nan_policy='unsupported')
    if (.not. all(ieee_is_nan(values))) then
        print *, 'FAIL lmoment unsupported nan policy'
        error stop 1
    end if

    rank15 = 0.0_dp
    rank15(1, :, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = &
        [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    values15 = lmoment(rank15, axis=2)
    shape15 = shape(values15)
    if (shape15(1) /= 4 .or. any(shape15(2:) /= 1)) then
        print *, 'FAIL lmoment rank15 shape', shape(values15)
        error stop 1
    end if
    call check_array('lmoment rank15 values', reshape(values15, [size(values15)]), &
        [2.5_dp, 0.8333333333333333_dp, 0.0_dp, 0.0_dp])
    keep15 = lmoment_keepdims(rank15, 2, axis=2)
    if (any(shape(keep15) /= 1)) then
        print *, 'FAIL lmoment rank15 scalar keepdims shape', shape(keep15)
        error stop 1
    end if
    call check_close('lmoment rank15 scalar keepdims', &
        keep15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        0.8333333333333333_dp)

    print *, 'descriptive lmoment tests passed'

contains

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value

        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(label, actual)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual !! computed value

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', trim(label), actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ' size', size(actual), size(expected)
            error stop 1
        end if
        if (any(abs(actual - expected) > tol * max(1.0_dp, abs(expected)))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_array_nan(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values including NaNs

        integer :: j

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ' size', size(actual), size(expected)
            error stop 1
        end if
        do j = 1, size(actual)
            if (ieee_is_nan(expected(j))) then
                if (.not. ieee_is_nan(actual(j))) then
                    print *, 'FAIL ', trim(label), ' NaN at', j, actual(j)
                    error stop 1
                end if
            else if (abs(actual(j) - expected(j)) > tol * max(1.0_dp, abs(expected(j)))) then
                print *, 'FAIL ', trim(label), ' at', j, actual(j), expected(j)
                error stop 1
            end if
        end do
    end subroutine check_array_nan

end program test_descriptive_lmoment
