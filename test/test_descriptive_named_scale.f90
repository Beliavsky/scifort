! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0.

program test_descriptive_named_scale
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_stats, only : dp, iqr, iqr_keepdims, median_abs_deviation, &
        median_abs_deviation_keepdims
    implicit none

    real(dp), parameter :: tol = 5.0e-14_dp
    real(dp), parameter :: iqr_normal_two = 0.7413011092528009_dp
    real(dp), parameter :: mad_normal_two = 1.482602218505602_dp
    real(dp) :: x1(2)
    real(dp) :: x2(2, 1)
    real(dp) :: x3(2, 1, 1)
    real(dp) :: x4(2, 1, 1, 1)
    real(dp) :: x5(2, 1, 1, 1, 1)
    real(dp) :: x6(2, 1, 1, 1, 1, 1)
    real(dp) :: x7(2, 1, 1, 1, 1, 1, 1)
    real(dp) :: x8(2, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x9(2, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x10(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: xnan(3)
    real(dp) :: matrix(3, 4)
    real(dp) :: nan
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: kept(:, :)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    x1 = reshape([1.0_dp, 3.0_dp], shape(x1))
    x2 = reshape([1.0_dp, 3.0_dp], shape(x2))
    x3 = reshape([1.0_dp, 3.0_dp], shape(x3))
    x4 = reshape([1.0_dp, 3.0_dp], shape(x4))
    x5 = reshape([1.0_dp, 3.0_dp], shape(x5))
    x6 = reshape([1.0_dp, 3.0_dp], shape(x6))
    x7 = reshape([1.0_dp, 3.0_dp], shape(x7))
    x8 = reshape([1.0_dp, 3.0_dp], shape(x8))
    x9 = reshape([1.0_dp, 3.0_dp], shape(x9))
    x10 = reshape([1.0_dp, 3.0_dp], shape(x10))
    x11 = reshape([1.0_dp, 3.0_dp], shape(x11))
    x12 = reshape([1.0_dp, 3.0_dp], shape(x12))
    x13 = reshape([1.0_dp, 3.0_dp], shape(x13))
    x14 = reshape([1.0_dp, 3.0_dp], shape(x14))
    x15 = reshape([1.0_dp, 3.0_dp], shape(x15))
    xnan = [1.0_dp, nan, 3.0_dp]
    matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 8.0_dp, 12.0_dp, &
        9.0_dp, 18.0_dp, 27.0_dp, 16.0_dp, 32.0_dp, 48.0_dp], [3, 4])

    call check_close(iqr(x1, scale='normal'), iqr_normal_two, 'iqr rank 1 normal')
    call check_close(median_abs_deviation(x1, scale='NORMAL'), mad_normal_two, &
        'mad rank 1 NORMAL')
    call check_close(sum(iqr_keepdims(x1, scale='normal')), iqr_normal_two, &
        'iqr rank 1 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x1, scale='normal')), &
        mad_normal_two, 'mad rank 1 keepdims normal')

    call check_close(sum(iqr(x2, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 2 axis normal')
    call check_close(iqr(x2, scale='normal'), iqr_normal_two, &
        'iqr rank 2 all normal')
    call check_close(sum(median_abs_deviation(x2, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 2 normal')
    call check_close(sum(iqr_keepdims(x2, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 2 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x2, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 2 keepdims normal')

    call check_close(sum(iqr(x3, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 3 axis normal')
    call check_close(iqr(x3, scale='normal'), iqr_normal_two, &
        'iqr rank 3 all normal')
    call check_close(sum(median_abs_deviation(x3, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 3 normal')
    call check_close(sum(iqr_keepdims(x3, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 3 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x3, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 3 keepdims normal')

    call check_close(sum(iqr(x4, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 4 axis normal')
    call check_close(iqr(x4, scale='normal'), iqr_normal_two, &
        'iqr rank 4 all normal')
    call check_close(sum(median_abs_deviation(x4, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 4 normal')
    call check_close(sum(iqr_keepdims(x4, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 4 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x4, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 4 keepdims normal')

    call check_close(sum(iqr(x5, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 5 axis normal')
    call check_close(iqr(x5, scale='normal'), iqr_normal_two, &
        'iqr rank 5 all normal')
    call check_close(sum(median_abs_deviation(x5, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 5 normal')
    call check_close(sum(iqr_keepdims(x5, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 5 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x5, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 5 keepdims normal')

    call check_close(sum(iqr(x6, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 6 axis normal')
    call check_close(iqr(x6, scale='normal'), iqr_normal_two, &
        'iqr rank 6 all normal')
    call check_close(sum(median_abs_deviation(x6, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 6 normal')
    call check_close(sum(iqr_keepdims(x6, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 6 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x6, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 6 keepdims normal')

    call check_close(sum(iqr(x7, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 7 axis normal')
    call check_close(iqr(x7, scale='normal'), iqr_normal_two, &
        'iqr rank 7 all normal')
    call check_close(sum(median_abs_deviation(x7, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 7 normal')
    call check_close(sum(iqr_keepdims(x7, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 7 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x7, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 7 keepdims normal')

    call check_close(sum(iqr(x8, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 8 axis normal')
    call check_close(iqr(x8, scale='normal'), iqr_normal_two, &
        'iqr rank 8 all normal')
    call check_close(sum(median_abs_deviation(x8, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 8 normal')
    call check_close(sum(iqr_keepdims(x8, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 8 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x8, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 8 keepdims normal')

    call check_close(sum(iqr(x9, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 9 axis normal')
    call check_close(iqr(x9, scale='normal'), iqr_normal_two, &
        'iqr rank 9 all normal')
    call check_close(sum(median_abs_deviation(x9, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 9 normal')
    call check_close(sum(iqr_keepdims(x9, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 9 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x9, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 9 keepdims normal')

    call check_close(sum(iqr(x10, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 10 axis normal')
    call check_close(iqr(x10, scale='normal'), iqr_normal_two, &
        'iqr rank 10 all normal')
    call check_close(sum(median_abs_deviation(x10, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 10 normal')
    call check_close(sum(iqr_keepdims(x10, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 10 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x10, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 10 keepdims normal')

    call check_close(sum(iqr(x11, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 11 axis normal')
    call check_close(iqr(x11, scale='normal'), iqr_normal_two, &
        'iqr rank 11 all normal')
    call check_close(sum(median_abs_deviation(x11, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 11 normal')
    call check_close(sum(iqr_keepdims(x11, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 11 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x11, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 11 keepdims normal')

    call check_close(sum(iqr(x12, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 12 axis normal')
    call check_close(iqr(x12, scale='normal'), iqr_normal_two, &
        'iqr rank 12 all normal')
    call check_close(sum(median_abs_deviation(x12, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 12 normal')
    call check_close(sum(iqr_keepdims(x12, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 12 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x12, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 12 keepdims normal')

    call check_close(sum(iqr(x13, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 13 axis normal')
    call check_close(iqr(x13, scale='normal'), iqr_normal_two, &
        'iqr rank 13 all normal')
    call check_close(sum(median_abs_deviation(x13, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 13 normal')
    call check_close(sum(iqr_keepdims(x13, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 13 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x13, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 13 keepdims normal')

    call check_close(sum(iqr(x14, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 14 axis normal')
    call check_close(iqr(x14, scale='normal'), iqr_normal_two, &
        'iqr rank 14 all normal')
    call check_close(sum(median_abs_deviation(x14, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 14 normal')
    call check_close(sum(iqr_keepdims(x14, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 14 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x14, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 14 keepdims normal')

    call check_close(sum(iqr(x15, axis=1, scale='normal')), iqr_normal_two, &
        'iqr rank 15 axis normal')
    call check_close(iqr(x15, scale='normal'), iqr_normal_two, &
        'iqr rank 15 all normal')
    call check_close(sum(median_abs_deviation(x15, axis=1, scale='normal')), &
        mad_normal_two, 'mad rank 15 normal')
    call check_close(sum(iqr_keepdims(x15, axis=1, scale='normal')), &
        iqr_normal_two, 'iqr rank 15 keepdims normal')
    call check_close(sum(median_abs_deviation_keepdims(x15, axis=1, &
        scale='normal')), mad_normal_two, 'mad rank 15 keepdims normal')

    call check_close(iqr(xnan, scale='normal', nan_policy='omit'), iqr_normal_two, &
        'iqr normal omit')
    call check_close(median_abs_deviation(xnan, scale='normal', nan_policy='omit'), &
        mad_normal_two, 'mad normal omit')
    call check_true(ieee_is_nan(iqr(xnan, scale='normal', nan_policy='propagate')), &
        'iqr normal propagate')
    call check_true(ieee_is_nan(median_abs_deviation(xnan, scale='normal', &
        nan_policy='propagate')), 'mad normal propagate')
    call check_true(ieee_is_nan(iqr(x1, scale='bad')), 'iqr invalid named scale')
    call check_true(ieee_is_nan(median_abs_deviation(x1, scale='bad')), &
        'mad invalid named scale')

    call check_close(iqr([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp], &
        rng=[10.0_dp, 90.0_dp], scale='normal'), 8.450832645481928_dp, &
        'iqr custom range normal')
    values = iqr(matrix, axis=1, scale='normal')
    call check_array(values, [0.7413011092528009_dp, 2.9652044370112036_dp, &
        6.671709983275209_dp, 11.860817748044814_dp], 'matrix iqr normal axis 1')
    values = median_abs_deviation(matrix, axis=2, scale='normal')
    call check_array(values, [5.930408874022408_dp, 11.860817748044816_dp, &
        17.791226622067225_dp], 'matrix mad normal axis 2')
    call check_close(iqr(matrix, scale='normal'), 12.231468302671214_dp, &
        'matrix iqr all normal')
    kept = iqr_keepdims(matrix, axis=2, scale='normal')
    call check_array(reshape(kept, [3]), [5.559758319396007_dp, &
        11.119516638792014_dp, 16.67927495818802_dp], &
        'matrix iqr keepdims normal axis 2')

    print *, 'descriptive named-scale tests passed'

contains

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy reference value
        character(len=*), intent(in) :: label !! test label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tol) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_array(actual, expected, label)
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        character(len=*), intent(in) :: label !! test label

        if (size(actual) /= size(expected)) then
            print *, 'FAIL size ', label
            error stop 1
        end if
        if (any(ieee_is_nan(actual))) then
            print *, 'FAIL NaN ', label
            error stop 1
        end if
        if (any(abs(actual - expected) > tol)) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_true(condition, label)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: label !! test label

        if (.not. condition) then
            print *, 'FAIL ', label
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_named_scale
