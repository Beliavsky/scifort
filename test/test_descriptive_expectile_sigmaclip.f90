! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_expectile_sigmaclip
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : expectile, sigmaclip
    implicit none

    real(dp) :: x(4)
    real(dp) :: matrix(3, 4)
    real(dp) :: weights4(4)
    real(dp) :: weights31(3, 1)
    real(dp) :: qnan
    real(dp), allocatable :: clipped(:)
    real(dp) :: lower
    real(dp) :: upper
    real(dp) :: clip_example(36)
    real(dp) :: rank15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    integer :: i

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, 4.0_dp, 2.0_dp, -1.0_dp]
    call check_close('expectile alpha 0', expectile(x, 0.0_dp), -1.0_dp)
    call check_close('expectile alpha 0.2', expectile(x, 0.2_dp), &
        0.4285714285714286_dp)
    call check_close('expectile default', expectile(x), 1.5_dp)
    call check_close('expectile alpha 0.8', expectile(x, 0.8_dp), &
        2.5714285714285716_dp)
    call check_close('expectile alpha 1', expectile(x, 1.0_dp), 4.0_dp)
    call check_close('expectile weighted', &
        expectile(x, 0.8_dp, [1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp]), &
        3.3333333333333335_dp)

    matrix = reshape([(real(i, dp), i = 1, 12)], shape(matrix))
    weights4 = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    weights31(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp]
    call check_close('expectile matrix', expectile(matrix, 0.3_dp), 5.25_dp)
    call check_close('expectile matrix vector weights', &
        expectile(matrix, 0.3_dp, weights4), 6.857142857142858_dp)
    call check_close('expectile matrix singleton weights', &
        expectile(matrix, 0.3_dp, reshape(weights4, [1, 4])), &
        6.857142857142858_dp)
    call check_close('expectile matrix column weights', &
        expectile(matrix, 0.3_dp, weights31), 5.693548387096775_dp)
    call check_close('expectile scalar weight', expectile(matrix, 0.3_dp, 2.0_dp), &
        5.25_dp)
    call check_close('expectile negative weights', &
        expectile([1.0_dp, 2.0_dp, 3.0_dp], 0.2_dp, &
            [1.0_dp, -1.0_dp, 1.0_dp]), 1.25_dp)
    call check_nan('expectile invalid alpha', expectile(x, -0.1_dp))
    call check_nan('expectile invalid broadcast', expectile(matrix, 0.3_dp, [1.0_dp, 2.0_dp]))
    call check_nan('expectile zero weight sum', &
        expectile(x, 0.5_dp, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]))
    call check_nan('expectile NaN input', expectile([1.0_dp, qnan, 3.0_dp], 0.3_dp))

    rank15 = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], shape(rank15))
    call check_close('expectile rank15', expectile(rank15, 0.25_dp), &
        2.7_dp)

    call sigmaclip([1.0_dp, 2.0_dp, 3.0_dp, 100.0_dp], clipped, lower, upper, &
        low=1.5_dp, high=1.5_dp)
    call check_array('sigmaclip simple values', clipped, [1.0_dp, 2.0_dp, 3.0_dp])
    call check_close('sigmaclip simple lower', lower, 0.7752551286084111_dp)
    call check_close('sigmaclip simple upper', upper, 3.224744871391589_dp)

    do i = 1, 31
        clip_example(i) = 9.5_dp + real(i - 1, dp) / 30.0_dp
    end do
    clip_example(32:36) = [0.0_dp, 5.0_dp, 10.0_dp, 15.0_dp, 20.0_dp]
    call sigmaclip(clip_example, clipped, lower, upper, low=1.5_dp, high=1.5_dp)
    call check_array('sigmaclip SciPy example', clipped, &
        [9.966666666666667_dp, 10.0_dp, 10.033333333333333_dp, 10.0_dp])
    call check_close('sigmaclip example lower', lower, 9.964644660940673_dp)
    call check_close('sigmaclip example upper', upper, 10.035355339059327_dp)

    call sigmaclip(reshape([1.0_dp, 2.0_dp, 3.0_dp, 100.0_dp], [2, 2]), &
        clipped, lower, upper, low=1.5_dp, high=1.5_dp)
    call check_array('sigmaclip matrix flatten', clipped, [1.0_dp, 2.0_dp, 3.0_dp])

    call sigmaclip([1.0_dp, 2.0_dp, 3.0_dp, 100.0_dp], clipped, lower, upper, &
        low=0.0_dp, high=0.0_dp)
    call check_true('sigmaclip zero factors empty', size(clipped) == 0)
    call check_nan('sigmaclip zero factors lower', lower)
    call check_nan('sigmaclip zero factors upper', upper)

    call sigmaclip([1.0_dp, qnan, 3.0_dp], clipped, lower, upper)
    call check_true('sigmaclip NaN empty', size(clipped) == 0)
    call check_nan('sigmaclip NaN lower', lower)
    call check_nan('sigmaclip NaN upper', upper)

    print *, 'descriptive expectile/sigmaclip tests passed'

contains

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), size(actual), size(expected)
            error stop 1
        end if
        if (any(abs(actual - expected) > 3.0e-12_dp)) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (abs(actual - expected) > 3.0e-12_dp * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(label, actual)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', trim(label), actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_expectile_sigmaclip
