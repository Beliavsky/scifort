! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_summary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : describe, describe_all, mode, mode_all, mode_keepdims
    implicit none

    real(dp), allocatable :: mode0
    integer, allocatable :: count0
    real(dp), allocatable :: modes1(:)
    integer, allocatable :: counts1(:)
    real(dp), allocatable :: modes2(:, :)
    integer, allocatable :: counts2(:, :)
    integer, allocatable :: nobs0
    real(dp), allocatable :: minimum0
    real(dp), allocatable :: maximum0
    real(dp), allocatable :: mean0
    real(dp), allocatable :: variance0
    real(dp), allocatable :: skewness0
    real(dp), allocatable :: kurtosis0
    integer, allocatable :: nobs1(:)
    real(dp), allocatable :: minimum1(:)
    real(dp), allocatable :: maximum1(:)
    real(dp), allocatable :: mean1(:)
    real(dp), allocatable :: variance1(:)
    real(dp), allocatable :: skewness1(:)
    real(dp), allocatable :: kurtosis1(:)
    real(dp) :: mode_scalar
    integer :: count_scalar
    integer :: nobs_scalar
    real(dp) :: minimum_scalar
    real(dp) :: maximum_scalar
    real(dp) :: mean_scalar
    real(dp) :: variance_scalar
    real(dp) :: skewness_scalar
    real(dp) :: kurtosis_scalar
    real(dp) :: qnan
    real(dp) :: pinf
    real(dp) :: matrix(3, 4)
    real(dp) :: masked_matrix(2, 3)
    real(dp) :: vector(4)
    real(dp) :: empty(0)
    real(dp) :: rank15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: modes15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    integer, allocatable :: counts15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    integer, allocatable :: nobs15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: minimum15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: maximum15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: mean15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: variance15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: skewness15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: kurtosis15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    pinf = ieee_value(0.0_dp, ieee_positive_inf)

    vector = [2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
    call mode(vector, mode0, count0)
    call check_close('mode vector tie', mode0, 1.0_dp)
    call check_int('mode vector count', count0, 2)

    call mode_all([qnan, qnan, 1.0_dp], mode_scalar, count_scalar, &
        nan_policy='propagate')
    call check_nan('mode NaN category', mode_scalar)
    call check_int('mode NaN category count', count_scalar, 2)
    call mode_all([qnan, qnan, 1.0_dp], mode_scalar, count_scalar, nan_policy='omit')
    call check_close('mode omit NaN', mode_scalar, 1.0_dp)
    call check_int('mode omit NaN count', count_scalar, 1)
    call mode_all(empty, mode_scalar, count_scalar)
    call check_nan('mode empty', mode_scalar)
    call check_int('mode empty count', count_scalar, 0)

    matrix = reshape([1.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, &
        5.0_dp, 5.0_dp, 6.0_dp, qnan, qnan, 7.0_dp], shape(matrix))
    call mode(matrix, modes1, counts1, axis=1, nan_policy='propagate')
    call check_array('mode matrix values', modes1(1:3), [1.0_dp, 2.0_dp, 5.0_dp])
    call check_nan('mode matrix NaN', modes1(4))
    call check_int_array('mode matrix counts', counts1, [2, 2, 2, 2])
    call mode_keepdims(matrix, modes2, counts2, axis=1, nan_policy='omit')
    call check_shape2('mode keepdims shape', shape(modes2), [1, 4])
    call check_shape2('mode keepdims count shape', shape(counts2), [1, 4])
    call check_close('mode keepdims omit', modes2(1, 4), 7.0_dp)
    call check_int('mode keepdims omit count', counts2(1, 4), 1)

    vector = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call describe(vector, nobs0, minimum0, maximum0, mean0, variance0, skewness0, kurtosis0)
    call check_int('describe vector nobs', nobs0, 4)
    call check_close('describe vector min', minimum0, 0.0_dp)
    call check_close('describe vector max', maximum0, 3.0_dp)
    call check_close('describe vector mean', mean0, 1.5_dp)
    call check_close('describe vector variance', variance0, 1.6666666666666665_dp)
    call check_close('describe vector skewness', skewness0, 0.0_dp)
    call check_close('describe vector kurtosis', kurtosis0, -1.36_dp)

    call describe_all(vector, nobs_scalar, minimum_scalar, maximum_scalar, mean_scalar, &
        variance_scalar, skewness_scalar, kurtosis_scalar, bias=.false.)
    call check_close('describe unbiased kurtosis', kurtosis_scalar, -1.2_dp)
    call describe_all(vector, nobs_scalar, minimum_scalar, maximum_scalar, mean_scalar, &
        variance_scalar, skewness_scalar, kurtosis_scalar, ddof=4)
    call check_true('describe ddof n infinity', .not. ieee_is_finite(variance_scalar) .and. &
        variance_scalar > 0.0_dp)
    call describe_all(vector, nobs_scalar, minimum_scalar, maximum_scalar, mean_scalar, &
        variance_scalar, skewness_scalar, kurtosis_scalar, ddof=5)
    call check_close('describe ddof above n', variance_scalar, -5.0_dp)

    masked_matrix = reshape([1.0_dp, 3.0_dp, qnan, 5.0_dp, pinf, 7.0_dp], &
        shape(masked_matrix))
    call describe(masked_matrix, nobs1, minimum1, maximum1, mean1, variance1, &
        skewness1, kurtosis1, axis=1, nan_policy='omit')
    call check_int_array('describe omit nobs', nobs1, [2, 1, 1])
    call check_array('describe omit min', minimum1, [1.0_dp, 5.0_dp, 7.0_dp])
    call check_array('describe omit max', maximum1, [3.0_dp, 5.0_dp, 7.0_dp])
    call check_array('describe omit mean', mean1, [2.0_dp, 5.0_dp, 7.0_dp])
    call check_close('describe omit variance first', variance1(1), 2.0_dp)
    call check_nan('describe omit variance second', variance1(2))
    call check_nan('describe omit variance third', variance1(3))
    call check_array('describe omit skewness', skewness1, [0.0_dp, 0.0_dp, 0.0_dp])
    call check_array('describe omit kurtosis', kurtosis1, [-2.0_dp, -3.0_dp, -3.0_dp])

    rank15 = reshape([1.0_dp, 4.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], shape(rank15))
    call mode(rank15, modes15, counts15, axis=2)
    call check_true('mode rank15 shape', size(modes15) == 2 .and. size(counts15) == 2)
    call check_array('mode rank15 values', reshape(modes15, [size(modes15)]), &
        [2.0_dp, 2.0_dp])
    call check_int_array('mode rank15 counts', reshape(counts15, [size(counts15)]), &
        [2, 1])
    call describe(rank15, nobs15, minimum15, maximum15, mean15, variance15, &
        skewness15, kurtosis15, axis=2)
    call check_true('describe rank15 shape', size(nobs15) == 2 .and. all(nobs15 == 3))
    call check_array('describe rank15 mean', reshape(mean15, [size(mean15)]), &
        [5.0_dp / 3.0_dp, 3.0_dp])

    print *, 'descriptive summary tests passed'

contains

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        integer :: i

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), size(actual), size(expected)
            error stop 1
        end if
        do i = 1, size(actual)
            if (ieee_is_nan(expected(i))) then
                if (.not. ieee_is_nan(actual(i))) then
                    print *, 'FAIL ', trim(label), i, actual(i), expected(i)
                    error stop 1
                end if
            else if (abs(actual(i) - expected(i)) > 2.0e-12_dp) then
                print *, 'FAIL ', trim(label), i, actual(i), expected(i)
                error stop 1
            end if
        end do
    end subroutine check_array

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (ieee_is_nan(expected)) then
            if (.not. ieee_is_nan(actual)) then
                print *, 'FAIL ', trim(label), actual, expected
                error stop 1
            end if
        else if (abs(actual - expected) > 2.0e-12_dp) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_int(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual
        integer, intent(in) :: expected

        if (actual /= expected) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_int

    subroutine check_int_array(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual(:)
        integer, intent(in) :: expected(:)

        if (size(actual) /= size(expected) .or. any(actual /= expected)) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_int_array

    subroutine check_nan(label, actual)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', trim(label), actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_shape2(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual(2)
        integer, intent(in) :: expected(2)

        if (any(actual /= expected)) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_shape2

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_summary
