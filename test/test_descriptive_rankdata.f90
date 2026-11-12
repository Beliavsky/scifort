! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py. Whole-array multidimensional cases
! use Fortran-order flattening, the natural analogue of NumPy raveling.

program test_descriptive_rankdata
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : rankdata, tiecorrect
    implicit none

    real(dp) :: matrix(2, 3)
    real(dp) :: qnan
    real(dp), allocatable :: flat_ranks(:)
    real(dp), allocatable :: ranked(:, :)
    real(dp) :: x15(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    matrix = reshape([3.0_dp, 2.0_dp, 1.0_dp, qnan, 1.0_dp, 2.0_dp], shape(matrix))

    flat_ranks = rankdata(matrix)
    call check_true('whole propagate', all(ieee_is_nan(flat_ranks)))
    flat_ranks = rankdata(matrix, nan_policy='omit')
    call check_array('whole average omit', flat_ranks, &
        [5.0_dp, 3.5_dp, 1.5_dp, qnan, 1.5_dp, 3.5_dp])
    flat_ranks = rankdata(matrix, 'min', nan_policy='omit')
    call check_array('whole min omit', flat_ranks, &
        [5.0_dp, 3.0_dp, 1.0_dp, qnan, 1.0_dp, 3.0_dp])
    flat_ranks = rankdata(matrix, 'max', nan_policy='omit')
    call check_array('whole max omit', flat_ranks, &
        [5.0_dp, 4.0_dp, 2.0_dp, qnan, 2.0_dp, 4.0_dp])
    flat_ranks = rankdata(matrix, 'dense', nan_policy='omit')
    call check_array('whole dense omit', flat_ranks, &
        [3.0_dp, 2.0_dp, 1.0_dp, qnan, 1.0_dp, 2.0_dp])
    flat_ranks = rankdata(matrix, 'ordinal', nan_policy='omit')
    call check_array('whole ordinal omit', flat_ranks, &
        [5.0_dp, 3.0_dp, 1.0_dp, qnan, 2.0_dp, 4.0_dp])

    ranked = rankdata(matrix, axis=1, nan_policy='omit')
    call check_array('axis1 omit', reshape(ranked, [size(ranked)]), &
        [2.0_dp, 1.0_dp, 1.0_dp, qnan, 1.0_dp, 2.0_dp])
    ranked = rankdata(matrix, axis=2, nan_policy='omit')
    call check_array('axis2 omit', reshape(ranked, [size(ranked)]), &
        [3.0_dp, 1.5_dp, 1.5_dp, qnan, 1.5_dp, 1.5_dp])
    ranked = rankdata(matrix, axis=2, method='ordinal', nan_policy='omit')
    call check_array('axis2 ordinal omit', reshape(ranked, [size(ranked)]), &
        [3.0_dp, 1.0_dp, 1.0_dp, qnan, 2.0_dp, 2.0_dp])

    flat_ranks = rankdata([1.0_dp, qnan, 2.0_dp], nan_policy='omit')
    call check_array('vector omit positions', flat_ranks, [1.0_dp, qnan, 2.0_dp])
    flat_ranks = rankdata([real(dp) ::], nan_policy='omit')
    if (size(flat_ranks) /= 0) then
        print *, 'FAIL empty rankdata'
        error stop 1
    end if
    flat_ranks = rankdata([2.0_dp, 1.0_dp], 'bad')
    call check_true('legacy invalid method NaNs', all(ieee_is_nan(flat_ranks)))

    call check_close('tiecorrect example', tiecorrect([1.0_dp, 2.5_dp, 2.5_dp, 4.0_dp]), 0.9_dp)
    call check_close('tiecorrect repeated ranks', &
        tiecorrect(rankdata([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 7.0_dp, 2.0_dp, &
            8.0_dp, 4.0_dp])), 0.983333333333333333_dp)
    call check_close('tiecorrect empty', tiecorrect([real(dp) ::]), 1.0_dp)
    call check_close('tiecorrect nan distinct', tiecorrect([1.0_dp, qnan, 1.0_dp]), 0.75_dp)

    x15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = [4.0_dp, 1.0_dp, 1.0_dp]
    r15 = rankdata(x15, axis=1)
    call check_array('rank15 axis', reshape(r15, [size(r15)]), [3.0_dp, 1.5_dp, 1.5_dp])
    flat_ranks = rankdata(x15)
    call check_array('rank15 all', flat_ranks, [3.0_dp, 1.5_dp, 1.5_dp])

    print *, 'descriptive rankdata tests passed'

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
            else if (abs(actual(i) - expected(i)) > 5.0e-13_dp) then
                print *, 'FAIL ', trim(label), i, actual(i), expected(i)
                error stop 1
            end if
        end do
    end subroutine check_array

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (abs(actual - expected) > 5.0e-13_dp) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_rankdata
