! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive_trimmed
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : tmax, tmax_all, tmax_keepdims, tmean, tmean_all, &
        tmean_keepdims, tmin, tmin_all, tmin_keepdims, tsem, tsem_all, tsem_keepdims, &
        tstd, tstd_all, tstd_keepdims, tvar, tvar_all, tvar_keepdims
    implicit none

    real(dp) :: matrix(3, 4)
    real(dp) :: nan
    real(dp) :: nan_matrix(3, 4)
    real(dp) :: rank15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r1(:)
    real(dp), allocatable :: r2(:, :)
    real(dp), allocatable :: r14(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: x(20)
    integer :: i

    x = [(real(i - 1, dp), i=1,20)]
    nan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_close(tmean(x, lower_limit=3.0_dp, upper_limit=17.0_dp), 10.0_dp, &
        'tmean vector')
    call check_close(tvar(x, lower_limit=3.0_dp, upper_limit=17.0_dp), 20.0_dp, &
        'tvar vector')
    call check_close(tstd(x, lower_limit=3.0_dp, upper_limit=17.0_dp), &
        4.4721359549995796_dp, 'tstd vector')
    call check_close(tsem(x, lower_limit=3.0_dp, upper_limit=17.0_dp), &
        1.1547005383792515_dp, 'tsem vector')
    call check_close(tmean(x, lower_limit=3.0_dp, upper_limit=17.0_dp, &
        inclusive=[.true., .false.]), 9.5_dp, 'tmean upper exclusive')
    call check_close(tmean(x, lower_limit=3.0_dp, upper_limit=17.0_dp, &
        inclusive=[.false., .true.]), 10.5_dp, 'tmean lower exclusive')
    call check_close(tmin(x, lower_limit=13.0_dp), 13.0_dp, 'tmin inclusive')
    call check_close(tmin(x, lower_limit=13.0_dp, inclusive=.false.), 14.0_dp, &
        'tmin exclusive')
    call check_close(tmax(x, upper_limit=13.0_dp), 13.0_dp, 'tmax inclusive')
    call check_close(tmax(x, upper_limit=13.0_dp, inclusive=.false.), 12.0_dp, &
        'tmax exclusive')

    matrix = reshape([(real(i, dp), i=1,12)], shape(matrix))
    r1 = tmean(matrix, axis=1, lower_limit=3.0_dp, upper_limit=10.0_dp)
    call check_vector(r1, [3.0_dp, 5.0_dp, 8.0_dp, 10.0_dp], 'tmean matrix axis1')
    r1 = tmean(matrix, axis=2, lower_limit=3.0_dp, upper_limit=10.0_dp)
    call check_vector(r1, [7.0_dp, 6.5_dp, 6.0_dp], 'tmean matrix axis2')
    r1 = tvar(matrix, axis=1, lower_limit=3.0_dp, upper_limit=10.0_dp)
    call check_vector_nan(r1, [nan, 1.0_dp, 1.0_dp, nan], 'tvar matrix axis1')
    r1 = tsem(matrix, axis=2, lower_limit=3.0_dp, upper_limit=10.0_dp)
    call check_vector(r1, [1.7320508075688772_dp, 1.5_dp, 1.7320508075688772_dp], &
        'tsem matrix axis2')
    r1 = tmin(matrix, axis=2, lower_limit=4.0_dp)
    call check_vector(r1, [4.0_dp, 5.0_dp, 6.0_dp], 'tmin matrix axis2')
    r1 = tmax(matrix, axis=2, upper_limit=9.0_dp)
    call check_vector(r1, [7.0_dp, 8.0_dp, 9.0_dp], 'tmax matrix axis2')

    call check_close(tmean(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp), &
        6.5_dp, 'tmean matrix default all')
    call check_close(tmean_all(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp), &
        6.5_dp, 'tmean_all matrix')
    call check_close(tvar_all(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp), &
        6.0_dp, 'tvar_all matrix')
    call check_close(tstd_all(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp), &
        2.449489742783178_dp, 'tstd_all matrix')
    call check_close(tsem_all(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp), &
        0.8660254037844386_dp, 'tsem_all matrix')
    call check_close(tmin_all(matrix, lower_limit=4.0_dp), 4.0_dp, 'tmin_all matrix')
    call check_close(tmax_all(matrix, upper_limit=9.0_dp), 9.0_dp, 'tmax_all matrix')

    r2 = tmean_keepdims(matrix, lower_limit=3.0_dp, upper_limit=10.0_dp)
    if (any(shape(r2) /= [1, 1])) error stop 'tmean keepdims all shape'
    call check_close(r2(1, 1), 6.5_dp, 'tmean keepdims all')
    r2 = tvar_keepdims(matrix, axis=1, lower_limit=3.0_dp, upper_limit=10.0_dp)
    if (any(shape(r2) /= [1, 4])) error stop 'tvar keepdims shape'
    call check_vector_nan(reshape(r2, [size(r2)]), [nan, 1.0_dp, 1.0_dp, nan], &
        'tvar keepdims values')
    r2 = tmin_keepdims(matrix, axis=2, lower_limit=4.0_dp)
    if (any(shape(r2) /= [3, 1])) error stop 'tmin keepdims shape'
    call check_vector(reshape(r2, [size(r2)]), [4.0_dp, 5.0_dp, 6.0_dp], &
        'tmin keepdims values')
    r2 = tmax_keepdims(matrix, axis=2, upper_limit=9.0_dp)
    if (any(shape(r2) /= [3, 1])) error stop 'tmax keepdims shape'

    nan_matrix = matrix
    nan_matrix(2, 2) = nan
    r1 = tmean(nan_matrix, axis=1, lower_limit=3.0_dp, upper_limit=10.0_dp, &
        nan_policy='propagate')
    call check_vector_nan(r1, [3.0_dp, nan, 8.0_dp, 10.0_dp], &
        'tmean nan propagate')
    r1 = tmean(nan_matrix, axis=1, lower_limit=3.0_dp, upper_limit=10.0_dp, &
        nan_policy='omit')
    call check_vector(r1, [3.0_dp, 5.0_dp, 8.0_dp, 10.0_dp], 'tmean nan omit')

    rank15 = 0.0_dp
    rank15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = 1.0_dp
    rank15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = 3.0_dp
    r14 = tvar(rank15, axis=1)
    call check_close(sum(r14), 2.0_dp, 'tvar rank15')
    call check_close(tmean(rank15), 2.0_dp, 'tmean rank15 all')

    if (.not. ieee_is_nan(tmean([1.0_dp, 2.0_dp], lower_limit=10.0_dp))) then
        error stop 'empty trimmed slice should return NaN'
    end if

    print *, 'descriptive trimmed statistics tests passed'

contains

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        character(len=*), intent(in) :: label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > 5.0e-13_dp) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_vector(actual, expected, label)
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        character(len=*), intent(in) :: label
        integer :: i

        if (size(actual) /= size(expected)) error stop 'vector size mismatch'
        do i = 1, size(actual)
            call check_close(actual(i), expected(i), label)
        end do
    end subroutine check_vector

    subroutine check_vector_nan(actual, expected, label)
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        character(len=*), intent(in) :: label
        integer :: i

        if (size(actual) /= size(expected)) error stop 'vector size mismatch'
        do i = 1, size(actual)
            if (ieee_is_nan(expected(i))) then
                if (.not. ieee_is_nan(actual(i))) then
                    print *, 'FAIL ', label, ' at ', i
                    error stop 1
                end if
            else
                call check_close(actual(i), expected(i), label)
            end if
        end do
    end subroutine check_vector_nan

end program test_descriptive_trimmed
