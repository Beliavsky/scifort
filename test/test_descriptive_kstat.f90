! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive_kstat
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_stats, only : kstat, kstat_all, kstat_keepdims, kstatvar, &
        kstatvar_all, kstatvar_keepdims
    implicit none

    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp) :: x(6)
    real(dp) :: xnan(5)
    real(dp) :: x2(3, 4)
    real(dp), allocatable :: a1(:)
    real(dp), allocatable :: a2(:)
    real(dp), allocatable :: kkeep(:, :)
    real(dp), allocatable :: kvkeep(:, :)
    real(dp), allocatable :: keep15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: x15(4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: qnan
    integer :: i

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    x2 = reshape([(real(i, dp), i=1,12)], shape(x2))

    call check_scalar('kstat n1', kstat(x, n=1), 3.5_dp)
    call check_scalar('kstat n2', kstat(x, n=2), 3.5_dp)
    call check_scalar('kstat default', kstat(x), 3.5_dp)
    call check_scalar('kstat n3', kstat(x, n=3), 0.0_dp)
    call check_scalar('kstat n4', kstat(x, n=4), -14.7_dp)
    call check_scalar('kstat positional n4', kstat(x, 4), -14.7_dp)
    call check_scalar('kstat_all n4', kstat_all(x, n=4), -14.7_dp)
    call check_scalar('kstatvar n1', kstatvar(x, n=1), 7.0_dp / 12.0_dp)
    call check_scalar('kstatvar n2', kstatvar(x, n=2), 1.75_dp)
    call check_scalar('kstatvar default', kstatvar(x), 1.75_dp)
    call check_scalar('kstatvar_all n1', kstatvar_all(x, n=1), 7.0_dp / 12.0_dp)

    a1 = kstat(x2, axis=1)
    call check_array('kstat default matrix axis1', a1, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
    a1 = kstat(x2, n=2, axis=1)
    call check_array('kstat matrix axis1', a1, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
    a2 = kstat(x2, n=2, axis=2)
    call check_array('kstat matrix axis2', a2, [15.0_dp, 15.0_dp, 15.0_dp])
    a1 = kstat(x2, n=3, axis=1)
    call check_array('kstat n3 matrix axis1', a1, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    a2 = kstat(x2, n=4, axis=2)
    call check_array('kstat n4 matrix axis2', a2, [-270.0_dp, -270.0_dp, -270.0_dp])
    a1 = kstatvar(x2, axis=2)
    call check_array('kstatvar default matrix axis2', a1, [49.5_dp, 49.5_dp, 49.5_dp])
    a1 = kstatvar(x2, n=1, axis=1)
    call check_array('kstatvar n1 matrix axis1', a1, &
        [1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp])
    a2 = kstatvar(x2, n=2, axis=2)
    call check_array('kstatvar n2 matrix axis2', a2, [49.5_dp, 49.5_dp, 49.5_dp])

    a1 = kstat(x2, n=4, axis=1)
    if (.not. all(ieee_is_nan(a1))) then
        print *, 'FAIL kstat insufficient axis slice'
        error stop 1
    end if
    a1 = kstatvar(x2, n=2, axis=1)
    if (.not. all(ieee_is_nan(a1))) then
        print *, 'FAIL kstatvar insufficient axis slice'
        error stop 1
    end if

    kkeep = kstat_keepdims(x2, n=2, axis=1)
    if (any(shape(kkeep) /= [1, 4])) then
        print *, 'FAIL kstat keepdims shape', shape(kkeep)
        error stop 1
    end if
    call check_array('kstat keepdims values', reshape(kkeep, [size(kkeep)]), &
        [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
    kvkeep = kstatvar_keepdims(x2, n=1, axis=2)
    if (any(shape(kvkeep) /= [3, 1])) then
        print *, 'FAIL kstatvar keepdims shape', shape(kvkeep)
        error stop 1
    end if
    call check_array('kstatvar keepdims values', reshape(kvkeep, [size(kvkeep)]), &
        [3.75_dp, 3.75_dp, 3.75_dp])

    xnan = [1.0_dp, 2.0_dp, qnan, 4.0_dp, 5.0_dp]
    call check_scalar('kstat omit n1', kstat(xnan, n=1, nan_policy='omit'), 3.0_dp)
    call check_scalar('kstat omit n2', kstat(xnan, n=2, nan_policy='omit'), &
        10.0_dp / 3.0_dp)
    call check_scalar('kstat omit n4', kstat(xnan, n=4, nan_policy='omit'), &
        -110.0_dp / 3.0_dp)
    call check_scalar('kstatvar omit n1', kstatvar(xnan, n=1, nan_policy='omit'), &
        5.0_dp / 6.0_dp)
    call check_scalar('kstatvar omit n2', kstatvar(xnan, n=2, nan_policy='omit'), &
        -19.0_dp / 18.0_dp)
    if (.not. ieee_is_nan(kstat(xnan, n=2, nan_policy='propagate'))) then
        print *, 'FAIL kstat propagate'
        error stop 1
    end if
    if (.not. ieee_is_nan(kstatvar(xnan, n=1, nan_policy='propagate'))) then
        print *, 'FAIL kstatvar propagate'
        error stop 1
    end if
    if (.not. ieee_is_nan(kstat(xnan, n=2, nan_policy='unsupported'))) then
        print *, 'FAIL kstat unsupported policy'
        error stop 1
    end if

    if (.not. ieee_is_nan(kstat([real(dp) ::], n=1))) then
        print *, 'FAIL kstat empty sample'
        error stop 1
    end if
    if (.not. ieee_is_nan(kstatvar([1.0_dp], n=1))) then
        print *, 'FAIL kstatvar insufficient sample'
        error stop 1
    end if

    x15 = 0.0_dp
    x15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = &
        [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call check_scalar('kstat rank15 all', kstat_all(x15, n=4), -10.0_dp / 3.0_dp)
    keep15 = kstat_keepdims(x15, n=2, axis=1)
    if (any(shape(keep15) /= 1)) then
        print *, 'FAIL kstat rank15 keepdims shape', shape(keep15)
        error stop 1
    end if
    call check_scalar('kstat rank15 keepdims', &
        keep15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5.0_dp / 3.0_dp)

    print *, 'descriptive kstat tests passed'

contains

    subroutine check_scalar(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value

        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_scalar

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

end program test_descriptive_kstat
