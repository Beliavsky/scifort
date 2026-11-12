! SPDX-License-Identifier: MIT
! Tests for periodic waves and unit impulses.
program test_signal_periodic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_periodic, only : sawtooth, square, unit_impulse, &
        periodic_status_success, periodic_status_invalid_input
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: t(9) = 0.5_dp*pi*[-3.0_dp,-2.0_dp,-1.0_dp,0.0_dp, &
        0.5_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    real(dp), allocatable :: y(:), z(:,:)
    real(dp) :: nan
    integer :: status, i

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    call check(sawtooth(t), [-0.5_dp,0.0_dp,0.5_dp,-1.0_dp, &
        -0.75_dp,-0.5_dp,0.0_dp,0.5_dp,-1.0_dp])
    call check(sawtooth(t,width=0.0_dp), [0.5_dp,0.0_dp,-0.5_dp,1.0_dp, &
        0.75_dp,0.5_dp,0.0_dp,-0.5_dp,1.0_dp])
    call check(sawtooth(t,width=0.5_dp), [0.0_dp,1.0_dp,0.0_dp,-1.0_dp, &
        -0.5_dp,0.0_dp,1.0_dp,0.0_dp,-1.0_dp])
    call check(sawtooth([0.0_dp,pi,1.5_dp*pi,2.0_dp*pi], &
        width=[0.0_dp,0.5_dp,0.25_dp,1.0_dp]), [1.0_dp,1.0_dp,-1.0_dp/3.0_dp,-1.0_dp])
    call check(square(t), [1.0_dp,-1.0_dp,-1.0_dp,1.0_dp, &
        1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,1.0_dp])
    call check(square(t,duty=0.0_dp), spread(-1.0_dp,1,size(t)))
    call check(square(t,duty=1.0_dp), spread(1.0_dp,1,size(t)))
    call check(square([0.0_dp,0.5_dp*pi,pi,1.5_dp*pi], &
        duty=[0.0_dp,0.25_dp,0.5_dp,1.0_dp]), [-1.0_dp,-1.0_dp,-1.0_dp,1.0_dp])
    call check(sawtooth([real(dp) ::],status=status), [real(dp) ::])
    if (status /= periodic_status_success) error stop 'empty sawtooth status'
    call check(square([real(dp) ::],status=status), [real(dp) ::])
    if (status /= periodic_status_success) error stop 'empty square status'

    call reject(sawtooth(t,width=-0.1_dp,status=status),status)
    call reject(sawtooth(t,width=1.1_dp,status=status),status)
    call reject(sawtooth(t,width=nan,status=status),status)
    call reject(sawtooth(t,width=[0.5_dp],status=status),status)
    call reject(sawtooth(t,width=spread(1.1_dp,1,size(t)),status=status),status)
    call reject(square(t,duty=-0.1_dp,status=status),status)
    call reject(square(t,duty=nan,status=status),status)
    call reject(square(t,duty=[0.5_dp],status=status),status)
    call reject(square([nan],status=status),status)

    y = unit_impulse(8,status=status)
    if (status /= periodic_status_success) error stop 'unit impulse default status'
    call check(y,[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    y = unit_impulse(7,idx=3,status=status)
    if (status /= periodic_status_success) error stop 'unit impulse offset status'
    call check(y,[0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    call check(unit_impulse(8,centered=.true.), &
        [0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp])
    y = unit_impulse(0,status=status)
    if (size(y) /= 0 .or. status /= periodic_status_success) error stop 'empty impulse'
    y = unit_impulse(-3,status=status)
    if (size(y) /= 0 .or. status /= periodic_status_invalid_input) then
        error stop 'negative impulse size'
    end if
    y = unit_impulse(5,idx=0,status=status)
    call reject(y,status)
    y = unit_impulse(5,idx=6,status=status)
    call reject(y,status)

    z = unit_impulse([3,4], centered=.true., status=status)
    if (status /= periodic_status_success) error stop '2D impulse status'
    if (any(shape(z) /= [3,4])) error stop '2D impulse shape'
    if (abs(z(2,3)-1.0_dp) > 1.0e-15_dp .or. &
            abs(sum(z)-1.0_dp) > 1.0e-15_dp) error stop '2D impulse position'
    z = unit_impulse([2,5],idx=[2,4],status=status)
    if (status /= periodic_status_success) error stop '2D impulse explicit status'
    if (abs(z(2,4)-1.0_dp) > 1.0e-15_dp .or. &
            abs(sum(z)-1.0_dp) > 1.0e-15_dp) error stop '2D impulse explicit'
    z = unit_impulse([0,3],status=status)
    if (size(z) /= 0 .or. status /= periodic_status_success) error stop '2D empty impulse'
    z = unit_impulse([3,4],idx=[3,5],status=status)
    if (status /= periodic_status_invalid_input) error stop '2D bad index status'
    if (.not. all(ieee_is_nan(z))) error stop '2D bad index values'

    ! Property checks: exact period for reasonable phase magnitudes.
    do i = 1, 20
        y = [real(i,dp)/17.0_dp - pi]
        y = sawtooth(y)
        if (abs(y(1)) > 1.0_dp) error stop 'sawtooth bounded'
    end do

    print '(a)', 'test_signal_periodic passed'
contains
    subroutine check(actual, reference)
        real(dp), intent(in) :: actual(:) !! computed samples
        real(dp), intent(in) :: reference(:) !! trusted samples

        if (size(actual) /= size(reference)) error stop 'periodic size mismatch'
        if (any(abs(actual-reference) > 8.0e-14_dp)) then
            print *, 'actual   ',actual
            print *, 'expected ',reference
            error stop 'periodic value mismatch'
        end if
    end subroutine check

    subroutine reject(actual, ierr)
        real(dp), intent(in) :: actual(:) !! invalid result
        integer, intent(in) :: ierr !! nonzero status

        if (ierr /= periodic_status_invalid_input) error stop 'expected periodic status'
        if (.not. all(ieee_is_nan(actual))) error stop 'expected NaNs'
    end subroutine reject
end program test_signal_periodic
