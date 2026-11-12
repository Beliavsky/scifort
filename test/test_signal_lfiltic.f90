! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_lfiltic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfiltic, lfilter, &
        lfilter_status_success, lfilter_status_invalid_input
    implicit none

    real(dp), allocatable :: zi(:), outputs(:), continued(:), sequence(:)
    complex(dp), allocatable :: zic(:), outputsc(:), continuedc(:)
    real(dp), parameter :: a(3) = [2.0_dp, -0.7_dp, 0.14_dp]
    real(dp), parameter :: b(3) = [0.2_dp, 0.1_dp, 0.05_dp]
    complex(dp), parameter :: bc(3) = [cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(0.5_dp,-0.25_dp,dp), cmplx(0.125_dp,0.3_dp,dp)]
    complex(dp), parameter :: ac(2) = [cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(-0.4_dp,0.1_dp,dp)]
    complex(dp) :: samples(12)
    integer :: i, status

    zi = lfiltic([2.0_dp,4.0_dp,8.0_dp], [2.0_dp,-1.0_dp], &
        [3.0_dp,4.0_dp], [5.0_dp,6.0_dp], status)
    call close_real(zi,[35.5_dp,20.0_dp])
    if (status /= lfilter_status_success) error stop 'lfiltic real status'

    zi = lfiltic(b,a,[1.0_dp,2.0_dp],[3.0_dp,4.0_dp],status)
    call close_real(zi,[0.46_dp,0.005_dp])
    if (status /= lfilter_status_success) error stop 'lfiltic normalized real'

    zic = lfiltic(bc,ac,[cmplx(1.0_dp,2.0_dp,dp),cmplx(3.0_dp,4.0_dp,dp)], &
        [cmplx(2.0_dp,-1.0_dp,dp),cmplx(4.0_dp,3.0_dp,dp)],status)
    call close_complex(zic,[cmplx(0.475_dp,0.6375_dp,dp), &
        cmplx(0.275_dp,0.2375_dp,dp)])
    if (status /= lfilter_status_success) error stop 'lfiltic complex status'

    zi = lfiltic([1.0_dp, 2.0_dp, 3.0_dp],[1.0_dp, -0.5_dp], [6.0_dp])
    call close_real(zi, [3.0_dp, 0.0_dp])
    zi = lfiltic([1.0_dp],[1.0_dp, -0.5_dp, 0.25_dp], [6.0_dp])
    call close_real(zi, [3.0_dp, -1.5_dp])
    zi = lfiltic([1.0_dp,2.0_dp],[1.0_dp], [1.0_dp])
    call close_real(zi,[0.0_dp])
    zi = lfiltic([1.0_dp,2.0_dp],[1.0_dp], [1.0_dp], [3.0_dp])
    call close_real(zi,[6.0_dp])
    zi = lfiltic([1.0_dp],[2.0_dp], [1.0_dp], status=status)
    if (size(zi) /= 0 .or. status /= 0) error stop 'lfiltic order zero'

    sequence = [(0.2_dp*real(i,dp)+sin(real(i,dp)),i=1,25)]
    outputs = lfilter(b,a,sequence)
    zi = lfiltic(b,a,outputs(8:1:-1),sequence(8:1:-1))
    continued = lfilter(b,a,sequence(9:),zi=zi)
    call close_real(continued,outputs(9:))
    do i = 1, size(samples)
        samples(i) = cmplx(0.2_dp*real(i,dp),sin(real(i,dp)),dp)
    end do
    outputsc = lfilter(bc,ac,samples)
    zic = lfiltic(bc,ac,outputsc(7:1:-1),samples(7:1:-1))
    continuedc = lfilter(bc,ac,samples(8:),zi=zic)
    call close_complex(continuedc,outputsc(8:))

    zi = lfiltic([1.0_dp], [1.0_dp, 0.5_dp], [1.0_dp], &
        status=status)
    if (status /= lfilter_status_success) error stop 'lfiltic valid history'
    zi = lfiltic([1.0_dp], [0.0_dp, 0.5_dp], [1.0_dp], status=status)
    if (status /= lfilter_status_invalid_input) error stop 'lfiltic zero a0'
    if (.not. all(ieee_is_nan(zi))) error stop 'lfiltic invalid state'
    zi = lfiltic([1.0_dp,2.0_dp],[1.0_dp],[0.0_dp], &
        [ieee_value(0.0_dp,ieee_quiet_nan)],status)
    if (status /= lfilter_status_invalid_input) error stop 'lfiltic nan history'
    zic = lfiltic([cmplx(1.0_dp,0.0_dp,dp)], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(0.25_dp,0.0_dp,dp)], &
        [cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)],status=status)
    if (status /= lfilter_status_invalid_input) error stop 'lfiltic complex nan'
    print '(a)', 'test_signal_lfiltic passed'
contains
    subroutine close_real(actual, expected)
        real(dp), intent(in) :: actual(:) !! computed state or samples
        real(dp), intent(in) :: expected(:) !! expected state or samples
        if (size(actual) /= size(expected)) error stop 'lfiltic real size'
        if (any(abs(actual-expected) > 5.0e-13_dp)) error stop 'lfiltic real values'
    end subroutine close_real

    subroutine close_complex(actual, expected)
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: expected(:) !! expected complex values
        if (size(actual) /= size(expected)) error stop 'lfiltic complex size'
        if (any(abs(actual-expected) > 5.0e-13_dp)) error stop 'lfiltic complex values'
    end subroutine close_complex
end program test_signal_lfiltic
