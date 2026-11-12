! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_group_delay
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_signal_design, only : group_delay, group_delay_result, &
        design_status_success, design_status_invalid_input, design_status_singular
    use scifort_kinds, only : dp
    implicit none
    type(group_delay_result) :: gd
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: eps = 8.0e-12_dp
    real(dp), parameter :: w(5) = [0.0_dp,0.3_dp,1.0_dp,2.0_dp,3.0_dp]
    complex(dp), parameter :: bc(2) = [cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.2_dp,0.3_dp,dp)]
    complex(dp), parameter :: ac(2) = [cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(-0.3_dp,0.1_dp,dp)]
    real(dp), allocatable :: delay_old(:)

    gd = group_delay([0.0_dp,1.0_dp],worn=15)
    if (gd%status /= design_status_success) error stop 'pure delay status'
    if (size(gd%delay) /= 15) error stop 'pure delay grid'
    if (any(abs(gd%delay-1.0_dp)>eps)) error stop 'pure delay values'
    if (abs(gd%frequencies(15)-14.0_dp*pi/15.0_dp)>eps) &
        error stop 'pure delay grid range'
    gd = group_delay([0.0_dp,0.0_dp,1.0_dp],whole=.true.,worn=24)
    if (gd%status /= design_status_success) error stop 'whole circle'
    if (any(abs(gd%delay-2.0_dp)>eps)) error stop 'two sample delay'
    if (abs(gd%frequencies(24)-2.0_dp*pi*23.0_dp/24.0_dp)>eps) &
        error stop 'whole circle frequencies'

    gd = group_delay([1.0_dp],[1.0_dp,-0.5_dp],w=w)
    if (gd%status /= design_status_success) error stop 'IIR status'
    call close(gd%frequencies,w)
    call close(gd%delay,[1.0_dp,0.77263806396392_dp,0.02839399521894_dp, &
        -0.27492980103894_dp,-0.33258872493139_dp])
    gd = group_delay([1.0_dp,2.0_dp,-0.5_dp],[1.0_dp,-0.4_dp,0.13_dp],w=w)
    if (gd%status /= design_status_success) error stop 'IIR second order'
    call close(gd%delay,[0.59178082191781_dp,0.66182961533365_dp, &
        1.01876963242983_dp,0.91375720103445_dp,1.53297356479_dp])

    gd = group_delay(bc,ac,w=w)
    if (gd%status /= design_status_success) error stop 'complex delay status'
    call close(gd%delay,[0.6156862745098_dp,0.51056559456428_dp, &
        0.24163029659276_dp,0.00557294048472_dp,-0.27005564320005_dp])
    delay_old = gd%delay
    gd = group_delay(1.0e100_dp*bc,1.0e-100_dp*ac,w=w)
    if (gd%status /= design_status_success) error stop 'coefficient scaling status'
    call close(gd%delay,delay_old)
    gd = group_delay([1.0_dp],[1.0_dp,-0.5_dp], &
        fs=200.0_dp,w=[2.0_dp,12.0_dp,50.0_dp,90.0_dp])
    if (gd%status /= design_status_success) error stop 'sample rate status'
    call close(gd%delay,[0.9882530899443_dp,0.67105703820718_dp, &
        -0.2_dp,-0.32962727343721_dp])
    if (abs(gd%frequencies(4)-90.0_dp)>eps) error stop 'sample rate units'

    gd = group_delay([1.0_dp,1.0_dp],w=[pi])
    if (gd%status /= design_status_singular) error stop 'singular FIR status'
    if (.not. ieee_is_nan(gd%delay(1))) error stop 'singular FIR output'
    gd = group_delay([1.0_dp],[1.0_dp,-1.0_dp],w=[0.0_dp,0.3_dp])
    if (gd%status /= design_status_singular) error stop 'singular IIR status'
    if (.not. ieee_is_nan(gd%delay(1))) error stop 'singular IIR output'
    if (.not. abs(gd%delay(2)) < 5.0_dp) error stop 'nonsingular IIR output'

    gd = group_delay([0.0_dp,0.0_dp])
    if (gd%status /= design_status_invalid_input) error stop 'zero numerator'
    if (allocated(gd%delay)) error stop 'zero numerator output'
    gd = group_delay([1.0_dp],[0.0_dp])
    if (gd%status /= design_status_invalid_input) error stop 'zero denominator'
    gd = group_delay([1.0_dp],worn=0)
    if (gd%status /= design_status_invalid_input) error stop 'zero grid'
    gd = group_delay([1.0_dp],fs=0.0_dp)
    if (gd%status /= design_status_invalid_input) error stop 'bad fs'
    gd = group_delay([1.0_dp],w=[ieee_value(0.0_dp,ieee_quiet_nan)])
    if (gd%status /= design_status_invalid_input) error stop 'nan grid'
    print '(a)', 'test_signal_group_delay passed'
contains
    subroutine close(actual,expected)
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy or analytical reference
        if (size(actual) /= size(expected)) error stop 'group delay shape'
        if (any(abs(actual-expected)>eps)) error stop 'group delay values'
    end subroutine close
end program test_signal_group_delay
