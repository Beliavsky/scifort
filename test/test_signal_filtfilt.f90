! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_filtfilt
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, filtfilt, sosfiltfilt, &
        filtfilt_status_success, filtfilt_status_invalid_input, filtfilt_status_no_steady_state
    implicit none
    real(dp) :: b(2), a(2), sos(1,6), x(72)
    complex(dp) :: cb(2), ca(2), csos(1,6), cx(72)
    real(dp), allocatable :: y(:), y2(:)
    complex(dp), allocatable :: cy(:), cy2(:)
    character(len=8) :: mode(4)
    integer :: i, j, status
    b = [0.25_dp,0.25_dp]
    a = [1.0_dp,-0.5_dp]
    sos(1,:) = [b(1),b(2),0.0_dp,a(1),a(2),0.0_dp]
    x = [(sin(0.13_dp*i)+0.04_dp*i, i=1,72)]
    cb = cmplx(b,0.0_dp,dp)
    ca = cmplx(a,0.0_dp,dp)
    cb(2) = cmplx(0.25_dp,0.05_dp,dp)
    ca(2) = cmplx(-0.5_dp,0.03_dp,dp)
    csos(1,:) = [cb(1),cb(2),cmplx(0.0_dp,0.0_dp,dp), &
        ca(1),ca(2),cmplx(0.0_dp,0.0_dp,dp)]
    cx = cmplx(x,0.1_dp*sin(0.4_dp*x),dp)
    mode = ['odd     ', 'even    ', 'constant', 'none    ']
    do j=1,size(mode)
        y = filtfilt(b,a,x,padtype=trim(mode(j)),status=status)
        if (status /= filtfilt_status_success) error stop 'filtfilt status'
        y2 = sosfiltfilt(sos,x,padtype=trim(mode(j)),status=status)
        if (status /= filtfilt_status_success) error stop 'sosfiltfilt status'
        if (any(abs(y-y2)>1.0e-13_dp)) error stop 'sosfiltfilt vs filtfilt'
        cy = filtfilt(cb,ca,cx,padtype=trim(mode(j)),status=status)
        if (status /= filtfilt_status_success) error stop 'complex filtfilt status'
        cy2 = sosfiltfilt(csos,cx,padtype=trim(mode(j)),status=status)
        if (status /= filtfilt_status_success) error stop 'complex sosfiltfilt status'
        if (any(abs(cy-cy2)>1.0e-13_dp)) error stop 'complex sos vs filtfilt'
    end do
    y = filtfilt(b,a,x,padlen=0,status=status)
    if (status /= filtfilt_status_success) error stop 'no pad'
    y2 = filtfilt(b,a,x,padtype='none',status=status)
    if (any(abs(y-y2)>1.0e-13_dp)) error stop 'explicit unpadded'
    y = sosfiltfilt(sos,x,padlen=9,status=status)
    if (status /= filtfilt_status_success) error stop 'override padding'
    y = filtfilt(b,a,x(:4),status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'short input filtfilt'
    y = sosfiltfilt(sos,x(:4),status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'short input sosfiltfilt'
    y = filtfilt(b,a,x,padtype='bogus',status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'invalid padtype'
    y = sosfiltfilt(sos,x,padlen=-1,status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'invalid padlen'
    y = filtfilt([1.0_dp],[1.0_dp,-1.0_dp],x,status=status)
    if (status /= filtfilt_status_no_steady_state) error stop 'no steady state'
    if (.not. all(ieee_is_nan(y))) error stop 'NaN on failure'
    print '(a)', 'test_signal_filtfilt passed'
end program
