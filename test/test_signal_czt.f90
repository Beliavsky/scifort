! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_czt
    use scifort_signal, only : dp, czt, zoom_fft, czt_status_success, czt_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), allocatable :: y(:), expected(:)
    complex(dp) :: z(4), w, a
    real(dp) :: x(4), fr(2), freq
    integer :: k, j, st

    x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    y = czt(x,status=st)
    if (st /= czt_status_success .or. size(y) /= 4) error stop 'default czt'
    if (maxval(abs(y-[cmplx(10.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,2.0_dp,dp),cmplx(-2.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,-2.0_dp,dp)])) > 1.0e-13_dp) error stop 'czt vs DFT'

    y = zoom_fft(x,[0.0_dp,4.0_dp],m=4,fs=4.0_dp,status=st)
    if (st /= czt_status_success) error stop 'zoom status'
    if (maxval(abs(y-[cmplx(10.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,2.0_dp,dp),cmplx(-2.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,-2.0_dp,dp)])) > 1.0e-13_dp) error stop 'zoom vs DFT'

    w = cmplx(0.8_dp,-0.1_dp,dp)
    a = cmplx(1.1_dp,0.2_dp,dp)
    z = cmplx(x,0.0_dp,dp)
    y = czt(z,m=7,w=w,a=a,status=st)
    if (st /= czt_status_success .or. size(y) /= 7) error stop 'nonunit czt'
    do k = 0,6
        block
            complex(dp) :: actual
            actual = cmplx(0.0_dp,0.0_dp,dp)
            do j = 0,3
                actual = actual+z(j+1)*w**(j*k)/a**j
            end do
            if (abs(y(k+1)-actual) > 1.0e-12_dp) error stop 'nonunit spiral'
        end block
    end do

    fr = [1.0_dp,3.0_dp]
    y = zoom_fft(z,fr,m=5,fs=8.0_dp,endpoint=.true.,status=st)
    if (st /= czt_status_success .or. size(y) /= 5) error stop 'zoom endpoints'
    do k = 0,4
        freq = fr(1)+(fr(2)-fr(1))*real(k,dp)/4.0_dp
        block
            complex(dp) :: actual
            actual = cmplx(0.0_dp,0.0_dp,dp)
            do j = 0,3
                actual = actual+z(j+1)*exp(cmplx(0.0_dp,-2.0_dp*pi*freq*real(j,dp)/8.0_dp,dp))
            end do
            if (abs(y(k+1)-actual) > 1.0e-12_dp) error stop 'zoom frequency grid'
        end block
    end do

    y = zoom_fft(x,4.0_dp,m=4,fs=4.0_dp,status=st)
    if (st /= czt_status_success) error stop 'scalar zoom'
    if (abs(y(1)-sum(x)) > 1.0e-12_dp) error stop 'scalar zoom first bin'
    y = czt([1.0_dp],m=5,status=st)
    if (st /= czt_status_success .or. any(abs(y-1.0_dp)>0.0_dp)) error stop 'one-sample czt'

    y = czt(x,m=0,status=st)
    if (st /= czt_status_invalid_input .or. size(y) /= 0) error stop 'invalid czt count'
    y = czt(x,w=cmplx(0.0_dp,0.0_dp,dp),status=st)
    if (st /= czt_status_invalid_input .or. size(y) /= 0) error stop 'invalid w'
    y = czt(x,a=cmplx(0.0_dp,0.0_dp,dp),status=st)
    if (st /= czt_status_invalid_input .or. size(y) /= 0) error stop 'invalid a'
    y = zoom_fft(x,[0.0_dp,1.0_dp],fs=0.0_dp,status=st)
    if (st /= czt_status_invalid_input .or. size(y) /= 0) error stop 'invalid fs'
    print '(a)', 'PASS: chirp Z-transform and zoom FFT'
end program test_signal_czt
