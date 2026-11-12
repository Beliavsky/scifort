! SPDX-License-Identifier: MIT
! Public signal facade must expose the caller-output FFT routine.
program test_signal_fftconvolve_assumed_rank_public
    use scifort_signal, only : dp, fftconvolve_nd, convolution_status_success
    implicit none
    real(dp) :: a(3,2), b(2,1), y(4,2)
    integer :: code
    a(:,1) = [1.0_dp,2.0_dp,3.0_dp]
    a(:,2) = [4.0_dp,5.0_dp,6.0_dp]
    b(:,1) = [1.0_dp,-0.5_dp]
    call fftconvolve_nd(a,b,y,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'public fftconvolve_nd status'
    if (maxval(abs(y(:,1)-[1.0_dp,1.5_dp,2.0_dp,-1.5_dp])) &
        > 1.0e-12_dp) error stop 'public fftconvolve_nd results'
    print '(a)', 'test_signal_fftconvolve_assumed_rank_public passed'
end program test_signal_fftconvolve_assumed_rank_public
