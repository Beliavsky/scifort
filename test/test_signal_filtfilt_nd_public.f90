! SPDX-License-Identifier: MIT
! Verify the public scifort_signal facade exports both assumed-rank routines.
program test_signal_filtfilt_nd_public
    use scifort_signal, only : dp, filtfilt, filtfilt_nd, sosfiltfilt_nd, &
        filtfilt_status_success
    implicit none
    real(dp) :: x(30), y(30), sos(1,6)
    integer :: code, k
    x = [(sin(real(k,dp)),k=1,size(x))]
    sos(1,:) = [0.5_dp,0.5_dp,0.0_dp,1.0_dp,-0.2_dp,0.0_dp]
    call filtfilt_nd([0.5_dp,0.5_dp],[1.0_dp],x,y,status=code)
    if (code /= filtfilt_status_success) error stop 'public filtfilt_nd missing'
    call filtfilt_nd([0.5_dp,0.5_dp],[1.0_dp],x,y, &
        method='gust',irlen=5,status=code)
    if (code /= filtfilt_status_success) error stop 'public gust filtfilt_nd missing'
    if (any(abs(y-filtfilt([0.5_dp,0.5_dp],[1.0_dp],x, &
        method='gust',irlen=5)) > 1.e-12_dp)) error stop 'public gust values'
    call sosfiltfilt_nd(sos,x,y,status=code)
    if (code /= filtfilt_status_success) error stop 'public sosfiltfilt_nd missing'
    print '(a)', 'test_signal_filtfilt_nd_public passed'
end program test_signal_filtfilt_nd_public
