! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_dpss_public
    use scifort_signal
    implicit none
    type(dpss_result) :: taper
    type(spectral_result) :: estimate
    real(dp) :: x(32)
    integer :: i
    taper = dpss(32,2.5_dp,kmax=4,sym=.false.)
    if (taper%status /= dpss_status_success) error stop 'public dpss'
    if (any(shape(taper%tapers) /= [32,4])) error stop 'public dpss dimensions'
    x = [(sin(2.0_dp*acos(-1.0_dp)*real(i-1,dp)*5.0_dp/32.0_dp),i=1,32)]
    estimate = multitaper_psd(x,2.5_dp,kmax=4,fs=32.0_dp)
    if (estimate%status /= spectral_status_success) error stop 'public multitaper'
    if (size(estimate%power) /= 17) error stop 'public multitaper dimensions'
    if (estimate%power(6) <= estimate%power(1)) error stop 'public multitaper peak'
    print '(a)', 'Public DPSS and multitaper API passed'
end program test_signal_dpss_public
