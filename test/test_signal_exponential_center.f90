! SPDX-License-Identifier: MIT
program test_signal_exponential_center
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, spectral_named_window, spectral_status_success
    use scifort_signal_windows_extra, only : window_samples
    implicit none
    integer :: st, i
    real(dp), allocatable :: w(:), w2(:)
    real(dp), parameter :: tau=2.25_dp, center=-1.25_dp
    w = get_window('exponential', 12, parameter=tau, parameter2=center, status=st)
    if (st/=spectral_status_success .or. size(w)/=12) error stop 'exponential window status'
    do i=1,12
        if (abs(w(i)-exp(-abs(real(i-1,dp)-center)/tau))>1.0e-14_dp) &
            error stop 'exponential window values'
    end do
    w2=spectral_named_window('exponential',12,st,parameter=tau,parameter2=center)
    if (st/=spectral_status_success) error stop 'spectral_named_window status'
    if (maxval(abs(w2-w))>1.0e-14_dp) error stop 'spectral_named_window values'
    w2=window_samples('poisson',12,parameter=tau,parameter2=center,status=st)
    if (st/=0) error stop 'poisson window status'
    if (maxval(abs(w2-w))>1.0e-14_dp) error stop 'poisson alias'
    w = get_window('exponential', 12, parameter=tau, parameter2=center, fftbins=.false., status=st)
    if (st==spectral_status_success .or. size(w)/=0) error stop 'symmetric center must fail'
    w=get_window('exponential',0,parameter=tau,parameter2=center,status=st)
    if (st/=spectral_status_success .or. size(w)/=0) error stop 'empty window'
    w=get_window('exponential',1,parameter=tau,parameter2=center,status=st)
    if (st/=spectral_status_success .or. size(w)/=1) error stop 'singleton window'
    if (abs(w(1)-1.0_dp)>1.0e-14_dp) error stop 'singleton is one'
    w=get_window('exponential',7,parameter=tau,status=st)
    if (st/=spectral_status_success) error stop 'old tau-only API'
    if (abs(w(4)-exp(-0.5_dp/tau))>1.0e-14_dp) error stop 'old centered taper'
    print *, 'PASS: explicit-center exponential/poisson windows'
end program test_signal_exponential_center
