! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_windows_extra
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, spectral_status_success, spectral_status_invalid_input
    implicit none
    real(dp), allocatable :: w(:), v(:)
    integer :: st, n, i

    w = get_window('kaiser', 7, fftbins=.false., parameter=6.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'kaiser status'
    call close(w, [1.4873337104763207e-02_dp,2.5370608027132219e-01_dp,7.3189541539782077e-01_dp, &
                   1.0_dp,7.3189541539782077e-01_dp,2.5370608027132219e-01_dp, &
                   1.4873337104763207e-02_dp])
    w = get_window('gaussian', 7, fftbins=.false., parameter=2.0_dp)
    call close(w, [0.32465246735834974_dp, 0.6065306597126334_dp, &
                   0.8824969025845953_dp,1.0_dp,0.8824969025845953_dp, &
                   0.6065306597126334_dp,0.32465246735834974_dp])
    w = get_window('tukey', 7, fftbins=.false.)
    call close(w,[0.0_dp,0.75_dp,1.0_dp,1.0_dp,1.0_dp,0.75_dp,0.0_dp])
    w = get_window('general_gaussian',7,fftbins=.false.,parameter=2.0_dp,parameter2=3.0_dp)
    call close(w,[0.6065306597126334_dp,0.9059551911095096_dp,0.9938461733264412_dp, &
                  1.0_dp,0.9938461733264412_dp,0.9059551911095096_dp,0.6065306597126334_dp])
    w = get_window('triang', 7)
    call close(w,[0.125_dp,0.375_dp,0.625_dp,0.875_dp,0.875_dp,0.625_dp,0.375_dp])
    w = get_window('blackmanharris', 1, status=st)
    if (st /= spectral_status_success .or. size(w) /= 1) error stop 'unit window'
    if (abs(w(1)-1.0_dp)>0.0_dp) error stop 'single element'
    w = get_window('parzen', 0, status=st)
    if (st /= spectral_status_success .or. size(w) /= 0) error stop 'empty window'
    w = get_window('gaussian', 5, status=st)
    if (st /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'missing parameter'
    w = get_window('general_gaussian', 5, parameter=2.0_dp, status=st)
    if (st /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'missing width'
    w = get_window('kaiser', 5, parameter=ieee_value(0.0_dp,ieee_quiet_nan),status=st)
    if (st /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'nonfinite parameter'
    w = get_window('none', 5, status=st)
    if (st /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'invalid window name'
    w = get_window('hann', -1, status=st)
    if (st /= spectral_status_invalid_input .or. size(w) /= 0) error stop 'negative size'
    do n=2,21
        w = get_window('lanczos', n, fftbins=.false.,status=st)
        if (st /= spectral_status_success) error stop 'lanczos status'
        if (maxval(abs(w-w(n:1:-1)))>1.0e-14_dp) error stop 'lanczos symmetry'
        v = get_window('bohman', n, fftbins=.false.)
        if (maxval(abs(v-v(n:1:-1)))>1.0e-14_dp) error stop 'bohman symmetry'
        v = get_window('kaiser', n, fftbins=.false.,parameter=40.0_dp)
        if (maxval(abs(v-v(n:1:-1)))>1.0e-12_dp) error stop 'kaiser symmetry'
    end do
    print '(a)', 'test_signal_windows_extra passed'
contains
    subroutine close(actual, expected)
        real(dp), intent(in) :: actual(:) !! calculated window
        real(dp), intent(in) :: expected(:) !! SciPy reference
        if (size(actual) /= size(expected)) error stop 'window length'
        if (any(abs(actual-expected)>5.0e-12_dp)) error stop 'window values'
    end subroutine close
end program test_signal_windows_extra
