! SPDX-License-Identifier: MIT
! Check that new window functions are reexported by scifort_signal.
program test_signal_windows_advanced_public
    use scifort_kinds, only : dp
    use scifort_signal, only : get_window, chebwin, taylor, &
        kaiser_bessel_derived, spectral_status_success
    implicit none
    real(dp), allocatable :: w(:), v(:)
    integer :: st
    w = chebwin(17, 80.0_dp, status=st)
    if (st /= spectral_status_success .or. size(w) /= 17) error stop 'public chebwin'
    v = get_window('cheb', 17, parameter=80.0_dp, fftbins=.false.)
    if (any(abs(w-v)>2.0e-13_dp)) error stop 'public chebwin alias'
    w = taylor(17, nbar=5, sll=65.0_dp)
    v = get_window('taylorwin', 17, parameter=5.0_dp, parameter2=65.0_dp, fftbins=.false.)
    if (any(abs(w-v)>2.0e-13_dp)) error stop 'public Taylor'
    w = kaiser_bessel_derived(16, 8.0_dp)
    v = get_window('kaiser_bessel_derived', 16, parameter=8.0_dp, fftbins=.false.)
    if (any(abs(w-v)>2.0e-13_dp)) error stop 'public KBD'
    print '(a)', 'test_signal_windows_advanced_public passed'
end program test_signal_windows_advanced_public
