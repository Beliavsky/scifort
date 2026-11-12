! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_windows_public
    use scifort_signal, only : get_window, spectral_status_success
    use scifort_kinds, only : dp
    implicit none
    real(dp), allocatable :: w(:)
    integer :: st
    w = get_window('nuttall',17,status=st)
    if (st /= spectral_status_success .or. size(w)/=17) error stop 'public nuttall'
    w = get_window('kaiser',13,parameter=8.0_dp,status=st)
    if (st /= spectral_status_success .or. size(w)/=13) error stop 'public kaiser'
    print '(a)', 'test_signal_windows_public passed'
end program test_signal_windows_public
