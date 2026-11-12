! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Public scifort_signal dispatch coverage (not just the design submodule).
program test_signal_firwin_public
    use scifort_signal, only : dp, firwin, design_status_success, &
        design_status_invalid_input
    implicit none
    real(dp), allocatable :: r(:), b(:)
    integer :: st

    r = firwin(31, 0.20_dp, status=st)
    if (st /= design_status_success .or. size(r) /= 31) error stop 'scalar public firwin'
    b = firwin(31, [0.20_dp], pass_zero=.true., status=st)
    if (st /= design_status_success) error stop 'vector/logical public firwin'
    if (maxval(abs(r-b)) > 1.0e-12_dp) error stop 'public generics disagree'

    r = firwin(30, [0.11_dp,0.26_dp,0.47_dp], pass_zero='bandstop', &
        status=st)
    if (st /= design_status_success .or. size(r) /= 30) &
        error stop 'public multiband firwin'
    b = firwin(30, [0.11_dp,0.26_dp,0.47_dp], pass_zero=.false., &
        status=st)
    if (st /= design_status_invalid_input .or. size(b) /= 0) &
        error stop 'public invalid design'
    print '(a)', 'PASS: scifort_signal public firwin generic'
end program test_signal_firwin_public
