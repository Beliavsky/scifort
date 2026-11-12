! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_group_delay_zpk_public
    use scifort_signal, only : dp, group_delay, group_delay_zpk, &
        group_delay_result, design_status_success
    implicit none
    type(group_delay_result) :: result
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: w(3) = [0.1_dp,0.5_dp,1.3_dp]
    result = group_delay_zpk([0.4_dp],[0.7_dp],2.0_dp,w=w)
    if (result%status /= design_status_success) error stop 'public ZPK status'
    if (size(result%delay) /= 3) error stop 'public ZPK shape'
    result = group_delay([cmplx(0.4_dp,0.0_dp,dp)], &
        [cmplx(0.7_dp,0.0_dp,dp)],cmplx(1.0_dp,0.2_dp,dp),w=w)
    if (result%status /= design_status_success) error stop 'public generic status'
    if (any(abs(result%frequencies-w)>1.0e-13_dp)) error stop 'public generic grid'
    result = group_delay_zpk([complex(dp) ::], [cmplx(0.0_dp,0.0_dp,dp)], &
        (1.0_dp,0.0_dp), worn=16,whole=.true.)
    if (result%status /= design_status_success) error stop 'public whole-circle status'
    if (any(abs(result%delay-1.0_dp)>1.0e-12_dp)) error stop 'public delay'
    if (abs(result%frequencies(16)-2.0_dp*pi*15.0_dp/16.0_dp)>1.0e-13_dp) &
        error stop 'public whole-circle grid'
    print '(a)', 'test_signal_group_delay_zpk_public passed'
end program test_signal_group_delay_zpk_public
