! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_group_delay_sos_public
    use scifort_signal, only : dp, group_delay, group_delay_sos, &
        group_delay_result, design_status_success
    implicit none
    type(group_delay_result) :: gd
    real(dp) :: cascade(1,6)
    cascade(1,:) = [0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp]
    gd = group_delay_sos(cascade,worn=17)
    if (gd%status /= design_status_success) error stop 'SOS public status'
    if (size(gd%delay) /= 17) error stop 'SOS public grid'
    if (any(abs(gd%delay-1.0_dp)>2.0e-13_dp)) error stop 'SOS public values'
    gd = group_delay(cascade,w=[0.0_dp,0.3_dp,1.4_dp])
    if (gd%status /= design_status_success) error stop 'SOS generic public status'
    if (any(abs(gd%delay-1.0_dp)>2.0e-13_dp)) error stop 'SOS generic public values'
    print '(a)', 'test_signal_group_delay_sos_public passed'
end program test_signal_group_delay_sos_public
