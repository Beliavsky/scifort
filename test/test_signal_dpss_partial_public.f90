! SPDX-License-Identifier: MIT
! Smoke test the new solver option through the public scifort_signal module.
program test_signal_dpss_partial_public
    use scifort_signal, only : dp, dpss, dpss_result, dpss_status_success
    implicit none
    type(dpss_result) :: tapers
    tapers = dpss(4096, 3.5_dp, kmax=4, solver='partial')
    if (tapers%status /= dpss_status_success) error stop 'public partial DPSS'
    if (any(shape(tapers%tapers) /= [4096, 4])) error stop 'public DPSS dimensions'
    if (abs(sum(tapers%tapers(:,1)**2) - 1.0_dp) > 1.0e-10_dp) then
        error stop 'public DPSS norm'
    end if
    print '(a)', 'Public selected DPSS eigenpair test passed'
end program test_signal_dpss_partial_public
