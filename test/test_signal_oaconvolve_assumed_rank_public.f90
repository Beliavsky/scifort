! SPDX-License-Identifier: MIT
! Ensure the scifort_signal facade exports overlap-add assumed-rank filtering.
program test_signal_oaconvolve_assumed_rank_public
    use scifort_signal, only : dp, oaconvolve_nd, convolution_status_success
    implicit none
    real(dp) :: a(4,2), b(2,1), y(5,2)
    integer :: status

    a(:,1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    a(:,2) = [-1.0_dp, 3.0_dp, 2.0_dp, -2.0_dp]
    b(:,1) = [1.0_dp, -0.5_dp]
    call oaconvolve_nd(a,b,y,axes=[1],status=status)
    if (status /= convolution_status_success) error stop 'public OA status'
    if (maxval(abs(y(:,1)-[1.0_dp,1.5_dp,2.0_dp,2.5_dp,-2.0_dp])) &
        > 1.e-12_dp) error stop 'public OA values'
    print '(a)', 'test_signal_oaconvolve_assumed_rank_public passed'
end program test_signal_oaconvolve_assumed_rank_public
