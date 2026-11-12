! SPDX-License-Identifier: MIT
program test_signal_dlsim_sos_public
    use scifort_signal, only : dp, dlsim, dlsim_sos, discrete_simulation_result, &
        discrete_sim_success
    implicit none
    real(dp), parameter :: sos(1,6) = reshape([0.2_dp,0.3_dp,0.1_dp, &
        1.0_dp,-0.4_dp,0.08_dp],[1,6])
    real(dp) :: u(7)
    complex(dp) :: csos(1,6), cu(7)
    type(discrete_simulation_result) :: a, b
    integer :: k

    u = [(sin(real(k,dp)), k=1,7)]
    a = dlsim_sos(sos,u,dt=0.1_dp)
    b = dlsim(sos,u,dt=0.1_dp)
    if (a%status/=discrete_sim_success .or. b%status/=discrete_sim_success) &
        error stop 'real public simulation failed'
    if (maxval(abs(a%y-b%y))>1.0e-13_dp) error stop 'public overload mismatch'
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = cmplx(0.3_dp,0.02_dp,dp)
    cu = cmplx(u,0.2_dp,dp)
    a = dlsim_sos(csos,cu,dt=0.1_dp)
    if (a%status/=discrete_sim_success) error stop 'complex public simulation failed'
    if (any(shape(a%zf)/=[1,2])) error stop 'public final state mismatch'
    print '(a)', 'PASS public dlsim_sos overloads'
end program test_signal_dlsim_sos_public
