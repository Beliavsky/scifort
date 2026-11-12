! SPDX-License-Identifier: MIT
! Public scifort_signal export smoke test for SOS-native discrete responses.
program test_signal_sos_responses_public
    use scifort_signal, only : dp, dstep, dimpulse, dstep_sos, dimpulse_sos, &
        discrete_response_result, discrete_sim_success
    implicit none
    real(dp) :: rsos(1,6)
    complex(dp) :: csos(1,6)
    type(discrete_response_result) :: result, alias

    rsos(1,:) = [0.25_dp,0.0_dp,0.0_dp,1.0_dp,-0.5_dp,0.0_dp]
    result = dimpulse(sos=rsos,dt=0.2_dp,n=4)
    alias = dimpulse_sos(rsos,0.2_dp,n=4)
    if (result%status /= discrete_sim_success) error stop 'public real SOS impulse status'
    if (maxval(abs(result%y(:,1,1) - &
        cmplx([0.25_dp,0.125_dp,0.0625_dp,0.03125_dp],0.0_dp,dp))) > 1.0e-13_dp) &
        error stop 'public real SOS impulse values'
    if (maxval(abs(result%y-alias%y)) > 1.0e-13_dp) error stop 'public real SOS alias'

    csos = cmplx(rsos,0.0_dp,dp)
    csos(1,5) = cmplx(-0.5_dp,-0.25_dp,dp)
    result = dstep(sos=csos,dt=0.1_dp,n=5)
    alias = dstep_sos(csos,0.1_dp,n=5)
    if (result%status /= discrete_sim_success) error stop 'public complex SOS step status'
    if (maxval(abs(result%y-alias%y)) > 1.0e-13_dp) error stop 'public complex SOS alias'
    if (abs(aimag(result%y(2,1,1))) < 0.01_dp) error stop 'complex coefficients preserved'
    print '(a)', 'PASS public SOS response interfaces'
end program test_signal_sos_responses_public
