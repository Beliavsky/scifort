! SPDX-License-Identifier: MIT
program test_signal_dlsim_sos_bank_public
    use scifort_signal, only : dp, dlsim, dlsim_sos, dlsim_sos_bank, &
        discrete_sos_bank_result, discrete_sim_success
    implicit none
    real(dp) :: sos(2,6), u(11,3), bank(2,6,3), zi(2,2,3)
    complex(dp) :: csos(2,6,3), cu(11,3)
    type(discrete_sos_bank_result) :: a, b
    integer :: i, j

    sos(1,:) = [0.2_dp,0.3_dp,0.1_dp,1.0_dp,-0.4_dp,0.08_dp]
    sos(2,:) = [0.7_dp,0.2_dp,0.1_dp,1.0_dp,0.1_dp,0.05_dp]
    zi = 0.0_dp
    do j = 1, 3
        bank(:,:,j) = sos
        do i = 1, 11
            u(i,j) = sin(real(i+j,dp))
        end do
    end do
    a = dlsim_sos_bank(sos,u,dt=0.1_dp,zi=zi)
    b = dlsim(bank,u,dt=0.1_dp,zi=zi)
    if (a%status/=discrete_sim_success .or. b%status/=discrete_sim_success) &
        error stop 'public bank status'
    if (maxval(abs(a%y-b%y))>1.0e-13_dp) error stop 'public bank overload mismatch'
    csos = cmplx(bank,0.0_dp,dp)
    csos(1,2,2) = cmplx(0.3_dp,0.04_dp,dp)
    cu = cmplx(u,0.3_dp,dp)
    a = dlsim_sos(csos,cu,dt=0.1_dp)
    if (a%status/=discrete_sim_success) error stop 'public complex bank simulation'
    if (any(shape(a%zf)/=[2,2,3])) error stop 'public bank final-state shape'
    print '(a)', 'PASS public multichannel SOS simulation overloads'
end program test_signal_dlsim_sos_bank_public
