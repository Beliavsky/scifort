! SPDX-License-Identifier: MIT
program test_signal_dlsim_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dlsim, dlsim_sos, dstep_sos, &
        discrete_simulation_result, discrete_response_result, discrete_sim_success, &
        discrete_sim_invalid_input, discrete_sim_numerical_failure
    use scifort_signal_sos, only : sosfilt
    implicit none
    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp) :: sos(2,6), zi(2,2), values(41), t(41), bad(2,6)
    real(dp) :: noisy(82)
    complex(dp) :: csos(2,6), czi(2,2), inputs(41)
    complex(dp), allocatable :: ref(:)
    type(discrete_simulation_result) :: r, r1, r2
    type(discrete_response_result) :: step
    integer :: i

    sos(1,:) = [0.22_dp, 0.11_dp, -0.05_dp, 1.0_dp, -0.37_dp, 0.10_dp]
    sos(2,:) = [0.83_dp, 0.22_dp, 0.02_dp, 1.0_dp, 0.17_dp, 0.04_dp]
    zi(1,:) = [0.03_dp, -0.02_dp]
    zi(2,:) = [-0.01_dp, 0.08_dp]
    do i = 1, 41
        values(i) = sin(0.21_dp*real(i,dp)) + 0.4_dp*cos(0.37_dp*real(i,dp))
        inputs(i) = cmplx(values(i),cos(0.18_dp*real(i,dp)),dp)
        t(i) = 0.075_dp*real(i-1,dp)*(1.0_dp+0.003_dp*real(i,dp))
    end do

    r = dlsim_sos(sos, values, 0.075_dp)
    call check(r%status == discrete_sim_success,'real status')
    call check(all(shape(r%y)==[41,1]),'real output shape')
    call check(all(shape(r%x)==[41,4]),'state history shape')
    call check(all(shape(r%zf)==[2,2]),'final state shape')
    ref = cmplx(sosfilt(sos,values),0.0_dp,dp)
    call check(maxval(abs(r%y(:,1)-ref)) < tol,'real sosfilt equality')
    call check(maxval(abs(r%x(1,:)))==0.0_dp,'default initial states')
    r1 = dlsim(sos, values, 0.075_dp, zi=zi)
    r2 = dlsim_sos(sos, values, 0.075_dp, zi=zi)
    call check(r1%status == discrete_sim_success,'generic real alias')
    call check(maxval(abs(r1%y-r2%y)) < tol,'generic value')
    call check(maxval(abs(r2%y(:,1)-cmplx(sosfilt(sos,values,zi=zi),0.0_dp,dp))) &
        < tol,'real initial state')
    call check(abs(r2%x(1,1)-cmplx(zi(1,1),0.0_dp,dp)) < tol,'first delay stored')
    call check(abs(r2%x(1,2)-cmplx(zi(1,2),0.0_dp,dp)) < tol,'second delay stored')
    call check(abs(r2%x(1,3)-cmplx(zi(2,1),0.0_dp,dp)) < tol,'second section order')

    ! The saved final state permits streaming a long signal in smaller blocks.
    r1 = dlsim_sos(sos,values(:19),0.075_dp,zi=zi)
    r2 = dlsim_sos(sos,values(20:),0.075_dp,zi=real(r1%zf,dp))
    r = dlsim_sos(sos,values,0.075_dp,zi=zi)
    call check(maxval(abs(r%y(:19,1)-r1%y(:,1)))<tol,'first block')
    call check(maxval(abs(r%y(20:,1)-r2%y(:,1)))<tol,'second block')
    call check(maxval(abs(r%zf-r2%zf))<tol,'final state continuity')

    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = cmplx(0.11_dp,0.06_dp,dp)
    csos(2,5) = cmplx(0.17_dp,-0.03_dp,dp)
    czi = cmplx(zi,0.025_dp,dp)
    r = dlsim_sos(csos,inputs,0.25_dp,zi=czi)
    call check(r%status==discrete_sim_success,'complex status')
    ref = sosfilt(csos,inputs,zi=czi)
    call check(maxval(abs(r%y(:,1)-ref))<tol,'complex filter')
    r1 = dlsim(csos,values,0.25_dp,zi=czi)
    r2 = dlsim(csos,cmplx(values,0.0_dp,dp),0.25_dp,zi=czi)
    call check(maxval(abs(r1%y-r2%y))<tol,'complex SOS real inputs')
    r1 = dlsim_sos(sos,inputs,0.25_dp,zi=czi)
    ref = sosfilt(cmplx(sos,0.0_dp,dp),inputs,zi=czi)
    call check(maxval(abs(r1%y(:,1)-ref))<tol,'real SOS complex inputs')
    r2 = dlsim_sos(csos,inputs(:20),0.25_dp,zi=czi)
    r1 = dlsim_sos(csos,inputs(21:),0.25_dp,zi=r2%zf)
    call check(maxval(abs(r%y(21:,1)-r1%y(:,1)))<tol,'complex state continuity')

    ! The input times are on a different grid from the output times.
    r = dlsim_sos(sos,values,0.1_dp,t=t)
    call check(r%status==discrete_sim_success,'irregular interpolation')
    call check(size(r%t)==int(floor(t(41)/0.1_dp))+1,'interpolation length')
    call check(abs(r%t(size(r%t))-t(41))<tol,'interpolation end')
    call check(abs(r%t(1))<tol,'interpolation origin')
    do i = 1, size(r%t)
        call check(abs(r%y(i,1))<10.0_dp,'interpolation bounded')
    end do
    ! Accept noncontiguous 1-D views for SOS and input signals.
    noisy(1::2) = values
    noisy(2::2) = -3.0_dp
    r = dlsim_sos(sos(:,1:6),noisy(::2),0.05_dp)
    call check(r%status==discrete_sim_success,'strided input')
    call check(maxval(abs(r%y(:,1)-cmplx(sosfilt(sos,values),0.0_dp,dp))) &
        < tol,'strided input matches')

    step = dstep_sos(sos,0.025_dp,n=19,zi=zi)
    r = dlsim_sos(sos,[(1.0_dp,i=1,19)],0.025_dp,zi=zi)
    call check(maxval(abs(step%y(:,1,1)-r%y(:,1)))<tol,'step compatibility')

    r = dlsim_sos(sos,values,-1.0_dp)
    call invalid(r,'negative dt')
    r = dlsim_sos(sos,values,ieee_value(0.0_dp,ieee_quiet_nan))
    call invalid(r,'nonfinite dt')
    r = dlsim_sos(sos,values(:0),0.1_dp)
    call invalid(r,'empty input')
    r = dlsim_sos(sos(:0,:),values,0.1_dp)
    call invalid(r,'empty sos')
    r = dlsim_sos(sos(:,:5),values,0.1_dp)
    call invalid(r,'five SOS columns')
    bad = sos
    bad(2,4) = 2.0_dp
    r = dlsim_sos(bad,values,0.1_dp)
    call invalid(r,'SOS a0 normalization')
    bad = sos
    bad(1,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = dlsim_sos(bad,values,0.1_dp)
    call invalid(r,'nonfinite coefficient')
    r = dlsim_sos(sos,values,0.1_dp,zi=zi(:1,:))
    call invalid(r,'initial state shape')
    r = dlsim_sos(csos,inputs,0.1_dp,zi=czi(:,1:1))
    call invalid(r,'complex initial state shape')
    r = dlsim_sos(sos,values,0.1_dp,t=t(:39))
    call invalid(r,'input time shape')
    r = dlsim_sos(sos,values,0.1_dp,t=t(41:1:-1))
    call invalid(r,'input time ordering')
    r = dlsim_sos(sos,values,0.1_dp,t=-t)
    call invalid(r,'negative input endpoint')
    bad = sos
    bad(1,1) = huge(0.0_dp)
    r = dlsim_sos(bad,values,0.1_dp)
    call check(r%status==discrete_sim_numerical_failure,'overflow diagnosis')
    print '(a)', 'PASS SOS-native arbitrary-input simulation'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion outcome
        character(*), intent(in) :: message !! description
        if (.not. ok) then
            print '(2a)', 'FAIL: ', message
            error stop 1
        end if
    end subroutine check
    subroutine invalid(v,message)
        type(discrete_simulation_result), intent(in) :: v !! outcome to inspect
        character(*), intent(in) :: message !! invalid test name
        call check(v%status==discrete_sim_invalid_input,message)
        call check(.not. allocated(v%y), 'invalid: output remains unallocated')
        call check(.not. allocated(v%zf), 'invalid: final state remains unallocated')
    end subroutine invalid
end program test_signal_dlsim_sos
