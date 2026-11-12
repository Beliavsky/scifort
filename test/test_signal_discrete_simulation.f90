! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_discrete_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, dlsim, dstep, dimpulse, dlsim_zpk, dstep_zpk, dimpulse_zpk, &
        discrete_simulation_result, discrete_response_result, discrete_sim_success, discrete_sim_invalid_input
    implicit none
    type(discrete_simulation_result) :: sim
    type(discrete_response_result) :: res, pulse
    real(dp) :: a(1,1), b(1,1), c(1,1), d(1,1), u(3,1), nan_value
    real(dp) :: two_b(1,2), two_c(2,1), two_d(2,2), us(3,2)
    real(dp) :: zero_a(0,0), zero_b(0,2), zero_c(1,0), zero_d(1,2)
    complex(dp) :: az(1,1), bz(1,1), cz(1,1), dz(1,1), uz(2,1)
    real(dp), parameter :: tol = 3.0e-12_dp

    a = 0.5_dp
    b = 2.0_dp
    c = 3.0_dp
    d = 4.0_dp
    u(:,1) = [1.0_dp,2.0_dp,0.0_dp]
    sim = dlsim(a,b,c,d,u,0.25_dp,x0=[1.0_dp])
    call check(sim%status==discrete_sim_success,'SISO state simulation status')
    call check(size(sim%t)==3 .and. size(sim%x,2)==1,'SISO shapes')
    call check(maxval(abs(sim%y(:,1)-cmplx([7.0_dp,15.5_dp,15.75_dp],0.0_dp,dp)))<tol, &
        'SISO output')
    call check(maxval(abs(sim%x(:,1)-cmplx([1.0_dp,2.5_dp,5.25_dp],0.0_dp,dp)))<tol, &
        'initial state and recurrence')
    call check(maxval(abs(sim%t-[0.0_dp,0.25_dp,0.5_dp]))<tol,'SISO times')

    ! Unequally sampled inputs are interpolated onto scipy.signal.dlsim's
    ! linspace(0,t(end),floor(t(end)/dt)+1) output grid.
    a = 0.0_dp
    b = 0.0_dp
    c = 0.0_dp
    d = 1.0_dp
    sim = dlsim(a,b,c,d,reshape([0.0_dp,0.3_dp,0.9_dp],[3,1]), &
        0.4_dp,t=[0.0_dp,0.3_dp,0.9_dp])
    call check(sim%status==discrete_sim_success,'interpolated input status')
    call check(size(sim%t)==3,'off-grid length')
    call check(maxval(abs(sim%t-[0.0_dp,0.45_dp,0.9_dp]))<tol,'off-grid output times')
    call check(maxval(abs(sim%y(:,1)-cmplx(sim%t,0.0_dp,dp)))<tol,'interpolation values')

    ! For two inputs, each third-index slice is a separate unit-step or
    ! unit-impulse experiment; outputs are nt x nout x nin.
    a = 0.3_dp
    two_b(1,:) = [1.0_dp,2.0_dp]
    two_c(:,1) = [1.0_dp,3.0_dp]
    two_d(1,:) = [0.0_dp,1.0_dp]
    two_d(2,:) = [2.0_dp,0.0_dp]
    res = dstep(a,two_b,two_c,two_d,0.5_dp,n=5)
    pulse = dimpulse(a,two_b,two_c,two_d,0.5_dp,n=5)
    call check(res%status==discrete_sim_success .and. pulse%status==discrete_sim_success, &
        'MIMO response statuses')
    call check(all(shape(res%y)==[5,2,2]),'MIMO output shape')
    call check(maxval(abs(res%y(1,:,2)-cmplx([1.0_dp,0.0_dp],0.0_dp,dp)))<tol, &
        'step feedthrough input 2')
    call check(maxval(abs(res%y(2,:,2)-cmplx([3.0_dp,6.0_dp],0.0_dp,dp)))<tol, &
        'step state response input 2')
    call check(maxval(abs(pulse%y(3,:,2)-cmplx([0.6_dp,1.8_dp],0.0_dp,dp)))<tol, &
        'impulse propagation')
    us(:,1) = [1.0_dp,1.0_dp,1.0_dp]
    us(:,2) = 0.0_dp
    sim = dlsim(a,two_b,two_c,two_d,us,0.5_dp)
    call check(maxval(abs(sim%y(:,1)-res%y(1:3,1,1)))<tol,'step agrees with dlsim')

    ! Transfer functions use polynomials in descending powers of z.
    sim = dlsim([0.0_dp,1.0_dp],[1.0_dp,-0.5_dp],[1.0_dp,1.0_dp,1.0_dp,1.0_dp],1.0_dp)
    call check(sim%status==discrete_sim_success,'TF simulation')
    call check(maxval(abs(sim%y(:,1)-cmplx([0.0_dp,1.0_dp,1.5_dp,1.75_dp],0.0_dp,dp)))<tol, &
        'TF integrator samples')
    res = dstep([0.0_dp,1.0_dp],[1.0_dp,-0.5_dp],1.0_dp,n=4)
    pulse = dimpulse([0.0_dp,1.0_dp],[1.0_dp,-0.5_dp],1.0_dp,n=4)
    call check(res%status==discrete_sim_success .and. pulse%status==discrete_sim_success,'TF response')
    call check(maxval(abs(res%y(:,1,1)-cmplx([0.0_dp,1.0_dp,1.5_dp,1.75_dp],0.0_dp,dp)))<tol, &
        'TF step values')
    call check(maxval(abs(pulse%y(:,1,1)-cmplx([0.0_dp,1.0_dp,0.5_dp,0.25_dp],0.0_dp,dp)))<tol, &
        'TF impulse values')
    sim = dlsim_zpk([real(dp) :: ],[0.5_dp],1.0_dp,[1.0_dp,1.0_dp,1.0_dp],1.0_dp)
    call check(sim%status==discrete_sim_success,'ZPK simulation')
    call check(maxval(abs(sim%y(:,1)-cmplx([0.0_dp,1.0_dp,1.5_dp],0.0_dp,dp)))<tol, &
        'ZPK simulation values')
    res = dstep_zpk([real(dp) :: ],[0.5_dp],1.0_dp,1.0_dp,n=3)
    pulse = dimpulse_zpk([real(dp) :: ],[0.5_dp],1.0_dp,1.0_dp,n=3)
    call check(res%status==discrete_sim_success .and. pulse%status==discrete_sim_success, &
        'ZPK responses')
    call check(abs(res%y(3,1,1)-1.5_dp)<tol,'ZPK step sample')
    call check(abs(pulse%y(3,1,1)-0.5_dp)<tol,'ZPK impulse sample')

    ! Static gain, including the zero-dimensional state case.
    zero_d(1,:) = [2.0_dp,3.0_dp]
    sim = dlsim(zero_a,zero_b,zero_c,zero_d, &
        reshape([1.0_dp,5.0_dp,2.0_dp,6.0_dp],[2,2]),0.2_dp)
    call check(sim%status==discrete_sim_success,'empty state')
    call check(size(sim%x,2)==0 .and. size(sim%y,1)==2,'empty state shape')
    call check(abs(sim%y(1,1)-8.0_dp)<tol .and. abs(sim%y(2,1)-28.0_dp)<tol, &
        'static gain output')
    res = dstep([2.0_dp],[1.0_dp],0.25_dp,n=3)
    call check(res%status==discrete_sim_success,'static TF step')
    call check(maxval(abs(res%y(:,1,1)-2.0_dp))<tol,'static TF response')

    az(1,1) = cmplx(0.25_dp,0.1_dp,dp)
    bz(1,1) = cmplx(1.0_dp,-0.3_dp,dp)
    cz(1,1) = cmplx(0.5_dp,0.2_dp,dp)
    dz(1,1) = cmplx(0.1_dp,0.0_dp,dp)
    uz(:,1) = cmplx([1.0_dp,2.0_dp],[0.5_dp,-0.25_dp],dp)
    sim = dlsim(az,bz,cz,dz,uz,0.1_dp,x0=[cmplx(0.3_dp,0.2_dp,dp)])
    call check(sim%status==discrete_sim_success,'complex simulation')
    call check(abs(sim%y(1,1)-cz(1,1)*cmplx(0.3_dp,0.2_dp,dp)-dz(1,1)*uz(1,1))<tol, &
        'complex initial output')
    call check(abs(sim%x(2,1)-az(1,1)*sim%x(1,1)-bz(1,1)*uz(1,1))<tol, &
        'complex recurrence')
    sim = dlsim([1.0_dp], [1.0_dp,-0.4_dp],cmplx([1.0_dp,2.0_dp],0.2_dp,dp),0.5_dp)
    call check(sim%status==discrete_sim_success,'real TF with complex input')

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    sim = dlsim(a,b,c,d,u,0.0_dp)
    call check(sim%status==discrete_sim_invalid_input,'invalid dt')
    sim = dlsim(a,b,c,d,u,0.25_dp,t=[0.0_dp,2.0_dp])
    call check(sim%status==discrete_sim_invalid_input,'wrong time length')
    sim = dlsim(a,b,c,d,u,0.25_dp,t=[0.0_dp,0.2_dp,0.1_dp])
    call check(sim%status==discrete_sim_invalid_input,'nonmonotone times')
    sim = dlsim(a,b,c,d,u,0.25_dp,x0=[1.0_dp,2.0_dp])
    call check(sim%status==discrete_sim_invalid_input,'wrong state dimension')
    sim = dlsim(a,b,c,d,u,nan_value)
    call check(sim%status==discrete_sim_invalid_input,'nan dt')
    res = dimpulse(a,b,c,d,0.25_dp,n=0)
    call check(res%status==discrete_sim_invalid_input,'zero response length')
    res = dstep([1.0_dp,3.0_dp,4.0_dp],[1.0_dp,0.5_dp],0.25_dp,n=3)
    call check(res%status==discrete_sim_invalid_input,'improper TF rejection')

    print '(a)', 'PASS: discrete simulation and responses'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! diagnostic label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_discrete_simulation
