! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_continuous_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, lsim, step, impulse, lsim_zpk, step_zpk, impulse_zpk, &
        continuous_simulation_result, continuous_response_result, &
        continuous_sim_success, continuous_sim_invalid_input
    implicit none
    type(continuous_simulation_result) :: sim
    type(continuous_response_result) :: res, pulse
    real(dp) :: a(1,1), b(1,1), c(1,1), d(1,1), grid(5), u(5,1), nan_value
    real(dp) :: two_b(1,2), two_c(2,1), two_d(2,2), empty_a(0,0), empty_b(0,1)
    real(dp) :: empty_c(1,0), empty_d(1,1)
    complex(dp) :: ac(1,1), bc(1,1), cc(1,1), dc(1,1), uc(5,1)
    integer :: k
    real(dp), parameter :: tol=2.0e-11_dp

    a = -1.0_dp
    b = 1.0_dp
    c = 2.0_dp
    d = 0.5_dp
    grid = [0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp]
    u(:,1) = grid
    sim = lsim(a,b,c,d,u,grid)
    call check(sim%status==continuous_sim_success,'FOH status')
    call check(all(shape(sim%y)==[5,1]) .and. all(shape(sim%x)==[5,1]),'FOH shapes')
    call check(maxval(abs(sim%x(:,1)-(grid-1.0_dp+exp(-grid))))<tol,'FOH ramp state')
    call check(maxval(abs(sim%y(:,1)-(2.0_dp*(grid-1.0_dp+exp(-grid))+0.5_dp*grid)))<tol, &
        'FOH ramp output')
    u(:,1) = [1.0_dp,1.0_dp,2.0_dp,2.0_dp,0.0_dp]
    sim = lsim(a,b,c,d,u,grid,interp=.false.)
    call check(sim%status==continuous_sim_success,'ZOH status')
    call check(abs(sim%x(2,1)-(1.0_dp-exp(-0.25_dp)))<tol,'ZOH first state')
    call check(abs(sim%x(4,1)-(exp(-0.25_dp)*sim%x(3,1)+ &
        2.0_dp*(1.0_dp-exp(-0.25_dp))))<tol,'ZOH new level')

    sim = lsim(a,b,c,d,reshape([0.0_dp,0.0_dp],[2,1]),[0.5_dp,1.0_dp],x0=[3.0_dp])
    call check(sim%status==continuous_sim_success,'nonzero start status')
    call check(abs(sim%x(1,1)-3.0_dp*exp(-0.5_dp))<tol,'nonzero start state')
    sim = lsim(a,b,c,d,reshape([0.0_dp],[1,1]),[0.25_dp],x0=[3.0_dp])
    call check(sim%status==continuous_sim_success,'one sample')
    call check(abs(sim%x(1,1)-3.0_dp*exp(-0.25_dp))<tol,'one sample state')

    two_b(1,:) = [1.0_dp,2.0_dp]
    two_c(:,1) = [1.0_dp,3.0_dp]
    two_d(1,:) = [0.25_dp,-1.0_dp]
    two_d(2,:) = [0.5_dp,0.75_dp]
    res = step(a,two_b,two_c,two_d,grid)
    pulse = impulse(a,two_b,two_c,two_d,grid)
    call check(res%status==continuous_sim_success .and. pulse%status==continuous_sim_success, &
        'MIMO response statuses')
    call check(all(shape(res%y)==[5,2,2]),'MIMO response shape')
    do k=1,5
        call check(maxval(abs(res%y(k,:,2)-(two_c(:,1)*2.0_dp*(1.0_dp-exp(-grid(k))) &
            +two_d(:,2))))<tol,'MIMO step channel')
        call check(maxval(abs(pulse%y(k,:,2)-two_c(:,1)*2.0_dp*exp(-grid(k))))<tol, &
            'MIMO impulse channel')
    end do
    call check(abs(res%y(1,1,1)-two_d(1,1))<tol,'step feedthrough')
    call check(abs(pulse%y(1,1,1)-1.0_dp)<tol,'impulse ignores feedthrough delta')

    ! Irregular response times and nonzero initial state.
    res = step(a,b,c,d,[0.0_dp,0.1_dp,0.5_dp,1.75_dp],x0=[2.0_dp])
    pulse = impulse(a,b,c,d,[0.0_dp,0.1_dp,0.5_dp,1.75_dp],x0=[2.0_dp])
    call check(res%status==continuous_sim_success .and. pulse%status==continuous_sim_success, &
        'irregular response times')
    call check(abs(res%y(4,1,1)-(2.0_dp+2.0_dp*exp(-1.75_dp)+0.5_dp))<tol, &
        'irregular step plus initial condition')
    call check(abs(pulse%y(4,1,1)-6.0_dp*exp(-1.75_dp))<tol, &
        'irregular impulse plus initial condition')

    res = step([1.0_dp],[1.0_dp,1.0_dp],grid)
    pulse = impulse([1.0_dp],[1.0_dp,1.0_dp],grid)
    call check(res%status==continuous_sim_success .and. pulse%status==continuous_sim_success, &
        'transfer function response status')
    call check(maxval(abs(res%y(:,1,1)-(1.0_dp-exp(-grid))))<tol,'TF step')
    call check(maxval(abs(pulse%y(:,1,1)-exp(-grid)))<tol,'TF impulse')
    sim = lsim([1.0_dp],[1.0_dp,1.0_dp],grid,grid)
    call check(sim%status==continuous_sim_success,'transfer function simulation')
    call check(maxval(abs(sim%y(:,1)-(grid-1.0_dp+exp(-grid))))<tol,'TF ramp')
    res = step_zpk([real(dp) ::],[-1.0_dp],1.0_dp,grid)
    pulse = impulse_zpk([real(dp) ::],[-1.0_dp],1.0_dp,grid)
    call check(res%status==continuous_sim_success .and. pulse%status==continuous_sim_success, &
        'ZPK response status')
    call check(maxval(abs(res%y(:,1,1)-(1.0_dp-exp(-grid))))<tol,'ZPK step')
    call check(maxval(abs(pulse%y(:,1,1)-exp(-grid)))<tol,'ZPK impulse')
    sim = lsim_zpk([real(dp) ::],[-1.0_dp],1.0_dp,grid,grid)
    call check(sim%status==continuous_sim_success,'ZPK simulation')
    call check(maxval(abs(sim%y(:,1)-(grid-1.0_dp+exp(-grid))))<tol,'ZPK ramp')

    ! Complex arithmetic, including the impulse's complex state jump.
    ac = cmplx(-1.0_dp,0.25_dp,dp)
    bc = cmplx(0.7_dp,0.3_dp,dp)
    cc = cmplx(1.0_dp,-0.4_dp,dp)
    dc = cmplx(0.1_dp,0.1_dp,dp)
    uc(:,1) = cmplx(grid,0.5_dp*grid,dp)
    sim = lsim(ac,bc,cc,dc,uc,grid)
    call check(sim%status==continuous_sim_success,'complex FOH status')
    call check(maxval(abs(sim%x(:,1)-bc(1,1)*(1.0_dp+cmplx(0.0_dp,0.5_dp,dp))* &
        (exp(ac(1,1)*grid)-1.0_dp-ac(1,1)*grid)/(ac(1,1)**2)))<tol, &
        'complex FOH ramp state')
    pulse = impulse(ac,bc,cc,dc,grid)
    call check(pulse%status==continuous_sim_success,'complex impulse status')
    call check(maxval(abs(pulse%y(:,1,1)-cc(1,1)*bc(1,1)*exp(ac(1,1)*grid)))<tol, &
        'complex impulse result')

    ! Empty states are legal: continuous static gain is y = D*u.
    empty_d = 3.0_dp
    sim = lsim(empty_a,empty_b,empty_c,empty_d,reshape(grid,[5,1]),grid)
    call check(sim%status==continuous_sim_success,'zero-dimensional system')
    call check(size(sim%x,2)==0 .and. maxval(abs(sim%y(:,1)-3.0_dp*grid))<tol, &
        'zero-dimensional output')
    res = step(empty_a,empty_b,empty_c,empty_d,grid)
    pulse = impulse(empty_a,empty_b,empty_c,empty_d,grid)
    call check(res%status==continuous_sim_success .and. pulse%status==continuous_sim_success, &
        'zero-dimensional responses')
    call check(maxval(abs(res%y(:,1,1)-3.0_dp))<tol .and. &
        maxval(abs(pulse%y(:,1,1)))<tol,'static step and smooth impulse')

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    sim = lsim(a,b,c,d,reshape([0.0_dp,1.0_dp,2.0_dp],[3,1]), &
        [0.0_dp,0.2_dp,0.5_dp])
    call check(sim%status==continuous_sim_invalid_input,'nonuniform time rejected')
    sim = lsim(a,b,c,d,reshape([0.0_dp,1.0_dp],[2,1]),[-0.2_dp,0.0_dp])
    call check(sim%status==continuous_sim_invalid_input,'negative time rejected')
    sim = lsim(a,b,c,d,reshape([0.0_dp,1.0_dp],[2,1]),[0.0_dp,0.5_dp],x0=[1.0_dp,2.0_dp])
    call check(sim%status==continuous_sim_invalid_input,'wrong initial-state shape')
    sim = lsim(a,b,c,d,reshape([0.0_dp,1.0_dp],[2,1]),[0.0_dp,0.0_dp])
    call check(sim%status==continuous_sim_invalid_input,'duplicate times rejected')
    res = step(a,b,c,d,[0.0_dp,nan_value])
    call check(res%status==continuous_sim_invalid_input,'NaN response times rejected')
    res = impulse([1.0_dp,3.0_dp,4.0_dp],[1.0_dp,0.5_dp],grid)
    call check(res%status==continuous_sim_invalid_input,'improper TF rejected')

    print '(a)', 'PASS: continuous-time LTI simulation and responses'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion
        character(*), intent(in) :: label !! assertion description
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_continuous_simulation
