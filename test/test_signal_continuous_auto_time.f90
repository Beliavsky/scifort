! SPDX-License-Identifier: MIT
! Check SciPy-style automatic response grids while retaining explicit grids.
program test_signal_continuous_auto_time
    use scifort_signal, only : dp, step, impulse, step_zpk, impulse_zpk, &
        continuous_response_result, continuous_sim_success, continuous_sim_invalid_input
    implicit none
    real(dp), parameter :: tol = 3.0e-8_dp
    type(continuous_response_result) :: resp, imp, manual
    real(dp) :: a(2,2), b(2,1), c(1,2), d(1,1), t(4)
    real(dp) :: empty_a(0,0), empty_b(0,1), empty_c(1,0), empty_d(1,1)
    complex(dp) :: ac(2,2), bc(2,1), cc(1,2), dc(1,1)
    integer :: j

    ! First order transfer function, which has a time constant of 1/2.
    resp = step([1.0_dp],[1.0_dp,2.0_dp])
    imp = impulse([1.0_dp],[1.0_dp,2.0_dp],n=51)
    call check(resp%status==continuous_sim_success,'default TF status')
    call check(imp%status==continuous_sim_success,'custom-count impulse status')
    call check(size(resp%t)==100 .and. size(imp%t)==51,'automatic grid lengths')
    call check(abs(resp%t(1))<tol .and. abs(resp%t(100)-3.5_dp)<tol,'default TF horizon')
    call check(abs(imp%t(51)-3.5_dp)<tol,'custom TF horizon')
    call check(maxval(abs(resp%y(:,1,1)-0.5_dp*(1.0_dp-exp(-2.0_dp*resp%t))))<tol, &
        'default TF step response')
    call check(maxval(abs(imp%y(:,1,1)-exp(-2.0_dp*imp%t)))<tol, &
        'default TF impulse response')

    ! Non-triangular two-state realization uses characteristic roots.
    a(1,:) = [-2.0_dp, 1.0_dp]
    a(2,:) = [-2.0_dp,-3.0_dp]
    b(:,1) = [1.0_dp,2.0_dp]
    c(1,:) = [0.5_dp,1.0_dp]
    d = 0.25_dp
    resp = step(a,b,c,d,n=29)
    call check(resp%status==continuous_sim_success,'general state-space default')
    call check(size(resp%t)==29,'general state-space grid length')
    call check(abs(resp%t(29)-2.8_dp)<tol,'general state-space eigenvalue horizon')
    manual = step(a,b,c,d,resp%t)
    call check(manual%status==continuous_sim_success,'manual-grid comparison status')
    call check(maxval(abs(resp%y-manual%y))<tol,'automatic/manual response equality')

    ! Pure imaginary eigenvalues and exact zero roots select a 7-unit interval.
    a(1,:) = [0.0_dp,-2.0_dp]
    a(2,:) = [2.0_dp,0.0_dp]
    resp = step(a,b,c,d,n=33)
    call check(resp%status==continuous_sim_success,'oscillator status')
    call check(abs(resp%t(33)-7.0_dp)<tol,'oscillator time horizon')
    a = 0.0_dp
    a(2,2) = -4.0_dp
    imp = impulse(a,b,c,d,n=15)
    call check(imp%status==continuous_sim_success,'integrator status')
    call check(abs(imp%t(15)-7.0_dp)<tol,'integrator time horizon')

    ! A complex conjugate-pair matrix, passed through the complex interface.
    ac = cmplx(a,0.0_dp,dp)
    ac(1,1) = cmplx(-0.5_dp,2.0_dp,dp)
    ac(2,2) = cmplx(-2.0_dp,-1.0_dp,dp)
    bc = cmplx(b,0.0_dp,dp)
    cc = cmplx(c,0.0_dp,dp)
    dc = cmplx(d,0.0_dp,dp)
    resp = step(ac,bc,cc,dc,n=70)
    call check(resp%status==continuous_sim_success,'complex state-space status')
    call check(abs(resp%t(70)-14.0_dp)<tol,'complex slow pole horizon')

    ! Static systems and a one-point request remain well-defined.
    empty_d = 3.0_dp
    resp = step(empty_a,empty_b,empty_c,empty_d,n=1)
    call check(resp%status==continuous_sim_success,'static one-point status')
    call check(size(resp%t)==1 .and. abs(resp%t(1))<tol,'single sample at zero')
    call check(abs(resp%y(1,1,1)-3.0_dp)<tol,'static one-point output')
    resp = step([3.0_dp],[1.0_dp],n=12)
    call check(resp%status==continuous_sim_success,'constant transfer function')
    call check(abs(resp%t(12)-7.0_dp)<tol,'constant time horizon')

    ! Explicit irregular grids take precedence over optional n, like SciPy T.
    t = [0.0_dp,0.2_dp,0.8_dp,1.7_dp]
    resp = step([1.0_dp],[1.0_dp,2.0_dp],t,n=0)
    call check(resp%status==continuous_sim_success,'T overrides N')
    call check(all(resp%t==t),'T is returned unmodified')
    imp = impulse_zpk([real(dp) ::],[-2.0_dp],1.0_dp,n=13)
    call check(imp%status==continuous_sim_success,'ZPK auto response')
    call check(abs(imp%t(13)-3.5_dp)<tol,'ZPK time horizon')
    resp = step_zpk([real(dp) ::],[-2.0_dp],1.0_dp,n=13)
    call check(resp%status==continuous_sim_success,'ZPK step response')
    call check(maxval(abs(resp%y(:,1,1)-0.5_dp*(1.0_dp-exp(-2.0_dp*resp%t))))<tol, &
        'ZPK automatic response values')

    ! Incorrect N and T arguments return an invalid status, not a crash.
    resp = step([1.0_dp],[1.0_dp,2.0_dp],n=0)
    call check(resp%status==continuous_sim_invalid_input,'zero sample count')
    resp = impulse([1.0_dp],[1.0_dp,2.0_dp],n=-5)
    call check(resp%status==continuous_sim_invalid_input,'negative sample count')
    resp = step([1.0_dp],[1.0_dp,2.0_dp],t=[0.0_dp,0.0_dp],n=100)
    call check(resp%status==continuous_sim_invalid_input,'bad explicit grid')

    do j=1,2
        if (j==1) then
            resp = step([1.0_dp],[1.0_dp,2.0_dp],n=200)
        else
            resp = impulse([1.0_dp],[1.0_dp,2.0_dp],n=200)
        end if
        call check(resp%status==continuous_sim_success,'longer automatic grid status')
        call check(all(resp%t(2:)>resp%t(:199)),'automatic grid monotonicity')
    end do
    print '(a)', 'PASS: automatic continuous-time response grids'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: label !! assertion description
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_continuous_auto_time
