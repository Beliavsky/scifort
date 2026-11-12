! SPDX-License-Identifier: MIT
! Discrete feedback, block connections, Kalman matrices and Lyapunov Gramians.
program test_signal_dlti_interconnections
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_ss, dlti_tf, &
        dlti_gramian_result, dlti_success, dlti_invalid_input, dlti_numerical_failure
    use scifort_signal_discrete_simulation, only : discrete_simulation_result
    implicit none
    complex(dp) :: ap(1,1),bp(1,1),cp(1,1),dp1(1,1)
    complex(dp) :: ak(1,1),bk(1,1),ck(1,1),dk(1,1)
    complex(dp) :: xp, xk, u, y, next_xp, next_xk, denom
    complex(dp), allocatable :: ctr(:,:), obs(:,:), ur(:,:), values(:,:)
    real(dp) :: t(30)
    type(dlti_system) :: p,k,closed,positive,combined,selected,invalid,pole,zero
    type(discrete_simulation_result) :: output, o1, o2
    type(dlti_gramian_result) :: g
    integer :: i,s

    ap(1,1) = cmplx(0.55_dp,0.09_dp,dp)
    bp(1,1) = cmplx(0.8_dp,0.1_dp,dp)
    cp(1,1) = cmplx(1.1_dp,-0.1_dp,dp)
    dp1(1,1) = cmplx(0.2_dp,-0.15_dp,dp)
    ak(1,1) = cmplx(-0.15_dp,0.02_dp,dp)
    bk(1,1) = cmplx(0.5_dp,-0.03_dp,dp)
    ck(1,1) = cmplx(0.4_dp,0.12_dp,dp)
    dk(1,1) = cmplx(0.3_dp,-0.04_dp,dp)
    p = dlti_ss(ap,bp,cp,dp1,dt=0.2_dp)
    k = dlti_ss(ak,bk,ck,dk,dt=0.2_dp)
    call check(p%status == dlti_success .and. k%status == dlti_success,'construct')
    do s=-1,1,2
        closed = p%feedback(k,sign=s)
        call check(closed%status == dlti_success,'feedback construction')
        call check(closed%nstates() == 2,'feedback dimensions')
        allocate(ur(30,1))
        do i=1,30
            ur(i,1) = cmplx(sin(0.11_dp*i),cos(0.23_dp*i)*0.4_dp,dp)
            t(i) = real(i-1,dp)*p%dt
        end do
        output = closed%output(ur,t=t)
        call check(output%status == dlti_success,'feedback output status')
        xp = (0.0_dp,0.0_dp)
        xk = (0.0_dp,0.0_dp)
        denom = 1.0_dp - real(s,dp)*dk(1,1)*dp1(1,1)
        do i=1,size(ur,1)
            u = (ur(i,1) + real(s,dp)*(dk(1,1)*cp(1,1)*xp+ck(1,1)*xk))/denom
            y = cp(1,1)*xp+dp1(1,1)*u
            call check(abs(output%y(i,1)-y)<5.0e-12_dp,'closed-loop recurrence')
            next_xp = ap(1,1)*xp+bp(1,1)*u
            next_xk = ak(1,1)*xk+bk(1,1)*y
            xp = next_xp
            xk = next_xk
        end do
        deallocate(ur)
    end do
    invalid = p%feedback(k,sign=0)
    call check(invalid%status == dlti_invalid_input,'feedback bad sign')
    invalid = p%feedback(dlti_ss(ak,bk,ck,dk,dt=0.3_dp))
    call check(invalid%status == dlti_invalid_input,'feedback bad dt')
    dk(1,1) = (5.0_dp,0.0_dp)
    dp1(1,1) = (0.2_dp,0.0_dp)
    pole = dlti_ss(ap,bp,cp,dp1,dt=0.2_dp)
    invalid = pole%feedback(dlti_ss(ak,bk,ck,dk,dt=0.2_dp),sign=1)
    call check(invalid%status == dlti_numerical_failure,'singular algebraic loop')

    combined = p%append(k)
    call check(combined%status == dlti_success,'append status')
    call check(combined%ninputs() == 2 .and. combined%noutputs() == 2,'append channel counts')
    allocate(ur(30,2))
    do i=1,30
        ur(i,1) = cmplx(sin(0.11_dp*i),cos(0.18_dp*i),dp)
        ur(i,2) = cmplx(cos(0.32_dp*i),sin(0.07_dp*i),dp)
    end do
    output = combined%output(ur)
    o1 = p%output(ur(:,1))
    o2 = k%output(ur(:,2))
    call check(maxval(abs(output%y(:,1)-o1%y(:,1))) < 1.0e-12_dp,'append first output')
    call check(maxval(abs(output%y(:,2)-o2%y(:,1))) < 1.0e-12_dp,'append second output')
    selected = combined%select_io([2,1,2],[2,1])
    call check(selected%status == dlti_success,'select channels')
    call check(selected%ninputs() == 3 .and. selected%noutputs() == 2,'selected dimensions')
    call check(maxval(abs(selected%d-combined%d([2,1],[2,1,2]))) < 1.0e-15_dp,'select feedthrough')
    invalid = combined%select_io([0],[1])
    call check(invalid%status == dlti_invalid_input,'select zero input index')
    invalid = combined%select_io([1],[3])
    call check(invalid%status == dlti_invalid_input,'select high output index')
    invalid = combined%append(dlti_tf([1.0_dp],[1.0_dp,-0.5_dp],dt=1.0_dp))
    call check(invalid%status == dlti_invalid_input,'append sampling mismatch')

    ctr = combined%controllability()
    obs = combined%observability()
    call check(all(shape(ctr)==[2,4]),'controllability dimensions')
    call check(all(shape(obs)==[4,2]),'observability dimensions')
    call check(maxval(abs(ctr(:,1:2)-combined%b)) < 1.0e-14_dp,'reachability first block')
    call check(maxval(abs(ctr(:,3:4)-matmul(combined%a,combined%b))) < 1.0e-14_dp, &
        'reachability second block')
    call check(maxval(abs(obs(3:4,:)-matmul(combined%c,combined%a))) < 1.0e-14_dp, &
        'observability second block')
    g = combined%gramians()
    call check(g%status == dlti_success,'gramian status')
    call check(g%iterations > 0,'gramian iterations')
    values = g%wc-matmul(matmul(combined%a,g%wc),transpose(conjg(combined%a)))- &
        matmul(combined%b,transpose(conjg(combined%b)))
    call check(maxval(abs(values)) < 2.0e-12_dp,'controllability Lyapunov residual')
    values = g%wo-matmul(matmul(transpose(conjg(combined%a)),g%wo),combined%a)- &
        matmul(transpose(conjg(combined%c)),combined%c)
    call check(maxval(abs(values)) < 2.0e-12_dp,'observability Lyapunov residual')
    call check(abs(g%h2_norm) > 1.0_dp,'H2 norm positive')
    g = combined%gramians(tol=0.0_dp)
    call check(g%status == dlti_invalid_input,'invalid tolerance')
    g = combined%gramians(max_iter=1)
    call check(g%status == dlti_numerical_failure,'insufficient iterations')
    pole = dlti_ss(reshape([cmplx(1.1_dp,0.0_dp,dp)],[1,1]), &
        bp,cp,dp1)
    g = pole%gramians()
    call check(g%status == dlti_numerical_failure,'unstable gramian')
    zero = dlti_ss(reshape([cmplx(0.0_dp,0.0_dp,dp)],[0,0]), &
        reshape([cmplx(0.0_dp,0.0_dp,dp)],[0,1]), &
        reshape([cmplx(0.0_dp,0.0_dp,dp)],[1,0]), &
        reshape([cmplx(0.4_dp,0.2_dp,dp)],[1,1]))
    call check(zero%status == dlti_success,'zero-state construction')
    g = zero%gramians()
    call check(g%status == dlti_success,'zero-state Gramians')
    call check(abs(g%h2_norm-sqrt(0.2_dp)) < 1.0e-14_dp,'zero-state H2')
    print '(a)', 'PASS: discrete feedback, append, selection, Kalman matrices and Lyapunov Gramians'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! asserted condition
        character(*), intent(in) :: label !! failure description
        if (.not. ok) then
            print '(a,1x,a)', 'FAIL:',label
            error stop 1
        end if
    end subroutine check
end program test_signal_dlti_interconnections
