! SPDX-License-Identifier: MIT
program test_linalg_riccati
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_not_positive_definite, linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_riccati, only : solve_continuous_are, solve_discrete_are
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(1,1), b(1,1), q(1,1), r(1,1), x(1,1), expected, err
    real(dp) :: a2(2,2), b2(2,2), q2(2,2), r2(2,2), x2(2,2), e2(2,2), s2(2,2)
    real(dp) :: q_bad(2,2), x_bad(3,3), x_large(21,21), a_large(21,21)
    real(dp) :: a20(20,20), b20(20,20), q20(20,20), r20(20,20), x20(20,20)
    complex(dp) :: ac(2,2), bc(2,2), qc(2,2), rc(2,2), xc(2,2), ec(2,2), sc(2,2)
    integer :: status, steps, i

    a(1,1)=-1.0_dp
    b(1,1)=1.0_dp
    q(1,1)=1.0_dp
    r(1,1)=1.0_dp
    call solve_continuous_are(a,b,q,r,x,status,iterations=steps,residual=err)
    expected=sqrt(2.0_dp)-1.0_dp
    if (status/=linalg_status_success .or. steps<1 .or. err>1.0e-12_dp) &
        error stop 'stable scalar CARE status'
    if (abs(x(1,1)-expected)>1.0e-11_dp) error stop 'stable scalar CARE solution'
    a=1.0_dp
    call solve_continuous_are(a,b,q,r,x,status)
    if (status/=linalg_status_success .or. abs(x(1,1)-(1.0_dp+sqrt(2.0_dp)))>1.0e-11_dp) &
        error stop 'unstable scalar CARE'
    q=0.0_dp
    call solve_continuous_are(a,b,q,r,x,status)
    if (status/=linalg_status_success .or. abs(x(1,1)-2.0_dp)>1.0e-11_dp) &
        error stop 'singular Q: stable graph from top sign block'
    a=1.2_dp
    call solve_discrete_are(a,b,q,r,x,status)
    if (status/=linalg_status_no_convergence) error stop 'nonstabilizing zero-cost DARE accepted'
    q=1.0_dp
    a=0.5_dp
    call solve_discrete_are(a,b,q,r,x,status,iterations=steps,residual=err)
    expected=(0.25_dp+sqrt(4.0625_dp))/2.0_dp
    if (status/=linalg_status_success .or. steps<1 .or. err>1.0e-12_dp) &
        error stop 'stable scalar DARE status'
    if (abs(x(1,1)-expected)>1.0e-11_dp) error stop 'stable scalar DARE solution'
    a=1.25_dp
    call solve_discrete_are(a,b,q,r,x,status)
    if (status/=linalg_status_success) error stop 'unstable scalar DARE'
    if (abs(a(1,1)-x(1,1)/(r(1,1)+x(1,1))*a(1,1))>=1.0_dp) &
        error stop 'unstable scalar DARE closed loop'

    ! Explicit descriptor E and cross cost S, real overloads.
    a2=reshape([0.3_dp,0.2_dp,-0.1_dp,0.65_dp],[2,2])
    b2=reshape([1.0_dp,0.2_dp,0.4_dp,0.8_dp],[2,2])
    q2=0.0_dp
    q2(1,1)=1.5_dp
    q2(2,2)=2.0_dp
    r2=0.0_dp
    r2(1,1)=1.1_dp
    r2(2,2)=1.5_dp
    e2=reshape([1.1_dp,0.04_dp,-0.02_dp,0.9_dp],[2,2])
    s2=reshape([0.01_dp,0.02_dp,0.03_dp,-0.01_dp],[2,2])
    call solve_continuous_are(a2,b2,q2,r2,x2,status,e=e2,s=s2,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-10_dp) &
        error stop 'generalized real CARE'
    if (maxval(abs(x2-transpose(x2)))>1.0e-11_dp) error stop 'real CARE symmetry'
    call solve_discrete_are(a2,b2,q2,r2,x2,status,e=e2,s=s2,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-10_dp) &
        error stop 'generalized real DARE'

    ac=cmplx(a2,0.0_dp,dp)
    ac(1,2)=ac(1,2)+cmplx(0.0_dp,0.2_dp,dp)
    ac(2,1)=ac(2,1)-cmplx(0.0_dp,0.15_dp,dp)
    bc=cmplx(b2,0.0_dp,dp)
    bc(2,2)=bc(2,2)+cmplx(0.0_dp,0.25_dp,dp)
    qc=cmplx(q2,0.0_dp,dp)
    rc=cmplx(r2,0.0_dp,dp)
    ec=cmplx(e2,0.0_dp,dp)
    ec(1,2)=ec(1,2)+cmplx(0.0_dp,0.02_dp,dp)
    sc=cmplx(s2,0.0_dp,dp)
    sc(1,2)=sc(1,2)+cmplx(0.0_dp,0.01_dp,dp)
    call solve_continuous_are(ac,bc,qc,rc,xc,status,e=ec,s=sc,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-10_dp) &
        error stop 'generalized complex CARE'
    if (maxval(abs(xc-transpose(conjg(xc))))>1.0e-11_dp) error stop 'complex CARE Hermitian'
    call solve_discrete_are(ac,bc,qc,rc,xc,status,e=ec,s=sc,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-10_dp) &
        error stop 'generalized complex DARE'
    ! The residual diagnostic independently evaluates the original cross term.
    if (maxval(abs(xc-transpose(conjg(xc))))>1.0e-11_dp) error stop 'complex DARE Hermitian'

    q_bad=q2
    q_bad(1,2)=0.2_dp
    call solve_discrete_are(a2,b2,q_bad,r2,x2,status)
    if (status/=linalg_status_invalid_shape) error stop 'non-Hermitian Q accepted'
    r2(1,1)=-1.0_dp
    call solve_continuous_are(a2,b2,q2,r2,x2,status)
    if (status/=linalg_status_not_positive_definite) error stop 'indefinite R accepted'
    r2(1,1)=1.1_dp
    e2=0.0_dp
    call solve_discrete_are(a2,b2,q2,r2,x2,status,e=e2)
    if (status/=linalg_status_singular) error stop 'singular E accepted'
    call solve_continuous_are(a2,b2,q2,r2,x_bad,status)
    if (status/=linalg_status_invalid_shape) error stop 'bad output shape accepted'
    call solve_discrete_are(a2,b2,q2,r2,x2,status,max_iter=0)
    if (status/=linalg_status_invalid_shape) error stop 'zero iteration budget accepted'
    call solve_continuous_are(a2,b2,q2,r2,x2,status,max_iter=1)
    if (status/=linalg_status_no_convergence) error stop 'false convergence with one sign iteration'
    call solve_discrete_are(a2,b2,q2,r2,x2,status,max_iter=1)
    if (status/=linalg_status_no_convergence) error stop 'false convergence with one doubling iteration'
    a2(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call solve_discrete_are(a2,b2,q2,r2,x2,status)
    if (status/=linalg_status_no_convergence) error stop 'NaN accepted'

    a_large=0.0_dp
    do i=1,21
        a_large(i,i)=-1.0_dp
    end do
    call solve_continuous_are(a_large,a_large,a_large,a_large,x_large,status)
    if (status/=linalg_status_invalid_shape) error stop 'large ARE not rejected'

    a20=0.0_dp
    b20=0.0_dp
    q20=0.0_dp
    r20=0.0_dp
    do i=1,20
        a20(i,i)=-0.3_dp-0.05_dp*real(i,dp)
        b20(i,i)=1.0_dp
        q20(i,i)=1.0_dp
        r20(i,i)=1.0_dp
    end do
    call solve_continuous_are(a20,b20,q20,r20,x20,status,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-11_dp) &
        error stop '20-dimensional continuous ARE'
    do i=1,20
        a20(i,i)=0.55_dp+0.012_dp*real(i,dp)
    end do
    call solve_discrete_are(a20,b20,q20,r20,x20,status,residual=err)
    if (status/=linalg_status_success .or. err>1.0e-11_dp) &
        error stop '20-dimensional discrete ARE'

    print '(a)', 'dense Riccati native tests passed'
end program test_linalg_riccati
