! SPDX-License-Identifier: MIT
program test_interpolate_global
    use scifort_kinds, only : dp
    use scifort_interpolate, only : krogh_interpolator, krogh_interpolate, &
        floater_hormann_interpolator, krogh_interpolator_type, &
        complex_krogh_interpolator_type, floater_hormann_interpolator_type, &
        complex_floater_hormann_interpolator_type
    implicit none
    type(krogh_interpolator_type) :: k
    type(complex_krogh_interpolator_type) :: kc
    type(floater_hormann_interpolator_type) :: r
    type(complex_floater_hormann_interpolator_type) :: rc
    real(dp) :: x(7), h, y(7,2), testp, dtest, vv(2), q
    complex(dp) :: zy(7,2), zval
    real(dp), allocatable :: v(:), weights(:), derivatives(:,:,:)
    integer :: code, i
    ! f(x)=x^4-3x^2+2x+1 and its confluent Hermite constraints.
    ! Includes up to a third derivative at one node.
    x=[-1.0_dp,-1.0_dp,-1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp]
    y(:,1)=[-3.0_dp,4.0_dp,6.0_dp,1.0_dp,2.0_dp,-6.0_dp,1.0_dp]
    ! f(-1)=1-3-2+1=-3; f'(-1)=-4+6+2=4;
    ! f''(-1)=12-6=6; f(1)=1-3+2+1=1.
    y(:,2)=2.0_dp*y(:,1)
    k=krogh_interpolator(x,y,status=code)
    if (code/=0 .or. .not.k%valid() .or. k%channels()/=2) error stop 'Krogh construction'
    do i=0,35
        q=-1.2_dp+0.08_dp*real(i,dp)
        testp=q**4-3.0_dp*q**2+2.0_dp*q+1.0_dp
        dtest=4.0_dp*q**3-6.0_dp*q+2.0_dp
        if (abs(k%value(q)-testp)>1.0e-11_dp) error stop 'Krogh polynomial value'
        if (abs(k%derivative(q)-dtest)>1.0e-11_dp) error stop 'Krogh polynomial derivative'
        if (abs(k%derivative(q,nu=2)-(12.0_dp*q*q-6.0_dp))>1.0e-10_dp) &
            error stop 'Krogh second derivative'
        if (abs(k%derivative(q,nu=4)-24.0_dp)>1.0e-10_dp) error stop 'Krogh fourth derivative'
        if (abs(k%value(q,channel=2)-2.0_dp*testp)>1.0e-11_dp) error stop 'Krogh channels'
    end do
    if (k%derivative(0.2_dp,nu=7)/=0.0_dp) error stop 'Krogh zero high derivative'
    derivatives=k%derivatives([0.3_dp,0.6_dp],der=9)
    if (any(shape(derivatives)/=[9,2,2])) error stop 'Krogh all derivatives shape'
    if (abs(derivatives(1,1,1)-(0.3_dp**4-3.0_dp*0.3_dp**2+2*0.3_dp+1))>1.0e-11_dp) &
        error stop 'Krogh all derivatives order 0'
    if (abs(derivatives(3,2,2)-2*(12*0.6_dp**2-6))>1.0e-11_dp) &
        error stop 'Krogh all derivatives order 2'
    if (any(derivatives(8:,:,:)/=0.0_dp)) error stop 'Krogh all derivatives high-order zero'
    v=krogh_interpolate(x,y(:,1),[0.3_dp],nu=1,status=code)
    if (code/=0 .or. abs(v(1)-(4*0.3_dp**3-6*0.3_dp+2))>1.0e-11_dp) &
        error stop 'Krogh convenience'
    zy=cmplx(y,3.0_dp*y,dp)
    kc=krogh_interpolator(x,zy,status=code)
    if (code/=0) error stop 'complex Krogh construction'
    zval=kc%derivative(0.4_dp,nu=2,channel=2)
    block
        complex(dp), allocatable :: zd(:,:,:)
        zd=kc%derivatives([0.4_dp],der=4)
        if (abs(zd(3,1,2)-kc%derivative(0.4_dp,nu=2,channel=2))>1.0e-12_dp) &
            error stop 'complex Krogh all derivatives'
    end block
    if (abs(zval-cmplx(2*(12*0.4_dp**2-6),6*(12*0.4_dp**2-6),dp))>1.0e-11_dp) &
        error stop 'complex Krogh derivative'
    k=krogh_interpolator([0.0_dp,1.0_dp,0.0_dp],[1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code==0 .or. k%valid()) error stop 'non-contiguous duplicate nodes'
    ! Floater-Hormann with d=n-1 is exactly the polynomial interpolant.
    r=floater_hormann_interpolator([-1.0_dp,0.0_dp,1.0_dp], &
        [1.0_dp,0.0_dp,1.0_dp],d=2,status=code)
    if (code/=0 .or. r%degree()/=2) error stop 'Floater polynomial construction'
    do i=0,20
        q=-1.1_dp+0.11_dp*real(i,dp)
        if (abs(r%value(q)-q**2)>1.0e-11_dp) error stop 'Floater polynomial limit'
    end do
    r=floater_hormann_interpolator([3.0_dp,2.0_dp,1.0_dp], &
        [9.0_dp,4.0_dp,1.0_dp],d=1,status=code)
    if (code/=0 .or. abs(r%value(2.0_dp)-4.0_dp)>1.0e-13_dp) &
        error stop 'descending Floater exact node'
    weights=r%weights()
    if (size(weights)/=3 .or. abs(maxval(abs(weights))-1.0_dp)>1.0e-14_dp) &
        error stop 'Floater weights'
    rc=floater_hormann_interpolator([-1.0_dp,0.0_dp,1.0_dp], &
        cmplx([1.0_dp,0.0_dp,1.0_dp],[2.0_dp,0.0_dp,2.0_dp],dp),d=2,status=code)
    if (code/=0 .or. abs(rc%value(0.5_dp)-cmplx(0.25_dp,0.5_dp,dp))>1.0e-12_dp) &
        error stop 'complex Floater polynomial'
    r=floater_hormann_interpolator([0.0_dp,0.0_dp],[1.0_dp,2.0_dp],d=1,status=code)
    if (code==0 .or. r%valid()) error stop 'invalid Floater nodes accepted'
    r=floater_hormann_interpolator([0.0_dp,1.0_dp],[1.0_dp,2.0_dp],d=2,status=code)
    if (code==0 .or. r%valid()) error stop 'invalid Floater degree accepted'
    print '(a)', 'Krogh and Floater-Hormann native regression PASS'
end program test_interpolate_global
