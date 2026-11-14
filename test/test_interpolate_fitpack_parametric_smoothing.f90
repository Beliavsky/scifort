! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_fitpack_parametric_smoothing
    use scifort_kinds, only : dp
    use scifort_interpolate, only : splprep,splev,fitpack_parametric_tck_type, &
        fitpack_tck_invalid_input,fitpack_tck_singular,fitpack_tck_unsupported_smoothing
    use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: n=61
    real(dp) :: u(n),x(n,3),w(n),bad_u(n),bad_x(n,3),bad_w(n),many(n,12)
    real(dp) :: target,fp,nan,inf
    real(dp), allocatable :: vals(:,:),other(:,:),point(:),parameters(:)
    type(fitpack_parametric_tck_type) :: fit,reference
    integer :: i,j,k,status,checks
    checks=0
    do i=1,n
        u(i)=real(i-1,dp)/real(n-1,dp)
        x(i,:)=[sin(7*u(i))+0.12_dp*cos(79*u(i)), &
            cos(5*u(i))+0.14_dp*sin(67*u(i)),u(i)+0.09_dp*cos(53*u(i))]
        w(i)=0.6_dp+u(i)
    end do
    target=0.35_dp
    do k=1,5
        fit=splprep(x,u=u,w=w,k=k,s=target,status=status)
        call require(status==0 .and. fit%valid(),'weighted shared-knot fit')
        call require(size(fit%t)>2*(k+1),'automatic interior knots')
        call require(size(fit%c,1)==size(fit%t)-k-1,'unpadded coefficient count')
        call require(size(fit%c,2)==3 .and. fit%channels()==3,'channel count')
        call require(maxval(abs(fit%u-u))==0.0_dp,'stored explicit parameters')
        vals=splev(u,fit,status=status)
        fp=sum(spread(w**2,2,3)*(vals-x)**2)
        call require(status==0 .and. abs(fp-target)<0.0011_dp*target,'joint residual')
        vals=splev([0.31_dp,0.63_dp],fit)
        point=splev(0.31_dp,fit)
        call require(maxval(abs(vals(1,:)-point))<1e-13_dp,'scalar vector agreement')
        other=splev([0.31_dp,0.63_dp],fit,der=1,status=status)
        call require(status==0 .and. all(abs(other)<1e8_dp),'finite derivatives')
        point=splev(-1.0_dp,fit,ext=1,status=status)
        call require(status==0 .and. maxval(abs(point))==0.0_dp,'outside zero policy')
        fit=splprep(x,u=u,w=w,k=k,s=1000.0_dp,status=status)
        call require(status==0 .and. size(fit%t)==2*(k+1),'polynomial limit')
        fit=splprep(x(:k+1,:),u=u(:k+1),k=k,s=1.0_dp,status=status)
        call require(status==0 .and. fit%valid(),'minimal sample count')
        call require(maxval(abs(splev(u(:k+1),fit)-x(:k+1,:)))<1e-11_dp, &
            'minimal polynomial reproduction')
        fit=splprep(x,w=w,k=k,s=target,status=status)
        call require(status==0 .and. fit%valid(),'automatic chord parameters')
        parameters=fit%u
        reference=splprep(x,u=parameters,w=w,k=k,s=target,status=status)
        call require(status==0 .and. reference%valid(),'explicit chord parameters')
        vals=splev(parameters,fit)
        other=splev(parameters,reference)
        call require(maxval(abs(vals-other))<1e-12_dp,'parameterization agreement')
        fit=splprep(x,u=u,w=w,k=k,s=target,ub=-0.2_dp,ue=1.3_dp,status=status)
        call require(status==0 .and. fit%valid(),'enlarged parameter bounds')
        call require(fit%t(1)==-0.2_dp .and. fit%t(size(fit%t))==1.3_dp,'stored bounds')
        vals=splev(u,fit)
        fp=sum(spread(w**2,2,3)*(vals-x)**2)
        call require(abs(fp-target)<0.0011_dp*target,'enlarged bounds residual')
    end do
    fit=splprep(x*1e100_dp,u=u,w=w*1e-100_dp,s=target,status=status)
    call require(status==0 .and. fit%valid(),'balanced coordinate and weight scales')
    vals=splev(u,fit)/1e100_dp
    call require(abs(sum(spread(w**2,2,3)*(vals-x)**2)-target)<0.0011_dp*target, &
        'balanced scales residual')
    reference=splprep(x,u=u,w=w,s=target,status=status)
    fit=splprep(x*1e-100_dp,u=u,w=w,s=target*1e-200_dp,status=status)
    call require(status==0 .and. fit%valid(),'small coordinate scale')
    vals=splev(u,fit)/1e-100_dp
    other=splev(u,reference)
    call require(maxval(abs(vals-other))<1e-9_dp,'small scale agreement')
    fit=splprep(x*1e200_dp,w=w*1e-200_dp,s=target,status=status)
    call require(status==0 .and. fit%valid(),'large automatic chord lengths')
    vals=splev(fit%u,fit)/1e200_dp
    call require(abs(sum(spread(w**2,2,3)*(vals-x)**2)-target)<0.0011_dp*target, &
        'large automatic residual')
    do i=1,n
        bad_x(i,:)=[(-1.0_dp)**i,0.5_dp,0.0_dp]*1e308_dp
    end do
    fit=splprep(bad_x,w=spread(1e-308_dp,1,n),s=1000.0_dp,status=status)
    call require(status==0 .and. fit%valid(),'opposite huge coordinates avoid chord overflow')
    do i=1,n
        bad_x(i,:)=[1e16_dp+2.0_dp*real(i-1,dp), &
            1e16_dp+2.0_dp*real(mod(i,3),dp),0.0_dp]
    end do
    bad_u(1)=0.0_dp
    do i=2,n
        bad_u(i)=bad_u(i-1)+sqrt(sum((bad_x(i,:)-bad_x(i-1,:))**2))
    end do
    bad_u=bad_u/bad_u(n)
    fit=splprep(bad_x,s=1e5_dp,status=status)
    call require(status==0 .and. fit%valid(),'large-offset coordinate fit')
    call require(maxval(abs(fit%u-bad_u))<1e-13_dp,'large-offset chord geometry')
    fit=splprep(x,u=u*1e200_dp,w=w,s=target,status=status)
    call require(status==0 .and. fit%valid(),'large parameter range')
    bad_u=1e16_dp+2.0_dp*[(real(i-1,dp),i=1,n)]
    fit=splprep(x,u=bad_u,w=w,s=target,status=status)
    call require(status==0 .and. fit%valid(),'large-offset parameters')
    vals=splev(bad_u,fit)
    call require(abs(sum(spread(w**2,2,3)*(vals-x)**2)-target)<0.0011_dp*target, &
        'large-offset physical residual')
    do j=1,12
        many(:,j)=x(:,1)/real(j,dp)
    end do
    fit=splprep(many,u=u,s=target,status=status)
    call require(status==0 .and. fit%channels()==12,'preserved higher channel support')
    vals=splev(u,fit)
    call require(abs(sum((vals-many)**2)-target)<0.0011_dp*target,'higher channel residual')
    bad_x(:,1)=x(:,1)
    bad_x(:,2)=0.0_dp
    bad_x(:,3)=2.0_dp
    fit=splprep(bad_x,u=u,s=target,status=status)
    call require(status==0 .and. fit%valid(),'constant and zero channels')
    vals=splev(u,fit)
    call require(maxval(abs(vals(:,2)))==0.0_dp .and. &
        maxval(abs(vals(:,3)-2.0_dp))<1e-12_dp,'constant channel reproduction')
    fit=splprep(0*x,u=u,s=target,status=status)
    call require(status==0 .and. fit%valid(),'zero curve with explicit parameters')
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    fit=splprep(x,s=nan,status=status)
    call invalid('NaN target')
    fit=splprep(x,s=inf,status=status)
    call invalid('infinite target')
    fit=splprep(x,s=-1.0_dp,status=status)
    call invalid('negative target')
    bad_w=w
    bad_w(3)=nan
    fit=splprep(x,w=bad_w,s=target,status=status)
    call invalid('NaN weight')
    bad_w(3)=0
    fit=splprep(x,w=bad_w,s=target,status=status)
    call invalid('zero weight')
    bad_x=x
    bad_x(3,2)=inf
    fit=splprep(bad_x,s=target,status=status)
    call invalid('infinite coordinate')
    bad_u=u
    bad_u(3)=nan
    fit=splprep(x,u=bad_u,s=target,status=status)
    call invalid('NaN parameter')
    bad_u=u
    bad_u(3)=bad_u(2)
    fit=splprep(x,u=bad_u,s=target,status=status)
    call invalid('duplicate parameter')
    bad_x=x
    bad_x(3,:)=bad_x(2,:)
    fit=splprep(bad_x,s=target,status=status)
    call invalid('zero chord')
    fit=splprep(x,u=u,s=target,ub=nan,status=status)
    call invalid('NaN bound')
    fit=splprep(x,u=u,s=target,ue=0.9_dp,status=status)
    call invalid('bound excludes sites')
    fit=splprep(x,u=u,s=target,t=[0.4_dp],status=status)
    call invalid('knots with automatic task')
    fit=splprep(x,u=u,s=target,per=.true.,status=status)
    call require(status==fitpack_tck_unsupported_smoothing .and. .not.fit%valid(), &
        'periodic behavior preserved')
    fit=splprep(x,u=u,s=1e-100_dp,status=status)
    call numerical('unattainable residual')
    bad_w=w
    bad_w(2)=1e-100_dp
    fit=splprep(x,u=u,w=bad_w,s=target,status=status)
    call numerical('unresolved relative weight')
    bad_u=u
    bad_u(2)=1e-100_dp
    fit=splprep(x,u=bad_u,s=target,k=5,status=status)
    call numerical('unresolved parameter gap')
    fit=splprep(x,u=u*1e-320_dp,s=target,status=status)
    call numerical('unrepresentable physical knot intervals')
    fit=splprep(x,u=(2*u-1)*1e308_dp,s=target,status=status)
    call numerical('unrepresentable physical parameter span')
    print '(a,i0,a)','Adaptive parametric FITPACK smoothing: ',checks,' checks passed.'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Condition required by this regression.
        character(len=*), intent(in) :: label !! Description reported on failure.
        if (.not.ok) then
            print '(a,a)','FAIL: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine require
    subroutine invalid(label)
        implicit none
        character(len=*), intent(in) :: label !! Invalid-input case description.
        call require(status==fitpack_tck_invalid_input .and. .not.fit%valid(),label)
    end subroutine invalid
    subroutine numerical(label)
        implicit none
        character(len=*), intent(in) :: label !! Numerical-range case description.
        call require(status==fitpack_tck_singular .and. .not.fit%valid(),label)
    end subroutine numerical
end program test_interpolate_fitpack_parametric_smoothing
