! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_fitpack_closed
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splprep,splev,splrep, &
        fitpack_parametric_tck_type,fitpack_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: n=61,ndim=3
    real(dp) :: x(n,ndim),copy(n,ndim),bad(n),u(n),w(n),theta,fp,nan,inf
    real(dp) :: many(n,12),a(ndim),b(ndim),chords(n),gap,tolerance,constant(ndim)
    real(dp), allocatable :: vals(:,:),single(:,:),v(:)
    type(fitpack_parametric_tck_type) :: fit,other
    type(fitpack_tck_type) :: scalar
    integer :: i,j,k,d,status,checks
    checks=0
    do i=1,n
        u(i)=real(i-1,dp)/real(n-1,dp)
        theta=2*acos(-1.0_dp)*u(i)
        x(i,:)=[cos(theta)+0.06_dp*sin(13*theta),2*sin(theta),0.2_dp*cos(3*theta)]
        w(i)=0.5_dp+u(i)
    end do
    x(n,:)=1234.0_dp
    copy=x
    do k=1,5
        fit=splprep(x,u=u,w=w,k=k,per=.true.,status=status)
        call require(status==0 .and. fit%valid(),'periodic interpolation')
        call require(size(fit%c,1)==size(fit%t)-k-1,'unpadded coefficients')
        call require(fit%channels()==ndim,'channel count')
        vals=splev(u,fit)
        call require(maxval(abs(vals(:n-1,:)-x(:n-1,:)))<2e-12_dp,'interpolates active rows')
        call seam(fit)
        fit=splprep(x,u=u,w=w,k=k,s=0.05_dp,per=.true.,status=status)
        call require(status==0 .and. fit%valid(),'periodic smoothing')
        vals=splev(u,fit)
        fp=sum((spread(w(:n-1),2,ndim)*(vals(:n-1,:)-x(:n-1,:)))**2)
        call require(abs(fp-0.05_dp)<0.000055_dp,'one total weighted residual')
        call seam(fit)
        single=splev([0.231_dp],fit)
        v=splev(0.231_dp,fit)
        call require(maxval(abs(v-single(1,:)))<1e-14_dp,'scalar array agreement')
        other=splprep(x,w=w,k=k,s=0.05_dp,per=.true.,status=status)
        call require(status==0 .and. other%valid(),'automatic closed parameters')
        chords(1)=0.0_dp
        do i=2,n-1
            chords(i)=chords(i-1)+norm2(x(i,:)-x(i-1,:))
        end do
        chords(n)=chords(n-1)+norm2(x(1,:)-x(n-1,:))
        chords=chords/chords(n)
        call require(maxval(abs(chords-other%u))<1e-14_dp,'includes actual closing leg')
        call seam(other)
        fit=splprep(x,u=u,w=w,k=k,s=1000.0_dp,per=.true.,status=status)
        call require(status==0 .and. size(fit%t)==2*(k+1),'constant-vector limit')
        do j=1,ndim
            constant(j)=sum(w(:n-1)**2*x(:n-1,j))/sum(w(:n-1)**2)
        end do
        vals=splev(u,fit)
        call require(maxval(abs(vals-spread(constant,1,n)))<1e-13_dp,'weighted constant')
        fit=splprep(x(:k+1,:),u=u(:k+1),k=k,per=.true.,status=status)
        call require(status==0 .and. fit%valid(),'minimal periodic sample count')
    end do
    call require(all(x==copy),'input coordinates unchanged')
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    copy(n,:)=[nan,inf,-inf]
    bad=w
    bad(n)=nan
    fit=splprep(copy,w=bad,s=0.05_dp,per=.true.,status=status)
    call require(status==0 .and. fit%valid(),'ignored nonfinite closing row and weight')
    other=splprep(x,w=w,s=0.05_dp,per=.true.,status=status)
    call require(maxval(abs(fit%c-other%c))==0.0_dp,'ignored entries do not affect fit')
    call require(maxval(abs(fit%u-other%u))==0.0_dp,'ignored entries do not affect chords')
    copy(2,1)=nan
    fit=splprep(copy,u=u,per=.true.,status=status)
    call require(status==1 .and. .not.fit%valid(),'active NaN coordinate')
    bad(2)=inf
    fit=splprep(x,u=u,w=bad,per=.true.,status=status)
    call require(status==1,'active infinite weight')
    bad=w
    bad(2)=0.0_dp
    fit=splprep(x,u=u,w=bad,per=.true.,status=status)
    call require(status==1,'zero active weight')
    bad=u
    bad(n)=inf
    fit=splprep(x,u=bad,per=.true.,status=status)
    call require(status==1,'closing parameter still required finite')
    bad=u
    bad(3)=bad(2)
    fit=splprep(x,u=bad,per=.true.,status=status)
    call require(status==1,'repeated parameter')
    fit=splprep(x,u=u,s=nan,per=.true.,status=status)
    call require(status==1,'NaN target')
    fit=splprep(x,u=u,s=-1.0_dp,per=.true.,status=status)
    call require(status==1,'negative target')
    fit=splprep(x,u=u,k=6,per=.true.,status=status)
    call require(status==1,'invalid degree')
    fit=splprep(x,u=u,ub=-1.0_dp,per=.true.,status=status)
    call require(status==1,'periodic bounds must match parameter endpoints')
    fit=splprep(x,u=u,task=-1,t=[0.0_dp,1.0_dp],per=.true.,status=status)
    call require(status==2,'periodic fixed knots explicitly unsupported')
    fit=splprep(x,u=u,task=1,per=.true.,status=status)
    call require(status==1,'stateless continuation rejected')
    copy=x
    copy(2,:)=copy(1,:)
    fit=splprep(copy,per=.true.,status=status)
    call require(status==1,'zero interior chord rejected')
    copy=x
    copy(n-1,:)=copy(1,:)
    fit=splprep(copy,per=.true.,status=status)
    call require(status==1,'zero closing chord rejected')
    fit=splprep(0*x,u=u,per=.true.,status=status)
    call require(status==0 .and. maxval(abs(splev(u,fit)))==0.0_dp,'zero interpolation')
    fit=splprep(0*x,per=.true.,status=status)
    call require(status==1,'undefined zero-length automatic parameters')
    do j=1,12
        many(:,j)=x(:,mod(j-1,ndim)+1)
    end do
    fit=splprep(many,u=u,s=0.2_dp,per=.true.,status=status)
    call require(status==0 .and. fit%channels()==12,'higher channels extension')
    fit=splprep(x,u=u,w=w,s=0.05_dp,per=.true.,status=status)
    other=splprep(x*1e100_dp,u=u,w=w*1e-100_dp,s=0.05_dp,per=.true.,status=status)
    call require(status==0,'balanced extreme scale')
    vals=splev(u,other)/1e100_dp
    fp=sum((spread(w(:n-1),2,ndim)*(vals(:n-1,:)-x(:n-1,:)))**2)
    call require(abs(fp-0.05_dp)<0.000055_dp,'balanced scale residual')
    other=splprep(x*1e300_dp,u=u,w=w*1e-310_dp,s=0.05e-20_dp, &
        per=.true.,status=status)
    call require(status==0,'subnormal weight scale')
    fit=splprep(x,u=u*1e200_dp,s=0.05_dp,per=.true.,status=status)
    call require(status==0,'large parameter scale')
    fit=splprep(x,u=u*1e-200_dp,s=0.05_dp,per=.true.,status=status)
    call require(status==0,'small parameter scale')
    fit=splprep(x,u=1e16_dp+2*real([(i-1,i=1,n)],dp),s=0.05_dp, &
        per=.true.,status=status)
    call require(status==0,'large parameter offset')
    fit=splprep(x,u=u,s=1e-310_dp,per=.true.,status=status)
    call require(status==3 .and. .not.fit%valid(),'unrepresentable target fails clearly')
    fit=splprep(x(:,1:1),u=u,w=w,s=0.05_dp,per=.true.,status=status)
    scalar=splrep(u,x(:,1),w=w,s=0.05_dp,per=.true.,status=status)
    vals=splev(u,fit)
    call require(maxval(abs(vals(:,1)-splev(u,scalar)))<1e-11_dp,'one channel scalar identity')
    fit=splprep(x,u=u,s=0.05_dp,status=status)
    other=splprep(x,u=u,s=0.05_dp,per=.false.,status=status)
    call require(maxval(abs(fit%c-other%c))==0.0_dp,'nonperiodic behavior preserved')
    print '(a,i0,a)','closed splprep: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required condition.
        character(*), intent(in) :: label !! Failure diagnostic.
        if (.not.ok) then
            print *, 'FAIL: ',label
            stop 1
        end if
        checks=checks+1
    end subroutine require
    subroutine seam(tck)
        implicit none
        type(fitpack_parametric_tck_type), intent(in) :: tck !! Closed curve to check.
        do d=0,tck%k-1
            a=splev(tck%u(1),tck,der=d)
            b=splev(tck%u(size(tck%u)),tck,der=d)
            gap=maxval(abs(a-b)/(1+abs(a)+abs(b)))
            tolerance=1e-9_dp
            if (d==0) tolerance=1e-12_dp
            ! Fourth derivatives near zero amplify knot and coefficient rounding.
            if (d==4) tolerance=3e-8_dp
            call require(gap<tolerance,'periodic seam continuity')
        end do
    end subroutine seam
end program test_interpolate_fitpack_closed
