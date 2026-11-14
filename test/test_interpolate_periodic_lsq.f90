! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_periodic_lsq
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splrep,splprep,splev,bspline, &
        bspline_interpolator_type,fitpack_tck_type,fitpack_parametric_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: m=81
    real(dp) :: u(m),y(m),w(m),x(m,12),nan,inf,scale
    real(dp), allocatable :: full(:),input(:),interior(:),c(:),values(:,:),bad(:)
    type(bspline_interpolator_type) :: model
    type(fitpack_tck_type) :: fit,other
    type(fitpack_parametric_tck_type) :: curve,again
    integer :: i,j,k,ni,n,code,d,checks
    checks=0
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
    w=0.6_dp+u
    w(m)=inf
    do k=1,5
        do ni=2,12
            n=ni+2*k+1
            full=[(real(i-k-1,dp)/real(ni,dp),i=1,n)]
            c=[(sin(real(modulo(i-1,ni),dp)+0.3_dp),i=1,ni+k)]
            interior=full(k+2:n-k-1)
            model=bspline(full,c,k)
            y=model%evaluate(u)
            y(m)=nan
            do j=1,12
                x(:,j)=y*real(j,dp)
            end do
            fit=splrep(u,y,w=w,k=k,task=-1,t=interior,per=.true.,status=code)
            call require(code==0 .and. fit%valid(),'scalar periodic LSQ')
            call require(maxval(abs(fit%c(:ni+k)-c))<1e-11_dp,'recover known cyclic coefficients')
            call require(size(fit%c)==n .and. all(fit%c(ni+k+1:)==0.0_dp),'padded scalar tuple')
            call require(maxval(abs(fit%t-full))<2e-15_dp,'canonical periodic full knots')
            call require(maxval(abs(splev(u(:m-1),fit)-y(:m-1)))<1e-12_dp,'exact curve residual')
            do d=0,k-1
                scale=max(1.0_dp,maxval(abs(splev(u,fit,der=d))))
                call require(abs(splev(0.0_dp,fit,der=d)-splev(1.0_dp,fit,der=d)) < &
                    1e-10_dp*scale,'unwrapped endpoint derivatives')
            end do
            input=[(0.0_dp,i=1,k+1),interior,(1.0_dp,i=1,k+1)]
            curve=splprep(x,u=u,w=w,k=k,task=-1,t=input,per=.true.,status=code)
            call require(code==0 .and. curve%valid(),'shared periodic LSQ')
            call require(all(curve%u==u) .and. size(curve%c,1)==ni+k,'unpadded curve tuple')
            do j=1,12
                call require(maxval(abs(curve%c(:,j)-real(j,dp)*c))<2e-10_dp, &
                    'recover twelve-channel coefficients')
            end do
            input(:k+1)=-123.0_dp
            input(n-k:)=456.0_dp
            again=splprep(x,u=u,w=w,k=k,s=1000.0_dp,task=-1,t=input, &
                per=.true.,status=code)
            call require(code==0 .and. all(again%t==curve%t),'full exterior template canonicalized')
            call require(maxval(abs(again%c-curve%c))==0.0_dp,'s and exterior do not change LSQ')
            other=splrep(u,y,w=w,k=k,s=1000.0_dp,task=-1,t=interior,per=.true.,status=code)
            call require(code==0 .and. maxval(abs(other%c-fit%c))==0.0_dp,'s does not affect LSQ')
            call require(all(input(:k+1)==-123.0_dp),'input knots not mutated')
        end do
    end do
    y=abs(u-0.5_dp)
    x=spread(y,2,12)
    interior=[0.001_dp,0.002_dp,0.5_dp]
    fit=splrep(u,y,k=1,task=-1,t=interior,per=.true.,status=code)
    call require(code==0 .and. fit%valid(),'rank-deficient minimum-norm LSQ')
    call require(maxval(abs(splev(u,fit)-y))<1e-12_dp,'rank-deficient exact residual')
    call require(abs(fit%c(2))<1e-12_dp,'minimum-norm unobserved coefficient')
    bad=[0.5_dp,0.5_dp]
    fit=splrep(u,y,task=-1,t=bad,per=.true.,status=code)
    call require(code==1 .and. .not.fit%valid(),'repeated interior rejected')
    fit=splrep(u,y,task=-1,t=[0.0_dp],per=.true.,status=code)
    call require(code==1,'boundary interior knot rejected')
    fit=splrep(u,y,task=-1,t=[nan],per=.true.,status=code)
    call require(code==1,'nonfinite interior rejected')
    fit=splrep(u,y,task=-1,t=[real(dp)::],per=.true.,status=code)
    call require(code==1,'empty interior matches FITPACK rejection')
    fit=splrep(u,y,task=-1,per=.true.,status=code)
    call require(code==1,'missing fixed knots rejected')
    fit=splrep(u,y,t=[0.5_dp],per=.true.,status=code)
    call require(code==1,'task zero rejects supplied knots')
    input=[(0.0_dp,i=1,4),0.5_dp,(1.0_dp,i=1,4)]
    curve=splprep(x,u=u,task=-1,t=input,per=.true.,status=code)
    call require(code==0,'periodic fixed full knots accepted')
    curve=splprep(x,u=u,task=-1,t=input(:8),per=.true.,status=code)
    call require(code==1 .and. .not.curve%valid(),'full knots too short')
    input(5)=nan
    curve=splprep(x,u=u,task=-1,t=input,per=.true.,status=code)
    call require(code==1 .and. .not.allocated(curve%u),'failure clears parameters')
    do i=1,m
        x(i,1)=cos(2*acos(-1.0_dp)*u(i))
        x(i,2)=sin(2*acos(-1.0_dp)*u(i))
    end do
    input=[(0.0_dp,i=1,4),0.3_dp,0.6_dp,(1.0_dp,i=1,4)]
    curve=splprep(x(:,:2),task=-1,t=input,per=.true.,status=code)
    call require(code==0 .and. curve%valid(),'automatic closed chord parameters')
    again=splprep(x(:,:2),u=curve%u,task=-1,t=input,per=.true.,status=code)
    call require(code==0 .and. maxval(abs(curve%c-again%c))<1e-13_dp,'chord/supplied agreement')
    fit=splrep(u,y*1e150_dp,w=[(1e-150_dp,i=1,m)],k=1,task=-1, &
        t=interior,per=.true.,status=code)
    call require(code==0,'balanced coordinate/weight scales')
    call require(maxval(abs(splev(u,fit)/1e150_dp-y))<1e-12_dp,'scaled LSQ values')
    fit=splrep(u*1e-100_dp,y,k=1,task=-1,t=interior*1e-100_dp,per=.true.,status=code)
    call require(code==0,'tiny physical period')
    fit=splrep(u*1e100_dp,y,k=1,task=-1,t=interior*1e100_dp,per=.true.,status=code)
    call require(code==0,'huge physical period')
    block
        real(dp) :: many_u(300),many_y(300),many_knots(257),many_x(m,257)
        many_u=[(real(i-1,dp)/299,i=1,300)]
        many_y=cos(2*acos(-1.0_dp)*many_u)
        many_knots=[(real(i,dp)/258,i=1,257)]
        fit=splrep(many_u,many_y,k=3,task=-1,t=many_knots,per=.true.,status=code)
        call require(code==3 .and. .not.fit%valid(),'independent coefficient limit')
        many_x=0.0_dp
        curve=splprep(many_x,u=u,k=3,task=-1, &
            t=[(0.0_dp,i=1,4),0.5_dp,(1.0_dp,i=1,4)],per=.true.,status=code)
        call require(code==3 .and. .not.curve%valid(),'coordinate channel limit')
    end block
    print '(a,i0,a)','Periodic fixed-knot LSQ: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required assertion.
        character(*), intent(in) :: label !! Failure description.
        checks=checks+1
        if (.not.ok) then
            print '(a)','FAIL: '//label
            error stop 1
        end if
    end subroutine require
end program test_interpolate_periodic_lsq
