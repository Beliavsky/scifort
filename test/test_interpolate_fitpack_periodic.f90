! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_fitpack_periodic
    use scifort_kinds, only: dp
    use scifort_interpolate_fitpack_compat, only: splrep,splev,splder,splint, &
        fitpack_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: n=61
    real(dp) :: x(n),y(n),w(n),bad(n),vals(n),pi,fp,a,b,nan,inf
    type(fitpack_tck_type) :: fit,other,deriv
    integer :: i,k,d,code,checks
    checks=0
    pi=acos(-1.0_dp)
    do i=1,n
        x(i)=2*pi*real(i-1,dp)/real(n-1,dp)
        y(i)=cos(x(i))+sin(2*x(i))+0.08_dp*cos(13*x(i))
        w(i)=0.5_dp+real(i,dp)/n
    end do
    do k=1,5
        fit=splrep(x,y,w=w,k=k,per=.true.,status=code)
        call require(code==0 .and. fit%valid(),'periodic interpolation')
        vals=splev(x,fit)
        call require(maxval(abs(vals(:n-1)-y(:n-1)))<1e-12_dp,'interpolates samples')
        call seam(fit,k)
        fit=splrep(x,y,w=w,k=k,s=0.05_dp,per=.true.,status=code)
        call require(code==0 .and. fit%valid(),'periodic adaptive smoothing')
        call require(size(fit%c)==size(fit%t),'padded tuple')
        vals=splev(x,fit)
        fp=sum((w(:n-1)*(vals(:n-1)-y(:n-1)))**2)
        call require(abs(fp-0.05_dp)<0.000055_dp,'weighted target')
        call seam(fit,k)
        deriv=splder(fit,status=code)
        call require(code==0,'derivative construction')
        call require(abs(splev(0.7_dp,deriv)-splev(0.7_dp,fit,der=1))<1e-11_dp, &
            'derivative identity')
        call require(abs(splint(x(1),x(n),fit)-splint(x(1),2.0_dp,fit)- &
            splint(2.0_dp,x(n),fit))<1e-12_dp,'integral additivity')
        bad=y
        bad(n)=123456.0_dp
        other=splrep(x,bad,w=w,k=k,s=0.05_dp,per=.true.,status=code)
        call require(code==0,'ignored closing ordinate')
        call require(maxval(abs(other%c-fit%c))==0.0_dp,'closing ordinate has no effect')
        fit=splrep(x,y,w=w,k=k,s=1000.0_dp,per=.true.,status=code)
        call require(code==0 .and. size(fit%t)==2*(k+1),'constant limit')
        a=sum(w(:n-1)**2*y(:n-1))/sum(w(:n-1)**2)
        call require(maxval(abs(splev(x,fit)-a))<1e-13_dp,'weighted constant')
        fit=splrep(x(:k+1),2+0*y(:k+1),k=k,per=.true.,status=code)
        call require(code==0 .and. abs(splev(0.3_dp,fit)-2.0_dp)<1e-13_dp, &
            'minimal periodic sample count')
    end do
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    bad=y
    bad(n)=nan
    fit=splrep(x,bad,per=.true.,status=code)
    call require(code==0,'ignored NaN final ordinate')
    bad=w
    bad(n)=inf
    fit=splrep(x,y,w=bad,per=.true.,status=code)
    call require(code==0,'ignored infinite final weight')
    bad(2)=nan
    fit=splrep(x,y,w=bad,per=.true.,status=code)
    call require(code==1 .and. .not.fit%valid(),'active NaN weight rejected')
    bad=y
    bad(2)=inf
    fit=splrep(x,bad,per=.true.,status=code)
    call require(code==1,'active infinite ordinate rejected')
    fit=splrep(x,y,s=nan,per=.true.,status=code)
    call require(code==1,'NaN target rejected')
    fit=splrep(x,y,s=-1.0_dp,per=.true.,status=code)
    call require(code==1,'negative target rejected')
    fit=splrep(x,y,k=6,per=.true.,status=code)
    call require(code==1,'invalid degree')
    fit=splrep(x,y,bbox=[-1.0_dp,7.0_dp],per=.true.,status=code)
    call require(code==1,'period fixes the bounds')
    fit=splrep(x,y,task=-1,t=[1.0_dp],per=.true.,status=code)
    call require(code==2,'periodic fixed knots explicitly unsupported')
    fit=splrep(x,y,task=1,per=.true.,status=code)
    call require(code==1,'stateless continuation rejected')
    bad=x
    bad(3)=bad(2)
    fit=splrep(bad,y,per=.true.,status=code)
    call require(code==1,'repeated abscissa rejected')
    fit=splrep(x,y*1e100_dp,w=w*1e-100_dp,s=0.05_dp,per=.true.,status=code)
    call require(code==0,'balanced extreme scales')
    vals=splev(x,fit)/1e100_dp
    call require(abs(sum((w(:n-1)*(vals(:n-1)-y(:n-1)))**2)-0.05_dp)<0.000055_dp, &
        'extreme scale target')
    fit=splrep(x,y*1e300_dp,w=w*1e-310_dp,s=0.05e-20_dp,per=.true.,status=code)
    call require(code==0,'subnormal weight scale')
    fit=splrep(x*1e200_dp,y,s=0.05_dp,per=.true.,status=code)
    call require(code==0,'large coordinate scale')
    fit=splrep(x*1e-200_dp,y,s=0.05_dp,per=.true.,status=code)
    call require(code==0,'small coordinate scale')
    fit=splrep(1e16_dp+2*real([(i-1,i=1,n)],dp),y,s=0.05_dp,per=.true.,status=code)
    call require(code==0,'large common coordinate offset')
    fit=splrep(x,y,s=1e-310_dp,per=.true.,status=code)
    call require(code==3 .and. .not.fit%valid(),'unrepresentable target fails clearly')
    fit=splrep(x,0*y,per=.true.,status=code)
    call require(code==0 .and. maxval(abs(splev(x,fit)))==0.0_dp,'zero interpolation')
    fit=splrep(x,y,s=0.05_dp,status=code)
    other=splrep(x,y,s=0.05_dp,status=code,per=.false.)
    call require(maxval(abs(fit%c-other%c))==0.0_dp,'false leaves nonperiodic unchanged')
    print '(a,i0,a)','periodic splrep: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Condition required for success.
        character(*), intent(in) :: label !! Failure diagnostic.
        if (.not.ok) then
            print *, 'FAIL: ',label
            stop 1
        end if
        checks=checks+1
    end subroutine require
    subroutine seam(tck,degree)
        implicit none
        type(fitpack_tck_type), intent(in) :: tck !! Periodic representation to check.
        integer, intent(in) :: degree !! Degree whose lower derivatives must match.
        do d=0,degree-1
            a=splev(x(1),tck,der=d)
            b=splev(x(n),tck,der=d)
            call require(abs(a-b)<1e-9_dp*(1+abs(a)+abs(b)),'endpoint continuity')
        end do
    end subroutine seam
end program test_interpolate_fitpack_periodic
