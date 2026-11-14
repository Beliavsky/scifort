! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_fitpack_smoothing
    use scifort_kinds, only : dp
    use scifort_interpolate, only : splrep,splev,splder,splantider,splint, &
        fitpack_tck_type,fitpack_tck_invalid_input,fitpack_tck_singular
    use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: n=61
    real(dp) :: x(n),y(n),w(n),bad(n),vals(n),single(1),target,fp,q,nan,inf
    type(fitpack_tck_type) :: fit,deriv,integ,poly
    integer :: i,k,status,checks
    checks=0
    do i=1,n
        x(i)=real(i-1,dp)/real(n-1,dp)
        y(i)=sin(8*x(i))+0.15_dp*cos(79*x(i))
        w(i)=0.5_dp+x(i)
    end do
    do k=1,5
        target=0.15_dp
        fit=splrep(x,y,w=w,k=k,s=target,status=status)
        call require(status==0 .and. fit%valid(),'adaptive fit')
        call require(size(fit%t)>2*(k+1),'interior knots selected')
        call require(size(fit%c)==size(fit%t),'padded coefficient tuple')
        vals=splev(x,fit,status=status)
        fp=sum((w*(vals-y))**2)
        call require(status==0 .and. abs(fp-target)<0.0011_dp*target,'residual target')
        q=splev(0.317_dp,fit)
        single=splev([0.317_dp],fit)
        call require(abs(q-single(1))<1e-14_dp,'scalar array evaluation')
        deriv=splder(fit,status=status)
        call require(status==0,'derivative construction')
        call require(abs(splev(0.317_dp,deriv)-splev(0.317_dp,fit,der=1))<1e-10_dp, &
            'derivative identity')
        integ=splantider(fit,status=status)
        call require(status==0,'antiderivative construction')
        call require(abs(splev(0.8_dp,integ)-splev(0.1_dp,integ)- &
            splint(0.1_dp,0.8_dp,fit))<1e-11_dp,'integral identity')
        poly=splrep(x,y,w=w,k=k,s=100.0_dp,status=status)
        call require(status==0 .and. poly%valid(),'polynomial limit')
        call require(size(poly%t)==2*(k+1),'minimal polynomial knots')
        fit=splrep(x,y,k=k,s=1e-6_dp,status=status)
        call require(status==0 .and. fit%valid(),'near interpolation')
        call require(sum((splev(x,fit)-y)**2)<=1.0011e-6_dp,'small residual')
        fit=splrep(x,y,w=w,k=k,s=target,bbox=[-0.3_dp,1.4_dp],status=status)
        call require(status==0 .and. fit%valid(),'enlarged support')
        call require(fit%t(1)==-0.3_dp .and. fit%t(size(fit%t))==1.4_dp,'support bounds')
        call require(abs(sum((w*(splev(x,fit)-y))**2)-target)<0.0011_dp*target, &
            'enlarged support residual')
        fit=splrep(x(:k+1),y(:k+1),k=k,s=1.0_dp,status=status)
        call require(status==0 .and. fit%valid(),'minimal sample count')
    end do
    fit=splrep(x,0*y,s=1.0_dp,status=status)
    call require(status==0 .and. maxval(abs(splev(x,fit)))==0.0_dp,'zero observations')
    fit=splrep(x,y*1e100_dp,w=w*1e-100_dp,s=0.15_dp,status=status)
    call require(status==0 .and. fit%valid(),'balanced large scale')
    call require(abs(sum((w*(splev(x,fit)/1e100_dp-y))**2)-0.15_dp)<0.000165_dp, &
        'large scale residual')
    fit=splrep(x,y*1e-100_dp,s=0.15e-200_dp,status=status)
    call require(status==0 .and. fit%valid(),'small data scale')
    call require(abs(sum((splev(x,fit)/1e-100_dp-y)**2)-0.15_dp)<0.000165_dp, &
        'small scale residual')
    fit=splrep(x*1e200_dp,y,s=0.15_dp,status=status)
    call require(status==0 .and. fit%valid(),'large coordinate scale')
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    fit=splrep(x,y,s=-1.0_dp,status=status)
    call invalid('negative smoothing')
    fit=splrep(x,y,s=nan,status=status)
    call invalid('NaN smoothing')
    fit=splrep(x,y,s=inf,status=status)
    call invalid('infinite smoothing')
    bad=w
    bad(2)=0
    fit=splrep(x,y,w=bad,s=0.1_dp,status=status)
    call invalid('zero weight')
    bad=x
    bad(2)=bad(1)
    fit=splrep(bad,y,s=0.1_dp,status=status)
    call invalid('duplicate abscissa')
    bad=y
    bad(2)=nan
    fit=splrep(x,bad,s=0.1_dp,status=status)
    call invalid('NaN observation')
    fit=splrep(x,y,s=0.1_dp,bbox=[nan,1.0_dp],status=status)
    call invalid('NaN support')
    fit=splrep(x,y,s=0.1_dp,bbox=[0.1_dp,1.0_dp],status=status)
    call invalid('support excludes samples')
    fit=splrep(x,y,s=0.1_dp,bbox=[0.0_dp],status=status)
    call invalid('support shape')
    fit=splrep(x,y,s=0.1_dp,t=[0.5_dp],status=status)
    call invalid('knots with automatic task')
    fit=splrep(x,y,s=0.1_dp,task=1,status=status)
    call invalid('unsupported continuation retains invalid status')
    fit=splrep(x,y*1e100_dp,s=tiny(1.0_dp),status=status)
    call require(status==fitpack_tck_singular .and. .not.fit%valid(), &
        'unrepresentable scaled target fails explicitly')
    fit=splrep(x,y,s=1e-100_dp,status=status)
    call require(status==fitpack_tck_singular .and. .not.fit%valid(), &
        'target below achieved rounding error fails explicitly')
    bad=w
    bad(2)=1e-100_dp
    fit=splrep(x,y,w=bad,s=0.1_dp,status=status)
    call require(status==fitpack_tck_singular .and. .not.fit%valid(), &
        'extreme relative weights fail explicitly')
    bad=x
    bad(2)=1e-100_dp
    fit=splrep(bad,y,k=5,s=0.1_dp,status=status)
    call require(status==fitpack_tck_singular .and. .not.fit%valid(), &
        'unresolved derivative-jump scale fails explicitly')
    bad=w
    bad(2)=nan
    fit=splrep(x,y,w=bad,s=0.1_dp,status=status)
    call invalid('NaN weight')
    bad=y
    bad(2)=inf
    fit=splrep(x,bad,s=0.1_dp,status=status)
    call invalid('infinite observation')
    print '(a,i0,a)','Adaptive scalar FITPACK smoothing: ',checks,' checks passed.'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Condition required by this regression.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not.ok) then
            print '(a,a)','FAIL: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine require
    subroutine invalid(label)
        implicit none
        character(len=*), intent(in) :: label !! Invalid-input regression description.
        call require(status==fitpack_tck_invalid_input .and. .not.fit%valid(),label)
    end subroutine invalid
end program test_interpolate_fitpack_smoothing
