! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_make_splrep
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_interp_spline,make_lsq_spline, &
        bspline_interpolator_type,splrep,splev,fitpack_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: m=61
    real(dp) :: x(m),y(m),w(m),bad(m),values(m),fp,nan,inf
    real(dp), allocatable :: t(:),c(:)
    type(bspline_interpolator_type) :: fit,reference,lsq
    type(fitpack_tck_type) :: legacy
    integer :: i,k,code,checks
    checks=0
    do i=1,m
        x(i)=real(i-1,dp)/real(m-1,dp)
        y(i)=sin(8*x(i))+0.15_dp*cos(79*x(i))
        w(i)=0.5_dp+x(i)
    end do
    do k=1,5
        fit=make_splrep(x,y,k=k,status=code)
        call require(code==0 .and. fit%valid(),'default interpolation')
        call require(maxval(abs(fit%evaluate(x)-y))<2e-12_dp,'interpolates data')
        reference=make_interp_spline(x,y,k=k,status=code)
        call require(maxval(abs(fit%knots()-reference%knots()))==0.0_dp,'modern default knots')
        fit=make_splrep(x,y,w=w,k=k,s=0.15_dp,status=code)
        call require(code==0 .and. fit%valid(),'automatic smoothing')
        fp=sum((w*(fit%evaluate(x)-y))**2)
        call require(abs(fp-0.15_dp)<0.000165_dp,'weighted residual target')
        t=fit%knots()
        c=fit%coefficients()
        call require(size(c)==size(t)-k-1,'unpadded native coefficients')
        legacy=splrep(x,y,w=w,k=k,s=0.15_dp,status=code)
        call require(maxval(abs(fit%evaluate(x)-splev(x,legacy)))<1e-10_dp,'legacy agreement')
        reference=make_splrep(x,y,w=w,k=k,s=0.15_dp,t=t,status=code)
        call require(code==0 .and. reference%valid(),'supplied-knot smoothing')
        call require(all(reference%knots()==t),'supplied knots preserved')
        call require(abs(sum((w*(reference%evaluate(x)-y))**2)-0.15_dp)<0.000165_dp, &
            'fixed smoothing target')
        lsq=make_lsq_spline(x,y,t,k=k,w=w,status=code)
        call require(sum((w*(lsq%evaluate(x)-y))**2)<= &
            sum((w*(reference%evaluate(x)-y))**2)+1e-13_dp,'smoothing differs from LSQ')
        reference=make_splrep(x,y,w=w,k=k,s=100.0_dp,t=t,status=code)
        call require(code==0 .and. reference%valid(),'polynomial on supplied knots')
        call require(all(reference%knots()==t),'polynomial retains supplied knots')
        fit=make_splrep(x,y,w=w,k=k,s=100.0_dp,status=code)
        call require(maxval(abs(reference%evaluate(x)-fit%evaluate(x)))<1e-11_dp, &
            'zero-jump polynomial limit')
        fit=make_splrep(x,y,w=w,k=k,s=0.15_dp,nest=m+k+1,status=code)
        call require(code==0 .and. fit%valid(),'explicit full capacity')
        fit=make_splrep(x,y,w=w,k=k,s=100.0_dp,nest=2*(k+1),status=code)
        call require(code==0 .and. fit%valid(),'minimum capacity when feasible')
        fit=make_splrep(x,y,w=w,k=k,s=0.00001_dp,nest=2*(k+1),status=code)
        call require(code==3 .and. .not.fit%valid(),'capacity cannot satisfy target')
        fit=make_splrep(x,y,w=w,k=k,s=0.15_dp,xb=-0.3_dp,xe=1.4_dp,status=code)
        call require(code==0 .and. fit%valid(),'enlarged smoothing bounds')
        call require(abs(sum((w*(fit%evaluate(x)-y))**2)-0.15_dp)<0.000165_dp, &
            'enlarged bounds residual')
    end do
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    fit=make_splrep(x,y,w=w,status=code)
    call invalid('weights rejected for interpolation')
    fit=make_splrep(x,y,nest=m+4,status=code)
    call invalid('capacity rejected for interpolation')
    fit=make_splrep(x,y,t=t,status=code)
    call invalid('knots rejected for interpolation')
    fit=make_splrep(x,y,s=-1.0_dp,status=code)
    call invalid('negative target')
    fit=make_splrep(x,y,s=nan,status=code)
    call invalid('NaN target')
    fit=make_splrep(x,y,s=inf,status=code)
    call invalid('infinite target')
    fit=make_splrep(x,y,s=0.15_dp,xb=nan,status=code)
    call invalid('NaN bound')
    bad=y
    bad(3)=nan
    fit=make_splrep(x,bad,s=0.15_dp,status=code)
    call invalid('NaN ordinate')
    bad=x
    bad(3)=bad(2)
    fit=make_splrep(bad,y,s=0.15_dp,status=code)
    call invalid('duplicate abscissa')
    bad=w
    bad(3)=0.0_dp
    fit=make_splrep(x,y,w=bad,s=0.15_dp,status=code)
    call invalid('zero weight')
    fit=make_splrep(x,y,k=6,status=code)
    call invalid('unsupported degree')
    fit=make_splrep(x,y,s=0.15_dp,nest=7,status=code)
    call invalid('invalid capacity')
    fit=make_splrep(x,y,s=0.15_dp,t=t,nest=m+4,status=code)
    call invalid('knots and capacity mutually exclusive')
    fit=make_splrep(x,y,s=0.15_dp,t=[0.0_dp,1.0_dp],status=code)
    call invalid('short knot vector')
    fit=make_splrep(x,y,s=0.15_dp,bc_type='natural',status=code)
    call invalid('unsupported boundary name')
    fit=make_splrep(x,y,s=0.15_dp,bc_type='periodic',status=code)
    call invalid('periodic endpoints must match')
    fit=make_splrep(x,y,s=0.15_dp,bc_type='not-a-knot',status=code)
    call require(code==0 .and. fit%valid(),'explicit default boundary')
    t=[spread(0.0_dp,1,4),0.5_dp,spread(1.0_dp,1,4)]
    fit=make_splrep(x,y,s=1e-10_dp,t=t,status=code)
    call require(code==3 .and. .not.fit%valid(),'infeasible fixed-knot target')
    fit=make_splrep(x,y*1e100_dp,w=w*1e-100_dp,s=0.15_dp,status=code)
    call require(code==0 .and. fit%valid(),'balanced scales')
    values=fit%evaluate(x)/1e100_dp
    call require(abs(sum((w*(values-y))**2)-0.15_dp)<0.000165_dp,'scaled target')
    fit=make_splrep(1e16_dp+2*real([(i-1,i=1,m)],dp),y,s=0.15_dp,status=code)
    call require(code==0 .and. fit%valid(),'large offset coordinates')
    fit=make_splrep(x,y,s=1e-310_dp,status=code)
    call require(code==3 .and. .not.fit%valid(),'unrepresentable target')
    print '(a,i0,a)','make_splrep: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required numerical condition.
        character(*), intent(in) :: label !! Failure diagnostic.
        if (.not.ok) then
            print *, 'FAIL: ',label
            stop 1
        end if
        checks=checks+1
    end subroutine require
    subroutine invalid(label)
        implicit none
        character(*), intent(in) :: label !! Invalid-input case under test.
        call require(code==1 .and. .not.fit%valid(),label)
    end subroutine invalid
end program test_interpolate_make_splrep
