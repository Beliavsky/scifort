! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_make_splprep
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splprep,make_interp_spline, &
        bspline_array_interpolator_type,splprep,splev,fitpack_parametric_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: m=61
    real(dp) :: u(m),x(m,3),w(m),bad(m),badx(m,3),chords(m),fp,nan,inf
    real(dp), allocatable :: t(:),p(:),v(:,:),c(:,:)
    type(bspline_array_interpolator_type) :: fit,reference,poly
    type(fitpack_parametric_tck_type) :: legacy
    integer :: i,k,j,code,checks
    checks=0
    do i=1,m
        u(i)=real(i-1,dp)/real(m-1,dp)
        x(i,:)=[sin(8*u(i)),cos(6*u(i)),u(i)**2]+0.08_dp*cos(79*u(i))
        w(i)=0.5_dp+u(i)
    end do
    chords(1)=0.0_dp
    do i=2,m
        chords(i)=chords(i-1)+sqrt(sum((x(i,:)-x(i-1,:))**2))
    end do
    chords=chords/chords(m)
    do k=1,5
        fit=make_splprep(x,u=u,k=k,status=code,parameters=p)
        call require(code==0 .and. fit%valid(),'supplied-site interpolation')
        call require(all(p==u),'returned supplied sites')
        call require(maxval(abs(fit%evaluate(u)-x))<3e-12_dp,'interpolates all channels')
        reference=make_interp_spline(u,x,k=k,status=code)
        call require(all(fit%knots()==reference%knots()),'modern interpolation knots')
        fit=make_splprep(x,k=k,status=code,parameters=p)
        call require(code==0 .and. fit%valid(),'chord-site interpolation')
        call require(p(1)==0.0_dp .and. p(m)==1.0_dp,'normalized chord endpoints')
        call require(all(p(2:)>p(:m-1)),'increasing chord parameters')
        call require(maxval(abs(p-chords))<5e-15_dp,'Euclidean chord-length rule')
        call require(maxval(abs(fit%evaluate(p)-x))<3e-12_dp,'chord interpolation values')
        fit=make_splprep(x,u=u,w=w,k=k,s=0.15_dp,status=code,parameters=p)
        call require(code==0 .and. fit%valid(),'automatic smoothing')
        v=fit%evaluate(u)
        fp=sum((spread(w,2,3)*(v-x))**2)
        call require(abs(fp-0.15_dp)<0.000165_dp,'one shared weighted residual')
        t=fit%knots()
        c=fit%coefficients()
        call require(size(c,1)==size(t)-k-1 .and. size(c,2)==3,'native coefficient shape')
        legacy=splprep(x,u=u,w=w,k=k,s=0.15_dp,status=code)
        call require(maxval(abs(v-splev(u,legacy)))<1e-11_dp,'legacy automatic agreement')
        reference=make_splprep(x,u=u,w=w,k=k,s=0.15_dp,t=t,status=code)
        call require(code==0 .and. reference%valid(),'fixed-knot smoothing')
        call require(all(reference%knots()==t),'full supplied knots preserved')
        v=reference%evaluate(u)
        call require(abs(sum((spread(w,2,3)*(v-x))**2)-0.15_dp)<0.000165_dp, &
            'fixed-knot residual target')
        call require(maxval(abs(reference%evaluate(u(7))-v(7,:)))<2e-12_dp, &
            'scalar and array evaluation')
        call require(maxval(abs(reference%integrate(0.0_dp,1.0_dp)+ &
            reference%integrate(1.0_dp,0.0_dp)))<2e-12_dp,'reversed integral')
        reference=make_splprep(x,u=u,w=w,k=k,s=100.0_dp,t=t,status=code)
        call require(code==0 .and. reference%valid(),'polynomial on supplied knots')
        poly=make_splprep(x,u=u,w=w,k=k,s=100.0_dp,status=code)
        call require(all(reference%knots()==t),'polynomial retains knots')
        call require(maxval(abs(reference%evaluate(u)-poly%evaluate(u)))<2e-11_dp, &
            'zero-jump polynomial limit')
        fit=make_splprep(x,u=u,w=w,k=k,s=0.15_dp,nest=m+k+1,status=code)
        call require(code==0 .and. fit%valid(),'explicit full capacity')
        fit=make_splprep(x,u=u,w=w,k=k,s=100.0_dp,nest=2*(k+1),status=code)
        call require(code==0 .and. fit%valid(),'minimal feasible capacity')
        fit=make_splprep(x,u=u,w=w,k=k,s=1e-7_dp,nest=2*(k+1),status=code,parameters=p)
        call require(code==3 .and. .not.fit%valid(),'insufficient capacity')
        call require(.not.allocated(p),'no sites returned on failure')
        fit=make_splprep(x,u=u,w=w,k=k,s=0.15_dp,ub=-0.3_dp,ue=1.4_dp,status=code)
        call require(code==0 .and. fit%valid(),'enlarged smoothing bounds')
        fit=make_splprep(x,k=k,s=0.15_dp,status=code,parameters=p)
        call require(code==0 .and. fit%valid(),'automatic chords and smoothing')
        call require(abs(sum((fit%evaluate(p)-x)**2)-0.15_dp)<0.000165_dp, &
            'chord smoothing target')
    end do
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    fit=make_splprep(x,w=w,status=code)
    call invalid('interpolation weights')
    fit=make_splprep(x,t=t,status=code)
    call invalid('interpolation knots')
    fit=make_splprep(x,nest=m+4,status=code)
    call invalid('interpolation capacity')
    do j=1,3
        badx=x
        badx(3,2)=nan
        if (j==2) badx(3,2)=inf
        if (j==3) badx(3,2)=-inf
        fit=make_splprep(badx,s=0.15_dp,status=code)
        call invalid('nonfinite coordinates')
    end do
    bad=u
    bad(3)=bad(2)
    fit=make_splprep(x,u=bad,status=code)
    call invalid('duplicate sites')
    bad(3)=nan
    fit=make_splprep(x,u=bad,status=code)
    call invalid('NaN site')
    fit=make_splprep(x,u=u(:m-1),status=code)
    call invalid('short sites')
    badx=x
    badx(3,:)=badx(2,:)
    fit=make_splprep(badx,status=code)
    call invalid('zero chord')
    fit=make_splprep(badx,u=u,status=code)
    call require(code==0 .and. fit%valid(),'supplied sites allow repeated coordinates')
    fit=make_splprep(x,s=-1.0_dp,status=code)
    call invalid('negative target')
    fit=make_splprep(x,s=nan,status=code)
    call invalid('NaN target')
    fit=make_splprep(x,s=inf,status=code)
    call invalid('infinite target')
    bad=w
    bad(3)=0.0_dp
    fit=make_splprep(x,s=0.15_dp,w=bad,status=code)
    call invalid('zero weight')
    bad(3)=inf
    fit=make_splprep(x,s=0.15_dp,w=bad,status=code)
    call invalid('infinite weight')
    fit=make_splprep(x,s=0.15_dp,w=w(:m-1),status=code)
    call invalid('short weights')
    fit=make_splprep(x,k=0,status=code)
    call invalid('invalid degree')
    fit=make_splprep(x(:3,:),status=code)
    call invalid('too few samples')
    fit=make_splprep(x(:,:0),status=code)
    call invalid('no channels')
    fit=make_splprep(x,ub=nan,status=code)
    call invalid('NaN bound')
    fit=make_splprep(x,u=u,ub=-1.0_dp,status=code)
    call invalid('interpolation bound differs')
    fit=make_splprep(x,s=0.15_dp,nest=7,status=code)
    call invalid('small capacity')
    fit=make_splprep(x,s=0.15_dp,nest=m+5,status=code)
    call invalid('excess capacity')
    fit=make_splprep(x,s=0.15_dp,t=t,nest=m+4,status=code)
    call invalid('knots and capacity together')
    fit=make_splprep(x,s=0.15_dp,t=[0.0_dp,1.0_dp],status=code)
    call invalid('short knots')
    fit=make_splprep(x,u=u,s=0.15_dp, &
        t=[(0.0_dp,i=1,4),0.4_dp,0.4_dp,(1.0_dp,i=1,4)],status=code)
    call invalid('repeated interior knot')
    fit=make_splprep(x,u=u,s=0.15_dp, &
        t=[(0.0_dp,i=1,4),nan,(1.0_dp,i=1,4)],status=code)
    call invalid('nonfinite knot')
    fit=make_splprep(x,u=u,s=0.15_dp, &
        t=[-0.1_dp,(0.0_dp,i=1,3),(1.0_dp,i=1,4)],status=code)
    call invalid('nonclamped knot')
    fit=make_splprep(x,s=0.15_dp,bc_type='natural',status=code)
    call invalid('invalid boundary name')
    fit=make_splprep(x,bc_type='periodic',status=code)
    call require(code==2 .and. .not.fit%valid(),'periodic unsupported')
    fit=make_splprep(x,u=u,w=w,s=1e-12_dp,t=[(0.0_dp,i=1,4),(1.0_dp,i=1,4)],status=code)
    call require(code==3 .and. .not.fit%valid(),'fixed knots cannot satisfy target')
    fit=make_splprep(x,u=u,s=0.15_dp, &
        t=[(0.0_dp,i=1,4),0.00001_dp,0.00002_dp,(1.0_dp,i=1,4)],status=code)
    call require(code==3 .and. .not.fit%valid(),'rank-deficient fixed collocation')
    block
        real(dp) :: high(m,12)
        do j=1,12
            high(:,j)=real(j,dp)*x(:,mod(j,3)+1)
        end do
        fit=make_splprep(high,u=u,s=3.0_dp,status=code)
        call require(code==0 .and. fit%channels()==12,'more than ten channels')
    end block
    bad=u+1e12_dp
    fit=make_splprep(x,u=bad,s=0.15_dp,status=code)
    call require(code==0 .and. fit%valid(),'large parameter offset')
    fit=make_splprep(x*1e150_dp,u=u,w=w*1e-150_dp,s=0.15_dp,status=code)
    call require(code==0 .and. fit%valid(),'common coordinate and weight scaling')
    fit=make_splprep(x*1e150_dp,w=w*1e-150_dp,s=0.15_dp,status=code,parameters=p)
    call require(code==0 .and. fit%valid(),'large coordinate chord smoothing')
    call require(maxval(abs(p-chords))<5e-15_dp,'scaled chord rule')
    fit=make_splprep(x,u=u*1e-310_dp,status=code)
    call require(code==3 .and. .not.fit%valid(),'subnormal interpolation span')
    fit=make_splprep(x,u=u*3e-308_dp,status=code)
    call require(code==3 .and. .not.fit%valid(),'subnormal interpolation gaps')
    fit=make_splprep(x,u=u*huge(1.0_dp),status=code)
    call require(code==3 .and. .not.fit%valid(),'extreme interpolation parameters')
    bad=u
    bad(2)=1e-100_dp
    fit=make_splprep(x,u=bad,s=0.15_dp,status=code)
    call require(code==3 .and. .not.fit%valid(),'unresolved parameter gap')
    fit=make_splprep(x,u=u,s=0.15_dp, &
        t=[(0.0_dp,i=1,4),1e-100_dp,(1.0_dp,i=1,4)],status=code)
    call require(code==3 .and. .not.fit%valid(),'unresolved supplied knot gap')
    bad=w
    bad(3)=tiny(1.0_dp)
    fit=make_splprep(x,u=u,w=bad,s=0.15_dp,status=code)
    call require(code==3 .and. .not.fit%valid(),'relative weight guard')
    print '(a,i0,a)','make_splprep: ',checks,' checks passed'
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
    subroutine invalid(label)
        implicit none
        character(*), intent(in) :: label !! Invalid-input description.
        call require(code==1 .and. .not.fit%valid(),label)
    end subroutine invalid
end program test_interpolate_make_splprep
