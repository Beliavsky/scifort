! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_large_selection
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep,make_lsq_spline,bspline, &
        bspline_interpolator_type,bspline_array_interpolator_type
    use scifort_interpolate_curfit_core, only: fitpack_periodic_lsq
    implicit none
    type(bspline_interpolator_type) :: fit,poly,source,scaled
    type(bspline_array_interpolator_type) :: curve
    real(dp), allocatable :: u(:),y(:),w(:),x(:,:),t(:),p(:),v(:,:),c(:,:),known(:)
    real(dp) :: target,fp,pi,split
    integer, parameter :: sizes(4)=[301,1001,3001,10000]
    integer :: m,choice,k,i,code,checks,j,ni
    checks=0
    pi=acos(-1.0_dp)
    do choice=1,size(sizes)
        m=sizes(choice)
        u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
        y=u
        w=0.5_dp+u
        if (allocated(x)) deallocate(x)
        allocate(x(m,3))
        do k=1,5
            if (m>1001 .and. k/=3) cycle
            split=u((m+2)/2)
            y=100*max(u-split,0.0_dp)**k
            t=[(0.0_dp,i=1,k+1),(1.0_dp,i=1,k+1)]
            poly=make_lsq_spline(u,y,t,k=k,w=w,status=code)
            call require(code==0,'large polynomial reference')
            target=0.5_dp*sum((w*(poly%evaluate(u)-y))**2)
            fit=make_splrep(u,y,k=k,w=w,s=target,nest=2*(k+1),status=code)
            call require(code==0 .and. fit%valid(),'large scalar soft count')
            t=fit%knots()
            call require(size(t)==2*(k+1)+1,'large soft count exceeds target')
            call require(t(k+2)==split,'large middle-site selection')
            fp=sum((w*(fit%evaluate(u)-y))**2)
            call require(abs(fp-target)<0.0011_dp*target,'large soft residual')
            x(:,1)=y
            x(:,2)=2*y
            x(:,3)=-0.5_dp*y
            curve=make_splprep(x,u=u,k=k,w=w,s=5.25_dp*target,nest=2*(k+1), &
                status=code,parameters=p)
            call require(code==0 .and. curve%valid(),'large shared soft count')
            call require(all(p==u),'large parameter output')
            call require(size(curve%knots())==size(t),'large scalar/shared knot count')
            v=curve%evaluate(u)
            call require(abs(sum((spread(w,2,3)*(v-x))**2)-5.25_dp*target) < &
                0.0011_dp*5.25_dp*target,'large shared residual')

            y=sin(14*pi*u)+0.07_dp*cos(58*pi*u)
            target=0.001_dp*sum((w*y)**2)
            fit=make_splrep(u,y,k=k,w=w,s=target,status=code)
            call require(code==0 .and. fit%valid(),'large ordinary automatic refinement')
            call require(size(fit%knots())>3*k+2,'many ordinary interior knots')
            call require(abs(sum((w*(fit%evaluate(u)-y))**2)-target) < &
                0.0011_dp*target,'large ordinary refined residual')
            y(m)=y(1)
            target=0.001_dp*sum((w(:m-1)*y(:m-1))**2)
            fit=make_splrep(u,y,k=k,w=w,s=target,bc_type='periodic',status=code)
            call require(code==0 .and. fit%valid(),'large cyclic automatic refinement')
            t=fit%knots()
            ni=size(t)-2*k-1
            call require(ni>k,'large cyclic banded solve required')
            call require(maxval(abs(t(ni+1:)-t(:size(t)-ni)-1.0_dp))<2e-15_dp, &
                'large canonical periodic knots')
            fp=sum((w(:m-1)*(fit%evaluate(u(:m-1))-y(:m-1)))**2)
            call require(abs(fp-target)<0.0011_dp*target,'large cyclic residual')
            do j=0,k-1
                call require(abs(fit%derivative(u(1),j)-fit%derivative(u(m),j)) < &
                    1e-9_dp*max(1.0_dp,maxval(abs(fit%derivative(u,j)))),'large cyclic seam')
            end do
            x(:,1)=y
            x(:,2)=2*y
            x(:,3)=-0.5_dp*y
            curve=make_splprep(x,u=u,k=k,w=w,s=5.25_dp*target,bc_type='periodic', &
                status=code,parameters=p)
            call require(code==0 .and. curve%valid(),'large cyclic shared refinement')
            v=curve%evaluate(u)
            call require(abs(sum((spread(w(:m-1),2,3)*(v(:m-1,:)-x(:m-1,:)))**2) - &
                5.25_dp*target)<0.0011_dp*5.25_dp*target,'large cyclic shared residual')
            if (m==1001 .and. k==3) then
                scaled=make_splrep(u,y*1e100_dp,k=k,w=w*1e-100_dp,s=target, &
                    bc_type='periodic',status=code)
                call require(code==0 .and. scaled%valid(),'large balanced scaling')
                call require(abs(sum((w(:m-1)*(scaled%evaluate(u(:m-1))/1e100_dp - &
                    y(:m-1)))**2)-target)<0.0011_dp*target,'balanced large cyclic residual')
            end if
        end do
    end do

    ! Exercise the reused cyclic QR directly with known coefficients and missing columns.
    m=1001
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
    w=0.5_dp+u
    do k=1,5
        ni=32
        t=[(real(i-k-1,dp)/real(ni,dp),i=1,ni+2*k+1)]
        known=[(sin(2*pi*real(modulo(i-1,ni),dp)/real(ni,dp)),i=1,ni+k)]
        source=bspline(t,known,k,status=code)
        y=source%evaluate(u)
        call fitpack_periodic_lsq(u,reshape(y,[m,1]),w,k,t,c,code)
        call require(code==0,'banded cyclic coefficient solve')
        call require(maxval(abs(c(:,1)-known))<2e-12_dp,'known cyclic coefficient recovery')
        p=u
        p(:m-1)=p(:m-1)*1e-4_dp
        call fitpack_periodic_lsq(p,reshape(y,[m,1]),w,k,t,c,code)
        call require(code==3 .and. .not.allocated(c),'cyclic rank failure has no coefficients')
    end do
    m=10000
    k=3
    ni=m-1
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
    w=0.5_dp+u
    t=[(real(i-k-1,dp)/real(ni,dp),i=1,ni+2*k+1)]
    known=[(sin(2*pi*real(modulo(i-1,ni),dp)/real(ni,dp)),i=1,ni+k)]
    source=bspline(t,known,k,status=code)
    y=source%evaluate(u)
    call fitpack_periodic_lsq(u,reshape(y,[m,1]),w,k,t,c,code)
    call require(code==0,'near-ten-thousand coefficient cyclic solve')
    call require(maxval(abs(c(:,1)-known))<2e-11_dp,'large coefficient recovery without dense design')
    write(*,'(a,i0,a)') 'Large modern selection: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Whether the numerical contract holds.
        character(*), intent(in) :: label !! Failure description.
        checks=checks+1
        if (.not.ok) then
            write(*,'(a)') 'FAIL: '//label
            stop 1
        end if
    end subroutine require
end program test_interpolate_large_selection
