! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_modern_selection
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep,make_lsq_spline,bspline, &
        bspline_interpolator_type,bspline_array_interpolator_type,splrep,splprep,splev, &
        fitpack_tck_type,fitpack_parametric_tck_type
    implicit none
    integer, parameter :: m=41,large=301
    type(bspline_interpolator_type) :: scalar,polynomial,source
    type(bspline_array_interpolator_type) :: curve
    type(fitpack_tck_type) :: legacy
    type(fitpack_parametric_tck_type) :: legacy_curve
    real(dp) :: u(m),y(m),w(m),x(m,3),target,fp,mean
    real(dp) :: ul(large),yl(large),xl(large,3)
    real(dp), allocatable :: t(:),c(:),p(:),v(:,:)
    integer :: i,k,j,code,checks,nest
    checks=0
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
    w=0.5_dp+u
    do k=1,5
        nest=2*(k+1)
        y=100.0_dp*max(u-0.5_dp,0.0_dp)**k
        t=[(0.0_dp,i=1,k+1),(1.0_dp,i=1,k+1)]
        polynomial=make_lsq_spline(u,y,t,k=k,w=w,status=code)
        call require(code==0,'independent polynomial reference')
        target=0.5_dp*sum((w*(polynomial%evaluate(u)-y))**2)
        scalar=make_splrep(u,y,w=w,k=k,s=target,nest=nest,status=code)
        call require(code==0 .and. scalar%valid(),'feasible soft scalar count')
        t=scalar%knots()
        call require(size(t)==nest+1,'scalar count exceeds soft target by one')
        call require(t(k+2)==u((m+1)/2),'middle data site selected')
        fp=sum((w*(scalar%evaluate(u)-y))**2)
        call require(abs(fp-target)<0.0011_dp*target,'soft scalar residual satisfied')
        x(:,1)=y
        x(:,2)=2*y
        x(:,3)=-0.5_dp*y
        curve=make_splprep(x,u=u,w=w,k=k,s=5.25_dp*target,nest=nest,status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'feasible soft shared count')
        call require(all(p==u),'soft fit returns supplied parameters')
        call require(size(curve%knots())==nest+1,'shared count exceeds soft target')
        v=curve%evaluate(u)
        call require(abs(sum((spread(w,2,3)*(v-x))**2)-5.25_dp*target) < &
            0.0011_dp*5.25_dp*target,'one shared residual at soft stop')
        call require(maxval(abs(v(:,1)-scalar%evaluate(u)))<1e-9_dp,'scalar/shared agreement')

        ! A known two-coefficient periodic curve is feasible after the first insertion.
        t=[(real(i-k-1,dp)/2.0_dp,i=1,2*k+3)]
        c=[(real(modulo(i-1,2),dp),i=1,k+2)]
        source=bspline(t,c,k,status=code)
        y=source%evaluate(u)
        y(m)=y(1)
        mean=sum(w(:m-1)**2*y(:m-1))/sum(w(:m-1)**2)
        target=0.5_dp*sum((w(:m-1)*(y(:m-1)-mean))**2)
        scalar=make_splrep(u,y,w=w,k=k,s=target,nest=nest,bc_type='periodic',status=code)
        call require(code==0 .and. scalar%valid(),'periodic soft scalar count')
        t=scalar%knots()
        call require(size(t)==nest+1,'periodic first vector exceeds soft count')
        call require(maxval(abs(t(3:)-t(:size(t)-2)-1.0_dp))<2e-15_dp,'cyclic soft-stop knots')
        fp=sum((w(:m-1)*(scalar%evaluate(u(:m-1))-y(:m-1)))**2)
        call require(abs(fp-target)<0.0011_dp*target,'periodic scalar residual')
        do j=0,k-1
            call require(abs(scalar%derivative(u(1),j)-scalar%derivative(u(m),j)) < &
                1e-9_dp*max(1.0_dp,maxval(abs(scalar%derivative(u,j)))),'periodic seam')
        end do
        x(:,1)=y
        x(:,2)=2*y
        x(:,3)=-0.5_dp*y
        curve=make_splprep(x,u=u,w=w,k=k,s=5.25_dp*target,nest=nest, &
            bc_type='periodic',status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'periodic soft shared count')
        call require(all(p==u),'periodic parameters returned')
        call require(size(curve%knots())==nest+1,'periodic shared soft count')
        v=curve%evaluate(u)
        call require(abs(sum((spread(w(:m-1),2,3)*(v(:m-1,:)-x(:m-1,:)))**2) - &
            5.25_dp*target)<0.0011_dp*5.25_dp*target,'periodic shared residual')
        call require(maxval(abs(v(:,1)-scalar%evaluate(u)))<1e-10_dp,'periodic scalar/shared agreement')
        scalar=make_splrep(u,y,w=w,k=k,s=1e-15_dp,nest=nest,bc_type='periodic',status=code)
        ! This curve is exactly representable on the selected knots, even for a tiny target.
        call require(code==0 .and. scalar%valid(),'small feasible soft target')
    end do
    ul=[(real(i-1,dp)/real(large-1,dp),i=1,large)]
    yl=cos(2*acos(-1.0_dp)*ul)
    yl(large)=yl(1)
    xl(:,1)=yl
    xl(:,2)=2*yl
    xl(:,3)=-0.5_dp*yl
    do j=0,1
        if (j==0) then
            scalar=make_splrep(ul,yl,s=0.2_dp,status=code)
            legacy=splrep(ul,yl,s=0.2_dp)
            curve=make_splprep(xl,u=ul,s=0.2_dp,status=code,parameters=p)
            legacy_curve=splprep(xl,u=ul,s=0.2_dp)
        else
            scalar=make_splrep(ul,yl,s=1000.0_dp,bc_type='periodic',status=code)
            legacy=splrep(ul,yl,s=1000.0_dp,per=.true.)
            curve=make_splprep(xl,u=ul,s=1000.0_dp,bc_type='periodic',status=code,parameters=p)
            legacy_curve=splprep(xl,u=ul,s=1000.0_dp,per=.true.)
        end if
        call require(code==0 .and. curve%valid() .and. scalar%valid(),'large-input compatibility')
        call require(maxval(abs(scalar%evaluate(ul)-splev(ul,legacy)))<2e-12_dp,'large scalar legacy path')
        call require(maxval(abs(curve%evaluate(ul)-splev(ul,legacy_curve)))<2e-12_dp, &
            'large shared legacy path')
        call require(all(p==ul),'large supplied parameters')
    end do
    write(*,'(a,i0,a)') 'Modern constructor selection: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Whether the numerical contract holds.
        character(*), intent(in) :: label !! Diagnostic printed on failure.
        checks=checks+1
        if (.not.ok) then
            write(*,'(a)') 'FAIL: '//label
            stop 1
        end if
    end subroutine require
end program test_interpolate_modern_selection
