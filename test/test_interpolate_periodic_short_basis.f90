! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_periodic_short_basis
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splrep,splprep,splev,fitpack_tck_type, &
        fitpack_parametric_tck_type,fitpack_continuation_type
    implicit none
    integer, parameter :: m=61
    real(dp), parameter :: factors(6)=[0.03_dp,0.2_dp,0.001_dp,10.0_dp,0.0_dp,0.03_dp]
    real(dp) :: u(m),y(m),w(m),data(m,1),energy,target,fp,a,b,scale
    type(fitpack_tck_type) :: fit,fresh
    type(fitpack_parametric_tck_type) :: curve
    type(fitpack_continuation_type) :: state
    real(dp), allocatable :: values(:,:),v(:)
    integer :: i,k,j,d,code,n7,checks,mode
    checks=0
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]**1.1_dp
    y=cos(2*acos(-1.0_dp)*u)+0.08_dp*sin(14*acos(-1.0_dp)*u)
    y(m)=y(1)
    w=0.8_dp+u
    energy=sum((w(:m-1)*y(:m-1))**2)
    data(:,1)=y
    do k=1,5
        call state%reset()
        do j=1,size(factors)
            target=factors(j)*energy
            mode=1
            if (j==1) mode=0
            fit=splrep(u,y,w=w,k=k,s=target,per=.true.,task=mode,state=state,status=code)
            call require(code==0 .and. fit%valid(),'continued scalar fit')
            fresh=splrep(u,y,w=w,k=k,s=target,per=.true.,status=code)
            call require(code==0 .and. fresh%valid(),'fresh scalar fit')
            call check_fit(fit,target,k)
            call check_fit(fresh,target,k)
            curve=splprep(data,u=u,w=w,k=k,s=target,per=.true.,status=code)
            call require(code==0 .and. curve%valid(),'independent one-channel fit')
            values=splev(u,curve)
            v=splev(u,fresh)
            call require(maxval(abs(values(:,1)-v))<2e-10_dp,'scalar/parametric agreement')
            if (j==1 .and. (k==3 .or. k==5)) then
                n7=size(fresh%t)-2*k-1
                call require(n7==2 .and. n7<k,'regression uses multiple cyclic wraps')
            end if
        end do
    end do
    print '(a,i0,a)', 'test_interpolate_periodic_short_basis: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Condition required for the regression to pass.
        character(*), intent(in) :: label !! Description printed on failure.
        checks=checks+1
        if (.not.ok) then
            print '(a)',label
            error stop 1
        end if
    end subroutine require
    subroutine check_fit(tck,s,degree)
        implicit none
        type(fitpack_tck_type), intent(in) :: tck !! Successful scalar periodic fit.
        real(dp), intent(in) :: s !! Nonnegative weighted residual target.
        integer, intent(in) :: degree !! Degree whose derivatives must meet at the seam.
        v=splev(u,tck)
        fp=sum((w(:m-1)*(v(:m-1)-y(:m-1)))**2)
        if (s==0.0_dp) then
            call require(fp<1e-22_dp,'interpolation residual')
        else if (s<energy) then
            call require(abs(fp-s)<0.0011_dp*s,'active weighted residual')
        else
            call require(fp<=1.0011_dp*s,'constant residual bound')
        end if
        n7=size(tck%t)-2*degree-1
        do i=1,degree
            call require(tck%c(n7+i)==tck%c(i),'cyclic coefficient identity')
        end do
        do d=0,degree-1
            a=splev(u(1),tck,der=d)
            b=splev(u(m),tck,der=d)
            scale=max(1.0_dp,maxval(abs(splev(u,tck,der=d))))
            call require(abs(a-b)<2e-10_dp*scale,'periodic derivative seam')
        end do
    end subroutine check_fit
end program test_interpolate_periodic_short_basis
