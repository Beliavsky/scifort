! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_fixed_general
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep, &
        bspline_interpolator_type,bspline_array_interpolator_type
    implicit none
    integer, parameter :: m=81
    real(dp) :: u(m),y(m),x(m,3),w(m),target,fp,scale
    real(dp), allocatable :: t(:),v(:,:),p(:)
    type(bspline_interpolator_type) :: fit
    type(bspline_array_interpolator_type) :: curve
    integer :: i,k,d,form,code,checks
    checks=0
    do i=1,m
        u(i)=real(i-1,dp)/real(m-1,dp)
        y(i)=cos(2*acos(-1.0_dp)*u(i))
        w(i)=0.8_dp+u(i)
    end do
    y(m)=y(1)
    x(:,1)=y
    x(:,2)=2*y
    x(:,3)=-0.5_dp*y
    do k=1,5
        do form=1,3
            if (allocated(t)) deallocate(t)
            allocate(t(2*k+9))
            t(k+1:k+9)=[(real(i,dp)/8,i=0,8)]
            do i=1,k
                t(i)=-real(k+1-i,dp)/10
                t(k+9+i)=1+real(i,dp)/7
            end do
            if (form==2) t(:k)=0.0_dp
            if (form==3) t(k+10:)=1.0_dp
            target=0.5_dp*sum((w(:m-1)*y(:m-1))**2)
            fit=make_splrep(u,y,w=w,k=k,s=target,t=t,bc_type='periodic',status=code)
            call require(code==0 .and. fit%valid(),'periodic arbitrary exterior knots')
            call require(all(fit%knots()==t),'scalar exterior knots retained')
            fp=sum((w(:m-1)*(fit%evaluate(u(:m-1))-y(:m-1)))**2)
            call require(abs(fp-target)<0.0011_dp*target,'periodic active target')
            do d=0,k-1
                scale=max(1.0_dp,maxval(abs(fit%derivative(u,d))))
                call require(abs(fit%derivative(0.0_dp,d)-fit%derivative(1.0_dp,d)) < &
                    1e-8_dp*scale,'explicit endpoint constraints')
            end do
            curve=make_splprep(x,u=u,w=w,k=k,s=target,t=t, &
                bc_type='periodic',parameters=p,status=code)
            call require(code==0 .and. curve%valid(),'shared arbitrary periodic knots')
            call require(all(curve%knots()==t) .and. all(p==u),'shared knots and sites retained')
            v=curve%evaluate(u)
            fp=sum((spread(w(:m-1),2,3)*(v(:m-1,:)-x(:m-1,:)))**2)
            call require(abs(fp-target)<0.0011_dp*target,'shared arbitrary target')
            fit=make_splrep(u,y,w=w,k=k,s=target,t=t,status=code)
            call require(code==0 .and. fit%valid(),'nonclamped nonperiodic scalar')
            curve=make_splprep(x,u=u,w=w,k=k,s=target,t=t,status=code)
            call require(code==0 .and. curve%valid(),'nonclamped nonperiodic shared')
            fit=make_splrep(u,y,w=w,k=k,s=1000.0_dp,t=t,bc_type='periodic',status=code)
            call require(code==0 .and. fit%valid(),'arbitrary exterior constant limit')
            scale=sum(w(:m-1)**2*y(:m-1))/sum(w(:m-1)**2)
            call require(maxval(abs(fit%evaluate(u)-scale))<1e-10_dp,'weighted constant limit')
        end do
    end do
    y=abs(u-0.5_dp)
    x(:,1)=y
    x(:,2)=2*y
    x(:,3)=-0.5_dp*y
    t=[-0.5_dp,0.0_dp,0.001_dp,0.002_dp,0.5_dp,1.0_dp,1.001_dp]
    fit=make_splrep(u,y,k=1,s=0.01_dp,t=t,bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'regularized rank-deficient periodic scalar')
    curve=make_splprep(x,u=u,k=1,s=0.01_dp,t=t,bc_type='periodic',status=code)
    call require(code==0 .and. curve%valid(),'regularized rank-deficient periodic shared')
    t=[0.0_dp,0.0_dp,0.001_dp,0.002_dp,0.5_dp,1.0_dp,1.0_dp]
    fit=make_splrep(u,y,k=1,s=0.01_dp,t=t,status=code)
    call require(code==0 .and. fit%valid(),'regularized rank-deficient clamped scalar')
    curve=make_splprep(x,u=u,k=1,s=0.01_dp,t=t,status=code)
    call require(code==0 .and. curve%valid(),'regularized rank-deficient clamped shared')
    y(20)=y(20)+1.0_dp
    fit=make_splrep(u,y,k=1,s=1e-12_dp,t=t,status=code)
    call require(code==3 .and. .not.fit%valid(),'unattainable rank-deficient target')
    t=[-0.3_dp,-0.2_dp,-0.1_dp,0.0_dp,0.2_dp,0.4_dp,0.6_dp,0.8_dp, &
        1.0_dp,1.2_dp,1.4_dp,1.6_dp]
    fit=make_splrep(u*1e-100_dp,y,k=3,s=1000.0_dp,t=t*1e-100_dp, &
        bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'general tiny physical interval')
    fit=make_splrep(u*1e100_dp,y,k=3,s=1000.0_dp,t=t*1e100_dp, &
        bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'general huge physical interval')
    w(:m-1)=1e-150_dp
    w(m)=1e200_dp
    curve=make_splprep(x*1e150_dp,u=u,w=w,k=3,s=1000.0_dp,t=t, &
        bc_type='periodic',status=code)
    call require(code==0 .and. curve%valid(),'general common scale and ignored closing weight')
    t=[-2.0_dp,-1.7_dp,-1.3_dp,-1.0_dp,0.1_dp,0.6_dp,2.0_dp,2.3_dp,2.5_dp,2.8_dp]
    fit=make_splrep(u,y,k=3,s=1000.0_dp,t=t,xb=-1.0_dp,xe=2.0_dp,status=code)
    call require(code==0 .and. fit%valid(),'general enclosing base interval')
    t=[(0.0_dp,i=1,6),0.001_dp,0.002_dp,0.003_dp,0.1_dp,0.4_dp,0.6_dp,0.9_dp, &
        (1.0_dp,i=1,6)]
    fit=make_splrep(u,y,k=5,s=1000.0_dp,t=t,bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'clustered degree-five constant limit')
    scale=sum(y(:m-1))/real(m-1,dp)
    call require(maxval(abs(fit%coefficients()-scale))<1e-12_dp, &
        'constant limit uses exact polynomial reproduction')
    print '(a,i0,a)','General fixed-knot smoothing: ',checks,' checks passed'
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
end program test_interpolate_fixed_general
