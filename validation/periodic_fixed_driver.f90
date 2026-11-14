! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program periodic_fixed_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep,bspline, &
        bspline_interpolator_type,bspline_array_interpolator_type
    implicit none
    integer :: m,ndim,k,scalar,provided,nq,capacity,nt,io
    do
        read(*,*,iostat=io) m,ndim,k,scalar,provided,nq,capacity,nt
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:,:),u(:),w(:),weights(:),q(:),p(:),t(:),fixed(:),c(:,:),v(:,:)
            integer, allocatable :: cap
            real(dp) :: target
            type(bspline_interpolator_type) :: fit
            type(bspline_array_interpolator_type) :: curve
            integer :: i,j,d,code
            allocate(x(m,ndim),w(m),q(nq))
            read(*,*) target
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            if (provided/=0) then
                allocate(u(m))
                read(*,*) u
            end if
            read(*,*) q
            allocate(fixed(nt))
            read(*,*) fixed
            if (target>0.0_dp) weights=w
            if (capacity>0) cap=capacity
            if (scalar/=0) then
                fit=make_splrep(u,x(:,1),w=weights,k=k,s=target,nest=cap, &
                    t=fixed,bc_type='periodic',status=code)
                if (code==0) then
                    t=fit%knots()
                    c=reshape(fit%coefficients(),[size(t)-k-1,1])
                    curve=bspline(t,c,k,periodic=target==0.0_dp .and. k>1)
                    p=u
                end if
            else
                curve=make_splprep(x,u=u,w=weights,k=k,s=target,nest=cap, &
                    t=fixed,bc_type='periodic',status=code,parameters=p)
            end if
            print '(i0)',code
            if (code/=0) cycle
            t=curve%knots()
            c=curve%coefficients()
            print '(i0)',size(t)
            print '(*(es26.17e3,1x))',p
            print '(*(es26.17e3,1x))',t
            do j=1,ndim
                print '(*(es26.17e3,1x))',c(:,j)
            end do
            v=curve%evaluate(p)
            do j=1,ndim
                print '(*(es26.17e3,1x))',v(:,j)
            end do
            do d=0,min(3,k)
                v=curve%derivative(q,d)
                do j=1,ndim
                    print '(*(es26.17e3,1x))',v(:,j)
                end do
            end do
            print '(*(es26.17e3,1x))',curve%integrate(p(1),p(m))
        end block
    end do
end program periodic_fixed_driver
