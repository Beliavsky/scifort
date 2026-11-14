! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program make_splprep_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splprep,bspline_array_interpolator_type
    implicit none
    integer :: m,ndim,k,provided,nq,mode,nt,io
    do
        read(*,*,iostat=io) m,ndim,k,provided,nq,mode,nt
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:,:),u(:),w(:),q(:),p(:),t(:),c(:,:),v(:,:)
            real(dp) :: target,bounds(2)
            type(bspline_array_interpolator_type) :: fit
            integer :: i,j,d,code
            allocate(x(m,ndim),w(m),q(nq),t(nt))
            read(*,*) target,bounds
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            if (provided/=0) then
                allocate(u(m))
                read(*,*) u
            end if
            if (nt>0) read(*,*) t
            read(*,*) q
            select case(mode)
            case(0)
                fit=make_splprep(x,u=u,k=k,status=code,parameters=p)
            case(1)
                fit=make_splprep(x,u=u,w=w,k=k,s=target,ub=bounds(1),ue=bounds(2), &
                    status=code,parameters=p)
            case(2)
                fit=make_splprep(x,u=u,w=w,k=k,s=target,ub=bounds(1),ue=bounds(2),t=t, &
                    status=code,parameters=p)
            case(3)
                fit=make_splprep(x,u=u,w=w,k=k,s=target,ub=bounds(1),ue=bounds(2), &
                    nest=m+k+1,status=code,parameters=p)
            end select
            print '(i0)',code
            if (code/=0) cycle
            t=fit%knots()
            c=fit%coefficients()
            print '(i0)',size(t)
            print '(*(es26.17e3,1x))',p
            print '(*(es26.17e3,1x))',t
            do j=1,ndim
                print '(*(es26.17e3,1x))',c(:,j)
            end do
            v=fit%evaluate(p)
            do j=1,ndim
                print '(*(es26.17e3,1x))',v(:,j)
            end do
            do d=0,min(3,k)
                v=fit%derivative(q,d)
                do j=1,ndim
                    print '(*(es26.17e3,1x))',v(:,j)
                end do
            end do
            print '(*(es26.17e3,1x))',fit%integrate(bounds(1),bounds(2))
        end block
    end do
end program make_splprep_driver
