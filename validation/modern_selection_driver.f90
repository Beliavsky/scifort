! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program modern_selection_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep, &
        bspline_interpolator_type,bspline_array_interpolator_type
    implicit none
    integer :: m,ndim,k,periodic,scalar,provided,capacity,nq,io
    do
        read(*,*,iostat=io) m,ndim,k,periodic,scalar,provided,capacity,nq
        if (io/=0) exit
        block
            type(bspline_interpolator_type) :: fit
            type(bspline_array_interpolator_type) :: curve
            real(dp), allocatable :: x(:,:),w(:),u(:),q(:),p(:),t(:),c(:,:),v(:,:)
            real(dp) :: s,bounds(2)
            integer :: i,j,d,code
            character(10) :: boundary
            logical :: valid
            allocate(x(m,ndim),w(m),q(nq))
            read(*,*) s,bounds
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            if (provided/=0) then
                allocate(u(m))
                read(*,*) u
            end if
            read(*,*) q
            boundary='not-a-knot'
            if (periodic/=0) boundary='periodic'
            if (scalar/=0) then
                if (s==0.0_dp) then
                    fit=make_splrep(u,x(:,1),k=k,bc_type=boundary,status=code)
                else if (capacity==0) then
                    fit=make_splrep(u,x(:,1),w=w,k=k,s=s,xb=bounds(1),xe=bounds(2), &
                        bc_type=boundary,status=code)
                else
                    fit=make_splrep(u,x(:,1),w=w,k=k,s=s,xb=bounds(1),xe=bounds(2), &
                        nest=capacity,bc_type=boundary,status=code)
                end if
                valid=fit%valid()
                if (code==0) then
                    p=u
                    t=fit%knots()
                    c=reshape(fit%coefficients(),[size(t)-k-1,1])
                end if
            else
                if (s==0.0_dp) then
                    curve=make_splprep(x,u=u,k=k,bc_type=boundary,status=code,parameters=p)
                else if (capacity==0) then
                    curve=make_splprep(x,u=u,w=w,k=k,s=s,ub=bounds(1),ue=bounds(2), &
                        bc_type=boundary,status=code,parameters=p)
                else
                    curve=make_splprep(x,u=u,w=w,k=k,s=s,ub=bounds(1),ue=bounds(2), &
                        nest=capacity,bc_type=boundary,status=code,parameters=p)
                end if
                valid=curve%valid()
                if (code==0) then
                    t=curve%knots()
                    c=curve%coefficients()
                end if
            end if
            write(*,'(3(i0,1x))') code,merge(1,0,valid),merge(1,0,allocated(p))
            if (code/=0) cycle
            write(*,'(i0)') size(t)
            write(*,'(*(es26.17e3,1x))') p
            write(*,'(*(es26.17e3,1x))') t
            do j=1,ndim
                write(*,'(*(es26.17e3,1x))') c(:,j)
            end do
            if (scalar/=0) then
                v=reshape(fit%evaluate(p),[m,1])
            else
                v=curve%evaluate(p)
            end if
            do j=1,ndim
                write(*,'(*(es26.17e3,1x))') v(:,j)
            end do
            do d=0,k
                if (scalar/=0) then
                    v=reshape(fit%derivative(q,d),[nq,1])
                else
                    v=curve%derivative(q,d)
                end if
                do j=1,ndim
                    write(*,'(*(es26.17e3,1x))') v(:,j)
                end do
            end do
        end block
    end do
end program modern_selection_driver
