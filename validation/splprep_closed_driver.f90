! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program splprep_closed_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : splprep,splev,fitpack_parametric_tck_type
    implicit none
    integer :: m,ndim,k,provided,nq,io
    do
        read(*,*,iostat=io) m,ndim,k,provided,nq
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:,:),u(:),w(:),q(:),values(:,:)
            real(dp) :: target,bounds(2)
            type(fitpack_parametric_tck_type) :: fit
            integer :: i,j,d,status
            allocate(x(m,ndim),u(m),w(m),q(nq))
            read(*,*) target,bounds
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            if (provided/=0) read(*,*) u
            read(*,*) q
            if (provided/=0) then
                fit=splprep(x,u=u,w=w,k=k,s=target,ub=bounds(1),ue=bounds(2),per=.true., &
                    status=status)
            else
                fit=splprep(x,w=w,k=k,s=target,ub=bounds(1),ue=bounds(2),per=.true.,status=status)
            end if
            print '(i0)',status
            if (status/=0) cycle
            print '(i0)',size(fit%t)
            print '(*(es26.17e3,1x))',fit%u
            print '(*(es26.17e3,1x))',fit%t
            do j=1,ndim
                print '(*(es26.17e3,1x))',fit%c(:,j)
            end do
            values=splev(fit%u,fit)
            do j=1,ndim
                print '(*(es26.17e3,1x))',values(:,j)
            end do
            do d=0,min(3,k)
                values=splev(q,fit,der=d)
                do j=1,ndim
                    print '(*(es26.17e3,1x))',values(:,j)
                end do
            end do
        end block
    end do
end program splprep_closed_driver
