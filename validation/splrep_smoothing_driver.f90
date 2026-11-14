! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program splrep_smoothing_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : splrep,splev,splint,fitpack_tck_type
    implicit none
    integer :: n,k,nq,io
    do
        read(*,*,iostat=io) n,k,nq
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:),y(:),w(:),q(:),values(:)
            real(dp) :: target,bounds(2)
            integer :: status,order
            type(fitpack_tck_type) :: fit
            allocate(x(n),y(n),w(n),q(nq))
            read(*,*) target,bounds
            read(*,*) x
            read(*,*) y
            read(*,*) w
            read(*,*) q
            fit=splrep(x,y,w=w,k=k,s=target,bbox=bounds,status=status)
            print '(i0)',status
            if (status/=0) cycle
            print '(i0)',size(fit%t)
            print '(*(es26.17e3,1x))',fit%t
            print '(*(es26.17e3,1x))',fit%c
            values=splev(x,fit)
            print '(*(es26.17e3,1x))',values
            do order=0,min(3,k)
                values=splev(q,fit,der=order)
                print '(*(es26.17e3,1x))',values
            end do
            print '(es26.17e3)',splint(bounds(1),bounds(2),fit)
        end block
    end do
end program splrep_smoothing_driver
