! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program make_splrep_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,bspline_interpolator_type
    implicit none
    integer :: m,k,mode,nt,capacity,nq,io
    do
        read(*,*,iostat=io) m,k,mode,nt,capacity,nq
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:),y(:),w(:),t(:),q(:),values(:)
            real(dp) :: target,bounds(2)
            type(bspline_interpolator_type) :: fit
            integer :: code,d
            allocate(x(m),y(m),w(m),t(nt),q(nq))
            read(*,*) target,bounds
            read(*,*) x
            read(*,*) y
            read(*,*) w
            if (nt>0) read(*,*) t
            read(*,*) q
            select case(mode)
            case(0)
                fit=make_splrep(x,y,k=k,status=code)
            case(1)
                fit=make_splrep(x,y,w=w,k=k,s=target,xb=bounds(1),xe=bounds(2),status=code)
            case(2)
                fit=make_splrep(x,y,w=w,k=k,s=target,xb=bounds(1),xe=bounds(2), &
                    t=t,status=code)
            case(3)
                fit=make_splrep(x,y,w=w,k=k,s=target,xb=bounds(1),xe=bounds(2), &
                    nest=capacity,status=code)
            end select
            print '(i0)',code
            if (code/=0) cycle
            t=fit%knots()
            values=fit%coefficients()
            print '(i0)',size(t)
            print '(*(es26.17e3,1x))',t
            print '(*(es26.17e3,1x))',values
            print '(*(es26.17e3,1x))',fit%evaluate(x)
            do d=0,min(3,k)
                print '(*(es26.17e3,1x))',fit%derivative(q,d)
            end do
            print '(es26.17e3)',fit%integrate(bounds(1),bounds(2))
        end block
    end do
end program make_splrep_driver
