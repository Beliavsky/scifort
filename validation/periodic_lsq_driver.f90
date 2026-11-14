! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program periodic_lsq_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splrep,splprep,splev,splint,bspline, &
        fitpack_tck_type,fitpack_parametric_tck_type,bspline_array_interpolator_type
    implicit none
    integer :: m,ndim,k,scalar,provided,nq,ni,io
    do
        read(*,*,iostat=io) m,ndim,k,scalar,provided,nq,ni
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:,:),u(:),w(:),q(:),interior(:),input(:),t(:),c(:,:)
            real(dp), allocatable :: p(:),v(:,:)
            type(bspline_array_interpolator_type) :: object
            type(fitpack_tck_type) :: fit
            type(fitpack_parametric_tck_type) :: curve
            integer :: i,j,d,code
            allocate(x(m,ndim),w(m),q(nq),interior(ni))
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            if (provided/=0) then
                allocate(u(m))
                read(*,*) u
            end if
            read(*,*) interior
            read(*,*) q
            if (scalar/=0) then
                fit=splrep(u,x(:,1),w=w,k=k,task=-1,t=interior,per=.true.,status=code)
                if (code==0) then
                    t=fit%t
                    c=reshape(fit%c(:size(t)-k-1),[size(t)-k-1,1])
                    p=u
                end if
            else
                input=[(-13.0_dp,i=1,k+1),interior,(29.0_dp,i=1,k+1)]
                curve=splprep(x,u=u,w=w,k=k,task=-1,t=input,per=.true.,status=code)
                if (code==0) then
                    t=curve%t
                    c=curve%c
                    p=curve%u
                end if
            end if
            print '(i0)',code
            if (code/=0) cycle
            print '(i0)',size(t)
            print '(*(es26.17e3,1x))',p
            print '(*(es26.17e3,1x))',t
            do j=1,ndim
                print '(*(es26.17e3,1x))',c(:,j)
            end do
            if (scalar/=0) then
                print '(*(es26.17e3,1x))',splev(p,fit)
                do d=0,k
                    print '(*(es26.17e3,1x))',splev(q,fit,der=d)
                end do
                print '(*(es26.17e3,1x))',splint(p(1),p(m),fit)
            else
                v=splev(p,curve)
                do j=1,ndim
                    print '(*(es26.17e3,1x))',v(:,j)
                end do
                do d=0,k
                    v=splev(q,curve,der=d)
                    do j=1,ndim
                        print '(*(es26.17e3,1x))',v(:,j)
                    end do
                end do
                object=bspline(curve%t,curve%c,k)
                print '(*(es26.17e3,1x))',object%integrate(p(1),p(m))
            end if
        end block
    end do
end program periodic_lsq_driver
