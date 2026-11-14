! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Standalone validation driver for the inspected, separately compiled perazz core.
program splrep_periodic_upstream_driver
    use fitpack_core, only : dp=>fp_real,fp_size,percur,splev,splder,splint
    implicit none
    integer(fp_size) :: m,k,nq,io
    do
        read(*,*,iostat=io) m,k,nq
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:),y(:),w(:),q(:),t(:),c(:),wrk(:),evalwrk(:)
            real(dp), allocatable :: values(:)
            integer(fp_size), allocatable :: iwrk(:)
            real(dp) :: target,bounds(2),fp
            integer(fp_size) :: n,nest,status,order,lwrk
            nest=m+2*k
            lwrk=m*(k+1)+nest*(8+5*k)
            allocate(x(m),y(m),w(m),q(nq),t(nest),c(nest),wrk(lwrk), &
                evalwrk(nest),iwrk(nest),values(max(m,nq)))
            read(*,*) target,bounds
            read(*,*) x
            read(*,*) y
            read(*,*) w
            read(*,*) q
            t=0
            c=0
            wrk=0
            iwrk=0
            n=0
            call percur(0_fp_size,m,x,y,w,k,target,nest,n,t,c,fp, &
                wrk,lwrk,iwrk,status)
            print '(i0)',max(0_fp_size,status)
            if (status>0) cycle
            c(n-k:n)=0
            print '(i0)',n
            print '(*(es26.17e3,1x))',t(:n)
            print '(*(es26.17e3,1x))',c(:n)
            call splev(t,n,c,k,x,values,m,0_fp_size,status)
            if (status/=0) error stop 'upstream sample evaluation'
            print '(*(es26.17e3,1x))',values(:m)
            do order=0,min(3_fp_size,k)
                call splder(t,n,c,k,order,q,values,nq,0_fp_size,evalwrk,status)
                if (status/=0) error stop 'upstream query evaluation'
                print '(*(es26.17e3,1x))',values(:nq)
            end do
            print '(es26.17e3)',splint(t,n,c,k,bounds(1),bounds(2),evalwrk)
        end block
    end do
end program splrep_periodic_upstream_driver
