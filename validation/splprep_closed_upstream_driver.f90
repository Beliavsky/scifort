! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Standalone driver for the separately compiled pinned perazz/FITPACK core.
program splprep_closed_upstream_driver
    use fitpack_core, only : dp=>fp_real,fp_size,clocur,splev,splder
    implicit none
    integer(fp_size) :: m,ndim,k,provided,nq,io
    do
        read(*,*,iostat=io) m,ndim,k,provided,nq
        if (io/=0) exit
        block
            real(dp), allocatable :: x(:,:),u(:),w(:),q(:),values(:),t(:),c(:),wrk(:)
            real(dp), allocatable :: evalwrk(:)
            integer(fp_size), allocatable :: iwrk(:)
            real(dp) :: target,bounds(2),fp
            integer(fp_size) :: i,j,d,n,nest,ncc,lwrk,status,offset
            nest=m+2*k
            ncc=nest*ndim
            lwrk=m*(k+1)+nest*(7+ndim+5*k)
            allocate(x(ndim,m),u(m),w(m),q(nq),values(max(m,nq)),t(nest), &
                c(ncc),wrk(lwrk),evalwrk(nest),iwrk(nest))
            read(*,*) target,bounds
            do i=1,m
                read(*,*) x(:,i)
            end do
            read(*,*) w
            if (provided/=0) read(*,*) u
            read(*,*) q
            x(:,m)=x(:,1)
            t=0
            c=0
            wrk=0
            iwrk=0
            n=0
            call clocur(0_fp_size,provided,ndim,m,u,m*ndim,x,w, &
                k,target,nest,n,t,ncc,c,fp,wrk,lwrk,iwrk,status)
            print '(i0)',max(0_fp_size,status)
            if (status>0) cycle
            print '(i0)',n
            print '(*(es26.17e3,1x))',u
            print '(*(es26.17e3,1x))',t(:n)
            do j=1,ndim
                offset=(j-1)*n
                print '(*(es26.17e3,1x))',c(offset+1:offset+n-k-1)
            end do
            do j=1,ndim
                call splev(t,n,c((j-1)*n+1),k,u,values,m,0_fp_size,status)
                if (status/=0) error stop 'upstream sample evaluation'
                print '(*(es26.17e3,1x))',values(:m)
            end do
            do d=0,min(3_fp_size,k)
                do j=1,ndim
                    call splder(t,n,c((j-1)*n+1),k,d,q,values,nq,0_fp_size,evalwrk,status)
                    if (status/=0) error stop 'upstream query evaluation'
                    print '(*(es26.17e3,1x))',values(:nq)
                end do
            end do
        end block
    end do
end program splprep_closed_upstream_driver
