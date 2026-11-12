! SPDX-License-Identifier: MIT
! Compact text protocol for independent SciPy comparisons.
program linalg_rq_driver
    use scifort_kinds, only : dp
    use scifort_linalg_rq, only : rq, rq_r
    implicit none
    real(dp), allocatable :: a(:,:), r(:,:), q(:,:)
    complex(dp), allocatable :: z(:,:), rz(:,:), qz(:,:)
    real(dp) :: xr, xi
    integer :: cases, id, which, m, n, i, j, nq, status
    character(len=8) :: mode

    read(*,*) cases
    do id = 1, cases
        read(*,*) which, mode, m, n
        nq = n
        if (mode == 'economic') nq = min(m,n)
        if (which == 1) then
            allocate(a(m,n),r(m,nq),q(nq,n))
            do j = 1, n
                do i = 1, m
                    read(*,*) a(i,j)
                end do
            end do
            if (mode == 'r') then
                call rq_r(a,r,status)
            else
                call rq(a,r,q,status,mode)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') &
                'CASE',id,status,size(r,1),size(r,2),merge(0,size(q,1),mode=='r')
            do j = 1, size(r,2)
                do i = 1, size(r,1)
                    write(*,'(es26.17e3)') r(i,j)
                end do
            end do
            if (mode /= 'r') then
                do j = 1, size(q,2)
                    do i = 1, size(q,1)
                        write(*,'(es26.17e3)') q(i,j)
                    end do
                end do
            end if
            deallocate(a,r,q)
        else
            allocate(z(m,n),rz(m,nq),qz(nq,n))
            do j = 1, n
                do i = 1, m
                    read(*,*) xr,xi
                    z(i,j) = cmplx(xr,xi,dp)
                end do
            end do
            if (mode == 'r') then
                call rq_r(z,rz,status)
            else
                call rq(z,rz,qz,status,mode)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') &
                'CASE',id,status,size(rz,1),size(rz,2),merge(0,size(qz,1),mode=='r')
            do j = 1, size(rz,2)
                do i = 1, size(rz,1)
                    write(*,'(es26.17e3,1x,es26.17e3)') real(rz(i,j),dp),aimag(rz(i,j))
                end do
            end do
            if (mode /= 'r') then
                do j = 1, size(qz,2)
                    do i = 1, size(qz,1)
                        write(*,'(es26.17e3,1x,es26.17e3)') real(qz(i,j),dp),aimag(qz(i,j))
                    end do
                end do
            end if
            deallocate(z,rz,qz)
        end if
    end do
end program linalg_rq_driver
