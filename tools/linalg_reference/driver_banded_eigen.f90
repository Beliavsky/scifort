! SPDX-License-Identifier: MIT
! Batch protocol for the scipy.linalg eigenvalue comparison harness.
program driver_banded_eigen
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eigvalsh_tridiagonal, &
        eig_banded, eigvals_banded
    implicit none
    integer :: nc, icase, kind, n, kd, i, j, status, status2
    real(dp), allocatable :: d(:), e(:), w(:), w2(:), zr(:,:), br(:,:)
    complex(dp), allocatable :: zc(:,:), bc(:,:)
    real(dp) :: re, im
    logical :: lower

    read (*,*) nc
    do icase = 1, nc
        read (*,*) kind,n,kd
        allocate(w(n),w2(n))
        if (kind == 1) then
            allocate(d(n),e(max(0,n-1)),zr(n,n))
            do i=1,n
                read (*,*) re,im
                d(i)=re
            end do
            do i=1,n-1
                read (*,*) re,im
                e(i)=re
            end do
            call eigh_tridiagonal(d,e,w,status,zr)
            if (status == 0) then
                call eigvalsh_tridiagonal(d,e,w2,status2)
                if (status2 /= status .or. any(w /= w2)) status = -3
            end if
        else if (kind == 2 .or. kind == 3) then
            allocate(br(kd+1,n),zr(n,n))
            do j=1,n
                do i=1,kd+1
                    read (*,*) re,im
                    br(i,j)=re
                end do
            end do
            lower = kind == 3
            call eig_banded(br,w,status,zr,lower=lower)
            if (status == 0) then
                call eigvals_banded(br,w2,status2,lower=lower)
                if (status2 /= status .or. any(w /= w2)) status = -3
            end if
        else
            allocate(bc(kd+1,n),zc(n,n))
            do j=1,n
                do i=1,kd+1
                    read (*,*) re,im
                    bc(i,j)=cmplx(re,im,dp)
                end do
            end do
            lower = kind == 5
            call eig_banded(bc,w,status,zc,lower=lower)
            if (status == 0) then
                call eigvals_banded(bc,w2,status2,lower=lower)
                if (status2 /= status .or. any(w /= w2)) status = -3
            end if
        end if
        print '(i0)',status
        if (status == 0) then
            do i=1,n
                print '(es26.17e3)',w(i)
            end do
            do j=1,n
                do i=1,n
                    if (kind < 4) then
                        print '(2(es26.17e3,1x))',zr(i,j),0.0_dp
                    else
                        print '(2(es26.17e3,1x))',real(zc(i,j),dp),aimag(zc(i,j))
                    end if
                end do
            end do
        end if
        if (allocated(d)) deallocate(d,e)
        if (allocated(zr)) deallocate(zr)
        if (allocated(zc)) deallocate(zc)
        if (allocated(br)) deallocate(br)
        if (allocated(bc)) deallocate(bc)
        deallocate(w,w2)
    end do
end program driver_banded_eigen
