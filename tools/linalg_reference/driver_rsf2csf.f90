! SPDX-License-Identifier: MIT
! Batch Fortran interface for reference comparisons of rsf2csf.
program driver_rsf2csf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_schur, only : schur, rsf2csf
    implicit none
    integer :: cases, icase, n, mode, i, j, status
    real(dp), allocatable :: a(:,:), tr(:,:), zr(:,:)
    complex(dp), allocatable :: tc(:,:), zc(:,:)

    read(*,*) cases
    do icase = 1, cases
        read(*,*) n, mode
        allocate(a(n,n),tr(n,n),zr(n,n),tc(n,n),zc(n,n))
        if (mode == 0) then
            do j = 1, n
                do i = 1, n
                    read(*,*) tr(i,j)
                end do
            end do
            do j = 1, n
                do i = 1, n
                    read(*,*) zr(i,j)
                end do
            end do
            status = linalg_status_success
        else
            do j = 1, n
                do i = 1, n
                    read(*,*) a(i,j)
                end do
            end do
            call schur(a,tr,zr,status)
        end if
        if (status == linalg_status_success) call rsf2csf(tr,zr,tc,zc,status)
        write(*,'(i0)') status
        if (status == linalg_status_success) then
            do j = 1,n
                do i = 1,n
                    write(*,'(4(es31.17e3,1x))') real(tc(i,j),dp), aimag(tc(i,j)), &
                        real(zc(i,j),dp), aimag(zc(i,j))
                end do
            end do
        end if
        deallocate(a,tr,zr,tc,zc)
    end do
end program driver_rsf2csf
