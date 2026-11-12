! SPDX-License-Identifier: MIT
! Streaming protocol for independent SciPy comparison of LU, det, and inv.
program linalg_lu_det_inv_driver
    use scifort_kinds, only : dp
    use scifort_linalg_lu, only : lu_indices
    use scifort_linalg_det_inv, only : det, inv
    implicit none
    real(dp), allocatable :: a(:, :), l(:, :), u(:, :), ai(:, :)
    complex(dp), allocatable :: z(:, :), lz(:, :), uz(:, :), zi(:, :)
    integer, allocatable :: perm(:)
    integer :: cases, idx, kind_id, m, n, k, i, j, status_lu, status_det, status_inv
    real(dp) :: value, xr, xi
    complex(dp) :: zvalue

    read(*,*) cases
    do idx = 1, cases
        read(*,*) kind_id, m, n
        k = min(m,n)
        allocate(perm(m))
        if (kind_id == 1) then
            allocate(a(m,n),l(m,k),u(k,n),ai(m,m))
            do j = 1, n
                do i = 1, m
                    read(*,*) a(i,j)
                end do
            end do
            call lu_indices(a,perm,l,u,status_lu)
            status_det = -1
            status_inv = -1
            value = 0.0_dp
            if (m == n) then
                call det(a,value,status_det)
                call inv(a,ai,status_inv)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0)') &
                'CASE',idx,status_lu,status_det,status_inv
            do i = 1, m
                write(*,'(i0)') perm(i)
            end do
            do j = 1, k
                do i = 1, m
                    write(*,'(es26.17e3)') l(i,j)
                end do
            end do
            do j = 1, n
                do i = 1, k
                    write(*,'(es26.17e3)') u(i,j)
                end do
            end do
            if (m == n) then
                write(*,'(es26.17e3)') value
                do j = 1, n
                    do i = 1, m
                        write(*,'(es26.17e3)') ai(i,j)
                    end do
                end do
            end if
            deallocate(a,l,u,ai)
        else
            allocate(z(m,n),lz(m,k),uz(k,n),zi(m,m))
            do j = 1, n
                do i = 1, m
                    read(*,*) xr, xi
                    z(i,j) = cmplx(xr,xi,dp)
                end do
            end do
            call lu_indices(z,perm,lz,uz,status_lu)
            status_det = -1
            status_inv = -1
            zvalue = cmplx(0.0_dp,0.0_dp,dp)
            if (m == n) then
                call det(z,zvalue,status_det)
                call inv(z,zi,status_inv)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0)') &
                'CASE',idx,status_lu,status_det,status_inv
            do i = 1, m
                write(*,'(i0)') perm(i)
            end do
            do j = 1, k
                do i = 1, m
                    write(*,'(es26.17e3,1x,es26.17e3)') real(lz(i,j),dp),aimag(lz(i,j))
                end do
            end do
            do j = 1, n
                do i = 1, k
                    write(*,'(es26.17e3,1x,es26.17e3)') real(uz(i,j),dp),aimag(uz(i,j))
                end do
            end do
            if (m == n) then
                write(*,'(es26.17e3,1x,es26.17e3)') real(zvalue,dp),aimag(zvalue)
                do j = 1, n
                    do i = 1, m
                        write(*,'(es26.17e3,1x,es26.17e3)') real(zi(i,j),dp),aimag(zi(i,j))
                    end do
                end do
            end if
            deallocate(z,lz,uz,zi)
        end if
        deallocate(perm)
    end do
end program linalg_lu_det_inv_driver
