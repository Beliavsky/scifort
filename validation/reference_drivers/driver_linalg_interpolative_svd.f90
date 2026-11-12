! SPDX-License-Identifier: MIT
! Text driver for independent SciPy comparisons of deterministic ID-based SVD.
program driver_linalg_interpolative_svd
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative_svd, only : interpolative_svd
    implicit none
    real(dp), allocatable :: a(:, :), u(:, :), s(:), v(:, :)
    complex(dp), allocatable :: z(:, :), uz(:, :), vz(:, :)
    integer :: cases, testcase, dtype, m, n, mode, rank, k, status, i, j
    real(dp) :: x, y, eps

    read(*,*) cases
    do testcase = 1,cases
        read(*,*) dtype, m, n, mode, rank, eps
        allocate(a(m,n),z(m,n))
        do i = 1,m
            do j = 1,n
                read(*,*) x,y
                a(i,j) = x
                z(i,j) = cmplx(x,y,dp)
            end do
        end do
        if (dtype == 0) then
            if (mode == 0) then
                k = rank
                allocate(u(m,k),s(k),v(n,k))
                call interpolative_svd(a,k,u,s,v,status)
            else
                call interpolative_svd(a,eps,k,u,s,v,status)
            end if
        else
            if (mode == 0) then
                k = rank
                allocate(uz(m,k),s(k),vz(n,k))
                call interpolative_svd(z,k,uz,s,vz,status)
            else
                call interpolative_svd(z,eps,k,uz,s,vz,status)
            end if
        end if
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',testcase,status,k
        if (status == 0) then
            write(*,'(*(es26.17e3,1x))') s
            if (dtype == 0) then
                do i = 1,m
                    write(*,'(*(es26.17e3,1x))') (u(i,j),j=1,k)
                end do
                do i = 1,n
                    write(*,'(*(es26.17e3,1x))') (v(i,j),j=1,k)
                end do
            else
                do i = 1,m
                    write(*,'(*(es26.17e3,1x))') (real(uz(i,j),dp),aimag(uz(i,j)),j=1,k)
                end do
                do i = 1,n
                    write(*,'(*(es26.17e3,1x))') (real(vz(i,j),dp),aimag(vz(i,j)),j=1,k)
                end do
            end if
        end if
        if (allocated(u)) deallocate(u)
        if (allocated(v)) deallocate(v)
        if (allocated(uz)) deallocate(uz)
        if (allocated(vz)) deallocate(vz)
        if (allocated(s)) deallocate(s)
        deallocate(a,z)
    end do
end program driver_linalg_interpolative_svd
