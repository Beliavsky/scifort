! SPDX-License-Identifier: MIT
! Line-oriented driver for NumPy/SciPy randomized interpolative SVD tests.
program driver_linalg_interpolative_svd_random
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative_svd_random, only : randomized_interpolative_svd
    implicit none
    real(dp), allocatable :: a(:, :), u(:, :), v(:, :), s(:)
    complex(dp), allocatable :: z(:, :), uz(:, :), vz(:, :)
    integer(int64) :: seed
    integer :: cases, testcase, dtype, m, n, mode, target, k, status, i, j, over
    real(dp) :: x, y, eps

    read(*,*) cases
    do testcase = 1, cases
        read(*,*) dtype, m, n, mode, target, eps, over, seed
        allocate(a(m,n),z(m,n))
        do i = 1, m
            do j = 1, n
                read(*,*) x, y
                a(i,j) = x
                z(i,j) = cmplx(x,y,dp)
            end do
        end do
        if (dtype == 0) then
            if (mode == 0) then
                k = target
                allocate(u(m,k),s(k),v(n,k))
                call randomized_interpolative_svd(a,k,u,s,v,status,seed=seed,oversampling=over)
            else
                call randomized_interpolative_svd(a,eps,k,u,s,v,status,seed=seed,oversampling=over)
            end if
        else
            if (mode == 0) then
                k = target
                allocate(uz(m,k),s(k),vz(n,k))
                call randomized_interpolative_svd(z,k,uz,s,vz,status,seed=seed,oversampling=over)
            else
                call randomized_interpolative_svd(z,eps,k,uz,s,vz,status,seed=seed,oversampling=over)
            end if
        end if
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',testcase,status,k
        if (status == 0) then
            write(*,'(*(es26.17e3,1x))') s
            if (dtype == 0) then
                do i = 1, m
                    write(*,'(*(es26.17e3,1x))') (u(i,j),j=1,k)
                end do
                do i = 1, n
                    write(*,'(*(es26.17e3,1x))') (v(i,j),j=1,k)
                end do
            else
                do i = 1, m
                    write(*,'(*(es26.17e3,1x))') (real(uz(i,j),dp),aimag(uz(i,j)),j=1,k)
                end do
                do i = 1, n
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
end program driver_linalg_interpolative_svd_random
