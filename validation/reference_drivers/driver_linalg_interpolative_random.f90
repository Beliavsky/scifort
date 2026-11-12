! SPDX-License-Identifier: MIT
! Line-based NumPy/SciPy regression driver for randomized IDs.
program driver_linalg_interpolative_random
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative_random, only : randomized_interp_decomp
    implicit none
    real(dp), allocatable :: a(:, :), proj(:, :)
    complex(dp), allocatable :: z(:, :), projz(:, :)
    integer, allocatable :: idx(:)
    integer(int64) :: seed
    integer :: cases, item, dtype, m, n, mode, target, k, status, i, j, over
    real(dp) :: x, y, eps

    read(*,*) cases
    do item=1,cases
        read(*,*) dtype,m,n,mode,target,eps,over,seed
        allocate(a(m,n),z(m,n))
        do i=1,m
            do j=1,n
                read(*,*) x,y
                a(i,j) = x
                z(i,j) = cmplx(x,y,dp)
            end do
        end do
        if (dtype == 0) then
            if (mode == 0) then
                k = target
                call randomized_interp_decomp(a,k,idx,proj,status,seed=seed,oversampling=over)
            else
                call randomized_interp_decomp(a,eps,k,idx,proj,status,seed=seed,oversampling=over)
            end if
        else
            if (mode == 0) then
                k = target
                call randomized_interp_decomp(z,k,idx,projz,status,seed=seed,oversampling=over)
            else
                call randomized_interp_decomp(z,eps,k,idx,projz,status,seed=seed,oversampling=over)
            end if
        end if
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',item,status,k
        if (status == 0) then
            write(*,'(*(i0,1x))') idx
            if (dtype == 0) then
                do i=1,k
                    write(*,'(*(es26.17e3,1x))') (proj(i,j),j=1,n-k)
                end do
            else
                do i=1,k
                    write(*,'(*(es26.17e3,1x))') (real(projz(i,j),dp),aimag(projz(i,j)),j=1,n-k)
                end do
            end if
        end if
        if (allocated(idx)) deallocate(idx)
        if (allocated(proj)) deallocate(proj)
        if (allocated(projz)) deallocate(projz)
        deallocate(a,z)
    end do
end program driver_linalg_interpolative_random
