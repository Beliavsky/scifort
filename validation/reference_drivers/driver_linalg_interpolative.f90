! SPDX-License-Identifier: MIT
! Batch protocol for test/compare_linalg_interpolative_scipy.py.
program driver_linalg_interpolative
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative, only : interp_decomp, id_to_svd
    implicit none
    real(dp), allocatable :: a(:, :), proj(:, :)
    complex(dp), allocatable :: z(:, :), projz(:, :), uz(:, :), vz(:, :)
    real(dp), allocatable :: u(:, :), v(:, :), s(:)
    integer, allocatable :: idx(:)
    real(dp) :: x, y, eps
    integer :: number, testcase, dtype, m, n, mode, k, rank, i, j, status

    read(*,*) number
    do testcase = 1, number
        read(*,*) dtype, m, n, mode, rank, eps
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
                call interp_decomp(a,rank,idx,proj,status)
                k = rank
            else
                call interp_decomp(a,eps,k,idx,proj,status)
            end if
        else
            if (mode == 0) then
                call interp_decomp(z,rank,idx,projz,status)
                k = rank
            else
                call interp_decomp(z,eps,k,idx,projz,status)
            end if
        end if
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',testcase,status,k
        if (status == 0) then
            write(*,'(*(i0,1x))') idx
            do j = 1, n-k
                do i = 1, k
                    if (dtype == 0) then
                        write(*,'(2(es26.17e3,1x))') proj(i,j),0.0_dp
                    else
                        write(*,'(2(es26.17e3,1x))') real(projz(i,j),dp),aimag(projz(i,j))
                    end if
                end do
            end do
        end if
        if (status == 0 .and. k > 0 .and. mod(testcase,7) == 0) then
            allocate(s(k))
            if (dtype == 0) then
                allocate(u(m,k),v(n,k))
                call id_to_svd(a(:,idx(1:k)),idx,proj,u,s,v,status)
                deallocate(u,v)
            else
                allocate(uz(m,k),vz(n,k))
                call id_to_svd(z(:,idx(1:k)),idx,projz,uz,s,vz,status)
                deallocate(uz,vz)
            end if
            write(*,'(a,1x,i0)') 'SVD',status
            if (status == 0) write(*,'(*(es26.17e3,1x))') s
            deallocate(s)
        end if
        if (allocated(idx)) deallocate(idx)
        if (allocated(proj)) deallocate(proj)
        if (allocated(projz)) deallocate(projz)
        deallocate(a,z)
    end do
end program driver_linalg_interpolative
