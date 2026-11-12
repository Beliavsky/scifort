! SPDX-License-Identifier: MIT
! Batch protocol: mode n first last low high; d(1:n); e(1:n-1).
! Modes 1/2 return selected eigenvalues, 3/4 include selected eigenvectors.
! Odd modes use inclusive one-based indices; even modes use (low,high].
program driver_tridiagonal_selection
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eigvalsh_tridiagonal
    implicit none

    integer :: nc, icase, mode, n, first, last, i, j, status, count
    real(dp) :: low, high
    real(dp), allocatable :: d(:), e(:), w(:), z(:, :)

    read (*,*) nc
    do icase = 1, nc
        read (*,*) mode, n, first, last, low, high
        allocate(d(n), e(max(0,n-1)), w(n))
        read (*,*) d
        if (n > 1) read (*,*) e
        if (mode >= 3) allocate(z(n,n))
        select case(mode)
        case(1)
            call eigvalsh_tridiagonal(d,e,w,status,subset_by_index=[first,last], &
                num_selected=count)
        case(2)
            call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[low,high], &
                num_selected=count)
        case(3)
            call eigh_tridiagonal(d,e,w,status,z,subset_by_index=[first,last], &
                num_selected=count)
        case(4)
            call eigh_tridiagonal(d,e,w,status,z,subset_by_value=[low,high], &
                num_selected=count)
        case default
            error stop 'invalid selection comparison mode'
        end select
        print '(i0,1x,i0)',status,count
        if (status == 0) then
            do i = 1, count
                print '(es26.17e3)',w(i)
            end do
            if (mode >= 3) then
                do j = 1, count
                    do i = 1, n
                        print '(es26.17e3)',z(i,j)
                    end do
                end do
            end if
        end if
        if (allocated(z)) deallocate(z)
        deallocate(d,e,w)
    end do
end program driver_tridiagonal_selection
