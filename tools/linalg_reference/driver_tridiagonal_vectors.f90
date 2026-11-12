! SPDX-License-Identifier: MIT
! Batch protocol: number of cases; per case n first last, diagonal, offdiagonal.
! Output: status/count followed by count eigenvalues and n*count eigenvector elements.
program driver_tridiagonal_vectors
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal
    implicit none

    integer :: nc, icase, n, first, last, i, j, status, count, k
    real(dp), allocatable :: d(:), e(:), w(:), z(:, :)

    read (*,*) nc
    do icase = 1, nc
        read (*,*) n, first, last
        k = last-first+1
        allocate(d(n),e(n-1),w(k),z(n,k))
        read (*,*) d
        if (n > 1) read (*,*) e
        call eigh_tridiagonal(d,e,w,status,z,subset_by_index=[first,last], &
            num_selected=count,eigenvector_method='inverse')
        print '(i0,1x,i0)', status, count
        if (status == 0) then
            do i = 1, count
                print '(es26.17e3)',w(i)
            end do
            do j = 1, count
                do i = 1, n
                    print '(es26.17e3)',z(i,j)
                end do
            end do
        end if
        deallocate(d,e,w,z)
    end do
end program driver_tridiagonal_vectors
