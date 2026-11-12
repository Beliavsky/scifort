! SPDX-License-Identifier: MIT
! Test-input/output driver for independent scipy.sparse.coo_array comparisons.
program sparse_coo_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_to_dense, csc_to_dense
    use scifort_sparse_coo, only : coo_matrix, coo_from_triplets, coo_from_csr, &
        coo_from_csc, coo_to_csr, coo_to_csc, coo_to_dense, &
        coo_matvec, coo_sum_duplicates, coo_eliminate_zeros, sparse_find
    implicit none
    type(coo_matrix) :: a, b, c
    type(csr_matrix) :: csr
    type(csc_matrix) :: csc
    integer, allocatable :: rr(:), cc(:), rfind(:), cfind(:)
    real(dp), allocatable :: values(:), x(:), y(:), dense(:,:), found(:)
    integer :: trials, icase, m, n, nnz, k, i, j, info

    read(*,*) trials
    do icase = 1, trials
        read(*,*) m, n, nnz
        allocate(rr(nnz),cc(nnz),values(nnz),x(n))
        do k = 1, nnz
            read(*,*) rr(k),cc(k),values(k)
        end do
        if (n > 0) read(*,*) x
        a = coo_from_triplets(m,n,rr,cc,values,info)
        call require(info == 0, 'coo from triplets')
        dense = coo_to_dense(a,info)
        call require(info == 0, 'coo to dense')
        csr = coo_to_csr(a,info)
        call require(info == 0, 'coo to CSR')
        csc = coo_to_csc(a,info)
        call require(info == 0, 'coo to CSC')
        call require(all(csr_to_dense(csr) == dense), 'CSR numerical equality')
        call require(all(csc_to_dense(csc) == dense), 'CSC numerical equality')
        b = coo_from_csr(csr,info)
        call require(info == 0 .and. all(coo_to_dense(b) == dense), 'CSR round trip')
        c = coo_from_csc(csc,info)
        call require(info == 0 .and. all(coo_to_dense(c) == dense), 'CSC round trip')
        y = coo_matvec(a,x,info)
        call require(info == 0, 'coo matvec')
        call sparse_find(a,rfind,cfind,found,info)
        call require(info == 0, 'coo find')
        b = a
        call coo_sum_duplicates(b,info)
        call require(info == 0, 'coo sum duplicates')
        write(*,'(A,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
            'CASE',icase,m,n,nnz,size(b%data),size(found)
        do k = 1, size(b%data)
            write(*,'(A,1X,I0,1X,I0,1X,ES26.17E3)') &
                'SUM',b%row(k),b%col(k),b%data(k)
        end do
        do k = 1, size(found)
            write(*,'(A,1X,I0,1X,I0,1X,ES26.17E3)') &
                'FIND',rfind(k),cfind(k),found(k)
        end do
        do i = 1,m
            do j = 1,n
                write(*,'(A,1X,I0,1X,I0,1X,ES26.17E3)') 'D',i,j,dense(i,j)
            end do
            write(*,'(A,1X,I0,1X,ES26.17E3)') 'Y',i,y(i)
        end do
        call coo_eliminate_zeros(b,info)
        call require(info == 0, 'coo eliminate zeros')
        write(*,'(A,1X,I0)') 'NONZERO',size(b%data)
        deallocate(rr,cc,values,x,rfind,cfind,found,y,dense)
    end do
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program sparse_coo_reference_driver
