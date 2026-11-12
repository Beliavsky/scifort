! SPDX-License-Identifier: MIT
! Text protocol used by compare_sparse_scipy.py for independent SciPy checks.
program sparse_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_to_dense, &
        csr_transpose, csr_matmul, csr_add, csr_elementwise_multiply, csr_matvec, &
        sparse_success
    implicit none
    type(csr_matrix) :: a, b, c, product, summed, hadamard, trans
    integer :: cases, trial, m, k, n, status, i
    real(dp), allocatable :: x(:), y(:)

    read (*,*) cases
    do trial = 1, cases
        read (*,*) m, k, n
        call read_matrix(m,k,a)
        call read_matrix(k,n,b)
        call read_matrix(m,k,c)
        trans = csr_transpose(a,status)
        call ensure(status == sparse_success, 'transpose')
        product = csr_matmul(a,b,status)
        call ensure(status == sparse_success, 'matmul')
        summed = csr_add(a,c,status,alpha=2.0_dp,beta=-3.0_dp)
        call ensure(status == sparse_success, 'add')
        hadamard = csr_elementwise_multiply(a,c,status)
        call ensure(status == sparse_success, 'hadamard')
        allocate(x(k))
        do i = 1, k
            x(i) = real(mod(i,7)-3,dp)/3.0_dp
        end do
        y = csr_matvec(a,x,status)
        call ensure(status == sparse_success, 'matvec')
        print '(a,1x,i0)', 'CASE', trial
        call emit('A',csr_to_dense(a))
        call emit('B',csr_to_dense(b))
        call emit('C',csr_to_dense(c))
        call emit('T',csr_to_dense(trans))
        call emit('P',csr_to_dense(product))
        call emit('S',csr_to_dense(summed))
        call emit('H',csr_to_dense(hadamard))
        call emit('V',reshape(y,[size(y),1]))
        deallocate(x)
    end do
contains
    subroutine ensure(ok,message)
        logical, intent(in) :: ok !! Whether a checked operation succeeded.
        character(len=*), intent(in) :: message !! Operation name on failure.
        if (.not. ok) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine ensure

    subroutine read_matrix(n_rows,n_cols,mat)
        integer, intent(in) :: n_rows !! Nonnegative row dimension.
        integer, intent(in) :: n_cols !! Nonnegative column dimension.
        type(csr_matrix), intent(out) :: mat !! Canonical CSR constructed from triplets.
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: nnz, j, code

        read (*,*) nnz
        allocate(rows(nnz),cols(nnz),vals(nnz))
        do j = 1, nnz
            read (*,*) rows(j),cols(j),vals(j)
        end do
        mat = csr_from_coo(n_rows,n_cols,rows,cols,vals,code)
        call ensure(code == sparse_success, 'COO input')
    end subroutine read_matrix

    subroutine emit(label, values)
        character(len=*), intent(in) :: label !! Output matrix identifier.
        real(dp), intent(in) :: values(:,:) !! Dense result for independent comparison.
        integer :: j

        write (*,'(a,1x,i0,1x,i0)') label, size(values,1), size(values,2)
        do j = 1, size(values,1)
            write (*,'(*(es26.17e3,1x))') values(j,:)
        end do
    end subroutine emit
end program sparse_reference_driver
