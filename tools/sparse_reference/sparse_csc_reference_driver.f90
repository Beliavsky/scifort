! SPDX-License-Identifier: MIT
! Reads randomized sparse matrices and diagonal arguments for SciPy comparison.
program sparse_csc_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, csr_to_dense, &
        csr_to_csc, csc_from_coo, csc_to_csr, csc_to_dense, csc_transpose, &
        csr_diags, csr_diagonal, csr_setdiag, csr_triu, csr_tril, sparse_success
    implicit none
    type(csr_matrix) :: a, b, d
    type(csc_matrix) :: c, t
    integer, allocatable :: rows(:), cols(:), offs(:), ptr(:)
    real(dp), allocatable :: values(:), diag_values(:), replacement(:), diagonal(:)
    integer :: ncases, caseid, m, n, nnz, n_diag, nvals, k, i, ierr, k_cut, k_set, k_diag

    read(*,*) ncases
    do caseid = 1, ncases
        read(*,*) m, n, nnz
        allocate(rows(nnz),cols(nnz),values(nnz))
        do i = 1, nnz
            read(*,*) rows(i), cols(i), values(i)
        end do
        read(*,*) n_diag
        allocate(offs(n_diag),ptr(n_diag+1))
        if (n_diag > 0) read(*,*) offs
        read(*,*) ptr
        nvals = ptr(n_diag+1)-1
        allocate(diag_values(nvals))
        do i = 1, nvals
            read(*,*) diag_values(i)
        end do
        read(*,*) k_cut, k_set, k_diag
        read(*,*) nvals
        allocate(replacement(nvals))
        do i = 1, nvals
            read(*,*) replacement(i)
        end do

        a = csr_from_coo(m,n,rows,cols,values,ierr)
        call ok(ierr, 'csr_from_coo')
        c = csr_to_csc(a,ierr)
        call ok(ierr, 'csr_to_csc')
        b = csc_to_csr(c,ierr)
        call ok(ierr, 'csc_to_csr')
        t = csc_transpose(c,ierr)
        call ok(ierr, 'csc_transpose')
        d = csr_diags(m,n,offs,ptr,diag_values,ierr)
        call ok(ierr, 'csr_diags')
        write(*,'(a,1x,i0)') 'CASE',caseid
        call output_matrix('A',csr_to_dense(a))
        call output_matrix('C',csc_to_dense(c))
        call output_matrix('R',csr_to_dense(b))
        call output_matrix('T',csc_to_dense(t))
        call output_matrix('B',csr_to_dense(d))
        b = csr_triu(a,k_cut,ierr)
        call ok(ierr, 'csr_triu')
        call output_matrix('U',csr_to_dense(b))
        b = csr_tril(a,k_cut,ierr)
        call ok(ierr, 'csr_tril')
        call output_matrix('L',csr_to_dense(b))
        diagonal = csr_diagonal(a,k_diag,ierr)
        call ok(ierr, 'csr_diagonal')
        call output_matrix('D',reshape(diagonal,[1,size(diagonal)]))
        call csr_setdiag(a,replacement,k_set,ierr)
        call ok(ierr, 'csr_setdiag')
        call output_matrix('S',csr_to_dense(a))
        deallocate(rows,cols,values,offs,ptr,diag_values,replacement,diagonal)
    end do
contains
    subroutine ok(code, title)
        integer, intent(in) :: code !! Status code.
        character(len=*), intent(in) :: title !! Operation being checked.
        if (code /= sparse_success) then
            write(*,'(a,1x,a,1x,i0)') 'FAIL',title,code
            error stop 'sparse_csc_reference_driver failure'
        end if
    end subroutine ok

    subroutine output_matrix(label, a)
        character(len=*), intent(in) :: label !! Short operation identifier.
        real(dp), intent(in) :: a(:,:) !! Matrix printed row-by-row.
        integer :: row
        write(*,'(a,1x,i0,1x,i0)') label,size(a,1),size(a,2)
        do row = 1, size(a,1)
            write(*,'(*(es24.16e3,1x))') a(row,:)
        end do
    end subroutine output_matrix
end program sparse_csc_reference_driver
