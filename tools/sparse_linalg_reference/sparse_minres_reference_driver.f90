! SPDX-License-Identifier: MIT
! Read zero-based? No; the test stream uses one-based COO entries.
program sparse_minres_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : minres
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: values(:), b(:), x0(:), x(:)
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: shift, rtol, residual
    integer :: cases, n, nnz, cap, jacobi, status, case_id, k, steps
    read(*,*) cases
    do case_id = 1, cases
        read(*,*) n, nnz, shift, rtol, cap, jacobi
        allocate(rows(nnz),cols(nnz),values(nnz),b(n),x0(n))
        do k = 1,nnz
            read(*,*) rows(k),cols(k),values(k)
        end do
        do k = 1,n
            read(*,*) b(k)
        end do
        do k = 1,n
            read(*,*) x0(k)
        end do
        a = csr_from_coo(n,n,rows,cols,values)
        call minres(a,b,x,status,rtol=rtol,atol=1.0e-12_dp,maxiter=cap, &
            shift=shift,jacobi=jacobi/=0,x0=x0,iterations=steps,residual_norm=residual)
        print '(a,1x,i0)', 'CASE', case_id
        print '(a,1x,i0,1x,i0,1x,es24.16)', 'MINRES', status, steps, residual
        write(*,'(a)',advance='no') 'X'
        do k = 1,size(x)
            write(*,'(1x,es24.16)',advance='no') x(k)
        end do
        write(*,*)
        deallocate(rows,cols,values,b,x0,x)
    end do
end program sparse_minres_reference_driver
