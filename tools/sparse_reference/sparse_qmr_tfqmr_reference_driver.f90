! SPDX-License-Identifier: MIT
! Stream CSR systems to the two new Krylov solvers for SciPy cross-checks.
program sparse_qmr_tfqmr_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, sparse_success
    use scifort_sparse_linalg, only : qmr, tfqmr, ilu0, ilu0_factor
    implicit none
    character(len=1) :: code
    type(csr_matrix) :: a
    type(ilu0_factor) :: factor
    integer, allocatable :: row(:), col(:)
    real(dp), allocatable :: va(:), rhs(:), x0(:), x(:)
    integer :: ios, n, nnz, pre, status, i, steps
    real(dp) :: residual

    do
        read(*,*,iostat=ios) code
        if (ios /= 0) exit
        read(*,*) n, nnz, pre
        allocate(row(nnz),col(nnz),va(nnz),rhs(n),x0(n))
        do i = 1, nnz
            read(*,*) row(i),col(i),va(i)
        end do
        read(*,*) rhs
        read(*,*) x0
        a = csr_from_coo(n,n,row,col,va,status)
        if (status /= sparse_success) error stop 'CSR assembly'
        if (code == 'Q') then
            call qmr(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                maxiter=350,x0=x0,jacobi=pre==1, &
                iterations=steps,residual_norm=residual)
        else if (code == 'T') then
            if (pre == 2) then
                call ilu0(a,factor,status)
                if (status /= 0) error stop 'ILU0 factorization'
                call tfqmr(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=350,x0=x0,ilu=factor, &
                    iterations=steps,residual_norm=residual)
            else
                call tfqmr(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=350,x0=x0,jacobi=pre==1, &
                    iterations=steps,residual_norm=residual)
            end if
        else
            error stop 'unknown solver code'
        end if
        print '(a,2(1x,i0),1x,es26.17e3)', 'SOLVE',status,steps,residual
        print '(*(es26.17e3,1x))', x
        deallocate(row,col,va,rhs,x0,x)
    end do
end program sparse_qmr_tfqmr_reference_driver
