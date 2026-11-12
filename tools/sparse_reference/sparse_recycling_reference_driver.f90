! SPDX-License-Identifier: MIT
! Stream sparse nonsymmetric systems and run two RHS solves with reused directions.
program sparse_recycling_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, sparse_success
    use scifort_sparse_linalg, only : lgmres, gcrotmk, ilu0, ilu0_factor
    implicit none
    character(len=1) :: code
    type(csr_matrix) :: a
    type(ilu0_factor) :: factor
    integer, allocatable :: row(:), col(:)
    real(dp), allocatable :: va(:), rhs(:), rhs2(:), x0(:), x02(:), x(:)
    real(dp), allocatable :: saved(:,:)
    integer :: ios, n, nnz, pre, status, i, steps, m, k, pass
    real(dp) :: residual

    do
        read(*,*,iostat=ios) code
        if (ios /= 0) exit
        read(*,*) n, nnz, pre, m, k
        allocate(row(nnz),col(nnz),va(nnz),rhs(n),rhs2(n),x0(n),x02(n))
        do i = 1, nnz
            read(*,*) row(i),col(i),va(i)
        end do
        read(*,*) rhs
        read(*,*) rhs2
        read(*,*) x0
        read(*,*) x02
        a = csr_from_coo(n,n,row,col,va,status)
        if (status /= sparse_success) error stop 'CSR assembly'
        if (pre == 2) then
            call ilu0(a,factor,status)
            if (status /= 0) error stop 'ILU0 factorization'
        end if
        do pass = 1, 2
            if (pass == 2) then
                rhs = rhs2
                x0 = x02
            end if
            if (code == 'L') then
                if (pre == 2) then
                    call lgmres(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=180,inner_m=m,outer_k=k,x0=x0,ilu=factor, &
                        iterations=steps,residual_norm=residual,outer_v=saved)
                else
                    call lgmres(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=180,inner_m=m,outer_k=k,x0=x0,jacobi=pre==1, &
                        iterations=steps,residual_norm=residual,outer_v=saved)
                end if
            else if (code == 'G') then
                if (pre == 2) then
                    call gcrotmk(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=180,m=m,k=k,x0=x0,ilu=factor, &
                        iterations=steps,residual_norm=residual,recycle_u=saved)
                else
                    call gcrotmk(a,rhs,x,status,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=180,m=m,k=k,x0=x0,jacobi=pre==1, &
                        iterations=steps,residual_norm=residual,recycle_u=saved)
                end if
            else
                error stop 'unknown solver code'
            end if
            print '(a,2(1x,i0),1x,es26.17e3,1x,i0)', 'SOLVE',status,steps,residual,size(saved,2)
            print '(*(es26.17e3,1x))', x
        end do
        deallocate(row,col,va,rhs,rhs2,x0,x02,x,saved)
    end do
end program sparse_recycling_reference_driver
