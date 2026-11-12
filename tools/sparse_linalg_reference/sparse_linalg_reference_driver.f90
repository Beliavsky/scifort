! SPDX-License-Identifier: MIT
! Driver for reproducible reference comparisons against scipy.sparse.linalg.
program sparse_linalg_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, sparse_success
    use scifort_sparse_linalg, only : cg, bicgstab, gmres, lsmr
    implicit none
    type(csr_matrix) :: a
    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: values(:), b(:), x(:)
    integer :: tests, case_id, kind, m, n, nnz, i, info, steps, status
    real(dp) :: rn

    read(*,*) tests
    do case_id = 1, tests
        read(*,*) kind, m, n, nnz
        allocate(rows(nnz), cols(nnz), values(nnz), b(m))
        do i = 1, nnz
            read(*,*) rows(i), cols(i), values(i)
        end do
        do i = 1, m
            read(*,*) b(i)
        end do
        a = csr_from_coo(m,n,rows,cols,values,status)
        if (status /= sparse_success) error stop 'invalid CSR input to reference driver'
        write(*,'(a,1x,i0)') 'CASE', case_id
        if (kind == 1) then
            call cg(a,b,x,info,rtol=1.0e-10_dp,maxiter=10*n, &
                iterations=steps,residual_norm=rn)
            call emit('C',info,steps,rn,x)
            call cg(a,b,x,info,rtol=1.0e-10_dp,maxiter=10*n, &
                jacobi=.true.,iterations=steps,residual_norm=rn)
            call emit('CJ',info,steps,rn,x)
        end if
        if (kind /= 3) then
            call bicgstab(a,b,x,info,rtol=1.0e-10_dp,maxiter=10*n, &
                iterations=steps,residual_norm=rn)
            call emit('B',info,steps,rn,x)
            call bicgstab(a,b,x,info,rtol=1.0e-10_dp,maxiter=10*n, &
                jacobi=.true.,iterations=steps,residual_norm=rn)
            call emit('BJ',info,steps,rn,x)
            call gmres(a,b,x,info,rtol=1.0e-10_dp,restart=min(8,n), &
                maxiter=100,iterations=steps,residual_norm=rn)
            call emit('G',info,steps,rn,x)
            call gmres(a,b,x,info,rtol=1.0e-10_dp,restart=min(8,n), &
                maxiter=100,jacobi=.true.,iterations=steps,residual_norm=rn)
            call emit('GJ',info,steps,rn,x)
        end if
        call lsmr(a,b,x,info,itn=steps,normr=rn,atol=1.0e-12_dp, &
            btol=1.0e-12_dp,maxiter=10*max(m,n))
        call emit('L',info,steps,rn,x)
        if (kind == 3) then
            call lsmr(a,b,x,info,itn=steps,normr=rn,damp=0.15_dp, &
                atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=10*max(m,n))
            call emit('LD',info,steps,rn,x)
        end if
        deallocate(rows,cols,values,b)
    end do
contains
    subroutine emit(label, code, it, norm, vec)
        character(len=*), intent(in) :: label !! One/two letter method identifier.
        integer, intent(in) :: code !! Solver status.
        integer, intent(in) :: it !! Number of solver steps.
        real(dp), intent(in) :: norm !! True residual norm.
        real(dp), intent(in) :: vec(:) !! Solution vector.
        integer :: j
        write(*,'(a,1x,i0,1x,i0,1x,es25.16e3)',advance='no') label, code, it, norm
        do j = 1, size(vec)
            write(*,'(1x,es25.16e3)',advance='no') vec(j)
        end do
        print *
    end subroutine emit
end program sparse_linalg_reference_driver
