! SPDX-License-Identifier: MIT
! Machine-readable original SciFort ILUT factorization/solve reference driver.
program sparse_spilu_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, csr_to_csc, csr_to_coo
    use scifort_sparse_linalg, only : ilu0_factor, spilu, spilu_solve, &
        gmres
    implicit none
    type(csr_matrix) :: a
    type(csc_matrix) :: ac
    type(ilu0_factor) :: factor
    integer :: nc, icase, n, nnz, nrhs, storage, k, j, status, solve_status, gmres_info, steps
    integer, allocatable :: rows(:), columns(:), fr(:), fc(:)
    real(dp), allocatable :: values(:), fv(:), b(:,:), x(:,:), solution(:)
    real(dp) :: tol, fill, shift, residual

    read(*,*) nc
    do icase = 1,nc
        read(*,*) n, nnz, nrhs, storage, tol, fill, shift
        allocate(rows(nnz),columns(nnz),values(nnz),b(n,nrhs))
        do k = 1,nnz
            read(*,*) rows(k),columns(k),values(k)
        end do
        do j = 1,nrhs
            do k = 1,n
                read(*,*) b(k,j)
            end do
        end do
        a = csr_from_coo(n,n,rows,columns,values)
        print '(a,1x,i0)', 'CASE',icase
        if (storage == 1) then
            ac = csr_to_csc(a)
            call spilu(ac,factor,status,drop_tol=tol,fill_factor=fill,diagonal_shift=shift)
        else
            call spilu(a,factor,status,drop_tol=tol,fill_factor=fill,diagonal_shift=shift)
        end if
        print '(a,1x,i0)', 'STATUS',status
        if (status == 0) then
            call csr_to_coo(factor%lu,fr,fc,fv)
            print '(a,1x,i0)', 'NNZ',size(fr)
            do k = 1,size(fr)
                write(*,'(a,1x,i0,1x,i0,1x,es26.18)') 'LU',fr(k),fc(k),fv(k)
            end do
            call spilu_solve(factor,b,x,solve_status)
            print '(a,1x,i0)', 'SOLVE',solve_status
            if (solve_status == 0) then
                do j = 1,nrhs
                    write(*,'(a,*(1x,es26.18))') 'X',x(:,j)
                end do
            end if
            call gmres(a,b(:,1),solution,gmres_info,rtol=1.0e-10_dp,atol=0.0_dp, &
                maxiter=100,restart=min(20,n),ilu=factor,iterations=steps,residual_norm=residual)
            print '(a,1x,i0,1x,i0,1x,es26.18)', 'GMRES',gmres_info,steps,residual
        end if
        deallocate(rows,columns,values,b)
    end do
end program sparse_spilu_reference_driver
