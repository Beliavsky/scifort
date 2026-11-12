! SPDX-License-Identifier: MIT
! Batch numerical reference driver for SciPy comparisons.
program sparse_lsqr_ilu0_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : lsqr, ilu0_factor, ilu0, ilu0_solve, &
        gmres, bicgstab
    implicit none
    type(csr_matrix) :: a
    type(ilu0_factor) :: fac
    integer :: nc, icase, m, n, nnz, cap, k, info, iterations
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: damping, r1, r2, rn, gradnorm
    real(dp), allocatable :: values(:), x0(:), b(:), x(:), variance(:)

    read(*,*) nc
    do icase = 1, nc
        read(*,*) m, n, nnz, damping, cap
        allocate(rows(nnz),cols(nnz),values(nnz),x0(n),b(m))
        do k = 1, nnz
            read(*,*) rows(k), cols(k), values(k)
        end do
        do k = 1, m
            read(*,*) b(k)
        end do
        do k = 1, n
            read(*,*) x0(k)
        end do
        a = csr_from_coo(m,n,rows,cols,values)
        call lsqr(a,b,x,info,itn=iterations,r1norm=r1,r2norm=r2,arnorm=gradnorm, &
            damp=damping,atol=1.0e-12_dp,btol=1.0e-12_dp, &
            conlim=1.0e14_dp,iter_lim=cap,x0=x0,var=variance,calc_var=.true.)
        print '(a,1x,i0)', 'CASE', icase
        print '(a,2(1x,i0),3(1x,es24.16))', 'LSQR', info, iterations, r1, r2, gradnorm
        write(*,'(a,*(1x,es24.16))') 'X',x
        write(*,'(a,*(1x,es24.16))') 'VAR',variance
        if (m == n) then
            call ilu0(a,fac,info)
            print '(a,1x,i0)', 'ILU',info
            if (info == 0) then
                write(*,'(a,*(1x,es24.16))') 'DATA',fac%lu%data
                call ilu0_solve(fac,b,x,info)
                print '(a,1x,i0)', 'ISOLVE',info
                write(*,'(a,*(1x,es24.16))') 'IX',x
                call bicgstab(a,b,x,info,ilu=fac,rtol=1.0e-11_dp,maxiter=10*n, &
                    iterations=iterations,residual_norm=rn)
                print '(a,2(1x,i0),1x,es24.16)', 'B',info,iterations,rn
                write(*,'(a,*(1x,es24.16))') 'BX',x
                call gmres(a,b,x,info,ilu=fac,rtol=1.0e-11_dp,restart=min(n,12), &
                    maxiter=80,iterations=iterations,residual_norm=rn)
                print '(a,2(1x,i0),1x,es24.16)', 'G',info,iterations,rn
                write(*,'(a,*(1x,es24.16))') 'GX',x
            end if
        end if
        deallocate(rows,cols,values,x0,b)
    end do
end program sparse_lsqr_ilu0_reference_driver
