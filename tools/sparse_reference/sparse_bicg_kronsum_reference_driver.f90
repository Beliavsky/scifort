! SPDX-License-Identifier: MIT
! Line-oriented driver for independent Python/SciPy numerical comparisons.
program sparse_bicg_kronsum_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_to_coo, &
        csr_kronsum, csr_matvec, sparse_success
    use scifort_sparse_linalg, only : bicg, cgs
    implicit none
    character(len=1) :: kind, method
    type(csr_matrix) :: a, b, c
    integer, allocatable :: row(:), col(:), rb(:), cb(:)
    real(dp), allocatable :: va(:), vb(:), rhs(:), x0(:), x(:)
    integer :: ios, n, m, na, nb, stat, stat2, i, iters, jac
    real(dp) :: residual

    do
        read(*,*,iostat=ios) kind
        if (ios /= 0) exit
        if (kind == 'K') then
            read(*,*) m,n,na,nb
            allocate(row(na),col(na),va(na),rb(nb),cb(nb),vb(nb))
            do i = 1, na
                read(*,*) row(i),col(i),va(i)
            end do
            do i = 1, nb
                read(*,*) rb(i),cb(i),vb(i)
            end do
            a = csr_from_coo(m,m,row,col,va,stat)
            b = csr_from_coo(n,n,rb,cb,vb,stat2)
            if (stat == sparse_success .and. stat2 == sparse_success) then
                c = csr_kronsum(a,b,stat)
            else
                stat = 1
            end if
            if (stat /= sparse_success) then
                print '(a)', 'ERR'
            else
                call csr_to_coo(c,row,col,va)
                print '(a,3(1x,i0))', 'OK',c%n_rows,c%n_cols,size(va)
                do i = 1, size(va)
                    print '(2(i0,1x),es26.17e3)', row(i),col(i),va(i)
                end do
            end if
            deallocate(row,col,va,rb,cb,vb)
        else if (kind == 'S') then
            read(*,*) n,na,method,jac
            allocate(row(na),col(na),va(na),rhs(n),x0(n))
            do i = 1, na
                read(*,*) row(i),col(i),va(i)
            end do
            read(*,*) rhs
            read(*,*) x0
            a = csr_from_coo(n,n,row,col,va,stat)
            if (stat == sparse_success) then
                if (method == 'B') then
                    call bicg(a,rhs,x,stat,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=300,x0=x0,jacobi=jac/=0,iterations=iters, &
                        residual_norm=residual)
                else
                    call cgs(a,rhs,x,stat,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                        maxiter=300,x0=x0,jacobi=jac/=0,iterations=iters, &
                        residual_norm=residual)
                end if
                print '(a,2(1x,i0),1x,es26.17e3)', 'SOLVE',stat,iters,residual
                print '(*(es26.17e3,1x))', x
            else
                print '(a)', 'ERR'
            end if
            deallocate(row,col,va,rhs,x0)
            if (allocated(x)) deallocate(x)
        else
            error stop 'unknown input case'
        end if
    end do
end program sparse_bicg_kronsum_reference_driver
