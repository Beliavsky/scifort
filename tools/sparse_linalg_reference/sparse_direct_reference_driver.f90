! SPDX-License-Identifier: MIT
! Batched SciFort sparse direct and triangular solver numerical reference driver.
program sparse_direct_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, csr_to_csc, csr_to_coo
    use scifort_sparse_linalg, only : splu, splu_factor, spsolve, splu_solve, &
        spsolve_triangular
    implicit none
    type(csr_matrix) :: a
    type(csc_matrix) :: ac
    type(splu_factor) :: fac
    integer :: nc, icase, n, nnz, nrhs, kind, storage, k, j, status, status2
    integer, allocatable :: rows(:), columns(:), fr(:), fc(:)
    real(dp), allocatable :: values(:), fv(:), b(:,:), x(:,:), y(:,:)

    read(*,*) nc
    do icase=1,nc
        read(*,*) n, nnz, nrhs, kind, storage
        allocate(rows(nnz), columns(nnz), values(nnz), b(n,nrhs))
        do k=1,nnz
            read(*,*) rows(k), columns(k), values(k)
        end do
        do j=1,nrhs
            do k=1,n
                read(*,*) b(k,j)
            end do
        end do
        a = csr_from_coo(n,n,rows,columns,values)
        ac = csr_to_csc(a)
        print '(a,1x,i0)', 'CASE', icase
        if (kind <= 1) then
            if (storage == 1) then
                call splu(ac,fac,status,pivot_tol=1.0e-13_dp)
            else
                call splu(a,fac,status,pivot_tol=1.0e-13_dp)
            end if
            print '(a,1x,i0)', 'LU', status
            if (status == 0) then
                print '(a,*(1x,i0))', 'PERM',fac%row_perm
                call csr_to_coo(fac%lu,fr,fc,fv)
                print '(a,1x,i0)', 'NNZ',size(fr)
                do k=1,size(fr)
                    print '(2(i0,1x),es26.18)', fr(k),fc(k),fv(k)
                end do
                call splu_solve(fac,b,x,status2)
                print '(a,1x,i0)', 'FSOLVE',status2
                if (status2 == 0) then
                    do j=1,nrhs
                        write(*,'(a,*(1x,es26.18))') 'X',x(:,j)
                    end do
                end if
            end if
            if (storage == 1) then
                call spsolve(ac,b,y,status2,pivot_tol=1.0e-13_dp)
            else
                call spsolve(a,b,y,status2,pivot_tol=1.0e-13_dp)
            end if
            print '(a,1x,i0)', 'SOLVE',status2
            if (status2 == 0) then
                do j=1,nrhs
                    write(*,'(a,*(1x,es26.18))') 'Y',y(:,j)
                end do
            end if
        else
            if (storage == 1) then
                call spsolve_triangular(ac,b,x,status,lower=(kind == 2 .or. kind == 4), &
                    unit_diagonal=(kind >= 4))
            else
                call spsolve_triangular(a,b,x,status,lower=(kind == 2 .or. kind == 4), &
                    unit_diagonal=(kind >= 4))
            end if
            print '(a,1x,i0)', 'TRI',status
            if (status == 0) then
                do j=1,nrhs
                    write(*,'(a,*(1x,es26.18))') 'X',x(:,j)
                end do
            end if
        end if
        deallocate(rows,columns,values,b)
    end do
end program sparse_direct_reference_driver
