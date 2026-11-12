! SPDX-License-Identifier: MIT
program sparse_eigs_shift_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : eigs
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: entries(:), residuals(:)
    complex(dp), allocatable :: values(:), vectors(:,:)
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: tol, sr, si
    integer :: cases, n, nnz, cap, nc, k, case_id, j, i, status, iters
    character(len=2) :: which
    read(*,*) cases
    do case_id=1,cases
        read(*,*) n, nnz, k, sr, si, tol, cap, nc, which
        allocate(rows(nnz),cols(nnz),entries(nnz))
        do i=1,nnz
            read(*,*) rows(i),cols(i),entries(i)
        end do
        a=csr_from_coo(n,n,rows,cols,entries)
        if (which == 'SM') then
            call eigs(a,k,values,vectors,status,which='SM',tol=tol,maxiter=cap, &
                ncv=nc,iterations=iters,residual_norms=residuals, &
                inner_rtol=1.0e-12_dp,inner_restart=n,inner_maxiter=max(100,4*n))
        else
            call eigs(a,k,values,vectors,status,which='LM',tol=tol,maxiter=cap, &
                ncv=nc,iterations=iters,residual_norms=residuals, &
                sigma=cmplx(sr,si,dp),inner_rtol=1.0e-12_dp, &
                inner_restart=n,inner_maxiter=max(100,4*n))
        end if
        print '(a,1x,i0)','CASE',case_id
        print '(a,1x,i0,1x,i0)','EIGS',status,iters
        write(*,'(a)',advance='no') 'LAMBDA'
        do j=1,size(values)
            write(*,'(2(1x,es24.16))',advance='no') real(values(j),dp),aimag(values(j))
        end do
        write(*,*)
        write(*,'(a)',advance='no') 'RESIDUAL'
        do j=1,size(residuals)
            write(*,'(1x,es24.16)',advance='no') residuals(j)
        end do
        write(*,*)
        do j=1,size(vectors,2)
            write(*,'(a,1x,i0)',advance='no') 'VECTOR',j
            do i=1,size(vectors,1)
                write(*,'(2(1x,es24.16))',advance='no') real(vectors(i,j),dp),aimag(vectors(i,j))
            end do
            write(*,*)
        end do
        deallocate(rows,cols,entries,values,vectors,residuals)
    end do
end program sparse_eigs_shift_reference_driver
