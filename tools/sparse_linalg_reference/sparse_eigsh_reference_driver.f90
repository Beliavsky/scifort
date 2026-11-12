! SPDX-License-Identifier: MIT
program sparse_eigsh_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : eigsh
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: entries(:), values(:), vectors(:,:), residuals(:)
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: tol
    integer :: cases, n, nnz, cap, nc, k, case_id, j, i, status, iters
    character(len=2) :: which
    read(*,*) cases
    do case_id = 1,cases
        read(*,*) n, nnz, k, which, tol, cap, nc
        allocate(rows(nnz),cols(nnz),entries(nnz))
        do i = 1,nnz
            read(*,*) rows(i),cols(i),entries(i)
        end do
        a = csr_from_coo(n,n,rows,cols,entries)
        call eigsh(a,k,values,vectors,status,which=which,tol=tol,maxiter=cap, &
            ncv=nc,iterations=iters,residual_norms=residuals)
        print '(a,1x,i0)', 'CASE', case_id
        print '(a,1x,i0,1x,i0)', 'EIGSH', status, iters
        write(*,'(a)',advance='no') 'LAMBDA'
        do j = 1,size(values)
            write(*,'(1x,es24.16)',advance='no') values(j)
        end do
        write(*,*)
        write(*,'(a)',advance='no') 'RESIDUAL'
        do j = 1,size(residuals)
            write(*,'(1x,es24.16)',advance='no') residuals(j)
        end do
        write(*,*)
        do j = 1,size(vectors,2)
            write(*,'(a,1x,i0)',advance='no') 'VECTOR', j
            do i = 1,size(vectors,1)
                write(*,'(1x,es24.16)',advance='no') vectors(i,j)
            end do
            write(*,*)
        end do
        deallocate(rows,cols,entries,values,vectors,residuals)
    end do
end program sparse_eigsh_reference_driver
