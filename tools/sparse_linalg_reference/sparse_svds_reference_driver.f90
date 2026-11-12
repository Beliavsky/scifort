! SPDX-License-Identifier: MIT
! Machine-readable SciPy reference driver for sparse singular value computations.
program sparse_svds_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : svds
    implicit none
    integer :: count_cases, ci, m, n, k, width, cap, nnz, j, info, iterations, status
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: tolerance
    real(dp), allocatable :: values(:), u(:,:), singular(:), vt(:,:), residuals(:)
    type(csr_matrix) :: a
    character(len=2) :: mode

    read(*,*) count_cases
    do ci = 1,count_cases
        read(*,*) m,n,k,width,cap,mode,tolerance,nnz
        allocate(rows(nnz),cols(nnz),values(nnz))
        do j = 1,nnz
            read(*,*) rows(j),cols(j),values(j)
        end do
        a = csr_from_coo(m,n,rows,cols,values,status)
        if (status /= 0) error stop 'invalid COO reference input'
        call svds(a,k,u,singular,vt,info,which=mode,ncv=width, &
            maxiter=cap,tol=tolerance,iterations=iterations,residual_norms=residuals)
        write(*,'(a,1x,i0)') 'CASE',ci
        write(*,'(2(i0,1x))') info,iterations
        if (size(singular) /= k) error stop 'missing singular values'
        do j = 1,k
            write(*,'(*(es26.17,1x))') singular(j),residuals(j),u(:,j),vt(j,:)
        end do
        deallocate(rows,cols,values)
    end do
end program sparse_svds_reference_driver
