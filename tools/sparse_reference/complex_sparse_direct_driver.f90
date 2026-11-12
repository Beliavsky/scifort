! SPDX-License-Identifier: MIT
! Streaming reference driver for complex sparse LU (input/output text pairs).
program complex_sparse_direct_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc
    use scifort_sparse_linalg_complex_direct, only : splu_complex_factor, &
        splu_complex, splu_complex_solve, spsolve_complex
    implicit none

    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    type(splu_complex_factor) :: factor
    integer, allocatable :: rows(:), columns(:)
    complex(dp), allocatable :: values(:), b(:,:), x(:,:)
    integer :: n, nnz, nrhs, format_code, api, operation, ierr, stat, i, j, read_status
    real(dp) :: real_part, imaginary_part, pivot_tol
    character(len=1) :: trans

    do
        read(*,*,iostat=read_status) n,nnz,nrhs,format_code,api,operation,pivot_tol
        if (read_status < 0) exit
        if (read_status /= 0) error stop 'invalid case header'
        if (n < 0 .or. nnz < 0 .or. nrhs < 0) error stop 'negative shape'
        allocate(rows(nnz),columns(nnz),values(nnz),b(n,nrhs))
        do i = 1, nnz
            read(*,*) rows(i),columns(i),real_part,imaginary_part
            values(i) = cmplx(real_part,imaginary_part,dp)
        end do
        do j = 1,nrhs
            do i = 1,n
                read(*,*) real_part,imaginary_part
                b(i,j) = cmplx(real_part,imaginary_part,dp)
            end do
        end do
        csr = csr_complex_from_coo(n,n,rows,columns,values,stat)
        if (stat /= 0) error stop 'invalid COO case'
        csc = csr_complex_to_csc(csr,stat)
        if (stat /= 0) error stop 'invalid CSC conversion'
        select case (operation)
        case (1)
            trans = 'N'
        case (2)
            trans = 'T'
        case (3)
            trans = 'H'
        case default
            error stop 'invalid trans code'
        end select
        if (api == 1) then
            if (format_code == 1) then
                call spsolve_complex(csr,b,x,ierr,pivot_tol)
            else
                call spsolve_complex(csc,b,x,ierr,pivot_tol)
            end if
        else
            if (format_code == 1) then
                call splu_complex(csr,factor,ierr,pivot_tol)
            else
                call splu_complex(csc,factor,ierr,pivot_tol)
            end if
            if (ierr == 0) call splu_complex_solve(factor,b,x,ierr,trans)
        end if
        write(*,'(i0)') ierr
        if (ierr == 0) then
            do j = 1,nrhs
                do i = 1,n
                    write(*,'(es25.16,1x,es25.16)') real(x(i,j),dp),aimag(x(i,j))
                end do
            end do
        end if
        deallocate(rows,columns,values,b)
    end do
end program complex_sparse_direct_driver
