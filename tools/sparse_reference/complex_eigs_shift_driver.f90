! SPDX-License-Identifier: MIT
! Streaming reference driver for scipy.sparse.linalg.eigs comparisons.
program complex_eigs_shift_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc, csr_complex_matvec
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: values(:), eigenvalues(:), eigenvectors(:,:), av(:)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: rows(:), cols(:)
    integer :: cases, case_id, n, nnz, k, width, cap, format, with_start, &
        i, j, status, cycles, products
    real(dp) :: re, im, tolerance, residual
    character(len=2) :: which
    complex(dp) :: shift
    real(dp) :: shift_real, shift_imag
    complex(dp), allocatable :: start(:)

    read(*,*) cases
    do case_id = 1,cases
        read(*,*) n, nnz, k, width, cap, which, tolerance, format, with_start, shift_real, shift_imag
        shift = cmplx(shift_real,shift_imag,dp)
        allocate(rows(nnz),cols(nnz),values(nnz))
        do j = 1,nnz
            read(*,*) rows(j),cols(j),re,im
            values(j) = cmplx(re,im,dp)
        end do
        if (with_start /= 0) then
            allocate(start(n))
            do j = 1,n
                read(*,*) re,im
                start(j) = cmplx(re,im,dp)
            end do
        end if
        csr = csr_complex_from_coo(n,n,rows,cols,values,status)
        if (status /= 0) error stop 'invalid generated COO matrix'
        csc = csr_complex_to_csc(csr,status)
        if (status /= 0) error stop 'CSC conversion failed'
        if (format == 0) then
            if (with_start == 0) then
                call eigs_complex(csr,k,eigenvalues,eigenvectors,status,which=which, &
                    tol=tolerance,ncv=width,maxiter=cap,iterations=cycles, &
                    matvecs=products,residual_norms=residuals,sigma=shift)
            else
                call eigs_complex(csr,k,eigenvalues,eigenvectors,status,which=which, &
                    tol=tolerance,ncv=width,maxiter=cap,v0=start,iterations=cycles, &
                    matvecs=products,residual_norms=residuals,sigma=shift)
            end if
        else
            if (with_start == 0) then
                call eigs_complex(csc,k,eigenvalues,eigenvectors,status,which=which, &
                    tol=tolerance,ncv=width,maxiter=cap,iterations=cycles, &
                    matvecs=products,residual_norms=residuals,sigma=shift)
            else
                call eigs_complex(csc,k,eigenvalues,eigenvectors,status,which=which, &
                    tol=tolerance,ncv=width,maxiter=cap,v0=start,iterations=cycles, &
                    matvecs=products,residual_norms=residuals,sigma=shift)
            end if
        end if
        write(*,'(a,1x,i0)') 'CASE',case_id
        write(*,'(4(i0,1x))') status,cycles,products,size(eigenvalues)
        do j = 1,size(eigenvalues)
            av = csr_complex_matvec(csr,eigenvectors(:,j),status)
            residual = norm2(abs(av-eigenvalues(j)*eigenvectors(:,j)))
            write(*,'(4(es26.17e3,1x))') real(eigenvalues(j),dp), aimag(eigenvalues(j)), &
                residuals(j), residual
        end do
        deallocate(rows,cols,values,eigenvalues,eigenvectors,residuals)
        if (allocated(start)) deallocate(start)
    end do
end program complex_eigs_shift_driver
