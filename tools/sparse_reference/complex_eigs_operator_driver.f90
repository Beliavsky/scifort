! SPDX-License-Identifier: MIT
! Streaming reference driver for matrix-free complex Arnoldi and cached shift solves.
program complex_eigs_operator_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_coo
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex_operator, eigs_complex, &
        eigs_complex_prepare_shift, complex_eigs_shift_cache
    implicit none
    type(csr_complex_matrix) :: csr
    type(complex_eigs_shift_cache) :: cache
    complex(dp), allocatable :: diagonal(:), upper(:), data(:), values(:), vectors(:,:), &
        cached_values(:), cached_vectors(:,:)
    real(dp), allocatable :: residuals(:), cached_residuals(:)
    integer, allocatable :: rows(:), cols(:)
    integer :: cases, case_id, n, k, shifted, status, cached_status, i, products, ierr
    real(dp) :: re, im
    complex(dp) :: sigma
    character(len=2) :: which

    read(*,*) cases
    do case_id = 1,cases
        read(*,*) n,k,shifted,which
        allocate(diagonal(n),upper(n-1),data(2*n-1),rows(2*n-1),cols(2*n-1))
        do i = 1,n
            read(*,*) re, im
            diagonal(i) = cmplx(re,im,dp)
            rows(i) = i
            cols(i) = i
        end do
        do i = 1,n-1
            read(*,*) re, im
            upper(i) = cmplx(re,im,dp)
            rows(n+i) = i
            cols(n+i) = i+1
        end do
        read(*,*) re, im
        sigma = cmplx(re,im,dp)
        data = [diagonal,upper]
        csr = csr_complex_from_coo(n,n,rows,cols,data,ierr)
        if (ierr /= 0) error stop 'operator driver: invalid random CSR'
        if (shifted == 0) then
            call eigs_complex_operator(n,k,apply_original,values,vectors,status, &
                which=which,tol=1.0e-9_dp,ncv=n,matvecs=products, &
                residual_norms=residuals,norm_estimate=real(n,dp))
            call eigs_complex(csr,k,cached_values,cached_vectors,cached_status, &
                which=which,tol=1.0e-9_dp,ncv=n,residual_norms=cached_residuals)
        else
            call eigs_complex_operator(n,k,apply_original,values,vectors,status, &
                which=which,tol=1.0e-9_dp,ncv=n,matvecs=products, &
                residual_norms=residuals,norm_estimate=real(n,dp), &
                sigma=sigma,inverse_matvec=apply_inverse)
            call eigs_complex_prepare_shift(csr,sigma,cache,ierr)
            if (ierr /= 0) error stop 'operator driver: shift factorization failed'
            call eigs_complex(csr,k,cached_values,cached_vectors,cached_status, &
                which=which,tol=1.0e-9_dp,ncv=n,residual_norms=cached_residuals, &
                sigma=sigma,shift_cache=cache)
        end if
        write(*,'(a,1x,i0)') 'CASE',case_id
        write(*,'(4(i0,1x))') status,cached_status,products,size(values)
        do i = 1,size(values)
            write(*,'(6(es25.16e3,1x))') real(values(i),dp),aimag(values(i)),residuals(i), &
                real(cached_values(i),dp),aimag(cached_values(i)),cached_residuals(i)
        end do
        deallocate(diagonal,upper,data,rows,cols,values,vectors,cached_values, &
            cached_vectors,residuals,cached_residuals)
    end do
contains
    subroutine apply_original(x,y,info)
        complex(dp), intent(in) :: x(:) !! Complex vector.
        complex(dp), intent(out) :: y(:) !! Upper bidiagonal operator times x.
        integer, intent(out) :: info !! Zero on success.
        y = diagonal*x
        y(:n-1) = y(:n-1)+upper*x(2:n)
        info = 0
    end subroutine apply_original

    subroutine apply_inverse(x,y,info)
        complex(dp), intent(in) :: x(:) !! Complex right-hand side.
        complex(dp), intent(out) :: y(:) !! Inverse shifted upper bidiagonal action.
        integer, intent(out) :: info !! Zero on success.
        integer :: j
        y(n) = x(n)/(diagonal(n)-sigma)
        do j = n-1,1,-1
            y(j) = (x(j)-upper(j)*y(j+1))/(diagonal(j)-sigma)
        end do
        info = 0
    end subroutine apply_inverse
end program complex_eigs_operator_driver
