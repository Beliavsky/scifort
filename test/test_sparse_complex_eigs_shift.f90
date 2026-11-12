! SPDX-License-Identifier: MIT
program test_sparse_complex_eigs_shift
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc, csr_complex_matvec
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex, &
        complex_eigs_success, complex_eigs_no_convergence, &
        complex_eigs_invalid_input, complex_eigs_singular_shift
    implicit none
    type(csr_complex_matrix) :: a
    type(csc_complex_matrix) :: ac
    complex(dp), allocatable :: lambda(:), x(:,:)
    integer, allocatable :: large_rows(:), large_cols(:)
    complex(dp), allocatable :: large_data(:)
    real(dp), allocatable :: errors(:)
    integer :: rows(13), cols(13), i, info, iterations, products, n
    complex(dp) :: data(13), shift
    real(dp), parameter :: tolerance = 2.0e-8_dp

    do i = 1,10
        rows(i) = i
        cols(i) = i
        data(i) = cmplx(real(i,dp),0.18_dp*real(i-5,dp),dp)
    end do
    rows(11:13) = [1,4,6]
    cols(11:13) = [3,7,10]
    data(11:13) = [cmplx(0.3_dp,-0.1_dp,dp), cmplx(0.2_dp,0.1_dp,dp), &
        cmplx(0.4_dp,0.3_dp,dp)]
    a = csr_complex_from_coo(10,10,rows,cols,data,info)
    if (info /= 0) error stop 'shift test: assembly failed'
    ac = csr_complex_to_csc(a,info)
    if (info /= 0) error stop 'shift test: CSC conversion failed'
    shift = cmplx(4.25_dp,0.09_dp,dp)
    call eigs_complex(a,2,lambda,x,info,which='LM',sigma=shift,ncv=10, &
        tol=1.0e-10_dp,iterations=iterations,matvecs=products,residual_norms=errors)
    if (info /= complex_eigs_success) error stop 'shift test: full CSR solve failed'
    if (abs(lambda(1)-data(4)) > tolerance .or. abs(lambda(2)-data(5)) > tolerance) &
        error stop 'shift test: wrong eigenvalues closest to sigma'
    if (iterations /= 1 .or. products < 10) error stop 'shift test: bad stats'
    call verify_eigenpairs(a,lambda,x,errors,tolerance)
    call eigs_complex(ac,2,lambda,x,info,which='LM',sigma=shift,ncv=10, &
        tol=1.0e-10_dp,residual_norms=errors)
    if (info /= complex_eigs_success) error stop 'shift test: CSC failed'
    call verify_eigenpairs(a,lambda,x,errors,tolerance)
    if (abs(lambda(1)-data(4)) > tolerance .or. abs(lambda(2)-data(5)) > tolerance) &
        error stop 'shift test: wrong CSC eigenvalues'

    ! A shift at an exact diagonal eigenvalue must not claim convergence.
    call eigs_complex(a,1,lambda,x,info,sigma=data(4),ncv=10)
    if (info /= complex_eigs_singular_shift) error stop 'shift test: singular shift undetected'
    if (size(lambda) /= 0 .or. size(x) /= 0) error stop 'shift test: nonempty singular output'
    call eigs_complex(a,1,lambda,x,info,sigma=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp))
    if (info /= complex_eigs_invalid_input) error stop 'shift test: NaN shift accepted'
    call eigs_complex(a,1,lambda,x,info,sigma=shift,ncv=1)
    if (info /= complex_eigs_invalid_input) error stop 'shift test: invalid ncv accepted'

    ! Nonconvergence must be reported, rather than incorrectly accepted.
    call eigs_complex(a,2,lambda,x,info,sigma=shift,ncv=3,maxiter=1, &
        tol=1.0e-14_dp,residual_norms=errors)
    if (info /= complex_eigs_no_convergence) error stop 'shift test: false convergence'
    if (size(lambda) /= 2 .or. size(errors) /= 2) error stop 'shift test: missing Ritz data'
    ! Interior eigenvalues of a larger nonsymmetric banded matrix, found with
    ! a genuinely truncated Krylov space (ncv << n), no dense factorization.
    n = 70
    allocate(large_rows(3*n-2),large_cols(3*n-2),large_data(3*n-2))
    do i = 1,n
        large_rows(i) = i
        large_cols(i) = i
        large_data(i) = cmplx(real(i,dp)/2.0_dp,0.03_dp*real(i,dp),dp)
    end do
    do i = 1,n-1
        large_rows(n+i) = i
        large_cols(n+i) = i+1
        large_data(n+i) = cmplx(0.02_dp,0.008_dp,dp)
        large_rows(2*n-1+i) = i+1
        large_cols(2*n-1+i) = i
        large_data(2*n-1+i) = cmplx(0.016_dp,-0.01_dp,dp)
    end do
    a = csr_complex_from_coo(n,n,large_rows,large_cols,large_data,info)
    if (info /= 0) error stop 'shift test: large CSR assembly failed'
    call eigs_complex(a,2,lambda,x,info,which='LM', &
        sigma=cmplx(17.3_dp,1.04_dp,dp),ncv=16,maxiter=45, &
        tol=1.0e-8_dp,iterations=iterations,matvecs=products,residual_norms=errors)
    if (info /= complex_eigs_success) error stop 'shift test: interior Arnoldi did not converge'
    if (maxval(abs(lambda-[large_data(35),large_data(34)])) > 1.0e-6_dp) &
        error stop 'shift test: incorrect interior eigenvalues'
    if (products < 16) error stop 'shift test: expected inverse applications'
    call verify_eigenpairs(a,lambda,x,errors,5.0e-7_dp)

    print '(a)', 'test_sparse_complex_eigs_shift: PASS'
contains
    subroutine verify_eigenpairs(matrix,values,vectors,residuals,threshold)
        type(csr_complex_matrix), intent(in) :: matrix !! Matrix in CSR format.
        complex(dp), intent(in) :: values(:) !! Returned eigenvalues.
        complex(dp), intent(in) :: vectors(:,:) !! Returned normalized right eigenvectors.
        real(dp), intent(in) :: residuals(:) !! Solver-reported true residual norms.
        real(dp), intent(in) :: threshold !! Maximum acceptable true residual.
        complex(dp), allocatable :: vector_result(:)
        integer :: j, status
        if (size(vectors,2) /= size(values) .or. size(residuals) /= size(values)) &
            error stop 'shift test: incorrect shapes'
        do j = 1,size(values)
            vector_result = csr_complex_matvec(matrix,vectors(:,j),status)
            if (status /= 0) error stop 'shift test: matvec failed'
            if (norm2(abs(vector_result-values(j)*vectors(:,j))) > threshold) &
                error stop 'shift test: inaccurate eigenvector residual'
            if (abs(norm2(abs(vectors(:,j)))-1.0_dp) > threshold) &
                error stop 'shift test: eigenvector not normalized'
            if (residuals(j) > threshold) error stop 'shift test: residual diagnostics incorrect'
        end do
    end subroutine verify_eigenpairs
end program test_sparse_complex_eigs_shift
