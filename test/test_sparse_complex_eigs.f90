! SPDX-License-Identifier: MIT
program test_sparse_complex_eigs
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc, csr_complex_matvec
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex, &
        complex_eigs_success, complex_eigs_invalid_input
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: vals(:),vecs(:,:),av(:)
    real(dp), allocatable :: residuals(:)
    integer :: rows(10),cols(10),i,status,cycles,products,n,nnz
    integer, allocatable :: sparse_rows(:),sparse_cols(:)
    complex(dp), allocatable :: sparse_values(:)
    complex(dp) :: v(10)
    real(dp) :: tolerance

    ! Triangular, genuinely nonsymmetric matrix with complex eigenvalues.
    do i = 1,8
        rows(i) = i
        cols(i) = i
        v(i) = cmplx(real(i,dp),0.21_dp*real(i-4,dp),dp)
    end do
    rows(9:10) = [1,4]
    cols(9:10) = [5,8]
    v(9) = cmplx(0.7_dp,1.4_dp,dp)
    v(10) = cmplx(-0.4_dp,0.3_dp,dp)
    csr = csr_complex_from_coo(8,8,rows,cols,v,status)
    if (status /= 0) error stop 'CSR assembly failed'
    csc = csr_complex_to_csc(csr,status)
    if (status /= 0) error stop 'CSC conversion failed'
    tolerance = 2.0e-9_dp
    call eigs_complex(csr,2,vals,vecs,status,which='LM',ncv=8, &
        iterations=cycles,matvecs=products,residual_norms=residuals)
    if (status /= complex_eigs_success) error stop 'CSR eigs did not converge'
    if (size(vals) /= 2 .or. size(vecs,1) /= 8 .or. size(vecs,2) /= 2) &
        error stop 'incorrect eigenpair output shapes'
    if (abs(vals(1)-v(8)) > tolerance .or. abs(vals(2)-v(7)) > tolerance) &
        error stop 'wrong leading eigenvalues'
    if (cycles /= 1 .or. products < 8) error stop 'bad iteration statistics'
    do i = 1,2
        av = csr_complex_matvec(csr,vecs(:,i),status)
        if (status /= 0) error stop 'matrix-vector multiply failed'
        if (norm2(abs(av-vals(i)*vecs(:,i))) > tolerance) error stop 'CSR true residual failed'
        if (abs(norm2(abs(vecs(:,i)))-1.0_dp) > tolerance) error stop 'vector not normalized'
        if (residuals(i) > tolerance) error stop 'reported CSR residual too large'
    end do
    call eigs_complex(csc,2,vals,vecs,status,which='LR',ncv=8,residual_norms=residuals)
    if (status /= complex_eigs_success) error stop 'CSC eigs did not converge'
    if (abs(vals(1)-v(8)) > tolerance .or. abs(vals(2)-v(7)) > tolerance) &
        error stop 'wrong largest-real eigenvalues'
    call eigs_complex(csr,1,vals,vecs,status,which='SI',ncv=8)
    if (status /= complex_eigs_success) error stop 'SI eigs failed'
    if (abs(vals(1)-v(1)) > tolerance) error stop 'wrong smallest-imag eigenvalue'
    call eigs_complex(csr,1,vals,vecs,status,which='SM',ncv=8)
    if (status /= complex_eigs_success) error stop 'SM eigs failed'
    if (abs(vals(1)-v(1)) > tolerance) error stop 'wrong smallest-magnitude eigenvalue'
    call eigs_complex(csr,1,vals,vecs,status,which='SM',ncv=4)
    if (status /= complex_eigs_invalid_input) error stop 'SM partial subspace not rejected'
    call eigs_complex(csr,1,vals,vecs,status,which='bad')
    if (status /= complex_eigs_invalid_input) error stop 'invalid which not rejected'
    call eigs_complex(csr,0,vals,vecs,status)
    if (status /= complex_eigs_invalid_input) error stop 'invalid k not rejected'
    ! Full-spectrum repeated dominant eigenvalue: do not report a spurious
    ! single duplicate vector as two independent eigenvectors.
    n = 10
    allocate(sparse_rows(n),sparse_cols(n),sparse_values(n))
    do i = 1,n
        sparse_rows(i) = i
        sparse_cols(i) = i
        sparse_values(i) = cmplx(real(i,dp),0.5_dp,dp)
    end do
    sparse_values(9:10) = cmplx(20.0_dp,-2.0_dp,dp)
    csr = csr_complex_from_coo(n,n,sparse_rows,sparse_cols,sparse_values,status)
    call eigs_complex(csr,2,vals,vecs,status,which='LM',ncv=n,residual_norms=residuals)
    if (status /= complex_eigs_success) error stop 'repeated eigenvalues failed'
    if (maxval(abs(vals-sparse_values(10))) > 1.0e-8_dp) error stop 'repeated values incorrect'
    if (abs(dot_product(vecs(:,1),vecs(:,2))) > 0.1_dp) &
        error stop 'repeated eigenvectors not independent'
    deallocate(sparse_rows,sparse_cols,sparse_values)

    ! Explicitly restarted search on a larger, complex nonsymmetric matrix.
    n = 48
    nnz = 3*n-2
    allocate(sparse_rows(nnz),sparse_cols(nnz),sparse_values(nnz))
    do i = 1,n
        sparse_rows(i) = i
        sparse_cols(i) = i
        sparse_values(i) = cmplx(real(n-i+1,dp)/3.0_dp,real(i-n/2,dp)/30.0_dp,dp)
    end do
    do i = 1,n-1
        sparse_rows(n+i) = i
        sparse_cols(n+i) = i+1
        sparse_values(n+i) = cmplx(0.02_dp,0.1_dp,dp)
        sparse_rows(2*n-1+i) = i+1
        sparse_cols(2*n-1+i) = i
        sparse_values(2*n-1+i) = cmplx(0.01_dp,-0.08_dp,dp)
    end do
    csr = csr_complex_from_coo(n,n,sparse_rows,sparse_cols,sparse_values,status)
    if (status /= 0) error stop 'large sparse CSR construction failed'
    call eigs_complex(csr,2,vals,vecs,status,which='LM',ncv=12,maxiter=30, &
        tol=1.0e-8_dp,iterations=cycles,matvecs=products,residual_norms=residuals)
    if (status /= complex_eigs_success) error stop 'restarted Arnoldi did not converge'
    if (cycles <= 1 .or. products <= 12) error stop 'restarted solver did not restart'
    do i = 1,2
        av = csr_complex_matvec(csr,vecs(:,i),status)
        if (norm2(abs(av-vals(i)*vecs(:,i))) > 2.0e-7_dp) &
            error stop 'restarted Arnoldi true residual too large'
    end do
    call eigs_complex(csr,2,vals,vecs,status,which='LM',ncv=6,maxiter=1, &
        tol=1.0e-14_dp,iterations=cycles,matvecs=products,residual_norms=residuals)
    if (status /= 1 .or. cycles /= 1) error stop 'budget exhaustion misreported'
    if (any(residuals < 1.0e-12_dp)) error stop 'unexpected false convergence'
    deallocate(sparse_rows,sparse_cols,sparse_values)
    print '(a)', 'test_sparse_complex_eigs: PASS'
end program test_sparse_complex_eigs
