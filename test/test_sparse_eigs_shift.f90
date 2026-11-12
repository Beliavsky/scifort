! SPDX-License-Identifier: MIT
program test_sparse_eigs_shift
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigs, sparse_linalg_invalid_input, sparse_linalg_breakdown
    implicit none
    integer, parameter :: n = 24
    type(csr_matrix) :: a
    real(dp) :: dense(n,n), err, dist
    complex(dp), allocatable :: values(:), vectors(:,:)
    real(dp), allocatable :: residuals(:)
    integer :: i, j, info, steps

    dense = 0.0_dp
    do i = 1,n
        dense(i,i) = -3.0_dp+0.25_dp*real(i,dp)
    end do
    do i = 1,n-1
        dense(i,i+1) = 0.008_dp
    end do
    a = csr_from_dense(dense)
    call eigs(a,2,values,vectors,info,sigma=cmplx(0.3_dp,0.0_dp,dp), &
        ncv=12,maxiter=120,tol=5.0e-8_dp,inner_rtol=1.0e-12_dp, &
        inner_restart=n,inner_maxiter=90,iterations=steps,residual_norms=residuals)
    print *, 'real shift status, steps:',info,steps
    if (info /= 0) error stop 'real shift convergence'
    if (maxval(residuals) > 1.0e-5_dp) error stop 'real shift residual'
    do j=1,2
        dist=abs(values(j)-cmplx(0.3_dp,0.0_dp,dp))
        if (dist > 0.23_dp) error stop 'real shift selection'
    end do

    dense = 0.0_dp
    dense(1,1)=1.0_dp
    dense(2,2)=1.0_dp
    dense(1,2)=-3.0_dp
    dense(2,1)=3.0_dp
    do i=3,n
        dense(i,i)=-5.0_dp+0.4_dp*real(i,dp)
        if (i<n) dense(i,i+1)=0.004_dp
    end do
    a=csr_from_dense(dense)
    call eigs(a,1,values,vectors,info,sigma=cmplx(1.1_dp,2.9_dp,dp), &
        ncv=10,maxiter=100,tol=1.0e-8_dp,inner_rtol=1.0e-12_dp, &
        inner_restart=n,inner_maxiter=90,iterations=steps,residual_norms=residuals)
    print *, 'complex shift status, steps, lambda:',info,steps,values(1)
    if (info /= 0) error stop 'complex shift convergence'
    if (abs(values(1)-cmplx(1.0_dp,3.0_dp,dp)) > 1.0e-5_dp) error stop 'complex shift eigenvalue'
    err=sqrt(sum(abs(matmul(dense,vectors(:,1))-values(1)*vectors(:,1))**2))
    if (err > 1.0e-6_dp) error stop 'complex shift residual'

    ! Reduced subspace SM now uses an implicit zero shift, not a full
    ! projected eigensystem.  The smallest modes are isolated and nonsingular.
    dense = 0.0_dp
    do i=1,n
        dense(i,i) = 2.0_dp+0.3_dp*real(i,dp)
    end do
    dense(1,1) = 0.18_dp
    dense(2,2) = -0.34_dp
    do i=1,n-1
        dense(i,i+1) = 0.004_dp
        dense(i+1,i) = -0.003_dp
    end do
    a=csr_from_dense(dense)
    call eigs(a,2,values,vectors,info,which='SM',ncv=10,maxiter=90, &
        tol=5.0e-8_dp,inner_restart=n,inner_maxiter=100, &
        inner_rtol=1.0e-12_dp,residual_norms=residuals)
    if (info/=0) error stop 'SM shift-invert convergence'
    if (minval(abs(values-cmplx(0.18_dp,0.0_dp,dp)))>0.01_dp .or. &
        minval(abs(values+cmplx(0.34_dp,0.0_dp,dp)))>0.01_dp) error stop 'SM selection'
    if (maxval(residuals)>1.0e-5_dp) error stop 'SM true residuals'

    ! Exact singular shifts are not spuriously reported as converged.
    dense = 0.0_dp
    do i=1,n
        dense(i,i) = real(i,dp)
    end do
    a=csr_from_dense(dense)
    call eigs(a,1,values,vectors,info,sigma=cmplx(1.0_dp,0.0_dp,dp), &
        ncv=8,maxiter=10,inner_restart=n,inner_maxiter=75)
    if (info /= sparse_linalg_breakdown) error stop 'singular exact shift'

    call eigs(a,1,values,vectors,info,sigma=cmplx(1.0_dp,0.0_dp,dp),which='LR')
    if (info /= sparse_linalg_invalid_input) error stop 'shift which validation'
    call eigs(a,1,values,vectors,info,sigma=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp))
    if (info /= sparse_linalg_invalid_input) error stop 'nonfinite sigma validation'
    call eigs(a,1,values,vectors,info,sigma=cmplx(0.0_dp,0.0_dp,dp),inner_rtol=0.0_dp)
    if (info /= sparse_linalg_invalid_input) error stop 'inner tolerance validation'
    call eigs(a,1,values,vectors,info,sigma=cmplx(0.0_dp,0.0_dp,dp),inner_restart=0)
    if (info /= sparse_linalg_invalid_input) error stop 'inner restart validation'
    call eigs(a,1,values,vectors,info,sigma=cmplx(0.0_dp,0.0_dp,dp),inner_maxiter=0)
    if (info /= sparse_linalg_invalid_input) error stop 'inner iterations validation'

    print '(a)','test_sparse_eigs_shift: PASS'
end program test_sparse_eigs_shift
