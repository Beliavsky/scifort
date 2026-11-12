! SPDX-License-Identifier: MIT
! Native regression: restarted Lanczos and Arnoldi shift-invert modes.
program test_sparse_restarted_shift
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigsh, eigsh_lanczos, eigs, &
        eigs_arnoldi, sparse_linalg_invalid_input
    implicit none
    integer, parameter :: n = 12
    real(dp) :: dense(n,n), diagonal(n), nan_value
    real(dp), allocatable :: rvalues(:), rvectors(:,:), residuals(:)
    complex(dp), allocatable :: cvalues(:), cvectors(:,:)
    type(csr_matrix) :: a
    integer :: i, info, cycles, locked, products

    dense = 0.0_dp
    do i = 1,n
        diagonal(i) = real(i,dp)-6.2_dp
        dense(i,i) = diagonal(i)
    end do
    a = csr_from_dense(dense)
    call eigsh_lanczos(a,2,rvalues,rvectors,info,sigma=0.0_dp,which='LM', &
        ncv=7,maxiter=70,tol=1.0e-8_dp,iterations=cycles,locked_count=locked, &
        matvecs=products,residual_norms=residuals)
    call require(info == 0 .and. cycles >= 1 .and. locked == 2, 'shifted symmetric convergence')
    call require(products >= 7 .and. maxval(residuals) < 1.0e-7_dp, 'shifted Lanczos diagnostics')
    call require(maxval(abs(rvalues-[-0.2_dp,0.8_dp])) < 1.0e-6_dp, 'interior symmetric eigenvalues')
    do i=1,2
        call require(abs(norm2(matmul(dense,rvectors(:,i))-rvalues(i)*rvectors(:,i)) &
            -residuals(i)) < 1.0e-10_dp, 'original matrix symmetric residual')
    end do
    call eigsh(a,2,rvalues,rvectors,info,method='lanczos',sigma=0.0_dp,which='LM', &
        ncv=7,maxiter=70,tol=1.0e-8_dp,residual_norms=residuals)
    call require(info == 0 .and. maxval(residuals) < 1.0e-7_dp, 'shifted Lanczos dispatch')

    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,which='LA', &
        ncv=7,maxiter=50,tol=1.0e-8_dp,residual_norms=residuals)
    call require(info == 0 .and. abs(rvalues(1)-0.8_dp)<1.0e-6_dp, 'inverse LA ordering')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,which='SA', &
        ncv=7,maxiter=50,tol=1.0e-8_dp,residual_norms=residuals)
    call require(info == 0 .and. abs(rvalues(1)+0.2_dp)<1.0e-6_dp, 'inverse SA ordering')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,which='SM', &
        ncv=n,maxiter=6,tol=1.0e-8_dp,residual_norms=residuals)
    call require(info == 0 .and. abs(rvalues(1)-5.8_dp)<1.0e-6_dp, 'inverse SM ordering')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,which='SM',ncv=7)
    call require(info == sparse_linalg_invalid_input, 'shifted SM reduced width invalid')
    nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=nan_value)
    call require(info == sparse_linalg_invalid_input, 'nonfinite real shift')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,inner_rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'negative MINRES tolerance')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=0.0_dp,inner_maxiter=0)
    call require(info == sparse_linalg_invalid_input, 'invalid MINRES cap')
    call eigsh_lanczos(a,1,rvalues,rvectors,info,sigma=diagonal(6),ncv=7, &
        inner_maxiter=30,maxiter=5)
    call require(info /= 0, 'singular real shift not claimed converged')

    ! A real nonsymmetric block with an isolated complex-conjugate pair.
    dense = 0.0_dp
    do i=1,n
        dense(i,i) = real(i,dp)-6.3_dp
    end do
    dense(1,1) = 0.5_dp
    dense(2,2) = 0.5_dp
    dense(1,2) = -1.0_dp
    dense(2,1) = 1.0_dp
    a=csr_from_dense(dense)
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(0.45_dp,0.8_dp,dp),ncv=7,maxiter=60, &
        tol=1.0e-8_dp,iterations=cycles,locked_count=locked, &
        matvecs=products,residual_norms=residuals)
    call require(info == 0 .and. cycles >= 1 .and. locked == 1, 'complex shifted Arnoldi')
    call require(abs(cvalues(1)-cmplx(0.5_dp,1.0_dp,dp))<1.0e-6_dp, 'complex interior eigenvalue')
    call require(products >= 7 .and. residuals(1)<1.0e-7_dp, 'Arnoldi shifted diagnostics')
    call require(abs(sqrt(sum(abs(matmul(dense,cvectors(:,1))-cvalues(1)*cvectors(:,1))**2)) &
        -residuals(1))<1.0e-10_dp, 'original matrix complex residual')
    call eigs(a,1,cvalues,cvectors,info,method='arnoldi',which='LM', &
        sigma=cmplx(0.45_dp,0.8_dp,dp),ncv=7,maxiter=60, &
        tol=1.0e-8_dp,residual_norms=residuals)
    call require(info == 0 .and. residuals(1)<1.0e-7_dp, 'shifted Arnoldi dispatch')
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(nan_value,0.0_dp,dp),ncv=7)
    call require(info == sparse_linalg_invalid_input, 'nonfinite complex shift')
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(0.0_dp,0.0_dp,dp),inner_restart=0)
    call require(info == sparse_linalg_invalid_input, 'invalid GMRES restart')
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(0.0_dp,0.0_dp,dp),inner_rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'invalid GMRES tolerance')
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(0.0_dp,0.0_dp,dp),inner_maxiter=0)
    call require(info == sparse_linalg_invalid_input, 'invalid GMRES limit')
    call eigs_arnoldi(a,1,cvalues,cvectors,info,which='LM', &
        sigma=cmplx(0.5_dp,1.0_dp,dp),ncv=7,inner_maxiter=25,maxiter=5)
    call require(info /= 0, 'singular complex shift not claimed converged')
    print '(a)', 'test_sparse_restarted_shift: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Condition under test.
        character(len=*), intent(in) :: message !! Failure description.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_restarted_shift
