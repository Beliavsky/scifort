! SPDX-License-Identifier: MIT
program test_sparse_eigsh
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : eigsh, sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a
    real(dp) :: matrix(6,6), values_diag(6)
    real(dp), allocatable :: values(:), vectors(:,:), residuals(:), y(:)
    integer :: info, i, j, steps, n
    character(len=2) :: mode

    values_diag = [-8.0_dp,-3.0_dp,0.0_dp,1.0_dp,4.0_dp,7.0_dp]
    matrix = 0.0_dp
    do i = 1,6
        matrix(i,i) = values_diag(i)
    end do
    a = csr_from_dense(matrix)
    do j = 1,4
        select case(j)
        case(1)
            mode = 'LM'
        case(2)
            mode = 'SM'
        case(3)
            mode = 'LA'
        case(4)
            mode = 'SA'
        end select
        call eigsh(a,2,values,vectors,info,which=mode, &
            tol=1.0e-11_dp,residual_norms=residuals,iterations=steps)
        call require(info == 0, 'full-rank small diagonal eigsh '//mode)
        call require(size(values) == 2 .and. size(vectors,1) == 6 .and. &
            size(vectors,2) == 2 .and. steps >= 1, 'correct output shapes')
        select case(j)
        case(1)
            call close(values,[-8.0_dp,7.0_dp],1.0e-9_dp,'LM values')
        case(2)
            call close(values,[0.0_dp,1.0_dp],1.0e-9_dp,'SM values')
        case(3)
            call close(values,[4.0_dp,7.0_dp],1.0e-9_dp,'LA values')
        case(4)
            call close(values,[-8.0_dp,-3.0_dp],1.0e-9_dp,'SA values')
        end select
        call require(maxval(residuals) < 1.0e-9_dp, 'eigen residuals')
        call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-10_dp, &
            'orthogonal eigenvectors')
    end do
    call eigsh(a,0,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'invalid k')
    call eigsh(a,6,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'k>=n')
    call eigsh(a,2,values,vectors,info,which='XX')
    call require(info == sparse_linalg_invalid_input, 'unknown which')
    call eigsh(a,2,values,vectors,info,maxiter=0)
    call require(info > 0 .and. size(values) == 2, 'zero update budget')
    call eigsh(a,2,values,vectors,info,ncv=2)
    call require(info == sparse_linalg_invalid_input, 'ncv<=k')
    call eigsh(a,2,values,vectors,info,tol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'negative tolerance')
    call eigsh(a,2,values,vectors,info,v0=[1.0_dp])
    call require(info == sparse_linalg_invalid_input, 'bad initial vector')
    matrix(1,2) = 0.5_dp
    a = csr_from_dense(matrix)
    call eigsh(a,2,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'non-symmetric matrix rejected')
    a = csr_matrix()
    call eigsh(a,2,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR rejected')

    ! Moderate n, small subspace (not full dense eigendecomposition).
    n = 100
    a = csr_from_coo(n,n,[(i,i=1,n), (i,i=1,n-1), (i,i=2,n)], &
        [(i,i=1,n), (i+1,i=1,n-1), (i-1,i=2,n)], &
        [ (2.0_dp+real(i,dp)/30.0_dp,i=1,n), &
          (-0.1_dp,i=1,n-1), (-0.1_dp,i=2,n)])
    call eigsh(a,2,values,vectors,info,which='LA',tol=1.0e-7_dp, &
        ncv=16,maxiter=2000,residual_norms=residuals,iterations=steps)
    call require(info == 0, 'large CSR LA eigsh')
    call require(maxval(residuals) < 1.0e-6_dp, 'large sparse eigenpair residuals')
    do j = 1,2
        y = csr_matvec(a,vectors(:,j))
        call require(norm2(y-values(j)*vectors(:,j)) < 1.0e-6_dp, &
            'large eigenvector true residual')
    end do
    ! Shift-invert around sigma: selection uses distance to sigma, not |lambda|.
    matrix=0.0_dp
    do i=1,6
        matrix(i,i)=real(i-3,dp)
    end do
    a=csr_from_dense(matrix)
    call eigsh(a,2,values,vectors,info,which='LM',sigma=0.4_dp, &
        tol=1.0e-10_dp,ncv=6,inner_rtol=1.0e-12_dp, &
        residual_norms=residuals,iterations=steps)
    call require(info == 0,'shift-invert full-subspace success')
    call close(values,[0.0_dp,1.0_dp],1.0e-8_dp,'shift-invert closest eigenvalues')
    call require(maxval(residuals) < 1.0e-8_dp,'shift-invert true residual')
    call eigsh(a,2,values,vectors,info,which='LM',sigma=0.0_dp, &
        inner_maxiter=60)
    call require(info /= 0,'singular exact shift not declared converged')
    call eigsh(a,2,values,vectors,info,which='LA',sigma=0.4_dp, &
        ncv=6,tol=1.0e-10_dp,residual_norms=residuals)
    call require(info == 0,'shift-invert largest-algebraic mode')
    call close(values,[1.0_dp,2.0_dp],1.0e-8_dp,'shift-invert LA transformed')
    call eigsh(a,2,values,vectors,info,which='SA',sigma=0.4_dp, &
        ncv=6,tol=1.0e-10_dp,residual_norms=residuals)
    call require(info == 0,'shift-invert smallest-algebraic mode')
    call close(values,[-1.0_dp,0.0_dp],1.0e-8_dp,'shift-invert SA transformed')
    call eigsh(a,2,values,vectors,info,which='SM',sigma=0.4_dp, &
        ncv=6,tol=1.0e-10_dp,residual_norms=residuals)
    call require(info == 0,'shift-invert smallest-magnitude mode')
    call close(values,[-2.0_dp,3.0_dp],1.0e-8_dp,'shift-invert SM transformed')
    call eigsh(a,2,values,vectors,info,which='SM',sigma=0.4_dp,ncv=4)
    call require(info == sparse_linalg_invalid_input, &
        'shift-invert SM reduced subspace intentionally unsupported')
    call eigsh(a,2,values,vectors,info,sigma=0.4_dp,inner_rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input,'negative inner tolerance')
    call eigsh(a,2,values,vectors,info,sigma=0.4_dp,inner_maxiter=-1)
    call require(info == sparse_linalg_invalid_input,'negative inner maxiter')

    n=60
    a=csr_from_coo(n,n,[(i,i=1,n), (i,i=1,n-1), (i,i=2,n)], &
        [(i,i=1,n), (i+1,i=1,n-1), (i-1,i=2,n)], &
        [(real(i-30,dp)/5.0_dp,i=1,n),(-0.015_dp,i=1,n-1), &
         (-0.015_dp,i=2,n)])
    call eigsh(a,2,values,vectors,info,which='LM',sigma=0.43_dp, &
        tol=1.0e-6_dp,ncv=18,maxiter=100,inner_maxiter=200, &
        residual_norms=residuals)
    call require(info == 0,'shift-invert reduced-subspace success')
    call close(values,[0.4_dp,0.6_dp],1.0e-3_dp,'shift-invert interior eigenvalues')
    call require(maxval(residuals) < 1.0e-5_dp,'shift-invert reduced true residuals')
    print '(a)', 'test_sparse_eigsh: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Required invariant.
        character(len=*), intent(in) :: message !! Assertion message.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
    subroutine close(actual,expected,tolerance,message)
        real(dp), intent(in) :: actual(:) !! Values obtained.
        real(dp), intent(in) :: expected(:) !! Expected values.
        real(dp), intent(in) :: tolerance !! Absolute tolerance.
        character(len=*), intent(in) :: message !! Assertion message.
        call require(size(actual) == size(expected),message//' size')
        call require(maxval(abs(actual-expected)) <= tolerance,message)
    end subroutine close
end program test_sparse_eigsh
