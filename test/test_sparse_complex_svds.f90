! SPDX-License-Identifier: MIT
program test_sparse_complex_svds
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc, csr_complex_from_coo
    use scifort_sparse_linalg_complex_svds, only : svds_complex, &
        complex_svds_success, complex_svds_invalid_input, complex_svds_no_convergence
    implicit none
    type(csr_complex_matrix) :: a, large
    type(csc_complex_matrix) :: ac
    complex(dp), allocatable :: matrix(:,:), u(:,:), vh(:,:), v0(:)
    real(dp), allocatable :: s(:), residuals(:)
    integer, allocatable :: rows(:), cols(:)
    complex(dp), allocatable :: values(:)
    integer :: info, j, itn
    real(dp) :: diag(3)

    diag = [4.0_dp, 2.0_dp, 0.6_dp]
    allocate(matrix(5,3))
    matrix = cmplx(0.0_dp,0.0_dp,dp)
    do j = 1,3
        matrix(j,j) = diag(j)*exp(cmplx(0.0_dp,real(j,dp)/3.0_dp,dp))
    end do
    a = csr_complex_from_dense(matrix,info)
    call require(info == 0, 'CSR construction')
    ac = csr_complex_to_csc(a,info)
    call require(info == 0, 'CSC conversion')
    call svds_complex(a,2,u,s,vh,info,which='LM',tol=1.0e-10_dp, &
        iterations=itn,residual_norms=residuals)
    call require(info == complex_svds_success, 'complex tall LM converged')
    call require(size(s) == 2 .and. size(u,1) == 5 .and. size(vh,2) == 3, 'LM dimensions')
    call require(maxval(abs(s-[2.0_dp,4.0_dp])) < 1.0e-9_dp, 'LM values')
    call require(itn >= 2 .and. itn <= 3, 'LM iterations')
    call check_triplets(matrix,u,s,vh,residuals,1.0e-8_dp)
    call svds_complex(ac,2,u,s,vh,info,which='SM',tol=1.0e-11_dp, &
        residual_norms=residuals)
    call require(info == complex_svds_success, 'CSC small values converged')
    call require(maxval(abs(s-[0.6_dp,2.0_dp])) < 1.0e-9_dp, 'SM values')
    call check_triplets(matrix,u,s,vh,residuals,1.0e-8_dp)

    ! Wide complex matrix: U and VH must be exchanged and *conjugated* correctly.
    matrix = conjg(transpose(matrix))
    a = csr_complex_from_dense(matrix,info)
    call require(info == 0, 'wide CSR construction')
    call svds_complex(a,2,u,s,vh,info,which='LM',residual_norms=residuals)
    call require(info == complex_svds_success, 'wide converged')
    call require(maxval(abs(s-[2.0_dp,4.0_dp])) < 1.0e-9_dp, 'wide singular values')
    call check_triplets(matrix,u,s,vh,residuals,1.0e-8_dp)

    ! Complex rank deficiency and repeated singular values, including zero.
    deallocate(matrix)
    allocate(matrix(4,4))
    matrix = cmplx(0.0_dp,0.0_dp,dp)
    matrix(1,1) = cmplx(2.0_dp,1.0_dp,dp)
    matrix(2,2) = cmplx(-1.0_dp,-2.0_dp,dp)
    matrix(3,3) = cmplx(0.3_dp,-0.8_dp,dp)
    a = csr_complex_from_dense(matrix,info)
    call svds_complex(a,2,u,s,vh,info,which='SM',residual_norms=residuals)
    call require(info == complex_svds_success,'rank-deficient SM converged')
    call require(s(1) < 1.0e-10_dp .and. abs(s(2)-sqrt(0.73_dp)) < 1.0e-9_dp, &
        'rank-deficient values')
    call check_triplets(matrix,u,s,vh,residuals,2.0e-8_dp)

    matrix = cmplx(0.0_dp,0.0_dp,dp)
    a = csr_complex_from_dense(matrix,info)
    call svds_complex(a,3,u,s,vh,info,which='SM',residual_norms=residuals)
    call require(info == 0 .and. all(s == 0.0_dp),'zero matrix')
    call check_triplets(matrix,u,s,vh,residuals,1.0e-10_dp)

    call svds_complex(a,5,u,s,vh,info)
    call require(info == complex_svds_invalid_input .and. size(s) == 0, &
        'invalid k')
    call svds_complex(a,2,u,s,vh,info,tol=-0.1_dp)
    call require(info == complex_svds_invalid_input,'invalid tolerance')
    call svds_complex(a,2,u,s,vh,info,ncv=1)
    call require(info == complex_svds_invalid_input,'invalid subspace')
    call svds_complex(a,2,u,s,vh,info,which='SM',ncv=3)
    call require(info == complex_svds_invalid_input,'SM requires full subspace')
    call svds_complex(a,2,u,s,vh,info,maxiter=1)
    call require(info == complex_svds_invalid_input,'invalid iteration cap')
    allocate(v0(4))
    v0 = cmplx(0.0_dp,0.0_dp,dp)
    call svds_complex(a,2,u,s,vh,info,v0=v0)
    call require(info == complex_svds_invalid_input,'zero initial vector')
    v0 = cmplx(1.0_dp,0.0_dp,dp)
    v0(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call svds_complex(a,2,u,s,vh,info,v0=v0)
    call require(info == complex_svds_invalid_input,'nonfinite initial vector')

    ! One-step Ritz vector cannot converge for the nontrivial diagonal problem.
    matrix = cmplx(0.0_dp,0.0_dp,dp)
    matrix(1,1) = cmplx(3.0_dp,0.5_dp,dp)
    matrix(2,2) = cmplx(2.0_dp,1.0_dp,dp)
    matrix(3,3) = cmplx(1.0_dp,-0.5_dp,dp)
    a = csr_complex_from_dense(matrix,info)
    call svds_complex(a,1,u,s,vh,info,ncv=1,maxiter=1, &
        residual_norms=residuals)
    call require(info == complex_svds_no_convergence, 'report unfinished partial SVD')
    call require(size(s) == 1 .and. residuals(1) > 1.0e-4_dp, &
        'nonconverged Ritz residual visible')

    ! Large, genuinely sparse rectangular matrix: no full dense A is built.
    allocate(rows(240),cols(240),values(240))
    do j = 1,120
        rows(j) = j
        cols(j) = j
        rows(120+j) = 120+j
        cols(120+j) = j
        if (j == 1) then
            values(j) = cmplx(18.0_dp,2.0_dp,dp)
        else if (j == 2) then
            values(j) = cmplx(10.0_dp,-3.0_dp,dp)
        else
            values(j) = cmplx(0.05_dp+0.007_dp*j,0.03_dp,dp)
        end if
        values(120+j) = cmplx(0.01_dp,-0.02_dp,dp)
    end do
    large = csr_complex_from_coo(4000,120,rows,cols,values,info)
    call require(info == 0, 'large sparse CSR')
    call svds_complex(large,2,u,s,vh,info,ncv=28,maxiter=28, &
        tol=1.0e-9_dp,iterations=itn,residual_norms=residuals)
    call require(info == complex_svds_success, 'large sparse partial convergence')
    call require(itn < 28 .and. maxval(residuals) < 1.0e-8_dp, &
        'large sparse residuals')
    call require(abs(s(1)-sqrt(109.0_dp+0.0005_dp)) < 1.0e-8_dp, &
        'large sparse second singular value')
    call require(abs(s(2)-sqrt(328.0_dp+0.0005_dp)) < 1.0e-8_dp, &
        'large sparse first singular value')

    print '(a)', 'PASS complex CSR/CSC sparse SVD native tests'
contains
    subroutine check_triplets(matrix,u,s,vh,residuals,tol)
        complex(dp), intent(in) :: matrix(:,:) !! Matrix A (m,n).
        complex(dp), intent(in) :: u(:,:) !! Left columns (m,k).
        real(dp), intent(in) :: s(:) !! Singular values (k).
        complex(dp), intent(in) :: vh(:,:) !! Right conjugate-transposed rows (k,n).
        real(dp), intent(in) :: residuals(:) !! Reported true residuals (k).
        real(dp), intent(in) :: tol !! Absolute error threshold.
        complex(dp), allocatable :: r1(:), r2(:), gram(:,:)
        integer :: k, j
        real(dp) :: actual
        k = size(s)
        call require(size(u,2) == k .and. size(vh,1) == k, 'triplet shapes')
        do j = 1,k
            r1 = matmul(matrix,conjg(vh(j,:)))-s(j)*u(:,j)
            r2 = matmul(conjg(transpose(matrix)),u(:,j))-s(j)*conjg(vh(j,:))
            actual = hypot(norm2(abs(r1)),norm2(abs(r2)))
            call require(actual < tol, 'SVD two-sided residual')
            call require(abs(actual-residuals(j)) < tol, 'reported residual')
        end do
        gram = matmul(conjg(transpose(u)),u)
        do j = 1,k
            gram(j,j) = gram(j,j)-cmplx(1.0_dp,0.0_dp,dp)
        end do
        call require(maxval(abs(gram)) < tol, 'left orthogonality')
        gram = matmul(vh,conjg(transpose(vh)))
        do j = 1,k
            gram(j,j) = gram(j,j)-cmplx(1.0_dp,0.0_dp,dp)
        end do
        call require(maxval(abs(gram)) < tol, 'right orthogonality')
    end subroutine check_triplets

    subroutine require(ok,message)
        logical, intent(in) :: ok !! Assertion predicate.
        character(len=*), intent(in) :: message !! Assertion label.
        if (.not. ok) then
            print '(a)', 'FAIL complex svds: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_complex_svds
