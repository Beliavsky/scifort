! SPDX-License-Identifier: MIT
! Partial sparse singular value decomposition: rectangular, repeated, rank-deficient.
program test_sparse_svds
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, &
        csr_matvec, csr_transpose
    use scifort_sparse_linalg, only : svds, sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a, at
    real(dp) :: x(8,6), nan_value
    real(dp), allocatable :: u(:,:), s(:), vt(:,:), res(:), v(:), tmp(:)
    real(dp), allocatable :: entries(:)
    integer, allocatable :: rows(:), cols(:)
    integer :: info, it, i, j

    x = 0.0_dp
    do i = 1,6
        x(i,i) = real(10-i,dp)
    end do
    a = csr_from_dense(x)
    call svds(a,2,u,s,vt,info,which='LM',ncv=6,tol=1.0e-9_dp, &
        iterations=it,residual_norms=res)
    call require(info == 0 .and. it <= 6, 'tall largest convergence')
    call close(s,[8.0_dp,9.0_dp],1.0e-9_dp,'tall largest values')
    call check_triplets(a,u,s,vt,res,2.0e-8_dp)

    at = csr_transpose(a)
    call svds(at,2,u,s,vt,info,which='LM',ncv=6,tol=1.0e-9_dp, &
        iterations=it,residual_norms=res)
    call require(info == 0, 'wide largest convergence')
    call close(s,[8.0_dp,9.0_dp],1.0e-9_dp,'wide largest values')
    call check_triplets(at,u,s,vt,res,2.0e-8_dp)
    call svds(at,2,u,s,vt,info,which='SM',ncv=6,tol=1.0e-9_dp, &
        residual_norms=res)
    call require(info == 0, 'wide smallest convergence')
    call close(s,[4.0_dp,5.0_dp],1.0e-9_dp,'wide smallest values')
    call check_triplets(at,u,s,vt,res,2.0e-8_dp)

    x = 0.0_dp
    x(1,1) = 5.0_dp
    x(2,2) = 5.0_dp
    x(3,3) = 2.0_dp
    x(4,4) = 1.0_dp
    a = csr_from_dense(x)
    call svds(a,2,u,s,vt,info,which='LM',ncv=6,tol=1.0e-9_dp, &
        residual_norms=res)
    call require(info == 0, 'repeated largest convergence')
    call close(s,[5.0_dp,5.0_dp],1.0e-9_dp,'repeated largest values')
    call check_triplets(a,u,s,vt,res,2.0e-8_dp)
    call svds(a,2,u,s,vt,info,which='SM',ncv=6,tol=1.0e-9_dp, &
        residual_norms=res)
    call require(info == 0, 'rank-deficient smallest convergence')
    call close(s,[0.0_dp,0.0_dp],1.0e-9_dp,'rank-deficient smallest values')
    call check_triplets(a,u,s,vt,res,2.0e-8_dp)

    x = 0.0_dp
    a = csr_from_dense(x)
    call svds(a,3,u,s,vt,info,which='SM',iterations=it,residual_norms=res)
    call require(info == 0 .and. it == 0 .and. all(s == 0.0_dp), 'zero matrix')
    call check_triplets(a,u,s,vt,res,1.0e-12_dp)

    ! Partial Golub--Kahan on a 1200-by-900 sparse matrix: no dense Gram.
    block
        integer, parameter :: m = 1200, n = 900
        allocate(rows(2*n-1),cols(2*n-1),entries(2*n-1))
        do i = 1,n
            rows(i) = i
            cols(i) = i
            entries(i) = 0.2_dp+0.8_dp*real(i,dp)/real(n,dp)
        end do
        entries(1) = 40.0_dp
        entries(2) = 25.0_dp
        do i = 1,n-1
            rows(n+i) = i
            cols(n+i) = i+1
            entries(n+i) = 0.002_dp
        end do
        a = csr_from_coo(m,n,rows,cols,entries)
    end block
    call svds(a,2,u,s,vt,info,which='LM',ncv=35,tol=1.0e-9_dp, &
        iterations=it,residual_norms=res)
    call require(info == 0 .and. it <= 35, 'large sparse partial convergence')
    call require(maxval(abs(s-[25.0_dp,40.0_dp])) < 1.0e-4_dp, 'large sparse singular values')
    call check_triplets(a,u,s,vt,res,1.0e-7_dp)

    call svds(a,2,u,s,vt,info,which='XYZ')
    call require(info == sparse_linalg_invalid_input, 'invalid selection')
    call svds(a,0,u,s,vt,info)
    call require(info == sparse_linalg_invalid_input, 'invalid k')
    call svds(a,901,u,s,vt,info)
    call require(info == sparse_linalg_invalid_input, 'k exceeds min dimension')
    call svds(a,2,u,s,vt,info,ncv=1)
    call require(info == sparse_linalg_invalid_input, 'invalid search width')
    call svds(a,2,u,s,vt,info,which='SM',ncv=15)
    call require(info == sparse_linalg_invalid_input, 'SM requires full search')
    call svds(a,2,u,s,vt,info,maxiter=1)
    call require(info == sparse_linalg_invalid_input, 'insufficient iteration budget')
    call svds(a,2,u,s,vt,info,v0=[1.0_dp,2.0_dp])
    call require(info == sparse_linalg_invalid_input, 'invalid v0 size')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call svds(a,2,u,s,vt,info,tol=nan_value)
    call require(info == sparse_linalg_invalid_input, 'NaN tolerance')

    ! An insufficient trial space must return nonconvergence rather than
    ! claiming convergence on unconverged projected singular values.
    block
        real(dp) :: dense(6,6)
        dense = reshape([ 1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
            2.0_dp,4.0_dp,1.0_dp,2.0_dp,7.0_dp,3.0_dp, &
            4.0_dp,1.0_dp,5.0_dp,3.0_dp,2.0_dp,8.0_dp, &
            6.0_dp,3.0_dp,2.0_dp,1.0_dp,4.0_dp,7.0_dp, &
            3.0_dp,8.0_dp,6.0_dp,1.0_dp,9.0_dp,2.0_dp, &
            1.0_dp,7.0_dp,3.0_dp,9.0_dp,2.0_dp,4.0_dp ],[6,6])
        a = csr_from_dense(dense)
    end block
    call svds(a,2,u,s,vt,info,ncv=2,tol=1.0e-12_dp, &
        iterations=it,residual_norms=res)
    call require(info > 0 .and. it == 2, 'true residual detects nonconvergence')
    call require(maxval(res) > 1.0e-9_dp, 'nonconverged true residual retained')

    print '(a)', 'test_sparse_svds: PASS'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Test assertion.
        character(len=*), intent(in) :: message !! Diagnostic description.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine close(actual,expected,threshold,message)
        real(dp), intent(in) :: actual(:) !! Numerical results.
        real(dp), intent(in) :: expected(:) !! Reference results.
        real(dp), intent(in) :: threshold !! Absolute tolerance.
        character(len=*), intent(in) :: message !! Diagnostic description.
        call require(size(actual)==size(expected), 'shape '//message)
        call require(maxval(abs(actual-expected)) < threshold, message)
    end subroutine close

    subroutine check_triplets(mat,left,singular,right_t,reported,threshold)
        type(csr_matrix), intent(in) :: mat !! CSR matrix under test.
        real(dp), intent(in) :: left(:,:) !! Left singular vectors.
        real(dp), intent(in) :: singular(:) !! Singular values.
        real(dp), intent(in) :: right_t(:,:) !! Transposed right singular vectors.
        real(dp), intent(in) :: reported(:) !! Returned two-sided residuals.
        real(dp), intent(in) :: threshold !! Required residual and orthogonality.
        type(csr_matrix) :: trans
        real(dp), allocatable :: a_v(:), at_u(:)
        real(dp) :: true_residual
        integer :: j, l

        trans = csr_transpose(mat)
        call require(size(left,1) == mat%n_rows, 'left vector shape')
        call require(size(right_t,2) == mat%n_cols, 'right vector shape')
        do j=1,size(singular)
            a_v = csr_matvec(mat,right_t(j,:))
            at_u = csr_matvec(trans,left(:,j))
            true_residual = hypot(norm2(a_v-singular(j)*left(:,j)), &
                norm2(at_u-singular(j)*right_t(j,:)))
            call require(true_residual < threshold, 'true singular residual')
            call require(abs(true_residual-reported(j)) < 5.0e-10_dp, 'reported residual')
            call require(abs(norm2(left(:,j))-1.0_dp) < threshold, 'left norm')
            call require(abs(norm2(right_t(j,:))-1.0_dp) < threshold, 'right norm')
            do l=1,j-1
                call require(abs(dot_product(left(:,j),left(:,l))) < threshold, &
                    'left orthogonality')
                call require(abs(dot_product(right_t(j,:),right_t(l,:))) < threshold, &
                    'right orthogonality')
            end do
        end do
    end subroutine check_triplets
end program test_sparse_svds
