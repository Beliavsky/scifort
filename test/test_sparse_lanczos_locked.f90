! SPDX-License-Identifier: MIT
! Hard-locking/orthogonal-deflation regression for symmetric sparse eigensolvers.
program test_sparse_lanczos_locked
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : eigsh, eigsh_lanczos, &
        sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a
    real(dp) :: dense(12,12), diagonal(12), nan_value
    real(dp), allocatable :: vals(:), vectors(:,:), errors(:), product(:)
    integer :: info, steps, locks, products, i, mode_index
    character(len=2) :: choice
    integer, parameter :: n = 80

    dense = 0.0_dp
    diagonal = [-8.0_dp,-4.0_dp,-2.0_dp,-1.0_dp,-0.2_dp, &
        0.8_dp,1.8_dp,3.0_dp,4.0_dp,5.0_dp,7.0_dp,9.0_dp]
    do i = 1,12
        dense(i,i) = diagonal(i)
    end do
    a = csr_from_dense(dense)
    do mode_index = 1,4
        select case (mode_index)
        case (1)
            choice = 'LA'
        case (2)
            choice = 'SA'
        case (3)
            choice = 'LM'
        case (4)
            choice = 'SM'
        end select
        call eigsh_lanczos(a,2,vals,vectors,info,which=choice, &
            ncv=12,tol=1.0e-10_dp,maxiter=10,hard_locking=.true., &
            iterations=steps,locked_count=locks,matvecs=products,residual_norms=errors)
        call require(info == 0 .and. locks == 2, 'locked full-subspace '//choice)
        call require(steps >= 1 .and. products > 0, 'locking diagnostics')
        call require(maxval(errors) < 1.0e-9_dp, 'locked residuals')
        select case (choice)
        case ('LA')
            call close(vals,[7.0_dp,9.0_dp],1.0e-9_dp,'LA')
        case ('SA')
            call close(vals,[-8.0_dp,-4.0_dp],1.0e-9_dp,'SA')
        case ('LM')
            call close(vals,[-8.0_dp,9.0_dp],1.0e-9_dp,'LM')
        case ('SM')
            call close(vals,[-0.2_dp,0.8_dp],1.0e-9_dp,'SM')
        end select
        call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-11_dp, &
            'locked vectors are orthogonal')
    end do

    ! The first extreme eigenvector converges in one cycle, while the
    ! second needs more cycles.  This exercises actual intermediate locking.
    block
        integer :: rows(3*n-4), cols(3*n-4)
        real(dp) :: entries(3*n-4)
        do i = 1,n
            rows(i) = i
            cols(i) = i
            entries(i) = 2.0_dp+real(i,dp)/30.0_dp
        end do
        entries(1) = 40.0_dp
        do i = 2,n-1
            rows(n+i-1) = i
            cols(n+i-1) = i+1
            entries(n+i-1) = -0.1_dp
            rows(2*n-2+i-1) = i+1
            cols(2*n-2+i-1) = i
            entries(2*n-2+i-1) = -0.1_dp
        end do
        a = csr_from_coo(n,n,rows,cols,entries)
    end block
    call eigsh_lanczos(a,2,vals,vectors,info,which='LA',ncv=12, &
        maxiter=1,tol=1.0e-9_dp,hard_locking=.true., &
        locked_count=locks,iterations=steps,residual_norms=errors)
    call require(info > 0 .and. steps == 1 .and. locks == 1, &
        'one mode hard-locked before nonconvergence')
    call eigsh_lanczos(a,2,vals,vectors,info,which='LA',ncv=12, &
        maxiter=80,tol=1.0e-9_dp,hard_locking=.true., &
        locked_count=locks,iterations=steps,residual_norms=errors)
    call require(info == 0 .and. locks == 2 .and. steps > 1, &
        'remaining eigenpair converges after deflation')
    call require(abs(vals(2)-40.0_dp) < 1.0e-8_dp, 'isolated largest eigenvalue')
    call require(maxval(errors) < 6.0e-8_dp, 'deflated eigenpair residuals')
    do i = 1,2
        product = csr_matvec(a,vectors(:,i))
        call require(abs(norm2(product-vals(i)*vectors(:,i))-errors(i)) < 1.0e-10_dp, &
            'true deflated residual')
    end do
    call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-10_dp, &
        'deflated Ritz-vector orthogonality')
    call eigsh(a,2,vals,vectors,info,method='lanczos_locked',which='LA', &
        ncv=12,maxiter=80,tol=1.0e-9_dp,residual_norms=errors)
    call require(info == 0 .and. maxval(errors)<6.0e-8_dp, 'locked eigsh method dispatch')

    ! Real shift-invert operates in a deflated inverse subspace, but
    ! validates eigenpairs using A rather than (A-sigma*I)**(-1).
    a = csr_from_dense(dense)
    call eigsh(a,2,vals,vectors,info,method='lanczos_locked', &
        sigma=0.0_dp,which='LM',ncv=12,maxiter=20,tol=1.0e-8_dp, &
        residual_norms=errors)
    call require(info == 0, 'deflated shift-invert dispatch')
    call close(vals,[-0.2_dp,0.8_dp],1.0e-8_dp,'shift-invert closest eigenvalues')
    call eigsh_lanczos(a,1,vals,vectors,info,which='LA',sigma=0.0_dp, &
        ncv=12,maxiter=15,tol=1.0e-8_dp,hard_locking=.true.)
    call require(info == 0 .and. abs(vals(1)-0.8_dp)<1.0e-8_dp, &
        'inverse largest algebraic')
    call eigsh_lanczos(a,1,vals,vectors,info,which='SA',sigma=0.0_dp, &
        ncv=12,maxiter=15,tol=1.0e-8_dp,hard_locking=.true.)
    call require(info == 0 .and. abs(vals(1)+0.2_dp)<1.0e-8_dp, &
        'inverse smallest algebraic')
    call eigsh_lanczos(a,1,vals,vectors,info,which='SM',sigma=0.0_dp, &
        ncv=12,maxiter=15,tol=1.0e-8_dp,hard_locking=.true.)
    call require(info == 0 .and. abs(vals(1)-9.0_dp)<1.0e-8_dp, &
        'inverse smallest magnitude is farthest from shift')

    ! A single-vector Krylov space can represent only one direction from
    ! an exactly repeated eigenvalue.  Two deterministic seed vectors
    ! recover the full double eigenspace, with genuine hard locking.
    block
        real(dp) :: repeated(40,40)
        repeated = 0.0_dp
        do i = 1,40
            repeated(i,i) = 3.0_dp/real(i,dp)
        end do
        repeated(1,1) = 10.0_dp
        repeated(2,2) = 10.0_dp
        repeated(3,3) = 9.0_dp
        a = csr_from_dense(repeated)
        call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true., &
            block_size=2,which='LA',ncv=12,maxiter=40,tol=1.0e-9_dp, &
            locked_count=locks,iterations=steps,residual_norms=errors)
        call require(info == 0 .and. locks == 2 .and. steps > 1, &
            'block-Krylov hard locking for repeated eigenvalues')
        call close(vals,[10.0_dp,10.0_dp],1.0e-7_dp,'repeated eigenvalue')
        call require(maxval(errors) < 1.0e-8_dp, 'repeated eigenpair residuals')
        call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-9_dp, &
            'independent repeated-eigenspace vectors')
    end block
    a = csr_from_dense(dense)
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true.,block_size=0)
    call require(info == sparse_linalg_invalid_input, 'invalid zero block size')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true.,block_size=3)
    call require(info == sparse_linalg_invalid_input, 'too many initial seeds')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.false.,block_size=2)
    call require(info == sparse_linalg_invalid_input, 'block size requires hard locking')

    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true.,ncv=2)
    call require(info == sparse_linalg_invalid_input, 'invalid locked subspace width')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true., &
        ncv=4,which='SM')
    call require(info == sparse_linalg_invalid_input, 'SM needs complete subspace')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true., &
        maxiter=0,iterations=steps)
    call require(info > 0 .and. steps == 0, 'zero locked cycle budget')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true.,tol=nan_value)
    call require(info == sparse_linalg_invalid_input, 'NaN locked tolerance')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true., &
        sigma=0.0_dp,inner_maxiter=0)
    call require(info == sparse_linalg_invalid_input, 'invalid inner MINRES limit')
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true., &
        sigma=0.8_dp,ncv=12,maxiter=4)
    call require(info /= 0, 'singular shift never claimed successful')
    dense(1,2) = 2.0_dp
    a = csr_from_dense(dense)
    call eigsh_lanczos(a,2,vals,vectors,info,hard_locking=.true.)
    call require(info == sparse_linalg_invalid_input, 'asymmetric matrix rejected')

    print '(a)', 'test_sparse_lanczos_locked: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Condition that must be satisfied.
        character(len=*), intent(in) :: message !! Error context.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine close(actual, expected, threshold, message)
        real(dp), intent(in) :: actual(:) !! Selected eigenvalues.
        real(dp), intent(in) :: expected(:) !! Target eigenvalues.
        real(dp), intent(in) :: threshold !! Permitted absolute error.
        character(len=*), intent(in) :: message !! Comparison label.
        call require(size(actual) == size(expected), 'shape '//message)
        call require(maxval(abs(actual-expected)) < threshold, message)
    end subroutine close
end program test_sparse_lanczos_locked
