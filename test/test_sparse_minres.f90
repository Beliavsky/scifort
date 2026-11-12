! SPDX-License-Identifier: MIT
program test_sparse_minres
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : minres, sparse_linalg_breakdown, &
        sparse_linalg_invalid_input, sparse_linalg_bad_diagonal
    implicit none
    type(csr_matrix) :: a
    real(dp) :: dense(4,4), rn
    real(dp), allocatable :: b(:), x(:), exact(:)
    integer :: info, steps, i, n

    dense = 0.0_dp
    dense(1,1) = -4.0_dp
    dense(2,2) = 3.0_dp
    dense(3,3) = 5.0_dp
    dense(4,4) = -6.0_dp
    do i = 1,3
        dense(i,i+1) = 0.35_dp
        dense(i+1,i) = 0.35_dp
    end do
    a = csr_from_dense(dense)
    exact = [1.0_dp,-2.0_dp,3.0_dp,-4.0_dp]
    b = matmul(dense,exact)
    call minres(a,b,x,info,rtol=1.0e-12_dp,iterations=steps,residual_norm=rn, &
        check_symmetry=.true.)
    call require(info == 0 .and. steps > 0, 'indefinite MINRES success')
    call require(maxval(abs(x-exact)) < 1.0e-10_dp .and. rn < 1.0e-10_dp, 'indefinite solution')
    call minres(a,b,x,info,rtol=1.0e-12_dp,x0=exact,iterations=steps)
    call require(info == 0 .and. steps == 0 .and. all(x == exact), 'exact starting solution')
    b = matmul(dense,exact)-1.25_dp*exact
    call minres(a,b,x,info,rtol=1.0e-12_dp,shift=1.25_dp)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp, 'shifted symmetric solve')
    call minres(a,b,x,info,rtol=1.0e-12_dp,shift=1.25_dp,maxiter=0)
    call require(info > 0 .and. maxval(abs(x)) == 0.0_dp, 'exhausted iteration cap')
    call minres(a,b,x,info,shift=ieee_value(0.0_dp,ieee_quiet_nan))
    call require(info == sparse_linalg_invalid_input, 'invalid shift')
    call minres(a,b,x,info,rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'invalid tolerance')
    call minres(a,b,x,info,jacobi=.true.)
    call require(info == sparse_linalg_bad_diagonal, 'indefinite Jacobi rejected')

    dense(1,2) = 1.0_dp
    a = csr_from_dense(dense)
    call minres(a,b,x,info,check_symmetry=.true.)
    call require(info == sparse_linalg_invalid_input, 'asymmetric input with optional check')

    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,2.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,0.0_dp,3.0_dp],[4,4]))
    call minres(a,[1.0_dp,0.0_dp,6.0_dp,9.0_dp],x,info,rtol=1.0e-12_dp)
    call require(info == 0, 'consistent singular system')
    call require(maxval(abs(x-[1.0_dp,0.0_dp,3.0_dp,3.0_dp])) < 1.0e-10_dp, &
        'singular solution')
    call minres(a,[1.0_dp,1.0_dp,6.0_dp,9.0_dp],x,info,rtol=1.0e-12_dp,maxiter=12)
    call require(info /= 0, 'inconsistent singular system not falsely converged')

    dense = 0.0_dp
    do i = 1,4
        dense(i,i) = real(i+3,dp)
    end do
    dense(1,2) = -1.0_dp
    dense(2,1) = -1.0_dp
    a = csr_from_dense(dense)
    exact = [-1.0_dp,2.0_dp,-3.0_dp,4.0_dp]
    b = matmul(dense,exact)
    call minres(a,b,x,info,rtol=1.0e-12_dp,jacobi=.true.,residual_norm=rn)
    call require(info == 0 .and. rn < 1.0e-10_dp, 'Jacobi preconditioned SPD MINRES')
    call require(maxval(abs(x-exact)) < 1.0e-10_dp, 'Jacobi solution accuracy')
    call minres(a,0.0_dp*b,x,info)
    call require(info == 0 .and. all(x == 0.0_dp), 'zero rhs')
    a = csr_matrix()
    call minres(a,[1.0_dp],x,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR')
    a = csr_from_coo(0,0,[integer::],[integer::],[real(dp)::])
    call minres(a,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0, 'empty matrix')
    a = csr_from_dense(reshape([0.0_dp],[1,1]))
    call minres(a,[1.0_dp],x,info)
    call require(info == sparse_linalg_breakdown, 'zero operator breakdown')

    ! An indefinite tridiagonal problem with O(n) sparse storage.
    n = 500
    a = csr_from_coo(n,n,[(i,i=1,n), (i,i=1,n-1), (i,i=2,n)], &
        [(i,i=1,n), (i+1,i=1,n-1), (i-1,i=2,n)], &
        [ (merge(-5.0_dp,5.0_dp,mod(i,2)==0),i=1,n), &
          (-0.25_dp,i=1,n-1), (-0.25_dp,i=2,n)])
    exact = [(sin(real(i,dp)*0.15_dp),i=1,n)]
    b = csr_matvec(a,exact)
    call minres(a,b,x,info,rtol=1.0e-11_dp,check_symmetry=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-8_dp, 'large indefinite CSR')

    print '(a)', 'test_sparse_minres: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Required assertion.
        character(len=*), intent(in) :: message !! Failure diagnostic.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_minres
