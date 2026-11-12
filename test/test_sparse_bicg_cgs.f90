! SPDX-License-Identifier: MIT
program test_sparse_bicg_cgs
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_eye
    use scifort_sparse_linalg, only : bicg, cgs, ilu0_factor, ilu0, &
        sparse_linalg_invalid_input, sparse_linalg_bad_diagonal, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a, empty
    type(ilu0_factor) :: factor
    real(dp) :: mat(4,4), rhs(4), expected(4), norm, bad
    real(dp), allocatable :: x(:)
    integer :: status, iterations

    mat = 0.0_dp
    mat(1,:) = [6.0_dp,2.0_dp,0.0_dp,-1.0_dp]
    mat(2,:) = [1.0_dp,5.0_dp,2.0_dp,0.0_dp]
    mat(3,:) = [0.0_dp,-1.0_dp,4.0_dp,1.0_dp]
    mat(4,:) = [1.0_dp,0.0_dp,2.0_dp,7.0_dp]
    a = csr_from_dense(mat)
    expected = [2.0_dp,-1.0_dp,3.0_dp,0.5_dp]
    rhs = matmul(mat,expected)
    call bicg(a,rhs,x,status,rtol=1.0e-12_dp,atol=1.0e-13_dp, &
        iterations=iterations,residual_norm=norm)
    call require(status == 0 .and. iterations >= 1 .and. iterations <= 4, 'BiCG status')
    call require(maxval(abs(x-expected)) < 1.0e-11_dp, 'BiCG solution')
    call require(norm < 1.0e-11_dp, 'BiCG true residual')
    call bicg(a,rhs,x,status,rtol=1.0e-12_dp,jacobi=.true.)
    call require(status == 0 .and. maxval(abs(x-expected)) < 1.0e-11_dp, 'BiCG Jacobi')
    call bicg(a,rhs,x,status,x0=expected,iterations=iterations,residual_norm=norm)
    call require(status == 0 .and. iterations == 0 .and. norm < 1.0e-13_dp, 'BiCG x0')

    call cgs(a,rhs,x,status,rtol=1.0e-12_dp,atol=1.0e-13_dp, &
        iterations=iterations,residual_norm=norm)
    call require(status == 0 .and. iterations >= 1, 'CGS status')
    call require(maxval(abs(x-expected)) < 1.0e-11_dp, 'CGS solution')
    call require(norm < 1.0e-11_dp, 'CGS true residual')
    call cgs(a,rhs,x,status,rtol=1.0e-12_dp,jacobi=.true.)
    call require(status == 0 .and. maxval(abs(x-expected)) < 1.0e-11_dp, 'CGS Jacobi')
    call ilu0(a,factor,status)
    call require(status == 0, 'ILU0 factorization')
    call cgs(a,rhs,x,status,rtol=1.0e-12_dp,ilu=factor)
    call require(status == 0 .and. maxval(abs(x-expected)) < 1.0e-11_dp, 'CGS ILU0')
    call cgs(a,rhs,x,status,x0=expected,iterations=iterations)
    call require(status == 0 .and. iterations == 0, 'CGS x0')

    empty = csr_eye(0)
    call bicg(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'BiCG empty')
    call cgs(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'CGS empty')
    call bicg(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'BiCG invalid rhs length')
    call cgs(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'CGS invalid rhs length')
    call bicg(a,rhs,x,status,maxiter=-1)
    call require(status == sparse_linalg_invalid_input, 'BiCG invalid cap')
    call cgs(a,rhs,x,status,rtol=-1.0_dp)
    call require(status == sparse_linalg_invalid_input, 'CGS invalid tolerance')
    call bicg(a,rhs,x,status,maxiter=0,iterations=iterations)
    call require(status > 0 .and. iterations == 0, 'BiCG zero-iteration cap')
    call cgs(a,rhs,x,status,maxiter=0,iterations=iterations)
    call require(status > 0 .and. iterations == 0, 'CGS zero-iteration cap')
    call bicg(a,rhs,x,status,maxiter=1,rtol=0.0_dp)
    call require(status > 0, 'BiCG nonconvergence')
    call cgs(a,rhs,x,status,maxiter=1,rtol=0.0_dp)
    call require(status > 0, 'CGS nonconvergence')

    ! Zero operator: nonzero RHS causes vanishing denominators.
    a = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]))
    call bicg(a,[1.0_dp,2.0_dp],x,status)
    call require(status == sparse_linalg_breakdown, 'BiCG zero operator breakdown')
    call cgs(a,[1.0_dp,2.0_dp],x,status)
    call require(status == sparse_linalg_breakdown, 'CGS zero operator breakdown')
    call bicg(a,[1.0_dp,2.0_dp],x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'BiCG Jacobi zero diagonal')
    call cgs(a,[1.0_dp,2.0_dp],x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'CGS Jacobi zero diagonal')

    bad = ieee_value(0.0_dp,ieee_quiet_nan)
    a = csr_eye(4)
    call bicg(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'BiCG NaN rhs')
    call cgs(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'CGS NaN rhs')
    call cgs(a,[0.0_dp,0.0_dp,0.0_dp,0.0_dp],x,status,iterations=iterations)
    call require(status == 0 .and. iterations == 0, 'CGS zero rhs')
    print '(a)', 'test_sparse_bicg_cgs: PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! Assertion condition.
        character(len=*), intent(in) :: message !! Failure context.
        if (.not. ok) then
            print '(a)', 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_bicg_cgs
