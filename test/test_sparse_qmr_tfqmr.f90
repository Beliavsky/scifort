! SPDX-License-Identifier: MIT
program test_sparse_qmr_tfqmr
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_eye
    use scifort_sparse_linalg, only : qmr, tfqmr, ilu0, ilu0_factor, &
        sparse_linalg_success, sparse_linalg_invalid_input, &
        sparse_linalg_bad_diagonal, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a, empty
    type(ilu0_factor) :: factor
    real(dp) :: mat(4,4), rhs(4), truth(4), res, bad, zrhs(2)
    real(dp) :: indef(3,3), rhs3(3), truth3(3)
    real(dp), allocatable :: x(:)
    integer :: status, steps

    mat = 0.0_dp
    mat(1,:) = [6.0_dp,2.0_dp,0.0_dp,-1.0_dp]
    mat(2,:) = [1.0_dp,5.0_dp,2.0_dp,0.0_dp]
    mat(3,:) = [0.0_dp,-1.0_dp,4.0_dp,1.0_dp]
    mat(4,:) = [1.0_dp,0.0_dp,2.0_dp,7.0_dp]
    a = csr_from_dense(mat)
    truth = [2.0_dp,-1.0_dp,3.0_dp,0.5_dp]
    rhs = matmul(mat,truth)
    call qmr(a,rhs,x,status,rtol=1.0e-12_dp,atol=1.0e-13_dp, &
        iterations=steps,residual_norm=res)
    call require(status == sparse_linalg_success, 'QMR convergence')
    call require(steps >= 1 .and. steps <= 10, 'QMR steps')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp, 'QMR solution')
    call require(res < 1.0e-11_dp, 'QMR true residual')
    call qmr(a,rhs,x,status,jacobi=.true.,rtol=1.0e-12_dp)
    call require(status == sparse_linalg_success, 'QMR Jacobi convergence')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp, 'QMR Jacobi solution')
    call qmr(a,rhs,x,status,x0=truth,iterations=steps,residual_norm=res)
    call require(status == 0 .and. steps == 0 .and. res < 1.0e-13_dp, 'QMR exact x0')

    call tfqmr(a,rhs,x,status,rtol=1.0e-12_dp,atol=1.0e-13_dp, &
        iterations=steps,residual_norm=res)
    call require(status == sparse_linalg_success, 'TFQMR convergence')
    call require(steps >= 1 .and. steps <= 20, 'TFQMR steps')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp, 'TFQMR solution')
    call require(res < 1.0e-11_dp, 'TFQMR true residual')
    call tfqmr(a,rhs,x,status,jacobi=.true.,rtol=1.0e-12_dp)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'TFQMR Jacobi')
    call ilu0(a,factor,status)
    call require(status == 0, 'ILU0 construction')
    call tfqmr(a,rhs,x,status,rtol=1.0e-12_dp,ilu=factor)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'TFQMR ILU0')
    call tfqmr(a,rhs,x,status,x0=truth,iterations=steps,residual_norm=res)
    call require(status == 0 .and. steps == 0 .and. res < 1.0e-13_dp, 'TFQMR exact x0')
    call tfqmr(a,rhs,x,status,ilu=factor,jacobi=.true.)
    call require(status == sparse_linalg_invalid_input, 'TFQMR two preconditioners')

    ! Indefinite and strongly nonsymmetric system from the SciPy QMR example.
    indef = reshape([3.0_dp,1.0_dp,0.0_dp, &
                     2.0_dp,-1.0_dp,5.0_dp, &
                     0.0_dp,0.0_dp,1.0_dp], [3,3])
    a = csr_from_dense(indef)
    truth3 = [1.5_dp,-2.0_dp,0.5_dp]
    rhs3 = matmul(indef,truth3)
    call qmr(a,rhs3,x,status,rtol=1.0e-12_dp)
    call require(status == 0 .and. maxval(abs(x-truth3)) < 1.0e-10_dp, &
        'QMR indefinite')
    call qmr(a,rhs3,x,status,rtol=1.0e-12_dp,jacobi=.true.)
    call require(status == 0 .and. maxval(abs(x-truth3)) < 1.0e-10_dp, &
        'QMR indefinite Jacobi')
    call tfqmr(a,rhs3,x,status,rtol=1.0e-12_dp)
    call require(status == 0 .and. maxval(abs(x-truth3)) < 1.0e-10_dp, &
        'TFQMR indefinite')
    call tfqmr(a,rhs3,x,status,rtol=1.0e-12_dp,jacobi=.true.)
    call require(status == 0 .and. maxval(abs(x-truth3)) < 1.0e-10_dp, &
        'TFQMR indefinite Jacobi')

    a = csr_from_dense(mat)
    empty = csr_eye(0)
    call qmr(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'QMR empty')
    call tfqmr(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'TFQMR empty')
    call qmr(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'QMR bad rhs')
    call tfqmr(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'TFQMR bad rhs')
    call qmr(a,rhs,x,status,maxiter=-1)
    call require(status == sparse_linalg_invalid_input, 'QMR bad cap')
    call tfqmr(a,rhs,x,status,rtol=-1.0_dp)
    call require(status == sparse_linalg_invalid_input, 'TFQMR bad rtol')
    call qmr(a,rhs,x,status,maxiter=0,iterations=steps)
    call require(status > 0 .and. steps == 0, 'QMR zero cap')
    call tfqmr(a,rhs,x,status,maxiter=0,iterations=steps)
    call require(status > 0 .and. steps == 0, 'TFQMR zero cap')
    call qmr(a,rhs,x,status,maxiter=1,rtol=0.0_dp)
    call require(status > 0, 'QMR nonconvergence')
    call tfqmr(a,rhs,x,status,maxiter=1,rtol=0.0_dp)
    call require(status > 0, 'TFQMR nonconvergence')

    a = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]))
    zrhs = [1.0_dp,2.0_dp]
    call qmr(a,zrhs,x,status)
    call require(status == sparse_linalg_breakdown, 'QMR zero operator breakdown')
    call tfqmr(a,zrhs,x,status)
    call require(status == sparse_linalg_breakdown, 'TFQMR zero operator breakdown')
    call qmr(a,zrhs,x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'QMR bad Jacobi')
    call tfqmr(a,zrhs,x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'TFQMR bad Jacobi')

    a = csr_eye(4)
    bad = ieee_value(0.0_dp,ieee_quiet_nan)
    call qmr(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'QMR NaN b')
    call tfqmr(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'TFQMR NaN b')
    call qmr(a,[0.0_dp,0.0_dp,0.0_dp,0.0_dp],x,status,iterations=steps)
    call require(status == 0 .and. steps == 0, 'QMR zero rhs')
    call tfqmr(a,[0.0_dp,0.0_dp,0.0_dp,0.0_dp],x,status,iterations=steps)
    call require(status == 0 .and. steps == 0, 'TFQMR zero rhs')
    print '(a)', 'test_sparse_qmr_tfqmr: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Assertion condition.
        character(len=*), intent(in) :: message !! Failure context.
        if (.not. ok) then
            print '(a)', 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_qmr_tfqmr
