! SPDX-License-Identifier: MIT
program test_sparse_recycling
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_eye
    use scifort_sparse_linalg, only : lgmres, gcrotmk, ilu0, ilu0_factor, &
        sparse_linalg_success, sparse_linalg_invalid_input, &
        sparse_linalg_bad_diagonal, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a, empty
    type(ilu0_factor) :: factor
    real(dp) :: mat(4,4), rhs(4), truth(4), shifted(4,4), res, bad
    real(dp) :: indef(5,5), b5(5), t5(5)
    real(dp), allocatable :: x(:), outer(:,:), recycled(:,:)
    integer :: status, steps

    mat = 0.0_dp
    mat(1,:) = [6.0_dp,2.0_dp,0.0_dp,-1.0_dp]
    mat(2,:) = [1.0_dp,5.0_dp,2.0_dp,0.0_dp]
    mat(3,:) = [0.0_dp,-1.0_dp,4.0_dp,1.0_dp]
    mat(4,:) = [1.0_dp,0.0_dp,2.0_dp,7.0_dp]
    a = csr_from_dense(mat)
    truth = [2.0_dp,-1.0_dp,3.0_dp,0.5_dp]
    rhs = matmul(mat,truth)
    call lgmres(a,rhs,x,status,rtol=1.0e-12_dp,inner_m=2,outer_k=2, &
        iterations=steps,residual_norm=res,outer_v=outer)
    call require(status == sparse_linalg_success, 'LGMRES converged')
    call require(steps >= 1 .and. steps < 30, 'LGMRES Arnoldi iterations')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp .and. res < 1.0e-11_dp, 'LGMRES solution')
    call require(allocated(outer) .and. size(outer,2) > 0 .and. size(outer,2) <= 2, 'LGMRES cache')
    call lgmres(a,rhs,x,status,rtol=1.0e-12_dp,inner_m=2,outer_k=2, &
        jacobi=.true.,outer_v=outer)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'LGMRES Jacobi with cache')
    call ilu0(a,factor,status)
    call require(status == 0, 'ILU0 construction')
    call lgmres(a,rhs,x,status,rtol=1.0e-12_dp,inner_m=1,outer_k=2,ilu=factor)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'LGMRES ILU0')
    call lgmres(a,rhs,x,status,rtol=1.0e-12_dp,inner_m=2,outer_k=0)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'LGMRES unaugmented')

    call gcrotmk(a,rhs,x,status,rtol=1.0e-12_dp,m=2,k=2, &
        iterations=steps,residual_norm=res,recycle_u=recycled)
    call require(status == 0 .and. steps > 0 .and. steps < 30, 'GCROT convergence')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp .and. res < 1.0e-11_dp, 'GCROT solution')
    call require(allocated(recycled) .and. size(recycled,2) > 0 .and. &
        size(recycled,2) <= 2, 'GCROT cache')
    call gcrotmk(a,rhs,x,status,rtol=1.0e-12_dp,m=1,k=2, &
        jacobi=.true.,recycle_u=recycled)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'GCROT Jacobi with cache')
    call gcrotmk(a,rhs,x,status,rtol=1.0e-12_dp,m=2,k=2,ilu=factor)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'GCROT ILU0')
    call gcrotmk(a,rhs,x,status,rtol=1.0e-12_dp,m=2,k=0)
    call require(status == 0 .and. maxval(abs(x-truth)) < 1.0e-11_dp, 'GCROT no recycling')

    ! Indefinite, substantially nonsymmetric nonnormal system. This also
    ! exercises the Arnoldi projection when the CG-type solvers are unsuitable.
    indef = 0.0_dp
    indef(1,:) = [-3.0_dp,2.0_dp,0.0_dp,0.0_dp,1.0_dp]
    indef(2,:) = [3.0_dp,2.0_dp,-2.0_dp,0.0_dp,0.0_dp]
    indef(3,:) = [0.0_dp,-1.0_dp,1.0_dp,2.0_dp,0.0_dp]
    indef(4,:) = [1.0_dp,0.0_dp,3.0_dp,-2.0_dp,3.0_dp]
    indef(5,:) = [0.0_dp,-1.0_dp,0.0_dp,-3.0_dp,4.0_dp]
    a = csr_from_dense(indef)
    t5 = [1.0_dp,-2.0_dp,0.5_dp,0.25_dp,1.75_dp]
    b5 = matmul(indef,t5)
    call lgmres(a,b5,x,status,rtol=1.0e-12_dp,inner_m=3,outer_k=2)
    call require(status == 0 .and. maxval(abs(x-t5)) < 1.0e-10_dp, 'LGMRES indefinite')
    call gcrotmk(a,b5,x,status,rtol=1.0e-12_dp,m=3,k=2)
    call require(status == 0 .and. maxval(abs(x-t5)) < 1.0e-10_dp, 'GCROT indefinite')
    call gcrotmk(a,b5,x,status,rtol=1.0e-12_dp,m=2,k=0,maxiter=4)
    call require(status > 0, 'GCROT indefinite capped stagnation')
    a = csr_from_dense(mat)

    ! An externally supplied exact correction must be used by LGMRES.
    outer = reshape(truth,[4,1])
    call lgmres(a,rhs,x,status,rtol=1.0e-12_dp,inner_m=1,outer_k=1, &
        iterations=steps,outer_v=outer)
    call require(status == 0 .and. steps == 2, 'LGMRES uses outer correction')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp, 'LGMRES outer correctness')
    ! Rebuild image basis when the operator changes, preserving the U space.
    recycled = reshape(truth,[4,1])
    shifted = mat
    shifted(1,1) = shifted(1,1)+1.0_dp
    a = csr_from_dense(shifted)
    rhs = matmul(shifted,truth)
    call gcrotmk(a,rhs,x,status,rtol=1.0e-12_dp,m=1,k=1, &
        iterations=steps,recycle_u=recycled)
    call require(status == 0 .and. steps == 0, 'GCROT uses externally supplied U for changed A')
    call require(maxval(abs(x-truth)) < 1.0e-11_dp, 'GCROT changed operator')
    a = csr_from_dense(mat)
    rhs = matmul(mat,truth)

    call lgmres(a,rhs,x,status,x0=truth,iterations=steps,residual_norm=res)
    call require(status == 0 .and. steps == 0 .and. res < 1.0e-13_dp, 'LGMRES exact x0')
    call gcrotmk(a,rhs,x,status,x0=truth,iterations=steps,residual_norm=res)
    call require(status == 0 .and. steps == 0 .and. res < 1.0e-13_dp, 'GCROT exact x0')
    empty = csr_eye(0)
    call lgmres(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'LGMRES empty')
    call gcrotmk(empty,[real(dp)::],x,status)
    call require(status == 0 .and. size(x) == 0, 'GCROT empty')
    call lgmres(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'LGMRES wrong rhs size')
    call gcrotmk(a,[real(dp)::],x,status)
    call require(status == sparse_linalg_invalid_input, 'GCROT wrong rhs size')
    call lgmres(a,rhs,x,status,maxiter=-1)
    call require(status == sparse_linalg_invalid_input, 'LGMRES negative maxiter')
    call gcrotmk(a,rhs,x,status,maxiter=-1)
    call require(status == sparse_linalg_invalid_input, 'GCROT negative maxiter')
    call lgmres(a,rhs,x,status,inner_m=0)
    call require(status == sparse_linalg_invalid_input, 'LGMRES zero inner_m')
    call gcrotmk(a,rhs,x,status,m=0)
    call require(status == sparse_linalg_invalid_input, 'GCROT zero m')
    call lgmres(a,rhs,x,status,outer_k=-1)
    call require(status == sparse_linalg_invalid_input, 'LGMRES negative outer_k')
    call gcrotmk(a,rhs,x,status,k=-1)
    call require(status == sparse_linalg_invalid_input, 'GCROT negative k')
    call lgmres(a,rhs,x,status,ilu=factor,jacobi=.true.)
    call require(status == sparse_linalg_invalid_input, 'LGMRES conflicting preconditioners')
    call gcrotmk(a,rhs,x,status,ilu=factor,jacobi=.true.)
    call require(status == sparse_linalg_invalid_input, 'GCROT conflicting preconditioners')
    call lgmres(a,rhs,x,status,maxiter=0,iterations=steps)
    call require(status > 0 .and. steps == 0, 'LGMRES cap zero')
    call gcrotmk(a,rhs,x,status,maxiter=0,iterations=steps)
    call require(status > 0 .and. steps == 0, 'GCROT cap zero')
    call lgmres(a,rhs,x,status,rtol=0.0_dp,atol=0.0_dp,maxiter=1,inner_m=1,outer_k=0)
    call require(status > 0, 'LGMRES nonconvergence')
    call gcrotmk(a,rhs,x,status,rtol=0.0_dp,atol=0.0_dp,maxiter=1,m=1,k=0)
    call require(status > 0, 'GCROT nonconvergence')

    a = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]))
    call lgmres(a,[1.0_dp,2.0_dp],x,status)
    call require(status == sparse_linalg_breakdown, 'LGMRES zero matrix')
    call gcrotmk(a,[1.0_dp,2.0_dp],x,status)
    call require(status == sparse_linalg_breakdown, 'GCROT zero matrix')
    call lgmres(a,[1.0_dp,2.0_dp],x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'LGMRES bad Jacobi')
    call gcrotmk(a,[1.0_dp,2.0_dp],x,status,jacobi=.true.)
    call require(status == sparse_linalg_bad_diagonal, 'GCROT bad Jacobi')

    a = csr_eye(4)
    bad = ieee_value(0.0_dp,ieee_quiet_nan)
    call lgmres(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'LGMRES NaN RHS')
    call gcrotmk(a,[bad,1.0_dp,2.0_dp,3.0_dp],x,status)
    call require(status == sparse_linalg_invalid_input, 'GCROT NaN RHS')
    outer = reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp],[4,1])
    outer(1,1) = bad
    call lgmres(a,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],x,status,outer_v=outer)
    call require(status == sparse_linalg_invalid_input, 'LGMRES NaN cache')
    recycled = outer
    call gcrotmk(a,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],x,status,recycle_u=recycled)
    call require(status == sparse_linalg_invalid_input, 'GCROT NaN cache')
    call lgmres(a,[0.0_dp,0.0_dp,0.0_dp,0.0_dp],x,status,iterations=steps)
    call require(status == 0 .and. steps == 0, 'LGMRES zero b')
    call gcrotmk(a,[0.0_dp,0.0_dp,0.0_dp,0.0_dp],x,status,iterations=steps)
    call require(status == 0 .and. steps == 0, 'GCROT zero b')
    print '(a)', 'test_sparse_recycling: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Assertion condition.
        character(len=*), intent(in) :: message !! Diagnostic failure context.
        if (.not. ok) then
            print '(a)', 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_recycling
