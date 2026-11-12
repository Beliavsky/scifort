! SPDX-License-Identifier: MIT
program test_sparse_linalg
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : cg, bicgstab, gmres, lsmr, &
        sparse_linalg_invalid_input, sparse_linalg_bad_diagonal
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: x(:), b(:), exact(:), nan_rhs(:)
    real(dp) :: matrix(5,5), rectangular(4,3), rn
    integer :: info, steps, i, n

    matrix = 0.0_dp
    do i = 1, 5
        matrix(i,i) = 4.0_dp
        if (i > 1) matrix(i,i-1) = -1.0_dp
        if (i < 5) matrix(i,i+1) = -1.0_dp
    end do
    a = csr_from_dense(matrix)
    exact = [1.0_dp,-2.0_dp,3.0_dp,-4.0_dp,5.0_dp]
    b = matmul(matrix,exact)
    call cg(a,b,x,info,rtol=1.0e-12_dp,iterations=steps,residual_norm=rn)
    call require(info == 0 .and. steps > 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp, &
        'CG SPD solution')
    call require(rn < 1.0e-10_dp, 'CG true residual')
    call cg(a,b,x,info,rtol=1.0e-12_dp,jacobi=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp, 'Jacobi CG')
    call cg(a,b,x,info,x0=exact)
    call require(info == 0 .and. all(x == exact), 'initial exact solution')
    call cg(a,b,x,info,maxiter=0)
    call require(info > 0, 'CG zero-iteration nonconvergence')
    call cg(a,b,x,info,rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'invalid tolerance')
    call cg(a,b,x,info,x0=[1.0_dp])
    call require(info == sparse_linalg_invalid_input, 'invalid initial guess')

    ! Nonsymmetric and nonsingular, known solution.
    matrix(1,3) = 0.75_dp
    matrix(2,4) = -0.2_dp
    matrix(3,5) = 1.1_dp
    a = csr_from_dense(matrix)
    b = matmul(matrix,exact)
    call bicgstab(a,b,x,info,rtol=1.0e-12_dp,residual_norm=rn)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-9_dp, 'BiCGSTAB nonsymmetric')
    call require(rn < 1.0e-9_dp, 'BiCGSTAB true residual')
    call bicgstab(a,b,x,info,rtol=1.0e-12_dp,jacobi=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-9_dp, 'Jacobi BiCGSTAB')
    call gmres(a,b,x,info,rtol=1.0e-12_dp,restart=5,iterations=steps,residual_norm=rn)
    call require(info == 0 .and. steps <= 5 .and. maxval(abs(x-exact)) < 1.0e-9_dp, &
        'full GMRES nonsymmetric')
    call require(rn < 1.0e-9_dp, 'GMRES true residual')
    call gmres(a,b,x,info,rtol=1.0e-10_dp,restart=2,maxiter=100,jacobi=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-8_dp, &
        'restarted Jacobi GMRES')
    call gmres(a,b,x,info,restart=0)
    call require(info == sparse_linalg_invalid_input, 'invalid restart')
    call gmres(a,b,x,info,maxiter=0)
    call require(info > 0, 'zero GMRES cycles')
    call bicgstab(a,b,x,info,maxiter=0)
    call require(info > 0, 'zero BiCGSTAB iterations')
    call cg(a,0.0_dp*b,x,info)
    call require(info == 0 .and. all(x == 0.0_dp), 'zero RHS')

    ! Explicit stored zeros and a missing Jacobi diagonal are different.
    a = csr_from_coo(2,2,[1,2,1],[1,2,2],[1.0_dp,0.0_dp,2.0_dp])
    call cg(a,[1.0_dp,1.0_dp],x,info,jacobi=.true.)
    call require(info == sparse_linalg_bad_diagonal, 'reject zero diagonal')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]))
    call cg(a,[1.0_dp,2.0_dp,3.0_dp],x,info)
    call require(info == sparse_linalg_invalid_input, 'RHS shape')
    nan_rhs = [ieee_value(0.0_dp,ieee_quiet_nan),2.0_dp]
    call bicgstab(a,nan_rhs,x,info)
    call require(info == sparse_linalg_invalid_input, 'nonfinite RHS')
    call gmres(a,nan_rhs,x,info)
    call require(info == sparse_linalg_invalid_input, 'nonfinite GMRES RHS')
    a = csr_from_dense(reshape([real(dp)::],[0,0]))
    call cg(a,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0, 'CG empty matrix')
    call gmres(a,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0, 'GMRES empty matrix')

    ! Rectangular full-column-rank least squares; reuse existing LSMR kernel.
    rectangular = reshape([ &
        1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
        2.0_dp,0.0_dp,1.0_dp,-1.0_dp, &
        0.0_dp,2.0_dp,1.0_dp,3.0_dp],[4,3])
    a = csr_from_dense(rectangular)
    exact = [2.0_dp,-1.0_dp,3.0_dp]
    b = matmul(rectangular,exact)
    call lsmr(a,b,x,info,itn=steps,normr=rn,atol=1.0e-12_dp,btol=1.0e-12_dp)
    call require(info >= 0 .and. maxval(abs(x-exact)) < 1.0e-8_dp, &
        'LSMR rectangular exact solution')
    call require(rn < 1.0e-8_dp, 'LSMR true residual')
    call lsmr(a,b,x,info,damp=0.2_dp,atol=1.0e-12_dp,btol=1.0e-12_dp)
    call require(info >= 0, 'damped LSMR status')
    call require(maxval(abs(matmul(transpose(rectangular),matmul(rectangular,x)) + &
        0.04_dp*x-matmul(transpose(rectangular),b))) < 1.0e-8_dp, &
        'damped LSMR normal equations')
    call lsmr(a,b,x,info,damp=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'LSMR invalid damping')
    call lsmr(a,b(:3),x,info)
    call require(info == sparse_linalg_invalid_input, 'LSMR invalid RHS')

    ! Consistent singular systems can still have valid Krylov solutions.
    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]))
    call gmres(a,[2.0_dp,0.0_dp],x,info,rtol=1.0e-12_dp)
    call require(info == 0 .and. maxval(abs(x-[2.0_dp,0.0_dp])) < 1.0e-12_dp, &
        'GMRES consistent singular system')
    ! Rank-deficient LSMR returns the minimum-norm least-squares solution.
    a = csr_from_dense(reshape([ &
        1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
        1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,1.0_dp,0.0_dp],[4,3]))
    call lsmr(a,[3.0_dp,6.0_dp,4.0_dp,0.0_dp],x,info, &
        atol=1.0e-12_dp,btol=1.0e-12_dp)
    call require(info >= 0 .and. maxval(abs(x-[1.5_dp,1.5_dp,4.0_dp])) < 1.0e-9_dp, &
        'rank-deficient minimum-norm LSMR')
    ! Breakdown is distinguishable from an exhausted iteration budget.
    a = csr_from_dense(reshape([0.0_dp],[1,1]))
    call cg(a,[1.0_dp],x,info)
    call require(info < 0, 'singular CG breakdown')
    call bicgstab(a,[1.0_dp],x,info)
    call require(info < 0, 'singular BiCGSTAB breakdown')
    call gmres(a,[1.0_dp],x,info)
    call require(info < 0, 'singular GMRES breakdown')
    call lsmr(a,[1.0_dp],x,info,normr=rn)
    call require(info >= 0 .and. size(x) == 1 .and. rn == 1.0_dp, &
        'LSMR zero operator')
    a = csr_matrix()
    call cg(a,[1.0_dp],x,info)
    call require(info == sparse_linalg_invalid_input, 'unallocated CSR invalid')

    ! Tridiagonal large-system check: sparse storage stays O(n).
    n = 800
    a = csr_from_coo(n,n,[(i,i=1,n), (i,i=1,n-1), (i,i=2,n)], &
        [(i,i=1,n), (i+1,i=1,n-1), (i-1,i=2,n)], &
        [ (4.0_dp,i=1,n), (-1.0_dp,i=1,n-1), (-1.0_dp,i=2,n)])
    exact = [(sin(real(i,dp)*0.15_dp),i=1,n)]
    b = csr_matvec(a,exact)
    call cg(a,b,x,info,rtol=1.0e-11_dp,jacobi=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-9_dp, &
        'large sparse Jacobi CG')
    print '(a)', 'test_sparse_linalg: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! True if the regression condition passed.
        character(len=*), intent(in) :: message !! Failure diagnostic.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_linalg
