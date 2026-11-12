! SPDX-License-Identifier: MIT
! Native regression tests for LSQR and zero-fill incomplete LU.
program test_sparse_lsqr_ilu0
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : lsqr, ilu0, ilu0_factor, ilu0_solve, valid_ilu0, &
        bicgstab, gmres, sparse_linalg_invalid_input, sparse_linalg_bad_diagonal
    implicit none
    type(csr_matrix) :: a
    type(ilu0_factor) :: fac
    real(dp), allocatable :: x(:), b(:), true_x(:), vv(:)
    real(dp) :: tall(5,3), square(4,4), rn, r2, gn, est, cond, xn, tol
    integer :: info, step, i, n_big

    tall(:,1) = [1.0_dp,2.0_dp,0.0_dp,4.0_dp,5.0_dp]
    tall(:,2) = [0.0_dp,3.0_dp,1.0_dp,-1.0_dp,2.0_dp]
    tall(:,3) = [2.0_dp,0.0_dp,5.0_dp,1.0_dp,-1.0_dp]
    a = csr_from_dense(tall)
    true_x = [1.0_dp,-3.0_dp,2.0_dp]
    b = matmul(tall,true_x)
    call lsqr(a,b,x,info,itn=step,r1norm=rn,r2norm=r2,anorm=est,acond=cond, &
        arnorm=gn,xnorm=xn,atol=1.0e-13_dp,btol=1.0e-13_dp)
    call check(info >= 0 .and. info <= 2, 'LSQR consistent stopping status')
    call check(step >= 1 .and. maxval(abs(x-true_x)) < 1.0e-10_dp, 'LSQR rectangular exact')
    call check(rn < 1.0e-10_dp .and. abs(rn-r2) < 1.0e-11_dp, 'LSQR residuals')
    call check(est > 0.0_dp .and. cond >= 1.0_dp .and. gn < 1.0e-9_dp, 'LSQR diagnostics')
    call check(abs(xn-norm2(x)) < 1.0e-12_dp, 'LSQR solution norm')

    b = b + [0.08_dp,-0.07_dp,0.02_dp,0.05_dp,0.12_dp]
    call lsqr(a,b,x,info,atol=1.0e-13_dp,btol=1.0e-13_dp,iter_lim=80, &
        r1norm=rn,arnorm=gn)
    call check(info >= 0, 'LSQR inconsistent least squares')
    call check(gn < 1.0e-8_dp, 'LSQR least-squares normal equations')
    call check(rn > 0.01_dp, 'LSQR residual for inconsistent system')
    call lsqr(a,b,x,info,damp=0.2_dp,x0=[1.0_dp,0.0_dp,-1.0_dp], &
        atol=1.0e-13_dp,btol=1.0e-13_dp,iter_lim=80,arnorm=gn,r1norm=rn,r2norm=r2)
    call check(info >= 0 .and. gn < 1.0e-8_dp, 'LSQR shifted damped normal equations')
    call check(r2 >= rn, 'LSQR damped residual norm')

    a = csr_from_dense(reshape([ &
        1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
        1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,1.0_dp,0.0_dp],[4,3]))
    call lsqr(a,[3.0_dp,6.0_dp,4.0_dp,0.0_dp],x,info,atol=1.0e-13_dp,btol=1.0e-13_dp)
    call check(info >= 0 .and. maxval(abs(x-[1.5_dp,1.5_dp,4.0_dp])) < 1.0e-9_dp, &
        'LSQR rank-deficient minimum-norm solution')
    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp], [2,3]))
    call lsqr(a,[2.0_dp,4.0_dp],x,info,atol=1.0e-13_dp,btol=1.0e-13_dp)
    call check(info >= 0 .and. norm2(csr_matvec(a,x)-[2.0_dp,4.0_dp]) < 1.0e-10_dp, &
        'LSQR underdetermined residual')
    call check(maxval(abs(x-[0.0_dp,2.0_dp,2.0_dp])) < 1.0e-9_dp, &
        'LSQR underdetermined minimum norm')

    a = csr_from_dense(reshape([2.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2]))
    call lsqr(a,[2.0_dp,8.0_dp],x,info,atol=1.0e-13_dp,btol=1.0e-13_dp, &
        calc_var=.true.,var=vv)
    call check(info >= 0 .and. maxval(abs(x-[1.0_dp,2.0_dp])) < 1.0e-10_dp, &
        'LSQR diagonal solve')
    call check(maxval(abs(vv-[0.25_dp,0.0625_dp])) < 1.0e-10_dp, &
        'LSQR diagonal variance estimate')
    call lsqr(a,[2.0_dp,8.0_dp],x,info,iter_lim=0, r1norm=rn)
    call check(info == 7 .and. rn > 1.0_dp, 'LSQR zero-iteration cap')
    call lsqr(a,[0.0_dp,0.0_dp],x,info)
    call check(info == 0 .and. all(x == 0.0_dp), 'LSQR zero RHS')
    call lsqr(a,[2.0_dp,8.0_dp],x,info,damp=-0.1_dp)
    call check(info == sparse_linalg_invalid_input, 'LSQR negative damp')
    call lsqr(a,[2.0_dp],x,info)
    call check(info == sparse_linalg_invalid_input, 'LSQR invalid RHS shape')
    call lsqr(a,[2.0_dp,8.0_dp],x,info,x0=[1.0_dp])
    call check(info == sparse_linalg_invalid_input, 'LSQR invalid initial guess shape')
    call lsqr(a,[ieee_value(0.0_dp,ieee_quiet_nan),8.0_dp],x,info)
    call check(info == sparse_linalg_invalid_input, 'LSQR invalid RHS nonfinite')
    a = csr_from_dense(reshape([real(dp) :: ],[0,3]))
    call lsqr(a,[real(dp) :: ],x,info)
    call check(info == 0 .and. size(x) == 3, 'LSQR empty rows')
    a = csr_from_dense(reshape([real(dp) :: ],[3,0]))
    call lsqr(a,[1.0_dp,2.0_dp,3.0_dp],x,info,r1norm=rn)
    call check(info == 0 .and. size(x) == 0 .and. rn > 0.0_dp, 'LSQR empty columns')
    a = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]))
    call lsqr(a,[1.0_dp,0.0_dp],x,info,r1norm=rn)
    call check(info == 0 .and. all(x == 0.0_dp) .and. rn == 1.0_dp, 'LSQR zero matrix')

    square = 0.0_dp
    do i = 1, 4
        square(i,i) = 4.0_dp
        if (i > 1) square(i,i-1) = -1.5_dp
        if (i < 4) square(i,i+1) = -0.5_dp
    end do
    a = csr_from_dense(square)
    true_x = [1.0_dp,2.0_dp,-1.0_dp,3.0_dp]
    b = matmul(square,true_x)
    call ilu0(a,fac,info)
    call check(info == 0 .and. valid_ilu0(fac), 'ILU0 tridiagonal valid')
    call ilu0_solve(fac,b,x,info)
    call check(info == 0 .and. maxval(abs(x-true_x)) < 1.0e-12_dp, &
        'ILU0 exact tridiagonal factor solve')
    call bicgstab(a,b,x,info,ilu=fac,iterations=step,rtol=1.0e-13_dp)
    call check(info == 0 .and. step <= 2 .and. maxval(abs(x-true_x)) < 1.0e-10_dp, &
        'ILU0 BiCGSTAB convergence')
    call gmres(a,b,x,info,ilu=fac,iterations=step,rtol=1.0e-13_dp)
    call check(info == 0 .and. step <= 2 .and. maxval(abs(x-true_x)) < 1.0e-10_dp, &
        'ILU0 GMRES convergence')
    call gmres(a,b,x,info,ilu=fac,jacobi=.true.)
    call check(info == sparse_linalg_invalid_input, 'ILU and Jacobi are mutually exclusive')
    call ilu0(a,fac,info,diagonal_shift=-1.0_dp)
    call check(info == sparse_linalg_invalid_input, 'ILU0 invalid diagonal shift')
    call ilu0(a,fac,info,diagonal_shift=0.1_dp)
    call check(info == 0 .and. valid_ilu0(fac), 'ILU0 shifted factorization')
    call ilu0_solve(fac,[1.0_dp],x,info)
    call check(info == sparse_linalg_invalid_input, 'ILU0 invalid RHS')
    a = csr_from_dense(reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp],[2,2]))
    call ilu0(a,fac,info)
    call check(info == sparse_linalg_bad_diagonal, 'ILU0 missing diagonal')
    a = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))
    call ilu0(a,fac,info)
    call check(info == sparse_linalg_bad_diagonal, 'ILU0 zero pivot')
    a = csr_from_dense(reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp],[2,2]))
    call ilu0(a,fac,info)
    call check(info == sparse_linalg_bad_diagonal, 'ILU0 zero Schur complement pivot')
    a = csr_from_coo(0,0,[integer :: ],[integer :: ],[real(dp) :: ])
    call ilu0(a,fac,info)
    call check(info == 0 .and. valid_ilu0(fac), 'ILU0 empty factorization')
    call ilu0_solve(fac,[real(dp) :: ],x,info)
    call check(info == 0 .and. size(x) == 0, 'ILU0 empty solve')
    ! ILU(0) on a large sparse nonsymmetric tridiagonal system has exact
    ! bandwidth-preserving factors, with no dense intermediate allocation.
    n_big = 1500
    a = csr_from_coo(n_big,n_big,[(i,i=1,n_big), (i,i=1,n_big-1), (i,i=2,n_big)], &
        [(i,i=1,n_big), (i+1,i=1,n_big-1), (i-1,i=2,n_big)], &
        [(3.0_dp,i=1,n_big),(-0.2_dp,i=1,n_big-1),(-0.1_dp,i=2,n_big)])
    true_x = [(cos(0.05_dp*real(i,dp)),i=1,n_big)]
    b = csr_matvec(a,true_x)
    call ilu0(a,fac,info)
    call check(info == 0 .and. size(fac%lu%data) == 3*n_big-2, &
        'ILU0 large sparse O(nnz) pattern')
    call ilu0_solve(fac,b,x,info)
    call check(info == 0 .and. maxval(abs(x-true_x)) < 1.0e-11_dp, &
        'ILU0 large sparse solve')
    call bicgstab(a,b,x,info,ilu=fac,rtol=1.0e-12_dp,iterations=step)
    call check(info == 0 .and. step <= 2 .and. maxval(abs(x-true_x)) < 1.0e-10_dp, &
        'ILU0 large sparse BiCGSTAB')
    print '(a)', 'test_sparse_lsqr_ilu0: PASS'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition !! Expected regression condition.
        character(len=*), intent(in) :: description !! Error message if check fails.
        if (.not. condition) then
            print '(a)', 'FAIL: '//description
            error stop 1
        end if
    end subroutine check
end program test_sparse_lsqr_ilu0
