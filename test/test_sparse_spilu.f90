! SPDX-License-Identifier: MIT
! Incomplete LU with threshold dropping, per-row fill, and existing Krylov solvers.
program test_sparse_spilu
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use scifort_kinds, only: dp
    use scifort_sparse, only: csr_matrix, csc_matrix, csr_from_dense, csr_from_coo, &
        csr_to_csc, csr_matvec, csr_to_dense
    use scifort_sparse_linalg, only: ilu0_factor, ilu0, spilu, spilu_solve, &
        valid_ilu0, bicgstab, gmres, sparse_linalg_bad_diagonal, &
        sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a, empty, missing
    type(csc_matrix) :: ac
    type(ilu0_factor) :: factor, factor2, ilu
    real(dp), allocatable :: x(:), y(:), b(:), dense(:,:), exact(:), block(:,:), xblock(:,:)
    real(dp) :: nanval
    integer :: info, n, i, k, nnz0

    dense = reshape([ &
        6.0_dp,-1.0_dp,0.0_dp,-0.5_dp, &
        -1.0_dp,6.0_dp,-1.0_dp,0.0_dp, &
        0.0_dp,-1.0_dp,6.0_dp,-1.0_dp, &
        -0.5_dp,0.0_dp,-1.0_dp,6.0_dp],[4,4])
    a = csr_from_dense(dense)
    exact = [1.0_dp,2.0_dp,-1.0_dp,3.0_dp]
    b = csr_matvec(a,exact)
    call ilu0(a,ilu,info)
    call require(info == 0,'baseline ILU0 succeeds')
    nnz0 = size(ilu%lu%data)
    call spilu(a,factor,info)
    call require(info == 0 .and. valid_ilu0(factor),'default ILUT')
    call spilu(a,factor,info,drop_tol=0.0_dp,fill_factor=10.0_dp)
    call require(info == 0 .and. valid_ilu0(factor),'exact-fill ILUT succeeds')
    call require(size(factor%lu%data) > nnz0,'extra fill is retained')
    call spilu_solve(factor,b,x,info)
    call require(info == 0,'ILUT apply succeeds')
    call require(maxval(abs(x-exact)) < 1.0e-11_dp,'unlimited ILUT exactly solves no-pivot system')
    call spilu_solve(factor,[1.0_dp],y,info)
    call require(info == sparse_linalg_invalid_input,'ILUT RHS dimension check')
    ac = csr_to_csc(a)
    call spilu(ac,factor2,info,drop_tol=0.0_dp,fill_factor=10.0_dp)
    call require(info == 0,'CSC adapter')
    call spilu_solve(factor2,b,y,info)
    call require(info == 0 .and. maxval(abs(y-x)) < 1.0e-12_dp,'CSC matching solve')
    block = spread(b,2,2)
    block(:,2) = 2.0_dp*b
    call spilu_solve(factor,block,xblock,info)
    call require(info == 0 .and. maxval(abs(xblock(:,2)-2.0_dp*exact)) < 1.0e-11_dp, &
        'multiple RHS solve')

    call spilu(a,factor,info,drop_tol=0.3_dp,fill_factor=1.0_dp)
    call require(info == 0 .and. valid_ilu0(factor),'dropped ILUT succeeds')
    call require(size(factor%lu%data) <= nnz0,'low fill has limited entries')
    call bicgstab(a,b,x,info,rtol=1.0e-12_dp,atol=0.0_dp,ilu=factor)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp,'ILUT BiCGSTAB')
    call gmres(a,b,x,info,rtol=1.0e-12_dp,atol=0.0_dp,ilu=factor)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp,'ILUT GMRES')

    call spilu(a,factor,info,drop_tol=-0.1_dp)
    call require(info == sparse_linalg_invalid_input,'negative tolerance rejected')
    call spilu(a,factor,info,drop_tol=1.1_dp)
    call require(info == sparse_linalg_invalid_input,'too-large tolerance rejected')
    call spilu(a,factor,info,fill_factor=0.5_dp)
    call require(info == sparse_linalg_invalid_input,'small fill rejected')
    call spilu(a,factor,info,diagonal_shift=-1.0_dp)
    call require(info == sparse_linalg_invalid_input,'negative shift rejected')
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    call spilu(a,factor,info,drop_tol=nanval)
    call require(info == sparse_linalg_invalid_input,'non-finite tolerance rejected')

    missing = csr_from_dense(reshape([0.0_dp,0.0_dp,0.0_dp,3.0_dp],[2,2]))
    call spilu(missing,factor,info)
    call require(info == sparse_linalg_bad_diagonal,'missing pivot rejected')
    call spilu(missing,factor,info,diagonal_shift=1.0_dp)
    call require(info == 0,'diagonal perturbation restores pivot')
    empty = csr_from_coo(0,0,[integer::],[integer::],[real(dp)::])
    call spilu(empty,factor,info)
    call require(info == 0 .and. valid_ilu0(factor),'empty factorization')
    call spilu_solve(factor,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0,'empty solve')

    n = 1200
    a = csr_from_coo(n,n,[(i,i=1,n),(i,i=2,n),(i,i=1,n-1)], &
        [(i,i=1,n),(i-1,i=2,n),(i+1,i=1,n-1)], &
        [(4.0_dp,i=1,n),(-1.0_dp,i=2,n),(-1.0_dp,i=1,n-1)])
    exact = [(cos(real(i,dp)*0.06_dp),i=1,n)]
    b = csr_matvec(a,exact)
    call spilu(a,factor,info,drop_tol=1.0e-4_dp,fill_factor=3.0_dp)
    call require(info == 0 .and. valid_ilu0(factor),'large sparse ILUT')
    call spilu_solve(factor,b,x,info)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-10_dp, &
        'large sparse ILUT accurate for tridiagonal matrix')
    print '(a)', 'test_sparse_spilu: PASS'
contains
    subroutine require(ok, desc)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: desc
        if (.not. ok) then
            print '(a)', 'FAIL: '//desc
            error stop 1
        end if
    end subroutine require
end program test_sparse_spilu
