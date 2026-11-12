! SPDX-License-Identifier: MIT
program test_sparse_direct
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_to_csc, csr_from_dense, &
        csr_from_coo, csr_to_dense, csr_matvec
    use scifort_sparse_linalg, only : spsolve_triangular, splu, splu_solve, &
        splu_factor, valid_splu, spsolve, sparse_linalg_invalid_input, &
        sparse_linalg_singular, sparse_linalg_bad_diagonal
    implicit none
    type(csr_matrix) :: a
    type(csc_matrix) :: ac
    type(splu_factor) :: factor
    real(dp), allocatable :: x(:), xblock(:,:), b(:), mat(:,:), exact(:), rhs(:,:)
    real(dp) :: matrix(4,4), lower(4,4), upper(4,4)
    integer :: info, n, i

    matrix = reshape([ &
        0.0_dp, 2.0_dp, 0.0_dp, -1.0_dp, &
        1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, &
        0.0_dp, -1.0_dp, 5.0_dp, 0.0_dp, &
        2.0_dp, 0.0_dp, 1.0_dp, 4.0_dp ],[4,4])
    a = csr_from_dense(matrix)
    exact = [2.0_dp,-3.0_dp,1.0_dp,4.0_dp]
    b = matmul(matrix,exact)
    call spsolve(a,b,x,info)
    call require(info == 0,'pivoted direct solve')
    call require(maxval(abs(x-exact)) < 1.0e-12_dp,'pivoted direct accuracy')
    call splu(a,factor,info)
    call require(info == 0 .and. valid_splu(factor),'valid pivoted factors')
    ac = csr_to_csc(a)
    call spsolve(ac,b,x,info)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'CSC direct solve')
    call splu(ac,factor,info)
    call require(info == 0 .and. valid_splu(factor),'CSC LU factor')
    call splu_solve(factor,b,x,info)
    call require(info == 0,'solve existing LU')
    call require(maxval(abs(x-exact)) < 1.0e-12_dp,'factor solve accuracy')
    rhs = reshape([b,2.0_dp*b,3.0_dp*b],[4,3])
    call spsolve(a,rhs,xblock,info)
    call require(info == 0,'multiple RHS direct solve')
    call require(maxval(abs(matmul(matrix,xblock)-rhs)) < 1.0e-11_dp,'multiple RHS residual')
    call splu_solve(factor,rhs,xblock,info)
    call require(info == 0,'multiple RHS factor solve')
    call require(maxval(abs(matmul(matrix,xblock)-rhs)) < 1.0e-11_dp,'reused factor accuracy')
    mat = csr_to_dense(factor%lu)
    block
        real(dp) :: l(4,4), u(4,4), product(4,4)
        integer :: j
        l = 0.0_dp
        u = 0.0_dp
        do i=1,4
            l(i,i) = 1.0_dp
            do j=1,4
                if (j < i) l(i,j) = mat(i,j)
                if (j >= i) u(i,j) = mat(i,j)
            end do
        end do
        product = matmul(l,u)
        do i=1,4
            call require(maxval(abs(product(i,:)-matrix(factor%row_perm(i),:))) &
                < 1.0e-12_dp,'sparse LU decomposition')
        end do
    end block

    lower = 0.0_dp
    do i=1,4
        lower(i,i) = real(i+1,dp)
        if (i>1) lower(i,i-1) = -0.5_dp
        if (i>2) lower(i,i-2) = 0.3_dp
    end do
    upper = transpose(lower)
    a = csr_from_dense(lower)
    b = matmul(lower,exact)
    call spsolve_triangular(a,b,x,info)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'lower triangular')
    ac = csr_to_csc(a)
    call spsolve_triangular(ac,b,x,info)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'lower CSC triangular')
    a = csr_from_dense(upper)
    b = matmul(upper,exact)
    call spsolve_triangular(a,b,x,info,lower=.false.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'upper triangular')
    ac = csr_to_csc(a)
    call spsolve_triangular(ac,b,x,info,lower=.false.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'upper CSC triangular')
    lower = 0.0_dp
    do i=2,4
        lower(i,i-1) = 0.2_dp
    end do
    a = csr_from_dense(lower)
    b = matmul(lower+eye(4),exact)
    call spsolve_triangular(a,b,x,info,unit_diagonal=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'implicit unit diagonal')
    ac = csr_to_csc(a)
    call spsolve_triangular(ac,b,x,info,unit_diagonal=.true.)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'CSC unit diagonal')
    call spsolve_triangular(a,b,x,info)
    call require(info == sparse_linalg_bad_diagonal,'missing diagonal rejected')

    a = csr_from_coo(3,3,[1,1,1,2,2,3,3], [1,1,2,2,3,2,3], &
        [1.0_dp,2.0_dp,0.0_dp,2.0_dp,1.0_dp,3.0_dp,4.0_dp])
    exact = [3.0_dp,-2.0_dp,5.0_dp]
    b = csr_matvec(a,exact)
    call spsolve(a,b,x,info)
    call require(info == 0 .and. maxval(abs(x-exact)) < 1.0e-12_dp,'duplicates and stored zeros')
    call splu(a,factor,info,pivot_tol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input,'invalid pivot tolerance')
    call splu(a,factor,info,pivot_tol=100.0_dp)
    call require(info == sparse_linalg_singular,'large minimum pivot rejected')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2]))
    call splu(a,factor,info)
    call require(info == sparse_linalg_singular,'singular factorization rejected')
    call spsolve(a,[1.0_dp,2.0_dp],x,info)
    call require(info == sparse_linalg_singular,'singular direct solve rejected')
    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))
    call spsolve(a,[1.0_dp],x,info)
    call require(info == sparse_linalg_invalid_input,'RHS shape rejected')
    call spsolve(a,[1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)],x,info)
    call require(info == sparse_linalg_invalid_input,'NaN RHS rejected')
    a = csr_from_dense(reshape([real(dp)::],[0,0]))
    call splu(a,factor,info)
    call require(info == 0 .and. valid_splu(factor),'empty LU')
    call spsolve(a,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0,'empty direct solve')
    call spsolve_triangular(a,[real(dp)::],x,info)
    call require(info == 0 .and. size(x) == 0,'empty triangular solve')

    n = 1500
    a = csr_from_coo(n,n,[(i,i=1,n),(i,i=2,n),(i,i=1,n-1)], &
        [(i,i=1,n),(i-1,i=2,n),(i+1,i=1,n-1)], &
        [(4.0_dp,i=1,n),(-1.0_dp,i=2,n),(-1.0_dp,i=1,n-1)])
    exact = [(sin(real(i,dp)/7.0_dp),i=1,n)]
    b = csr_matvec(a,exact)
    call spsolve(a,b,x,info)
    call require(info == 0,'large sparse LU')
    call require(maxval(abs(x-exact)) < 1.0e-11_dp,'large sparse LU accuracy')
    print '(a)', 'test_sparse_direct: PASS'
contains
    function eye(n) result(m)
        integer, intent(in) :: n
        real(dp) :: m(n,n)
        integer :: i
        m = 0.0_dp
        do i=1,n
            m(i,i) = 1.0_dp
        end do
    end function eye
    subroutine require(ok, desc)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: desc
        if (.not. ok) then
            print '(a)', 'FAIL: '//desc
            error stop 1
        end if
    end subroutine require
end program test_sparse_direct
