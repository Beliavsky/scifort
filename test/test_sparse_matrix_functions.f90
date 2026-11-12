! SPDX-License-Identifier: MIT
program test_sparse_matrix_functions
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_dense, csr_to_csc, &
        csc_to_csr, csr_eye, csr_from_coo, csr_to_dense, valid_csr, valid_csc
    use scifort_sparse_linalg, only : matrix_power, inv, sparse_linalg_success, &
        sparse_linalg_invalid_input, sparse_linalg_singular, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a, b, p, r, empty, bad
    type(csc_matrix) :: ac, pc, ic
    real(dp), allocatable :: d(:,:), dense(:,:), diagonal(:), inverse(:,:)
    real(dp) :: q(3,3), z(2,2), nan_value
    integer :: info, i, k, n

    q = reshape([0.0_dp, 2.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 4.0_dp], [3,3])
    a = csr_from_dense(q)
    ac = csr_to_csc(a)
    do k = 0,13
        p = matrix_power(a,k,info)
        call require(info == sparse_linalg_success .and. valid_csr(p), 'CSR power valid')
        pc = matrix_power(ac,k,info)
        call require(info == sparse_linalg_success .and. valid_csc(pc), 'CSC power valid')
        r = csc_to_csr(pc)
        call require(maxval(abs(csr_to_dense(p)-csr_to_dense(r))) < 1.0e-11_dp, 'CSC power')
        dense = csr_to_dense(p)
        d = eye(3)
        do i = 1,k
            d = matmul(d,q)
        end do
        call require(maxval(abs(dense-d)) < 1.0e-11_dp, 'matrix power accuracy')
    end do

    p = matrix_power(a,-1,info)
    call require(info == sparse_linalg_invalid_input .and. .not. valid_csr(p), 'negative power')
    p = matrix_power(a,-huge(0),info)
    call require(info == sparse_linalg_invalid_input, 'large negative power')
    b = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    p = matrix_power(b,0,info)
    call require(info == sparse_linalg_invalid_input, 'nonsquare zero power')
    ac = csr_to_csc(b)
    pc = matrix_power(ac,3,info)
    call require(info == sparse_linalg_invalid_input, 'nonsquare CSC power')
    p = inv(b,info)
    call require(info == sparse_linalg_invalid_input, 'nonsquare inverse')

    deallocate(d)
    allocate(d(0,0))
    empty = csr_from_dense(d)
    ac = csr_to_csc(empty)
    p = matrix_power(empty,0,info)
    call require(info == 0 .and. valid_csr(p) .and. p%n_rows == 0, 'empty identity')
    p = matrix_power(empty,huge(0),info)
    call require(info == 0 .and. valid_csr(p), 'empty huge power')
    pc = matrix_power(ac,0,info)
    call require(info == 0 .and. valid_csc(pc), 'empty CSC power')
    p = inv(empty,info)
    call require(info == 0 .and. valid_csr(p), 'empty inverse')
    ic = inv(ac,info)
    call require(info == 0 .and. valid_csc(ic), 'empty CSC inverse')

    ! A strictly upper triangular matrix has fourth power zero.
    deallocate(d)
    allocate(d(4,4))
    d = 0.0_dp
    do i = 1,3
        d(i,i+1) = 2.0_dp
    end do
    a = csr_from_dense(d)
    p = matrix_power(a,huge(0),info)
    call require(info == 0 .and. size(p%data) == 0, 'huge exponent nilpotent')

    ! A pivot is required at the first elimination step.
    a = csr_from_dense(q)
    p = inv(a,info,block_size=2)
    call require(info == 0 .and. valid_csr(p), 'pivoted inverse')
    inverse = csr_to_dense(p)
    call require(maxval(abs(matmul(q,inverse)-eye(3))) < 1.0e-12_dp, 'inverse residual')
    ac = csr_to_csc(a)
    ic = inv(ac,info,block_size=1)
    call require(info == 0 .and. valid_csc(ic), 'CSC inverse valid')
    r = csc_to_csr(ic)
    call require(maxval(abs(csr_to_dense(r)-inverse)) < 1.0e-12_dp, 'CSC inverse')

    z = reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2])
    a = csr_from_dense(z)
    p = inv(a,info)
    call require(info == sparse_linalg_singular .and. .not. valid_csr(p), 'singular inverse')
    ac = csr_to_csc(a)
    ic = inv(ac,info)
    call require(info == sparse_linalg_singular, 'singular CSC inverse')

    a = csr_from_coo(2,2,[1,1,2,1],[1,2,2,2], &
        [2.0_dp,0.0_dp,-4.0_dp,0.0_dp])
    call require(size(a%data) == 3, 'stored zero preserved')
    p = inv(a,info)
    call require(info == 0 .and. size(p%data) == 2, 'diagonal inverse shortcut')
    inverse = csr_to_dense(p)
    call require(abs(inverse(1,1)-0.5_dp) < 1.0e-15_dp, 'first reciprocal')
    call require(abs(inverse(2,2)+0.25_dp) < 1.0e-15_dp, 'second reciprocal')
    p = matrix_power(a,1,info)
    call require(info == 0 .and. size(p%data) == size(a%data), 'identity exponent preserves zeros')

    a = csr_from_dense(reshape([2.0_dp,0.0_dp,1.0e-8_dp,3.0_dp],[2,2]))
    p = inv(a,info,drop_tol=1.0e-8_dp)
    call require(info == 0 .and. size(p%data) == 2, 'inverse sparsification')
    p = inv(a,info,drop_tol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'invalid drop tolerance')
    p = inv(a,info,block_size=0)
    call require(info == sparse_linalg_invalid_input, 'invalid block size')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    p = inv(a,info,pivot_tol=nan_value)
    call require(info == sparse_linalg_invalid_input, 'invalid pivot threshold')
    p = inv(a,info,pivot_tol=10.0_dp)
    call require(info == sparse_linalg_singular, 'pivot rejection')

    a = csr_from_dense(reshape([1.0e-320_dp],[1,1]))
    p = inv(a,info)
    call require(info == sparse_linalg_breakdown, 'inverse overflow')
    a = csr_from_dense(reshape([1.0e200_dp],[1,1]))
    p = matrix_power(a,2,info)
    call require(info == sparse_linalg_breakdown, 'power overflow')

    ! Sparse diagonal problems must not create an n-by-n dense workspace.
    n = 12000
    allocate(diagonal(n))
    diagonal = 2.0_dp
    a = csr_eye(n)
    a%data = diagonal
    p = matrix_power(a,27,info)
    call require(info == 0 .and. size(p%data) == n, 'large sparse power structure')
    call require(abs(p%data(n)-2.0_dp**27) < 1.0e-5_dp, 'large sparse power values')
    p = inv(a,info)
    call require(info == 0 .and. size(p%data) == n, 'large diagonal inverse structure')
    call require(abs(p%data(n)-0.5_dp) < 1.0e-15_dp, 'large diagonal inverse')

    bad = csr_matrix()
    p = matrix_power(bad,2,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR power')
    p = inv(bad,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR inverse')

    print '(a)', 'test_sparse_matrix_functions: PASS'
contains
    function eye(size_) result(out)
        integer, intent(in) :: size_ !! Identity matrix order.
        real(dp) :: out(size_,size_)
        integer :: j
        out = 0.0_dp
        do j = 1,size_
            out(j,j) = 1.0_dp
        end do
    end function eye

    subroutine require(condition, message)
        logical, intent(in) :: condition !! Assertion condition.
        character(len=*), intent(in) :: message !! Assertion description.
        if (.not. condition) then
            print '(2a)', 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_matrix_functions
