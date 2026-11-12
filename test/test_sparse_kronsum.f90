! SPDX-License-Identifier: MIT
program test_sparse_kronsum
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, &
        csr_kronsum, csr_to_dense, csr_eye, valid_csr, &
        sparse_success, sparse_invalid_input
    implicit none
    type(csr_matrix) :: a, b, c, empty
    real(dp), allocatable :: matrix(:,:)
    real(dp) :: da(3,3), db(2,2), expected(6,6)
    integer :: i, j, k, l, status

    da = reshape([1.0_dp,2.0_dp,0.0_dp,0.0_dp,3.0_dp,-2.0_dp, &
                  1.0_dp,0.0_dp,4.0_dp],[3,3])
    db = reshape([2.0_dp,-3.0_dp,1.0_dp,5.0_dp],[2,2])
    a = csr_from_dense(da)
    b = csr_from_dense(db)
    c = csr_kronsum(a,b,status)
    call require(status == sparse_success .and. valid_csr(c), 'valid rectangular-factor kronsum')
    call require(c%n_rows == 6 .and. c%n_cols == 6, 'kronsum shape')
    expected = 0.0_dp
    do i = 1, 2
        do j = 1, 3
            do k = 1, 3
                expected((i-1)*3+j,(i-1)*3+k) = da(j,k)
            end do
            do l = 1, 2
                expected((i-1)*3+j,(l-1)*3+j) = &
                    expected((i-1)*3+j,(l-1)*3+j)+db(i,l)
            end do
        end do
    end do
    matrix = csr_to_dense(c)
    call require(maxval(abs(matrix-expected)) < 1.0e-15_dp, 'kronsum arithmetic')

    ! Preserve a stored zero produced by cancellation on the diagonal.
    a = csr_from_coo(2,2,[1,2],[1,2],[1.0_dp,0.0_dp])
    b = csr_from_coo(1,1,[1],[1],[-1.0_dp])
    c = csr_kronsum(a,b,status)
    call require(status == 0 .and. valid_csr(c), 'explicit zero cancellation valid')
    call require(size(c%data) == 2 .and. any(c%data == 0.0_dp), 'explicit cancelled zero retained')
    call require(all(c%data == [0.0_dp,-1.0_dp]), 'diagonal sums')
    empty = csr_eye(0)
    c = csr_kronsum(a,empty,status)
    call require(status == 0 .and. c%n_rows == 0 .and. valid_csr(c), 'zero-sized B')
    c = csr_kronsum(empty,a,status)
    call require(status == 0 .and. c%n_rows == 0 .and. valid_csr(c), 'zero-sized A')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    c = csr_kronsum(a,b,status)
    call require(status == sparse_invalid_input, 'nonsquare factor rejected')
    a = csr_eye(3)
    b%data(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    c = csr_kronsum(a,b,status)
    call require(status == sparse_invalid_input, 'nonfinite input rejected')

    ! Overlapping, individually finite diagonals can overflow on addition.
    a = csr_from_coo(1,1,[1],[1],[huge(1.0_dp)])
    b = csr_from_coo(1,1,[1],[1],[huge(1.0_dp)])
    c = csr_kronsum(a,b,status)
    call require(status == sparse_invalid_input, 'nonfinite overlapping sum rejected')

    ! One billion logical columns in a factor, no dense storage or costly
    ! work arrays: rejected on representability, not due to allocation.
    a = csr_from_coo(47000,47000,[1],[47000],[1.0_dp])
    b = csr_eye(47000)
    c = csr_kronsum(a,b,status)
    call require(status == sparse_invalid_input, 'dimension product overflow rejected')
    print '(a)', 'test_sparse_kronsum: PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! Assertion condition.
        character(len=*), intent(in) :: message !! Failure context.
        if (.not. ok) then
            print '(a)', 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_kronsum
