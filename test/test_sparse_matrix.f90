! SPDX-License-Identifier: MIT
program test_sparse_matrix
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, valid_csr, csr_from_coo, csr_from_dense, &
        csr_to_coo, csr_to_dense, csr_eye, csr_transpose, csr_matvec, csr_add, &
        csr_matmul, csr_elementwise_multiply, csr_scale, csr_eliminate_zeros, &
        sparse_success, sparse_invalid_input
    implicit none
    type(csr_matrix) :: a, b, c, t, bad
    real(dp), allocatable :: dense(:,:), lhs(:,:), rhs(:,:), y(:), vv(:)
    integer, allocatable :: ri(:), ci(:)
    integer :: status, n, i

    a = csr_from_coo(3,4,[3,1,2,1,3,1,2],[4,2,1,2,4,4,3], &
        [2.0_dp,3.0_dp,-1.0_dp,-2.0_dp,5.0_dp,0.0_dp,1.0_dp],status)
    call require(status == sparse_success .and. valid_csr(a), 'COO construction')
    call require(size(a%data) == 5, 'duplicate coalescing and stored zero')
    call require(all(a%indptr == [1,3,5,6]), 'row pointers')
    call require(all(a%indices == [2,4,1,3,4]), 'sorted canonical columns')
    call require(all(a%data == [1.0_dp,0.0_dp,-1.0_dp,1.0_dp,7.0_dp]), &
        'summation of COO duplicates')
    dense = csr_to_dense(a,status)
    call require(status == sparse_success, 'dense conversion')
    call require(size(dense,1) == 3 .and. size(dense,2) == 4, 'rectangular shape')
    call require(all(dense(1,:) == [0.0_dp,1.0_dp,0.0_dp,0.0_dp]), 'first row')
    call require(all(dense(2,:) == [-1.0_dp,0.0_dp,1.0_dp,0.0_dp]), 'second row')
    call require(all(dense(3,:) == [0.0_dp,0.0_dp,0.0_dp,7.0_dp]), 'third row')
    call csr_to_coo(a,ri,ci,vv,status)
    call require(status == sparse_success .and. all(ri == [1,1,2,2,3]), 'to COO rows')
    call require(all(ci == a%indices) .and. all(vv == a%data), 'to COO cols and values')

    t = csr_transpose(a,status)
    call require(status == sparse_success .and. valid_csr(t), 'CSR transpose')
    call require(all(csr_to_dense(t) == transpose(dense)), 'transpose values')
    c = csr_transpose(t,status)
    call require(all(c%indices == a%indices) .and. all(c%indptr == a%indptr), &
        'double transpose indices')
    call require(all(c%data == a%data), 'double transpose data')
    y = csr_matvec(a,[5.0_dp,6.0_dp,7.0_dp,8.0_dp],status)
    call require(status == sparse_success .and. all(y == [6.0_dp,2.0_dp,56.0_dp]), &
        'matrix-vector product')

    b = csr_from_dense(reshape([1.0_dp,3.0_dp,0.0_dp,2.0_dp, &
        -1.0_dp,0.0_dp,4.0_dp,5.0_dp], [4,2]),status)
    call require(status == sparse_success .and. valid_csr(b), 'dense constructor')
    c = csr_matmul(a,b,status)
    call require(status == sparse_success .and. valid_csr(c), 'matrix-matrix multiplication')
    lhs = csr_to_dense(c)
    rhs = matmul(dense,csr_to_dense(b))
    call require(all(abs(lhs-rhs) < 1.0e-12_dp), 'matrix product numerical values')
    call require(size(c%indices) <= 6, 'matrix multiply sparse output')

    c = csr_add(a,a,status,alpha=2.0_dp,beta=-2.0_dp)
    call require(status == sparse_success .and. valid_csr(c), 'CSR addition')
    call require(size(c%data) == 0 .and. all(csr_to_dense(c) == 0.0_dp), &
        'cancellation removes arithmetic zeros')
    c = csr_add(a,a,status)
    call require(all(csr_to_dense(c) == 2.0_dp*dense), 'CSR numeric addition')
    c = csr_elementwise_multiply(a,a,status)
    call require(status == sparse_success .and. all(csr_to_dense(c) == dense*dense), &
        'Hadamard product')
    c = csr_scale(a,-3.0_dp,status)
    call require(status == sparse_success .and. size(c%data) == 5, 'scale keeps zeros')
    call require(all(csr_to_dense(c) == -3.0_dp*dense), 'scalar multiplication')
    call csr_eliminate_zeros(c,status)
    call require(status == sparse_success .and. size(c%data) == 4, 'drop stored zeros')
    call require(all(csr_to_dense(c) == -3.0_dp*dense), 'zero elimination values')

    c = csr_eye(4,status=status)
    call require(status == sparse_success .and. size(c%data) == 4, 'identity matrix')
    call require(all(csr_to_dense(c) == csr_to_dense(csr_matmul(c,c))), &
        'identity multiplied by identity')
    c = csr_eye(3,5,offset=2,status=status)
    call require(size(c%data) == 3 .and. all(c%indices == [3,4,5]), 'upper eye')
    c = csr_eye(5,3,offset=-2,status=status)
    call require(size(c%data) == 3 .and. all(c%indices == [1,2,3]), 'lower eye')
    c = csr_eye(2,2,offset=5,status=status)
    call require(size(c%data) == 0 .and. all(c%indptr == 1), 'empty diagonal')

    c = csr_from_coo(0,4,[integer::],[integer::],[real(dp)::],status)
    call require(status == sparse_success .and. valid_csr(c), 'zero rows')
    call require(all(shape(csr_to_dense(c)) == [0,4]), 'zero rows preserve shape')
    t = csr_transpose(c,status)
    call require(status == sparse_success .and. valid_csr(t), 'transpose zero rows')
    call require(all(shape(csr_to_dense(t)) == [4,0]), 'zero columns preserved')
    c = csr_from_coo(4,0,[integer::],[integer::],[real(dp)::],status)
    call require(status == sparse_success .and. valid_csr(c), 'zero columns')
    t = csr_matmul(c,csr_from_coo(0,3,[integer::],[integer::],[real(dp)::]),status)
    call require(status == sparse_success .and. all(shape(csr_to_dense(t)) == [4,3]), &
        'rectangular empty matrix multiply')

    c = csr_from_coo(2,2,[1,1],[2,2],[4.0_dp,-4.0_dp],status)
    call require(status == sparse_success .and. size(c%data) == 1, &
        'exact cancellation of duplicate COO retains structural zero')
    call require(c%data(1) == 0.0_dp, 'summed duplicate value zero')
    call csr_eliminate_zeros(c,status)
    call require(size(c%data) == 0 .and. valid_csr(c), 'remove cancellation zero')
    c = csr_from_coo(2,2,[0],[1],[1.0_dp],status)
    call require(status == sparse_invalid_input .and. .not. valid_csr(c), 'bad COO row')
    c = csr_from_coo(2,2,[1],[3],[1.0_dp],status)
    call require(status == sparse_invalid_input, 'bad COO column')
    c = csr_from_coo(2,2,[1],[1,2],[1.0_dp],status)
    call require(status == sparse_invalid_input, 'mismatched coordinate lengths')
    c = csr_from_coo(2,2,[1],[1],[ieee_value(0.0_dp,ieee_quiet_nan)],status)
    call require(status == sparse_invalid_input, 'reject NaN')
    c = csr_from_coo(1,1,[1,1],[1,1],[huge(1.0_dp),huge(1.0_dp)],status)
    call require(status == sparse_invalid_input .and. .not. valid_csr(c), 'reject duplicate overflow')
    bad = a
    bad%indices(2) = 2
    call require(.not. valid_csr(bad), 'reject duplicate CSR indices')
    c = csr_add(bad,a,status)
    call require(status == sparse_invalid_input, 'reject invalid operand')
    y = csr_matvec(a,[1.0_dp],status)
    call require(status == sparse_invalid_input .and. size(y) == 0, 'reject wrong vector length')
    c = csr_matmul(a,a,status)
    call require(status == sparse_invalid_input, 'reject matrix shape mismatch')
    c = csr_add(a,b,status)
    call require(status == sparse_invalid_input, 'reject addition shape mismatch')
    c = csr_eye(-1,status=status)
    call require(status == sparse_invalid_input, 'reject negative size')
    c = csr_from_dense(reshape([huge(1.0_dp)],[1,1]))
    t = csr_matmul(c,c,status)
    call require(status == sparse_invalid_input, 'reject nonfinite product')

    ! O(n) nnz, not O(n**2) storage.  Test a sizeable rectangular product.
    n = 12000
    a = csr_from_coo(n,n+1,[(i,i=1,n)],[(i,i=1,n)], &
        [(2.0_dp,i=1,n)],status)
    b = csr_from_coo(n+1,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(3.0_dp,i=1,n)],status)
    c = csr_matmul(a,b,status)
    call require(status == sparse_success .and. valid_csr(c), 'large sparse multiply')
    call require(size(c%data) == n .and. all(c%data == 6.0_dp), &
        'large rectangular product stays sparse')
    call require(all(c%indices == [(i,i=1,n)]), 'large product structure')
    print '(a)', 'test_sparse_matrix: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! True if the regression condition passed.
        character(len=*), intent(in) :: message !! Test assertion diagnostic.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_matrix
