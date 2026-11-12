! SPDX-License-Identifier: MIT
program test_sparse_coo
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_to_dense, csc_to_dense, &
        sparse_invalid_input, sparse_success
    use scifort_sparse_coo, only : coo_matrix, valid_coo, coo_from_triplets, &
        coo_from_dense, coo_to_dense, coo_from_csr, coo_from_csc, &
        coo_to_csr, coo_to_csc, coo_matvec, coo_transpose, coo_sum_duplicates, &
        coo_eliminate_zeros, coo_get, sparse_find
    implicit none
    type(coo_matrix) :: a, b, c, invalid
    type(csr_matrix) :: csr
    type(csc_matrix) :: csc
    real(dp), allocatable :: dense(:,:), values(:), y(:)
    integer, allocatable :: rows(:), columns(:)
    integer :: info
    real(dp) :: val

    ! Unsorted coordinates, duplicates, cancellation to explicit zero.
    a = coo_from_triplets(3,4, [3,1,2,1,1,1,2,3], &
        [2,3,4,3,1,1,4,2], [3.0_dp,2.0_dp,5.0_dp,-2.0_dp, &
        1.0_dp,-1.0_dp,0.0_dp,4.0_dp],info)
    call require(info == sparse_success .and. valid_coo(a), 'coo construction')
    call require(size(a%data) == 8 .and. all(a%row == [3,1,2,1,1,1,2,3]), 'unsorted coordinates retained')
    dense = coo_to_dense(a,info)
    call require(info == 0 .and. all(shape(dense) == [3,4]), 'dense shape')
    call require(dense(3,2) == 7.0_dp .and. dense(2,4) == 5.0_dp, 'dense summation')
    call require(count(dense /= 0.0_dp) == 2, 'duplicate cancellation')
    val = coo_get(a,3,2,info)
    call require(info == 0 .and. val == 7.0_dp, 'get duplicate sum')
    val = coo_get(a,1,1,info)
    call require(info == 0 .and. val == 0.0_dp, 'get cancellation zero')
    val = coo_get(a,4,1,info)
    call require(info == sparse_invalid_input .and. val == 0.0_dp, 'get invalid coordinate')
    y = coo_matvec(a,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],info)
    call require(info == 0 .and. all(abs(y-[0.0_dp,20.0_dp,14.0_dp]) < 1.0e-12_dp), 'matvec')
    b = coo_transpose(a,info)
    call require(info == 0 .and. all(shape(coo_to_dense(b)) == [4,3]), 'transpose shape')
    call require(all(transpose(dense) == coo_to_dense(b)), 'transpose values')
    call require(all(b%row == a%col) .and. all(b%col == a%row), 'transpose stored order')

    csr = coo_to_csr(a,info)
    call require(info == 0, 'CSR conversion')
    call require(all(csr_to_dense(csr) == dense), 'CSR conversion values')
    call require(size(csr%data) == 4, 'CSR preserves cancellation zeros')
    csc = coo_to_csc(a,info)
    call require(info == 0 .and. size(csc%data) == 4, 'CSC preserves cancellation zeros')
    call require(all(csc_to_dense(csc) == dense), 'CSC conversion values')
    b = coo_from_csr(csr,info)
    call require(info == 0 .and. valid_coo(b) .and. size(b%data) == 4, 'CSR to COO')
    c = coo_from_csc(csc,info)
    call require(info == 0 .and. valid_coo(c) .and. size(c%data) == 4, 'CSC to COO')
    call require(all(coo_to_dense(c) == dense), 'CSC roundtrip values')

    call sparse_find(a,rows,columns,values,info)
    call require(info == 0 .and. size(values) == 2, 'find zero elimination')
    call require(all(rows == [2,3]) .and. all(columns == [4,2]), 'find row-major order')
    call require(all(values == [5.0_dp,7.0_dp]), 'find sums')
    call sparse_find(csr,rows,columns,values,info)
    call require(info == 0 .and. all(values == [5.0_dp,7.0_dp]), 'find CSR')
    call sparse_find(csc,rows,columns,values,info)
    call require(info == 0 .and. all(values == [5.0_dp,7.0_dp]), 'find CSC')

    b = a
    call coo_eliminate_zeros(b,info)
    call require(info == 0 .and. size(b%data) == 7, 'eliminate explicit zeros without summing')
    call coo_sum_duplicates(b,info)
    call require(info == 0 .and. size(b%data) == 4, 'duplicate summation keeps zeros')
    call require(all(b%row == [1,1,2,3]) .and. all(b%col == [1,3,4,2]), 'canonical order')
    call require(all(b%data == [0.0_dp,0.0_dp,5.0_dp,7.0_dp]), 'canonical values')
    call coo_eliminate_zeros(b,info)
    call require(info == 0 .and. size(b%data) == 2, 'canonical zero elimination')

    c = coo_from_dense(dense,info)
    call require(info == 0 .and. size(c%data) == 2, 'from dense nonzeros only')
    call require(all(coo_to_dense(c) == dense), 'dense roundtrip')
    a = coo_from_triplets(0,5,[integer::],[integer::],[real(dp)::],info)
    call require(info == 0 .and. valid_coo(a), 'zero rows')
    dense = coo_to_dense(a,info)
    call require(info == 0 .and. all(shape(dense) == [0,5]), 'zero rows dense shape')
    a = coo_from_triplets(3,0,[integer::],[integer::],[real(dp)::],info)
    call require(info == 0 .and. valid_coo(a), 'zero columns')
    csc = coo_to_csc(a,info)
    call require(info == 0 .and. size(csc%indptr) == 1, 'zero columns CSC')
    csr = coo_to_csr(a,info)
    call require(info == 0 .and. size(csr%indptr) == 4, 'zero columns CSR')

    a = coo_from_triplets(120000,120000,[1,120000],[120000,1],[2.0_dp,-3.0_dp],info)
    call require(info == 0 .and. valid_coo(a), 'large sparse COO')
    csr = coo_to_csr(a,info)
    call require(info == 0 .and. size(csr%data) == 2, 'large sparse conversion')
    y = coo_matvec(a,[real(dp)::],info) ! wrong length: returns empty result
    call require(info == sparse_invalid_input .and. size(y) == 0, 'matvec wrong size')

    invalid = coo_from_triplets(2,2,[1,3],[1,1],[1.0_dp,2.0_dp],info)
    call require(info == sparse_invalid_input .and. .not. valid_coo(invalid), 'invalid row rejected')
    invalid = coo_from_triplets(1,2,[1],[3],[2.0_dp],info)
    call require(info == sparse_invalid_input, 'invalid column rejected')
    invalid = coo_from_triplets(1,1,[1],[1],[ieee_value(0.0_dp,ieee_quiet_nan)],info)
    call require(info == sparse_invalid_input, 'NaN rejected')
    invalid = coo_from_triplets(1,1,[1],[1,1],[1.0_dp],info)
    call require(info == sparse_invalid_input, 'mismatched triplets rejected')
    a = coo_from_triplets(1,1,[1,1],[1,1],[huge(1.0_dp),huge(1.0_dp)],info)
    call require(info == sparse_success .and. valid_coo(a), 'finite entries can overflow on aggregation')
    dense = coo_to_dense(a,info)
    call require(info == sparse_invalid_input .and. size(dense) == 0, 'dense aggregate overflow rejected')
    csr = coo_to_csr(a,info)
    call require(info == sparse_invalid_input, 'CSR aggregate overflow rejected')
    call coo_sum_duplicates(a,info)
    call require(info == sparse_invalid_input .and. size(a%data) == 2, 'failed coalescing preserves input')

    print *, 'sparse COO and find: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo
