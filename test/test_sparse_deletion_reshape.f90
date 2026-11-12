! SPDX-License-Identifier: MIT
! Native checks for sparse structural deletion, assignment, and C-order reshape.
program test_sparse_deletion_reshape
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_from_dense, &
        csr_to_dense, csr_get, csr_assign, csr_delete_entries, csr_delete_rows, &
        csr_delete_columns, csr_reshape, valid_csr, sparse_success, sparse_invalid_input
    implicit none
    type(csr_matrix) :: a, b, src, unchanged
    integer :: status, i
    real(dp) :: expected(3,5)
    integer, allocatable :: big_rows(:), big_cols(:)

    a = csr_from_coo(3,5,[1,1,2,3,3],[2,5,4,1,5], &
        [0.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp])
    call csr_delete_entries(a,[1,1,2,2],[2,2,4,1],status)
    call require(status == 0 .and. valid_csr(a), 'delete stored coordinates')
    call require(size(a%data) == 3, 'stored zero removed')
    call require(a%indptr(2) == 2, 'row pointer after deletion')
    call csr_delete_entries(a,[3],[3],status)
    call require(status == 0 .and. size(a%data) == 3, 'delete missing entry')
    unchanged = a
    call csr_delete_entries(a,[0],[2],status)
    call require(status == sparse_invalid_input, 'out-of-bounds deletion rejected')
    call equal(a,csr_to_dense(unchanged),'invalid entry deletion is transactional')
    call csr_delete_entries(a,[1,2],[1],status)
    call require(status == sparse_invalid_input, 'mismatched delete lengths rejected')
    call equal(a,csr_to_dense(unchanged),'mismatch is transactional')

    a = csr_from_coo(4,5,[1,2,2,3,4],[1,1,4,2,5], &
        [0.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp])
    call csr_delete_rows(a,[2,2,4],status)
    call require(status == 0 .and. valid_csr(a), 'delete rows')
    call require(a%n_rows == 2 .and. a%n_cols == 5, 'row deletion shape')
    call require(size(a%data) == 2 .and. a%data(1) == 0.0_dp, 'stored zero preserved')
    call require(a%indptr(2) == 2 .and. a%indptr(3) == 3, 'row deletion pointers')
    unchanged = a
    call csr_delete_rows(a,[3],status)
    call require(status == sparse_invalid_input, 'invalid deleted row')
    call equal(a,csr_to_dense(unchanged),'invalid row deletion transactional')
    call csr_delete_rows(a,[1,2],status)
    call require(status == 0 .and. a%n_rows == 0 .and. size(a%data) == 0, 'all rows deleted')
    call require(valid_csr(a), 'zero-row CSR')

    a = csr_from_coo(3,6,[1,1,2,2,3],[1,6,2,4,3], &
        [1.0_dp,0.0_dp,2.0_dp,3.0_dp,4.0_dp])
    call csr_delete_columns(a,[2,2,4,5],status)
    call require(status == 0 .and. valid_csr(a), 'delete columns')
    call require(a%n_cols == 3 .and. a%n_rows == 3, 'column deletion shape')
    call require(a%indices(2) == 3 .and. a%data(2) == 0.0_dp, 'zero and remapped column')
    call require(size(a%data) == 3, 'deleted column values removed')
    call csr_delete_columns(a,[1,2,3],status)
    call require(status == 0 .and. a%n_cols == 0 .and. size(a%data) == 0, 'delete all columns')
    call require(valid_csr(a), 'zero-column CSR')

    a = csr_from_coo(2,6,[1,1,2,2],[1,6,2,4], &
        [0.0_dp,2.0_dp,3.0_dp,4.0_dp])
    b = csr_reshape(a,3,4,status)
    call require(status == 0 .and. valid_csr(b), 'C-order reshape')
    expected = 0.0_dp
    call require(size(b%data) == 4, 'reshape stored zeros')
    call require(b%indices(1) == 1 .and. b%data(1) == 0.0_dp, 'reshaped explicit zero')
    call require(b%indptr(3) == 4, 'reshaped row pointers')
    call require(csr_get(b,2,2) == 2.0_dp .and. csr_get(b,2,4) == 3.0_dp, 'reshaped positions')
    call require(csr_get(b,3,2) == 4.0_dp, 'reshaped last row')
    unchanged = a
    b = csr_reshape(a,5,3,status)
    call require(status == sparse_invalid_input, 'unequal reshape dimensions rejected')
    call equal(a,csr_to_dense(unchanged),'reshape never changes source')
    a = csr_from_coo(0,9,[integer::],[integer::],[real(dp)::])
    b = csr_reshape(a,7,0,status)
    call require(status == 0 .and. b%n_rows == 7 .and. b%n_cols == 0, 'empty reshape')
    call require(valid_csr(b), 'empty reshape canonical')

    ! Repeated row and column selectors: last assignment wins, including an
    ! implicit zero at the last source location; only explicitly stored zeros
    ! should consume output storage.
    a = csr_from_dense(reshape([(real(i,dp),i=1,15)],[3,5]))
    src = csr_from_coo(3,3,[1,2,3],[1,2,3],[4.0_dp,0.0_dp,9.0_dp])
    call csr_assign(a,[2,2,1],[5,1,5],src,status)
    call require(status == 0 .and. valid_csr(a), 'repeated sparse selection')
    expected = reshape([(real(i,dp),i=1,15)],[3,5])
    expected(2,5) = 0.0_dp
    expected(2,1) = 0.0_dp
    expected(1,5) = 9.0_dp
    expected(1,1) = 0.0_dp
    call equal(a,expected,'last sparse source occurrence')
    call require(.not. any(a%indices(a%indptr(2):a%indptr(3)-1) == 5), &
        'implicit zero does not consume storage')
    call require(any(a%indices(a%indptr(2):a%indptr(3)-1) == 1), &
        'explicit zero from source stored')

    ! A 50,000 by 50,000 replacement with just two source entries must
    ! use memory proportional to the stored entries, not 2.5 billion cells.
    a = csr_from_coo(50000,1000000000,[100,2],[100,999999999], &
        [77.0_dp,8.0_dp])
    src = csr_from_coo(50000,50000,[1,50000],[1,50000],[0.0_dp,99.0_dp])
    big_rows = [(i,i=1,50000)]
    big_cols = [(i,i=1,50000)]
    call csr_assign(a,big_rows,big_cols,src,status)
    call require(status == 0 .and. valid_csr(a), 'huge sparse-to-sparse overwrite')
    call require(size(a%data) == 3, 'huge selection keeps only explicit entries')
    call require(csr_get(a,100,100) == 0.0_dp .and. &
        csr_get(a,50000,50000) == 99.0_dp, 'huge selection overwrites correctly')
    call require(csr_get(a,2,999999999) == 8.0_dp, 'untouched sparse entry')

    ! Very wide matrices must avoid a column-sized work array.
    a = csr_from_coo(2,1000000000,[1,2],[999999999,1000000000], &
        [7.0_dp,8.0_dp])
    src = csr_from_coo(2,1000000000,[1],[1],[0.0_dp])
    call csr_assign(a,[2,1],[integer::],csr_from_coo(2,0, &
        [integer::],[integer::],[real(dp)::]),status)
    call require(status == 0 .and. size(a%data) == 2, 'zero-width selection')
    call csr_delete_columns(a,[1000000000,4,4],status)
    call require(status == 0 .and. a%n_cols == 999999998, 'huge logical column deletion')
    call require(size(a%data) == 1 .and. a%indices(1) == 999999998, &
        'huge column deletion scales with nnz')
    call require(csr_get(a,1,999999998) == 7.0_dp, 'huge column remap')
    call csr_delete_entries(a,[1],[999999998],status)
    call require(status == 0 .and. size(a%data) == 0, 'huge coordinate deletion')

    print '(a)', 'sparse deletion, sparse assignment and reshape native tests passed'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    subroutine equal(matrix, reference, message)
        type(csr_matrix), intent(in) :: matrix
        real(dp), intent(in) :: reference(:,:)
        character(len=*), intent(in) :: message
        real(dp), allocatable :: actual(:,:)
        call require(valid_csr(matrix),message//' structure')
        actual = csr_to_dense(matrix)
        call require(size(actual,1) == size(reference,1) .and. &
            size(actual,2) == size(reference,2),message//' shape')
        call require(all(abs(actual-reference) < 1.0e-12_dp), message//' entries')
    end subroutine equal
end program test_sparse_deletion_reshape
