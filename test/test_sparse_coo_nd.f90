! SPDX-License-Identifier: MIT
program test_sparse_coo_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo, only : coo_matrix, coo_from_triplets, coo_to_dense
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_from_coo, coo_nd_to_coo, &
        coo_nd_is_canonical, coo_nd_sum_duplicates, coo_nd_eliminate_zeros, &
        coo_nd_reshape, coo_nd_transpose, coo_nd_get, &
        coo_nd_from_dense, coo_nd_to_dense, &
        coo_nd_from_dense_flat, coo_nd_to_dense_flat
    implicit none
    type(coo_nd_array) :: a, b, c
    type(coo_matrix) :: coo, coo2
    integer :: ierr, k
    integer, allocatable :: many_coords(:,:)
    real(dp), allocatable :: many_values(:)
    real(dp) :: value
    real(dp), allocatable :: flat(:), orig_flat(:)
    real(dp) :: dense(2,3,4), filled(2,3,4)
    real(dp) :: dense7(2,1,1,1,1,1,2), filled7(2,1,1,1,1,1,2)
    real(dp) :: dense8(1,1,1,1,1,1,1,1)
    integer :: xyz(3,7)

    xyz(:,1) = [2,2,4]
    xyz(:,2) = [1,2,2]
    xyz(:,3) = [1,1,4]
    xyz(:,4) = [1,2,2]
    xyz(:,5) = [1,2,2]
    xyz(:,6) = [2,1,1]
    xyz(:,7) = [2,1,1]
    a = coo_nd_from_triplets([2,3,4],xyz, &
        [5.0_dp,1.0_dp,-1.0_dp,2.0_dp,-3.0_dp,0.0_dp,7.0_dp],ierr)
    call require(ierr == 0 .and. valid_coo_nd(a), '3-D triplets')
    call require(size(a%data) == 7 .and. .not. coo_nd_is_canonical(a), 'unsorted duplicates retained')
    value = coo_nd_get(a,[1,2,2],ierr)
    call require(ierr == 0 .and. value == 0.0_dp, 'duplicate cancellation')
    value = coo_nd_get(a,[2,1,1],ierr)
    call require(ierr == 0 .and. value == 7.0_dp, 'get with stored zero')
    value = coo_nd_get(a,[3,1,1],ierr)
    call require(ierr == sparse_invalid_input .and. value == 0.0_dp, 'invalid get')
    flat = coo_nd_to_dense_flat(a,ierr)
    call require(ierr == 0 .and. size(flat) == 24, 'C-order dense flatten')
    call require(flat(20) == 5.0_dp .and. flat(4) == -1.0_dp, 'C-order entries')
    call require(flat(13) == 7.0_dp .and. flat(6) == 0.0_dp, 'C-order duplicate sums')
    call require(count(flat /= 0.0_dp) == 3, 'C-order support count')
    orig_flat = flat
    call coo_nd_to_dense(a,filled,ierr)
    call require(ierr == 0 .and. all(filled(1,2,:) == 0.0_dp), 'assumed-rank dense output')
    call require(filled(2,2,4) == 5.0_dp .and. filled(2,1,1) == 7.0_dp, '3-D dense positions')
    b = coo_nd_transpose(a,[3,1,2],ierr)
    call require(ierr == 0 .and. all(b%shape == [4,2,3]), 'permuted transpose')
    call require(all(b%coords(:,1) == [4,2,2]), 'permuted transpose coordinates')
    call require(all(b%data == a%data), 'permuted transpose values')
    c = coo_nd_transpose(a,status=ierr)
    call require(ierr == 0 .and. all(c%shape == [4,3,2]), 'default reversed axes')
    call require(all(c%coords(:,2) == [2,2,1]), 'reversed coordinates')
    b = coo_nd_transpose(a,[1,1,3],ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_coo_nd(b), 'repeated axis rejected')
    b = coo_nd_transpose(a,[1,2],ierr)
    call require(ierr == sparse_invalid_input, 'incomplete axis list rejected')
    b = coo_nd_reshape(a,[4,3,2],status=ierr)
    call require(ierr == 0 .and. all(b%shape == [4,3,2]), 'C-order reshape')
    call require(all(b%coords(:,2) == [1,3,2]), 'C-order reshape coordinate')
    c = coo_nd_reshape(b,[2,3,4],status=ierr)
    call require(ierr == 0 .and. all(coo_nd_to_dense_flat(c) == orig_flat), 'C reshape roundtrip')
    b = coo_nd_reshape(a,[4,3,2],order='F',status=ierr)
    call require(ierr == 0 .and. all(b%coords(:,1) == [2,3,2]), 'F-order reshape coordinate')
    c = coo_nd_reshape(b,[2,3,4],order='F',status=ierr)
    call require(ierr == 0 .and. all(coo_nd_to_dense_flat(c) == orig_flat), 'F reshape roundtrip')
    b = coo_nd_reshape(a,[3,8],status=ierr)
    call require(ierr == 0 .and. all(b%shape == [3,8]), 'reshape changes rank')
    b = coo_nd_reshape(a,[2,2],status=ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_coo_nd(b), 'mismatched reshape rejected')
    b = coo_nd_reshape(a,[24],order='Q',status=ierr)
    call require(ierr == sparse_invalid_input, 'invalid layout rejected')

    b = a
    call coo_nd_eliminate_zeros(b,ierr)
    call require(ierr == 0 .and. size(b%data) == 6, 'stored zeros eliminated without coalescing')
    call coo_nd_sum_duplicates(b,ierr)
    call require(ierr == 0 .and. coo_nd_is_canonical(b), 'canonical sort')
    call require(size(b%data) == 4, 'duplicates combined retaining cancellation zeros')
    call require(all(b%coords(:,1) == [1,1,4]), 'canonical C-order first')
    call require(all(b%coords(:,2) == [1,2,2]), 'canonical C-order second')
    call require(all(b%coords(:,3) == [2,1,1]), 'canonical C-order third')
    call require(all(b%coords(:,4) == [2,2,4]), 'canonical C-order fourth')
    call require(all(b%data == [-1.0_dp,0.0_dp,7.0_dp,5.0_dp]), 'canonical values')
    call coo_nd_eliminate_zeros(b,ierr)
    call require(ierr == 0 .and. size(b%data) == 3, 'post-canonical zero elimination')

    dense = 0.0_dp
    dense(1,2,2) = 9.0_dp
    dense(2,1,3) = -6.0_dp
    b = coo_nd_from_dense(dense,ierr)
    call require(ierr == 0 .and. coo_nd_is_canonical(b), 'rank-3 dense construction')
    call require(size(b%data) == 2, 'dense constructor excludes zeros')
    filled = 42.0_dp
    call coo_nd_to_dense(b,filled,ierr)
    call require(ierr == 0 .and. all(filled == dense), 'rank-3 assumed-rank roundtrip')
    flat = coo_nd_to_dense_flat(b,ierr,'F')
    call require(ierr == 0 .and. all(flat == reshape(dense,[24])), 'F flat agrees with Fortran storage')
    c = coo_nd_from_dense_flat([2,3,4],flat,ierr,'F')
    call require(ierr == 0 .and. all(coo_nd_to_dense_flat(c) == coo_nd_to_dense_flat(b)), &
                 'F dense-flat constructor')
    flat = coo_nd_to_dense_flat(b,ierr,'C')
    c = coo_nd_from_dense_flat([2,3,4],flat,ierr,'C')
    call require(ierr == 0 .and. coo_nd_is_canonical(c), 'C dense-flat constructor')
    call require(all(coo_nd_to_dense_flat(c) == flat), 'C dense-flat roundtrip')

    dense7 = 0.0_dp
    dense7(1,1,1,1,1,1,2) = 3.0_dp
    dense7(2,1,1,1,1,1,1) = -4.0_dp
    a = coo_nd_from_dense(dense7,ierr)
    call require(ierr == 0 .and. all(a%shape == shape(dense7)), 'rank 7 dense constructor')
    call require(size(a%data) == 2, 'rank 7 sparse entries')
    call coo_nd_to_dense(a,filled7,ierr)
    call require(ierr == 0 .and. all(filled7 == dense7), 'rank 7 dense roundtrip')
    dense8 = 0.0_dp
    b = coo_nd_from_dense(dense8,ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_coo_nd(b), &
                 'unsupported assumed-rank 8 dense constructor rejected')

    coo = coo_from_triplets(2,3,[2,1,2],[3,2,3], &
                            [1.0_dp,0.0_dp,2.0_dp],ierr)
    a = coo_nd_from_coo(coo,ierr)
    call require(ierr == 0 .and. size(a%data) == 3, '2D COO to N-D retains entries')
    coo2 = coo_nd_to_coo(a,ierr)
    call require(ierr == 0 .and. all(coo2%row == coo%row) .and. all(coo2%col == coo%col), &
                 '2D COO conversion roundtrip coordinates')
    call require(all(coo_to_dense(coo) == coo_to_dense(coo2)), '2D COO conversion values')
    b = coo_nd_from_triplets([3],reshape([1,1,3],[1,3]),[1.0_dp,-1.0_dp,5.0_dp],ierr)
    call require(ierr == 0 .and. valid_coo_nd(b), 'rank 1 COO')
    coo2 = coo_nd_to_coo(b,ierr)
    call require(ierr == sparse_invalid_input, 'rank 1 rejected by 2D converter')
    call coo_nd_sum_duplicates(b,ierr)
    call require(ierr == 0 .and. size(b%data) == 2 .and. b%data(1) == 0.0_dp, &
                 'rank 1 duplicate zero preserved')
    a = coo_nd_from_triplets([0,3,2],reshape([integer::],[3,0]),[real(dp)::],ierr)
    call require(ierr == 0 .and. valid_coo_nd(a), 'empty zero-extent dimension')
    flat = coo_nd_to_dense_flat(a,ierr)
    call require(ierr == 0 .and. size(flat) == 0, 'empty dense output')
    b = coo_nd_reshape(a,[3,0,2,1],status=ierr)
    call require(ierr == 0 .and. all(b%shape == [3,0,2,1]), 'zero-element reshape')
    b = coo_nd_transpose(a,status=ierr)
    call require(ierr == 0 .and. all(b%shape == [2,3,0]), 'empty transpose')

    a = coo_nd_from_triplets([1000000,1000000,1000000], &
        reshape([1000000,1,1000000,1,1000000,2],[3,2]),[3.0_dp,-1.0_dp],ierr)
    call require(ierr == 0 .and. valid_coo_nd(a), 'sparse large shape without dense allocation')
    value = coo_nd_get(a,[1000000,1,1000000],ierr)
    call require(ierr == 0 .and. value == 3.0_dp, 'large index lookup')
    b = coo_nd_reshape(a,[1000000000,1000000000],status=ierr)
    call require(ierr == 0 .and. valid_coo_nd(b), 'int64 large-shape reshape')
    call require(all(b%coords(:,1) == [999999001,1000000]), '64-bit flattened index first')
    call require(all(b%coords(:,2) == [1000,999000002]), '64-bit flattened index second')
    flat = coo_nd_to_dense_flat(a,ierr)
    call require(ierr == sparse_invalid_input .and. size(flat) == 0, 'dense size cap')
    flat = coo_nd_to_dense_flat(a,ierr,max_elements=-1_int64)
    call require(ierr == sparse_invalid_input, 'invalid size cap')
    b = coo_nd_from_triplets([huge(0),huge(0),huge(0)],reshape([integer::],[3,0]),[real(dp)::],ierr)
    call require(ierr == sparse_success, 'overflow-size sparse allocation valid')
    c = coo_nd_reshape(b,[1,1,1],status=ierr)
    call require(ierr == sparse_invalid_input, 'reshape product overflow rejected')
    b = coo_nd_from_triplets([-1,2],reshape([integer::],[2,0]),[real(dp)::],ierr)
    call require(ierr == sparse_invalid_input, 'negative shape rejected')
    b = coo_nd_from_triplets([2,3],reshape([0,1],[2,1]),[1.0_dp],ierr)
    call require(ierr == sparse_invalid_input, 'zero coordinate rejected')
    b = coo_nd_from_triplets([2,3],reshape([3,1],[2,1]),[1.0_dp],ierr)
    call require(ierr == sparse_invalid_input, 'out of bounds coordinate rejected')
    b = coo_nd_from_triplets([2,3],reshape([1,1],[2,1]), &
                             [ieee_value(0.0_dp,ieee_quiet_nan)],ierr)
    call require(ierr == sparse_invalid_input, 'NaN rejected')
    b = coo_nd_from_triplets([2,3],reshape([1,1],[2,1]),[1.0_dp,2.0_dp],ierr)
    call require(ierr == sparse_invalid_input, 'wrong triplet length rejected')
    b = coo_nd_from_dense_flat([2,2],[1.0_dp,3.0_dp,5.0_dp],ierr)
    call require(ierr == sparse_invalid_input, 'incorrect dense-flat length rejected')
    b = coo_nd_from_triplets([2,2],reshape([1,1,1,1],[2,2]),[huge(1.0_dp),huge(1.0_dp)],ierr)
    call require(ierr == 0, 'finite duplicate components')
    call coo_nd_sum_duplicates(b,ierr)
    call require(ierr == sparse_invalid_input .and. size(b%data) == 2, 'overflow transactional')
    flat = coo_nd_to_dense_flat(b,ierr)
    call require(ierr == sparse_invalid_input .and. size(flat) == 0, 'dense overflow rejected')
    value = coo_nd_get(b,[1,1],ierr)
    call require(ierr == sparse_invalid_input .and. value == 0.0_dp, 'get overflow rejected')

    ! Stress the O(nnz log nnz) merge sort with many duplicates across axes.
    allocate(many_coords(4,5000),many_values(5000))
    do k = 1, 5000
        many_coords(:,k) = [1+mod(11*k,7),1+mod(13*k,5), &
            1+mod(17*k,4),1+mod(19*k,3)]
        many_values(k) = real(mod(k,7)-3,dp)
    end do
    a = coo_nd_from_triplets([7,5,4,3],many_coords,many_values,ierr)
    call require(ierr == 0, 'large COO stress input')
    orig_flat = coo_nd_to_dense_flat(a,ierr)
    call require(ierr == 0, 'large COO stress reference')
    call coo_nd_sum_duplicates(a,ierr)
    call require(ierr == 0 .and. coo_nd_is_canonical(a), 'large COO sorted')
    call require(size(a%data) == 420, 'large COO unique positions')
    flat = coo_nd_to_dense_flat(a,ierr)
    call require(ierr == 0 .and. all(flat == orig_flat), 'large COO aggregation values')

    print *, 'sparse N-D COO: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd
