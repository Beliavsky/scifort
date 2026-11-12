! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_filter
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_get, coo_nd_is_canonical
    use scifort_sparse_coo_nd_filter, only : coo_nd_compare_scalar, &
        coo_nd_count_nonzero, coo_nd_filter_sparse, coo_nd_filter_points, &
        coo_nd_filter_value
    implicit none
    type(coo_nd_array) :: a, mask, b, large, empty, bad
    integer :: status, k
    integer, allocatable :: stress_coords(:,:), stress_points(:,:)
    real(dp), allocatable :: stress_values(:)
    integer(int64) :: n

    a = coo_nd_from_triplets([3,4],reshape([ &
        1,2, 1,2, 2,4, 2,4, 3,1, 1,3, 3,3],[2,7]), &
        [2.0_dp,-2.0_dp,7.0_dp,0.0_dp,-5.0_dp,0.0_dp,4.0_dp],status)
    call require(status == sparse_success, 'COO source')
    n = coo_nd_count_nonzero(a,status)
    call require(status == sparse_success .and. n == 3_int64, 'count aggregate nonzeros')
    mask = coo_nd_compare_scalar(a,'>',0.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 2, 'positive predicate')
    call require(coo_nd_is_canonical(mask) .and. all(mask%data == 1.0_dp), 'canonical logical mask')
    call require(coo_nd_get(mask,[2,4]) == 1.0_dp .and. &
        coo_nd_get(mask,[3,3]) == 1.0_dp, 'positive mask coordinates')
    mask = coo_nd_compare_scalar(a,'<',0.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 1, 'negative predicate')
    mask = coo_nd_compare_scalar(a,'/=',0.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 3, 'nonzero predicate')
    mask = coo_nd_compare_scalar(a,'abs>',4.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 2, 'absolute predicate')
    mask = coo_nd_compare_scalar(a,'==',7.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 1, 'equality predicate')
    mask = coo_nd_compare_scalar(a,'>=',4.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 2, 'greater-equal predicate')
    mask = coo_nd_compare_scalar(a,'<=',-5.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 1, 'less-equal predicate')
    mask = coo_nd_compare_scalar(a,'==',0.0_dp,status)
    call require(status == sparse_invalid_input .and. .not. valid_coo_nd(mask), &
        'implicit-zero equality unsupported')
    mask = coo_nd_compare_scalar(a,'>',-1.0_dp,status)
    call require(status == sparse_invalid_input, 'implicit-zero comparison unsupported')
    mask = coo_nd_compare_scalar(a,'bad',0.0_dp,status)
    call require(status == sparse_invalid_input, 'invalid comparison operator')
    mask = coo_nd_compare_scalar(a,'>',0.0_dp,status,1_int64)
    call require(status == sparse_invalid_input, 'predicate limit')

    mask = coo_nd_from_triplets([3,4],reshape([ &
        1,2, 3,1, 2,4, 2,4, 3,1],[2,5]), &
        [1.0_dp,1.0_dp,5.0_dp,-5.0_dp,2.0_dp],status)
    call require(status == sparse_success, 'coordinate sparse mask')
    b = coo_nd_filter_sparse(a,mask,status)
    call require(status == sparse_success .and. size(b%data) == 3, &
        'filter retains duplicate source entries')
    call require(all(b%coords(:,1) == [1,2]) .and. &
        all(b%coords(:,3) == [3,1]), 'filter stable coordinate order')
    call require(coo_nd_get(b,[1,2]) == 0.0_dp .and. coo_nd_get(b,[3,1]) == -5.0_dp, &
        'mask combines duplicate coordinates by sum')
    b = coo_nd_filter_sparse(a,mask,status,complement=.true.)
    call require(status == sparse_success .and. size(b%data) == 4, 'complement filter')
    call require(any(b%data == 0.0_dp), 'complement retains stored zeros')
    b = coo_nd_filter_sparse(a,mask,status,max_nnz=2_int64)
    call require(status == sparse_invalid_input .and. .not. valid_coo_nd(b), 'filter nnz limit')
    b = coo_nd_filter_points(a,reshape([3,1,1,2,3,1],[2,3]),status)
    call require(status == sparse_success .and. size(b%data) == 3, 'point-set deduplication')
    b = coo_nd_filter_points(a,reshape([3,1],[2,1]),status,complement=.true.)
    call require(status == sparse_success .and. size(b%data) == 6, 'point complement')
    b = coo_nd_filter_points(a,reshape([4,1],[2,1]),status)
    call require(status == sparse_invalid_input, 'out-of-range point')
    b = coo_nd_filter_points(a,reshape([3,1,1],[3,1]),status)
    call require(status == sparse_invalid_input, 'point rank mismatch')
    b = coo_nd_filter_value(a,'>',0.0_dp,status)
    call require(status == sparse_success .and. size(b%data) == 3, 'aggregate-positive filter')
    call require(coo_nd_get(b,[2,4]) == 7.0_dp .and. &
        coo_nd_get(b,[3,3]) == 4.0_dp, 'filter values preserved')
    b = coo_nd_filter_value(a,'==',0.0_dp,status)
    call require(status == sparse_invalid_input, 'filter invalid implicit-zero predicate')

    large = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([1999999999,2000000000, 2000000000,1],[2,2]), &
        [2.0_dp,-3.0_dp],status)
    call require(status == sparse_success, 'huge sparse shape')
    b = coo_nd_filter_points(large,reshape([1999999999,2000000000],[2,1]),status)
    call require(status == sparse_success .and. size(b%data) == 1, 'huge point filter')
    mask = coo_nd_compare_scalar(large,'<',0.0_dp,status)
    call require(status == sparse_success .and. coo_nd_get(mask,[2000000000,1]) == 1.0_dp, &
        'huge sparse predicate')

    empty = coo_nd_from_triplets([0,3],reshape([integer ::],[2,0]), &
        [real(dp) ::],status)
    call require(status == sparse_success, 'empty shape input')
    mask = coo_nd_compare_scalar(empty,'>',0.0_dp,status)
    call require(status == sparse_success .and. size(mask%data) == 0, 'empty mask')
    b = coo_nd_filter_points(empty,reshape([integer ::],[2,0]),status)
    call require(status == sparse_success .and. size(b%data) == 0, 'empty filter')

    bad = coo_nd_from_triplets([1],reshape([1,1],[1,2]), &
        [huge(1.0_dp)*0.9_dp,huge(1.0_dp)*0.9_dp],status)
    call require(status == sparse_success, 'finite overflowing duplicate values')
    n = coo_nd_count_nonzero(bad,status)
    call require(status == sparse_invalid_input .and. n == -1_int64, 'count rejects overflow')
    mask = coo_nd_compare_scalar(bad,'>',0.0_dp,status)
    call require(status == sparse_invalid_input, 'predicate rejects duplicate overflow')

    ! A 12,000-entry source and 7,000 coordinate selectors: sparse-sized
    ! allocation only, including 1,000 duplicate selector positions.
    allocate(stress_coords(2,12000),stress_values(12000),stress_points(2,7000))
    do k = 1, 12000
        stress_coords(:,k) = [k,mod(k,10000)+1]
    end do
    stress_values = 1.0_dp
    do k = 1, 6000
        stress_points(:,k) = stress_coords(:,2*k)
    end do
    stress_points(:,6001:7000) = stress_points(:,1:1000)
    large = coo_nd_from_triplets([2000000000,2000000000],stress_coords, &
        stress_values,status)
    call require(status == sparse_success, 'large coordinate source')
    b = coo_nd_filter_points(large,stress_points,status)
    call require(status == sparse_success .and. size(b%data) == 6000, &
        'stable large coordinate-set join')
    n = coo_nd_count_nonzero(b,status)
    call require(status == sparse_success .and. n == 6000_int64, &
        'large filtered distinct count')

    print '(a)', 'PASS test_sparse_coo_nd_filter' 
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Assertion expression.
        character(*), intent(in) :: message !! Assertion label.
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_filter
