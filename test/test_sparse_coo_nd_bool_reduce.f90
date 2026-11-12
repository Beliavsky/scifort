! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_bool_reduce
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd_bool, only : coo_nd_bool_mask, coo_nd_bool_from_points, &
        coo_nd_bool_get, valid_coo_nd_bool_mask, coo_nd_bool_not
    use scifort_sparse_coo_nd_bool_reduce
    implicit none

    type(coo_nd_bool_mask) :: a, b, c, d, empty, huge_mask
    type(coo_nd_bool_count_array) :: counts
    integer :: info, i, j
    logical :: value
    integer(int64) :: n

    a = coo_nd_bool_from_points([2,1],reshape([1,1],[2,1]),info)
    call require(info == sparse_success,'construct source')
    c = coo_nd_bool_broadcast_to(a,[2,3],info)
    call require(info == sparse_success .and. valid_coo_nd_bool_mask(c),'broadcast to 2x3')
    call require(all(c%shape == [2,3]) .and. size(c%coords,2) == 3,'broadcast exceptions')
    do j = 1, 3
        call require(coo_nd_bool_get(c,[1,j]), 'broadcast first row')
        call require(.not. coo_nd_bool_get(c,[2,j]), 'broadcast second row')
    end do
    b = coo_nd_bool_from_points([1,3],reshape([1,2],[2,1]),info)
    c = coo_nd_bool_and_broadcast(a,b,info)
    call require(info == sparse_success .and. size(c%coords,2) == 1,'AND broadcast')
    call require(coo_nd_bool_get(c,[1,2]),'AND selected element')
    c = coo_nd_bool_or_broadcast(a,b,info)
    call require(info == sparse_success .and. size(c%coords,2) == 4,'OR broadcast')
    c = coo_nd_bool_xor_broadcast(a,b,info)
    call require(info == sparse_success .and. size(c%coords,2) == 3,'XOR broadcast')
    call require(.not. coo_nd_bool_get(c,[1,2]) .and. coo_nd_bool_get(c,[2,2]), &
        'XOR broadcast membership')
    c = coo_nd_bool_difference_broadcast(a,b,info)
    call require(info == sparse_success .and. size(c%coords,2) == 2,'difference broadcast')
    d = coo_nd_bool_not(b,info)
    c = coo_nd_bool_and_broadcast(a,d,info)
    call require(info == sparse_success .and. size(c%coords,2) == 2, &
        'complement-aware broadcast binary')

    c = coo_nd_bool_broadcast_to(b,[4,3],info)
    call require(info == sparse_success .and. size(c%coords,2) == 4, &
        'expand first singleton dimension')
    c = coo_nd_bool_broadcast_to(b,[2,4,3],info)
    call require(info == sparse_success .and. size(c%coords,2) == 8, &
        'expand missing leading dimension')
    c = coo_nd_bool_broadcast_to(b,[0,3],info)
    call require(info == sparse_success .and. size(c%coords,2) == 0, &
        'broadcast singleton to empty dimension')
    call require(.not. coo_nd_bool_any(c,info) .and. info == sparse_success, &
        'empty broadcast ANY false')
    call require(coo_nd_bool_all(c,info) .and. info == sparse_success, &
        'empty broadcast ALL true')

    a = coo_nd_bool_from_points([2,3,4],reshape([1,1,1, 2,1,1, 2,2,1, &
        2,2,2, 2,2,3, 2,2,4],[3,6]),info)
    counts = coo_nd_bool_count_axis(a,3,info)
    call require(info == sparse_success .and. valid_coo_nd_bool_count_array(counts), &
        'axis count canonical')
    call require(all(counts%shape == [2,3]) .and. counts%default_count == 0_int64, &
        'axis count shape and background')
    call require(coo_nd_bool_count_get(counts,[1,1],info) == 1_int64, &
        'axis count one exception')
    call require(coo_nd_bool_count_get(counts,[2,2],info) == 4_int64, &
        'axis count four exceptions')
    call require(coo_nd_bool_count_get(counts,[1,3],info) == 0_int64, &
        'axis count implicit zero')
    c = coo_nd_bool_any_axis(a,3,info)
    call require(info == sparse_success .and. size(c%coords,2) == 3, &
        'axis ANY with false background')
    c = coo_nd_bool_all_axis(a,3,info)
    call require(info == sparse_success .and. size(c%coords,2) == 1, &
        'axis ALL with false background')
    call require(coo_nd_bool_get(c,[2,2]),'axis ALL full coverage')
    d = coo_nd_bool_not(a,info)
    counts = coo_nd_bool_count_axis(d,3,info,keepdims=.true.)
    call require(info == sparse_success .and. all(counts%shape == [2,3,1]), &
        'axis count keepdims')
    call require(counts%default_count == 4_int64, 'axis count true background')
    call require(coo_nd_bool_count_get(counts,[2,2,1],info) == 0_int64, &
        'axis count complement full exception')
    c = coo_nd_bool_any_axis(d,3,info)
    call require(info == sparse_success .and. c%default_value .and. &
        size(c%coords,2) == 1, 'ANY complement with fully excluded row')
    c = coo_nd_bool_all_axis(d,3,info)
    call require(info == sparse_success .and. c%default_value .and. &
        size(c%coords,2) == 3, 'ALL complement with any exception')

    counts = coo_nd_bool_count_axis(a,1,info)
    call require(info == sparse_success .and. all(counts%shape == [3,4]), &
        'reduce leading axis re-sorts projected coordinate set')
    call require(coo_nd_bool_count_get(counts,[2,1]) == 1_int64, &
        'projected axis leading aggregation')
    counts = coo_nd_bool_count_axis(a,2,info)
    call require(info == sparse_success .and. all(counts%shape == [2,4]), &
        'reduce middle axis')
    call require(coo_nd_bool_count_get(counts,[2,1]) == 2_int64, &
        'projected middle-axis aggregation')

    b = coo_nd_bool_from_points([4],reshape([1,2,3,4],[1,4]),info)
    counts = coo_nd_bool_count_axis(b,1,info)
    call require(info == sparse_success .and. all(counts%shape == [1]), &
        'rank one reduction represented by length-one array')
    call require(coo_nd_bool_count_get(counts,[1]) == 4_int64,'rank-one count')
    call require(coo_nd_bool_all(b,info) .and. info == sparse_success, 'full ALL true')
    c = coo_nd_bool_not(b,info)
    call require(.not. coo_nd_bool_any(c,info) .and. info == sparse_success, &
        'full ANY false after complement')

    empty = coo_nd_bool_from_points([3,0,2],reshape([integer ::],[3,0]), &
        info,default_value=.true.)
    counts = coo_nd_bool_count_axis(empty,2,info)
    call require(info == sparse_success .and. counts%default_count == 0_int64, &
        'empty reduced axis count zero')
    c = coo_nd_bool_any_axis(empty,2,info)
    call require(info == sparse_success .and. .not. c%default_value, &
        'empty reduction ANY false')
    c = coo_nd_bool_all_axis(empty,2,info)
    call require(info == sparse_success .and. c%default_value, &
        'empty reduction ALL vacuously true')

    huge_mask = coo_nd_bool_from_points([2000000000,2000000000,2000000000], &
        reshape([1,1,1],[3,1]),info,default_value=.true.)
    counts = coo_nd_bool_count_axis(huge_mask,3,info)
    call require(info == sparse_success .and. counts%default_count == 2000000000_int64, &
        'large implicit-true baseline without total-volume multiplication')
    call require(coo_nd_bool_count_get(counts,[1,1]) == 1999999999_int64, &
        'large exact count value')
    c = coo_nd_bool_any_axis(huge_mask,3,info)
    call require(info == sparse_success .and. c%default_value .and. &
        size(c%coords,2) == 0, 'large ANY remains uniformly true')
    c = coo_nd_bool_all_axis(huge_mask,3,info)
    call require(info == sparse_success .and. c%default_value .and. &
        size(c%coords,2) == 1, 'large ALL has one exception')
    call require(coo_nd_bool_any(huge_mask,info) .and. info == sparse_success, &
        'huge full ANY avoids volume overflow')
    call require(.not. coo_nd_bool_all(huge_mask,info) .and. info == sparse_success, &
        'huge full ALL avoids volume overflow')

    c = coo_nd_bool_broadcast_to(b,[2,5],info)
    call require(info == sparse_invalid_input .and. .not. allocated(c%shape), &
        'reject incompatible broadcast')
    c = coo_nd_bool_broadcast_to(a,[2,3,4],info,max_exceptions=2_int64)
    call require(info == sparse_invalid_input .and. .not. allocated(c%shape), &
        'enforce broadcast exception cap')
    c = coo_nd_bool_any_axis(a,4,info)
    call require(info == sparse_invalid_input .and. .not. allocated(c%shape), &
        'reject invalid reduction axis')
    counts = coo_nd_bool_count_axis(a,0,info)
    call require(info == sparse_invalid_input .and. .not. allocated(counts%shape), &
        'reject invalid count axis')
    d = coo_nd_bool_from_points([5],reshape([1],[1,1]),info)
    c = coo_nd_bool_and_broadcast(a,d,info)
    call require(info == sparse_invalid_input, 'reject incompatible binary broadcast')
    c = coo_nd_bool_broadcast_to( &
        coo_nd_bool_from_points([1,1],reshape([1,1],[2,1])), &
        [2000000000,2000000000],info,max_exceptions=1000000_int64)
    call require(info == sparse_invalid_input .and. .not. allocated(c%shape), &
        'prevent enormous broadcast expansion')

    ! 12,000 points, count along a huge axis without dense masks.
    block
        integer, allocatable :: points(:,:)
        allocate(points(2,12000))
        do i = 1, 12000
            points(:,i) = [i,modulo(i,200)+1]
        end do
        a = coo_nd_bool_from_points([2000000000,2000000000],points,info)
        call require(info == sparse_success,'stress input')
        counts = coo_nd_bool_count_axis(a,1,info)
        call require(info == sparse_success .and. size(counts%values) == 200, &
            'large projection groups 12,000 exceptions')
        do j = 1, 200
            n = coo_nd_bool_count_get(counts,[j],info)
            call require(info == sparse_success .and. n == 60_int64, &
                'large projection exact count')
        end do
    end block

    print '(a)', 'PASS test_sparse_coo_nd_bool_reduce'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Assertion result.
        character(*), intent(in) :: message !! Assertion description.
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_bool_reduce
