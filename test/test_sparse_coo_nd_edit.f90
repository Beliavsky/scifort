! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_edit
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_get, coo_nd_is_canonical
    use scifort_sparse_coo_nd_edit, only : coo_nd_mask_select, coo_nd_mask_axis, &
        coo_nd_scatter, coo_nd_mask_assign
    implicit none
    type(coo_nd_array) :: a, b, before, empty, large
    integer :: info, i
    integer, allocatable :: coords(:,:), points(:,:)
    real(dp), allocatable :: vals(:)
    real(dp) :: infinity

    a = coo_nd_from_triplets([2,3],reshape([ &
        2,3, 1,2, 2,3, 1,1, 2,1, 1,2], [2,6]), &
        [3.0_dp,5.0_dp,-1.0_dp,0.0_dp,4.0_dp,-2.0_dp],info)
    call require(info == sparse_success, 'input COO')
    b = coo_nd_mask_select(a,[.true.,.true.,.false.,.false.,.false.,.true.],info)
    call require(info == 0 .and. valid_coo_nd(b), 'boolean gather')
    call require(all(b%shape == [3]) .and. size(b%data) == 5, 'boolean gather duplicates and zero')
    call require(coo_nd_get(b,[1]) == 0.0_dp .and. coo_nd_get(b,[2]) == 3.0_dp .and. &
        coo_nd_get(b,[3]) == 2.0_dp, 'boolean gather values')
    call require(any(b%data == 0.0_dp), 'explicit zero retained by gather')
    b = coo_nd_mask_select(a,[logical ::],info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'invalid flat mask length')
    b = coo_nd_mask_select(a,[.true.,.true.,.false.,.false.,.false.,.true.],info,4_int64)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'mask nnz limit')
    b = coo_nd_mask_select(a,[.false.,.false.,.false.,.false.,.false.,.false.],info)
    call require(info == 0 .and. all(b%shape == [0]) .and. size(b%data) == 0, 'all false mask')

    b = coo_nd_mask_axis(a,1,[.false.,.true.],info)
    call require(info == 0 .and. all(b%shape == [1,3]), 'boolean axis selection')
    call require(size(b%data) == 3 .and. coo_nd_get(b,[1,3]) == 2.0_dp, 'axis values')
    b = coo_nd_mask_axis(a,2,[.true.,.false.,.true.],info)
    call require(info == 0 .and. all(b%shape == [2,2]), 'second axis shape')
    call require(coo_nd_get(b,[2,2]) == 2.0_dp, 'second axis coordinates')
    b = coo_nd_mask_axis(a,3,[.true.],info)
    call require(info == sparse_invalid_input, 'invalid axis')
    b = coo_nd_mask_axis(a,1,[.true.],info)
    call require(info == sparse_invalid_input, 'invalid axis mask size')

    before = a
    call coo_nd_scatter(a,reshape([2,3, 1,2, 2,3, 2,2],[2,4]), &
        [9.0_dp,0.0_dp,-2.0_dp,7.0_dp],info)
    call require(info == 0 .and. coo_nd_is_canonical(a), 'replace scatter canonical output')
    call require(size(a%data) == 4, 'replace deletes zero and coalesces source')
    call require(coo_nd_get(a,[2,3]) == -2.0_dp .and. coo_nd_get(a,[2,2]) == 7.0_dp, &
        'last-write-wins and insertion')
    call require(coo_nd_get(a,[1,2]) == 0.0_dp .and. &
        any(a%data == 0.0_dp), 'replacement deletes touched zero, retains untouched zero')
    call coo_nd_scatter(a,reshape([2,3, 2,1],[2,2]), [2.0_dp,-4.0_dp],info, &
        accumulate=.true.)
    call require(info == 0 .and. size(a%data) == 2, 'additive scatter cancellation removes entries')
    call require(coo_nd_get(a,[2,3]) == 0.0_dp .and. coo_nd_get(a,[2,1]) == 0.0_dp, &
        'additive scatter cancellation values')

    a = before
    call coo_nd_scatter(a,reshape([1,2,1,2],[2,2]),[4.0_dp,-7.0_dp],info, &
        accumulate=.true.,keep_zeros=.true.)
    call require(info == 0 .and. size(a%data) == 4, 'keep zero after additive cancellation')
    call require(coo_nd_get(a,[1,2]) == 0.0_dp, 'additive duplicate update order')
    call require(count(a%data == 0.0_dp) == 2, 'two retained explicit zeros')

    a = before
    call coo_nd_mask_assign(a,[.false.,.true.,.false.,.false.,.false.,.true.], &
        0.0_dp,info)
    call require(info == 0 .and. size(a%data) == 2, 'masked zero assignment deletes selected entries')
    call require(coo_nd_get(a,[2,3]) == 0.0_dp .and. coo_nd_get(a,[2,1]) == 4.0_dp, &
        'masked zero assignment values')
    a = before
    call coo_nd_mask_assign(a,[.false.,.true.,.false.,.false.,.false.,.true.], &
        0.0_dp,info,keep_zeros=.true.)
    call require(info == 0 .and. size(a%data) == 4, 'masked assignment keep zero entries')
    call require(coo_nd_is_canonical(a), 'masked assignment canonical COO')
    a = before
    call coo_nd_mask_assign(a,[.false.,.true.,.false.,.false.,.false.,.true.], &
        2.0_dp,info,accumulate=.true.)
    call require(info == 0 .and. coo_nd_get(a,[1,2]) == 5.0_dp .and. &
        coo_nd_get(a,[2,3]) == 4.0_dp, 'masked accumulate')

    ! Every failed edit must preserve the source arrays byte-for-byte.
    before = a
    call coo_nd_scatter(a,reshape([3,1],[2,1]),[1.0_dp],info)
    call unchanged(info,'invalid scatter index')
    call coo_nd_scatter(a,reshape([1,1],[2,1]),[1.0_dp,2.0_dp],info)
    call unchanged(info,'scatter point/value mismatch')
    infinity = ieee_value(0.0_dp,ieee_positive_inf)
    call coo_nd_scatter(a,reshape([1,1],[2,1]),[infinity],info)
    call unchanged(info,'nonfinite scatter value')
    call coo_nd_scatter(a,reshape([1,3],[2,1]),[5.0_dp],info,max_nnz=2_int64)
    call unchanged(info,'scatter output nnz limit')
    call coo_nd_mask_assign(a,[.true.,.false.],4.0_dp,info)
    call unchanged(info,'invalid mask assignment size')
    call coo_nd_mask_assign(a,[.true.,.true.,.true.,.true.,.true.,.true.], &
        1.0_dp,info,max_nnz=4_int64)
    call unchanged(info,'mask update cap')
    call coo_nd_mask_assign(a,[.true.,.false.,.false.,.false.,.false.,.false.], &
        infinity,info)
    call unchanged(info,'nonfinite mask assignment')
    call coo_nd_mask_assign(a,[.false.,.false.,.false.,.false.,.false.,.false.], &
        1.0_dp,info)
    call require(info == 0 .and. all(a%coords == before%coords) &
        .and. all(a%data == before%data), 'empty mask assignment is no-op')

    a = coo_nd_from_triplets([1],reshape([1],[1,1]),[huge(0.0_dp)],info)
    before = a
    call coo_nd_scatter(a,reshape([1],[1,1]),[huge(0.0_dp)],info,accumulate=.true.)
    call unchanged(info,'overflow on additive update')

    empty = coo_nd_from_triplets([0,4],reshape([integer ::],[2,0]),[real(dp) ::],info)
    call require(info == 0, 'empty sparse input')
    b = coo_nd_mask_select(empty,[logical ::],info)
    call require(info == 0 .and. all(b%shape == [0]), 'mask on empty shape')
    b = coo_nd_mask_axis(empty,1,[logical ::],info)
    call require(info == 0 .and. all(b%shape == [0,4]), 'empty axis mask')
    call coo_nd_mask_assign(empty,[logical ::],3.0_dp,info)
    call require(info == 0 .and. size(empty%data) == 0, 'empty assignment')

    large = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([2000000000,2000000000],[2,1]),[4.0_dp],info)
    call require(info == 0, 'large sparse input')
    call coo_nd_scatter(large,reshape([1,1],[2,1]),[2.0_dp],info)
    call require(info == 0 .and. size(large%data) == 2 .and. &
        coo_nd_get(large,[2000000000,2000000000]) == 4.0_dp, 'huge dimensions no dense allocation')
    b = coo_nd_mask_select(large,[.true.],info)
    call require(info == sparse_invalid_input, 'huge shape rejected with tiny mask')

    ! Stress the stable merge sort and repeated overwrites without dense work.
    allocate(points(2,12000),vals(12000),coords(2,0))
    do i = 1,12000
        points(:,i) = [1+mod(i*11,61),1+mod(i*13,59)]
        vals(i) = real(mod(i,17)-8,dp)
    end do
    a = coo_nd_from_triplets([61,59],coords,[real(dp)::],info)
    call coo_nd_scatter(a,points,vals,info,max_nnz=4000_int64)
    call require(info == 0 .and. coo_nd_is_canonical(a), '12k scatter events canonicalized')
    call require(size(a%data) <= 3599, 'scatter unique key bound')

    print '(a)', 'test_sparse_coo_nd_edit: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Test condition.
        character(*), intent(in) :: message !! Failure label.
        if (.not. ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine require

    subroutine unchanged(code, message)
        integer, intent(in) :: code !! Actual operation status.
        character(*), intent(in) :: message !! Failure label.
        call require(code == sparse_invalid_input,message//': status')
        call require(all(a%shape == before%shape),message//': shape')
        call require(all(a%coords == before%coords),message//': coordinates')
        call require(all(a%data == before%data),message//': values')
    end subroutine unchanged
end program test_sparse_coo_nd_edit
