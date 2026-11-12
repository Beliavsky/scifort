! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_bool
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_get
    use scifort_sparse_coo_nd_bool
    implicit none

    type(coo_nd_bool_mask) :: a, b, c, d, x, y, z, huge_mask, empty_mask
    type(coo_nd_array) :: data, filtered, converted
    integer :: info, da, db, i, j, op, n, k
    integer, allocatable :: stress_points(:,:)
    logical :: av, bv, truth, result
    integer(int64) :: count_true

    ! Exhaust the four uniform-background combinations and four operations.
    do da = 0, 1
        do db = 0, 1
            a = coo_nd_bool_from_points([3,4], &
                reshape([1,2, 3,4, 1,2],[2,3]),info,default_value=da == 1)
            call require(info == sparse_success, 'make left mask')
            b = coo_nd_bool_from_points([3,4], &
                reshape([2,1, 3,4, 2,1],[2,3]),info,default_value=db == 1)
            call require(info == sparse_success, 'make right mask')
            call require(valid_coo_nd_bool_mask(a) .and. valid_coo_nd_bool_mask(b), 'canonical masks')
            do op = 1, 4
                select case (op)
                case (1)
                    c = coo_nd_bool_and(a,b,info)
                case (2)
                    c = coo_nd_bool_or(a,b,info)
                case (3)
                    c = coo_nd_bool_xor(a,b,info)
                case (4)
                    c = coo_nd_bool_difference(a,b,info)
                end select
                call require(info == sparse_success, 'binary status')
                call require(valid_coo_nd_bool_mask(c), 'binary canonical')
                n = 0
                do i = 1, 3
                    do j = 1, 4
                        av = coo_nd_bool_get(a,[i,j],info)
                        call require(info == sparse_success, 'get left status')
                        bv = coo_nd_bool_get(b,[i,j],info)
                        call require(info == sparse_success, 'get right status')
                        select case (op)
                        case (1)
                            truth = av .and. bv
                        case (2)
                            truth = av .or. bv
                        case (3)
                            truth = av .neqv. bv
                        case (4)
                            truth = av .and. .not. bv
                        end select
                        result = coo_nd_bool_get(c,[i,j],info)
                        call require(info == sparse_success .and. (result .eqv. truth), &
                            'Boolean truth table')
                        if (truth) n = n+1
                    end do
                end do
                count_true = coo_nd_bool_count_true(c,info)
                call require(info == sparse_success .and. count_true == int(n,int64), &
                    'count true matches truth table')
            end do
            x = coo_nd_bool_not(a,info)
            call require(info == sparse_success, 'complement status')
            y = coo_nd_bool_not(x,info)
            call require(info == sparse_success .and. (y%default_value .eqv. a%default_value), &
                'complement involution background')
            call require(all(y%coords == a%coords), 'complement involution exceptions')
            ! De Morgan: !(a & b) == !a | !b at every position.
            c = coo_nd_bool_and(a,b,info)
            x = coo_nd_bool_not(c,info)
            y = coo_nd_bool_not(a,info)
            z = coo_nd_bool_not(b,info)
            d = coo_nd_bool_or(y,z,info)
            call require(info == sparse_success, 'De Morgan combination')
            call require(x%default_value .eqv. d%default_value, 'De Morgan background')
            call require(all(x%coords == d%coords), 'De Morgan coords')
        end do
    end do

    ! Sparse-valued input: duplicate sums, cancellation, explicit stored zeros.
    data = coo_nd_from_triplets([2,3],reshape([1,1, 1,1, 2,3, 1,3, 2,2],[2,5]), &
        [2.0_dp,-2.0_dp,4.0_dp,0.0_dp,-3.0_dp],info)
    call require(info == sparse_success, 'real COO construction')
    a = coo_nd_bool_from_coo(data,info)
    call require(info == sparse_success, 'real COO conversion')
    call require(size(a%coords,2) == 2, 'duplicate cancel and explicit zeros ignored')
    call require(.not. coo_nd_bool_get(a,[1,1]), 'canceled coordinate is false')
    call require(coo_nd_bool_get(a,[2,3]) .and. coo_nd_bool_get(a,[2,2]), &
        'true nonzero coordinates')
    b = coo_nd_bool_from_coo(data,info,complement=.true.)
    call require(info == sparse_success, 'complement conversion')
    call require(coo_nd_bool_get(b,[1,1]) .and. .not. coo_nd_bool_get(b,[2,3]), &
        'complement conversion values')
    converted = coo_nd_bool_to_coo(a,info)
    call require(info == sparse_success .and. size(converted%data) == 2, 'materialize sparse false-bg mask')
    converted = coo_nd_bool_to_coo(b,info)
    call require(info == sparse_invalid_input, 'reject materializing dense complement')

    c = coo_nd_bool_compare_scalar(data,'==',0.0_dp,info)
    call require(info == sparse_success .and. c%default_value, &
        'equality-to-zero implicit true supported')
    call require(coo_nd_bool_get(c,[1,1]) .and. .not. coo_nd_bool_get(c,[2,3]), &
        'equality to zero exceptions')
    call require(coo_nd_bool_count_true(c,info) == 4_int64 .and. &
        info == sparse_success, 'equality to zero count')
    c = coo_nd_bool_compare_scalar(data,'/=',1.0_dp,info)
    call require(info == sparse_success .and. c%default_value, &
        'inequality-to-nonzero implicit true supported')
    c = coo_nd_bool_compare_scalar(data,'abs>',0.5_dp,info)
    call require(info == sparse_success .and. .not. c%default_value .and. &
        size(c%coords,2) == 2, 'abs predicate sparse true exceptions')

    ! Complement-aware filtering retains original COO duplicates and zeros.
    c = coo_nd_bool_compare_scalar(data,'==',0.0_dp,info)
    filtered = coo_nd_bool_filter(data,c,info)
    call require(info == sparse_success .and. size(filtered%data) == 3, &
        'filter zero-valued coordinates retains all original stored entries')
    call require(all(filtered%data == [2.0_dp,-2.0_dp,0.0_dp]), &
        'filter retained duplicates and explicit zero in source order')
    call require(coo_nd_get(filtered,[1,1]) == 0.0_dp, 'filtered duplicate sum')

    ! Input-validation and output cap errors leave result unallocated.
    c = coo_nd_bool_and(a,a,info,max_exceptions=0_int64)
    call require(info == sparse_invalid_input .and. .not. allocated(c%shape), &
        'binary cap failure unchanged output')
    c = coo_nd_bool_from_points([2],reshape([3],[1,1]),info)
    call require(info == sparse_invalid_input, 'reject invalid point')
    c = coo_nd_bool_compare_scalar(data,'invalid',0.0_dp,info)
    call require(info == sparse_invalid_input, 'reject invalid comparison')
    result = coo_nd_bool_get(a,[0,1],info)
    call require(info == sparse_invalid_input .and. .not. result, 'reject invalid lookup')
    filtered = coo_nd_bool_filter(data,coo_nd_bool_from_points([2,4], &
        reshape([1,1],[2,1])),info)
    call require(info == sparse_invalid_input, 'reject shape mismatch')
    c = coo_nd_bool_and(a,coo_nd_bool_from_points([2,4], &
        reshape([1,1],[2,1])),info)
    call require(info == sparse_invalid_input, 'reject mismatched binary shapes')

    empty_mask = coo_nd_bool_from_points([0,5],reshape([integer ::],[2,0]), &
        info,default_value=.true.)
    call require(info == sparse_success, 'empty shape true default')
    call require(coo_nd_bool_count_true(empty_mask,info) == 0_int64 .and. &
        info == sparse_success, 'empty count true')
    converted = coo_nd_bool_to_coo(empty_mask,info)
    call require(info == sparse_success .and. size(converted%data) == 0, &
        'empty shape complement materialization')

    huge_mask = coo_nd_bool_from_points([2000000000,2000000000], &
        reshape([2000000000,2000000000],[2,1]),info,default_value=.true.)
    count_true = coo_nd_bool_count_true(huge_mask,info)
    call require(info == sparse_success .and. count_true == 3999999999999999999_int64, &
        'large true background count within int64')
    call require(.not. coo_nd_bool_get(huge_mask,[2000000000,2000000000]) .and. &
        coo_nd_bool_get(huge_mask,[2000000000,1]), 'large coordinate lookup')
    x = coo_nd_bool_from_points([2000000000,2000000000,2000000000], &
        reshape([1,1,1],[3,1]),info,default_value=.true.)
    count_true = coo_nd_bool_count_true(x,info)
    call require(info == sparse_invalid_input .and. count_true == -1_int64, &
        'reject true count overflow')
    y = coo_nd_bool_not(x,info)
    count_true = coo_nd_bool_count_true(y,info)
    call require(info == sparse_success .and. count_true == 1_int64, &
        'false background count needs no dense-size product')

    ! 12,000 points on enormous tensor extents; no dense allocation.
    allocate(stress_points(2,12000))
    do k = 1, 12000
        stress_points(:,k) = [k,mod(k,10000)+1]
    end do
    a = coo_nd_bool_from_points([2000000000,2000000000],stress_points,info)
    call require(info == sparse_success .and. size(a%coords,2) == 12000, &
        'large sparse set construction')
    b = coo_nd_bool_not(a,info)
    c = coo_nd_bool_and(a,b,info)
    call require(info == sparse_success .and. .not. c%default_value .and. &
        size(c%coords,2) == 0, 'large set and its complement empty')
    c = coo_nd_bool_or(a,b,info)
    call require(info == sparse_success .and. c%default_value .and. &
        size(c%coords,2) == 0, 'large set or its complement universal')

    print '(a)', 'PASS test_sparse_coo_nd_bool'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Assertion result.
        character(*), intent(in) :: message !! Assertion description.

        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_bool
