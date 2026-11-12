! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_ops
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, &
        coo_nd_to_dense_flat, coo_nd_get, valid_coo_nd, coo_nd_is_canonical
    use scifort_sparse_coo_nd_ops, only : coo_nd_sum_all, coo_nd_sum_axes, &
        coo_nd_sum_axis, coo_nd_broadcast_to, coo_nd_add, coo_nd_multiply, coo_nd_tensordot
    implicit none
    type(coo_nd_array) :: a, b, c, d, v, t
    integer :: info, i, j, k, l
    real(dp) :: total, x, dense(2,3,2), matrix_b(2,4), reference(2,3,4)
    real(dp), allocatable :: flat(:)

    a = coo_nd_from_triplets([2,3,2],reshape([ &
        1,1,1, 1,1,1, 2,3,2, 2,3,1, 1,2,2],[3,5]), &
        [1.0_dp,2.0_dp,4.0_dp,-1.0_dp,0.0_dp],info)
    call require(info == 0,'construct test tensor')
    total = coo_nd_sum_all(a,info)
    call require(info == 0 .and. total == 6.0_dp,'sum all, duplicates and zeros')
    b = coo_nd_sum_axis(a,1,status=info)
    call require(info == 0 .and. all(b%shape == [3,2]),'single axis output shape')
    call require(coo_nd_is_canonical(b),'single axis canonical')
    call require(coo_nd_get(b,[1,1]) == 3.0_dp, 'single axis aggregate')
    call require(coo_nd_get(b,[3,2]) == 4.0_dp, 'single axis aggregate second')
    b = coo_nd_sum_axes(a,[1,3],status=info)
    call require(info == 0 .and. all(b%shape == [3]),'two axes output')
    call require(coo_nd_get(b,[1]) == 3.0_dp, 'two axes value')
    call require(coo_nd_get(b,[3]) == 3.0_dp, 'two axes cancellation')
    b = coo_nd_sum_axes(a,[3,1],keepdims=.true.,status=info)
    call require(info == 0 .and. all(b%shape == [1,3,1]), 'keepdims shape')
    call require(coo_nd_get(b,[1,1,1]) == 3.0_dp, 'keepdims value')
    b = coo_nd_sum_axes(a,[1,2,3],keepdims=.true.,status=info)
    call require(info == 0 .and. all(b%shape == [1,1,1]), 'full keepdims')
    call require(coo_nd_get(b,[1,1,1]) == total,'full keepdims value')
    b = coo_nd_sum_axes(a,[1,2,3],status=info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'scalar requires scalar function')
    b = coo_nd_sum_axes(a,[2,2],status=info)
    call require(info == sparse_invalid_input,'repeated reduced axis rejected')
    b = coo_nd_sum_axis(a,4,status=info)
    call require(info == sparse_invalid_input,'out of bounds reduced axis')
    b = coo_nd_sum_axes(a,[integer::],status=info)
    call require(info == 0 .and. size(b%data) == 4, 'identity reduction coalesces duplicates')

    v = coo_nd_from_triplets([1,3,1],reshape([1,2,1, 1,3,1],[3,2]), &
                            [2.0_dp,-1.0_dp],info)
    b = coo_nd_broadcast_to(v,[2,3,4],info)
    call require(info == 0 .and. all(b%shape == [2,3,4]),'broadcast target shape')
    call require(size(b%data) == 16, 'broadcast storage count')
    call require(coo_nd_get(b,[1,2,3]) == 2.0_dp .and. &
                 coo_nd_get(b,[2,3,4]) == -1.0_dp, 'broadcast values')
    call require(coo_nd_get(b,[2,1,4]) == 0.0_dp,'broadcast implicit zeros')
    b = coo_nd_broadcast_to(v,[2,3,4],info,max_nnz=15_int64)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'broadcast cap')
    b = coo_nd_broadcast_to(v,[2,2,4],info)
    call require(info == sparse_invalid_input, 'broadcast mismatch')
    b = coo_nd_broadcast_to(v,[1,3],info)
    call require(info == sparse_invalid_input,'broadcast cannot reduce rank')
    b = coo_nd_broadcast_to(v,[0,3,4],info)
    call require(info == 0 .and. size(b%data) == 0, 'zero dimension broadcasting')
    b = coo_nd_broadcast_to(v,[1000000000,3,1000000000],info)
    call require(info == sparse_invalid_input, 'guard huge replication')

    c = coo_nd_from_triplets([2,3],reshape([1,2, 1,2, 2,3, 1,2],[2,4]), &
        [1.0_dp,2.0_dp,4.0_dp,-1.0_dp],info)
    d = coo_nd_from_triplets([3],reshape([2,3],[1,2]),[3.0_dp,-2.0_dp],info)
    b = coo_nd_add(c,d,info)
    call require(info == 0 .and. all(b%shape == [2,3]), 'broadcast addition shape')
    call require(coo_nd_is_canonical(b), 'broadcast addition canonical')
    call require(coo_nd_get(b,[1,2]) == 5.0_dp .and. coo_nd_get(b,[2,2]) == 3.0_dp, &
                 'broadcast addition values')
    call require(coo_nd_get(b,[2,3]) == 2.0_dp .and. coo_nd_get(b,[1,3]) == -2.0_dp, &
                 'broadcast addition filled zeros')
    b = coo_nd_multiply(c,d,info)
    call require(info == 0 .and. all(b%shape == [2,3]),'multiply shape')
    call require(coo_nd_is_canonical(b) .and. size(b%data) == 2,'sparse intersection')
    call require(coo_nd_get(b,[1,2]) == 6.0_dp .and. coo_nd_get(b,[2,3]) == -8.0_dp, &
                 'broadcast multiply values')
    b = coo_nd_add(c,d,info,max_nnz=5_int64)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'binary add cap')
    b = coo_nd_multiply(c,d,info,max_nnz=3_int64)
    call require(info == sparse_invalid_input, 'binary multiply cap')

    ! Contract the last axis of a 3D tensor with the first axis of a matrix.
    matrix_b = 0.0_dp
    matrix_b(1,2) = 2.0_dp
    matrix_b(1,4) = -3.0_dp
    matrix_b(2,1) = 1.0_dp
    matrix_b(2,2) = 4.0_dp
    matrix_b(2,3) = -2.0_dp
    t = coo_nd_from_triplets([2,4],reshape([1,2, 1,4, 2,1, 2,2, 2,3],[2,5]), &
                            [2.0_dp,-3.0_dp,1.0_dp,4.0_dp,-2.0_dp],info)
    b = coo_nd_tensordot(a,t,[3],[1],info)
    call require(info == 0 .and. all(b%shape == [2,3,4]), 'tensordot 3D by 2D shape')
    call require(coo_nd_is_canonical(b), 'tensordot canonical result')
    dense = 0.0_dp
    dense(1,1,1) = 3.0_dp
    dense(2,3,2) = 4.0_dp
    dense(2,3,1) = -1.0_dp
    reference = 0.0_dp
    do i = 1,2
        do j = 1,3
            do k = 1,4
                do l = 1,2
                    reference(i,j,k) = reference(i,j,k) + dense(i,j,l)*matrix_b(l,k)
                end do
            end do
        end do
    end do
    do i = 1,2
        do j = 1,3
            do k = 1,4
                x = coo_nd_get(b,[i,j,k],info)
                call require(info == 0 .and. x == reference(i,j,k),'tensor contraction numeric')
            end do
        end do
    end do
    b = coo_nd_tensordot(c,d,[2],[1],info)
    call require(info == 0 .and. all(b%shape == [2]), 'full contraction one free axis')
    call require(coo_nd_get(b,[1]) == 6.0_dp .and. coo_nd_get(b,[2]) == -8.0_dp, &
                 'vector contraction numeric')
    b = coo_nd_tensordot(d,d,[1],[1],info)
    call require(info == 0 .and. all(b%shape == [1]), 'fully contracted scalar shape')
    call require(coo_nd_get(b,[1]) == 13.0_dp, 'fully contracted scalar value')
    b = coo_nd_tensordot(d,d,[integer::],[integer::],info)
    call require(info == 0 .and. all(b%shape == [3,3]) .and. size(b%data) == 4, &
                 'sparse outer product')
    b = coo_nd_tensordot(d,d,[integer::],[integer::],info,max_nnz=3_int64)
    call require(info == sparse_invalid_input, 'tensor contraction output cap')
    b = coo_nd_tensordot(a,t,[2],[1],info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), &
                 'contracted axis dimension mismatch')
    b = coo_nd_tensordot(a,t,[3,3],[1,1],info)
    call require(info == sparse_invalid_input, 'duplicate contracted axis')

    c = coo_nd_from_triplets([1],reshape([1,1],[1,2]),[huge(1.0_dp),huge(1.0_dp)],info)
    total = coo_nd_sum_all(c,info)
    call require(info == sparse_invalid_input .and. total == 0.0_dp, 'sum overflow detected')
    b = coo_nd_sum_axes(c,[1],.true.,info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'reduction overflow detected')
    d = coo_nd_from_triplets([1],reshape([1],[1,1]),[2.0_dp],info)
    b = coo_nd_multiply(c,d,info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'multiply overflow detected')
    b = coo_nd_tensordot(c,d,[integer::],[integer::],info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'contract duplicate overflow')

    c = coo_nd_from_triplets([0,2],reshape([integer::],[2,0]),[real(dp)::],info)
    b = coo_nd_sum_axes(c,[2],status=info)
    call require(info == 0 .and. all(b%shape == [0]) .and. size(b%data) == 0, 'empty reduction')
    b = coo_nd_add(c,c,info)
    call require(info == 0 .and. all(b%shape == [0,2]) .and. size(b%data) == 0, 'empty sum')
    b = coo_nd_multiply(c,c,info)
    call require(info == 0 .and. size(b%data) == 0, 'empty multiply')
    b = coo_nd_tensordot(c,c,[2],[2],info)
    call require(info == 0 .and. all(b%shape == [0,0]) .and. size(b%data) == 0, 'empty contraction')
    flat = coo_nd_to_dense_flat(b,info)
    call require(info == 0 .and. size(flat) == 0, 'empty dense conversion')

    ! Tensor contractions on enormous sparse shapes must not compute dense
    ! shape products or allocate storage proportional to the dense output.
    c = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([2000000000,1],[2,1]),[5.0_dp],info)
    d = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([1,2000000000],[2,1]),[6.0_dp],info)
    b = coo_nd_tensordot(c,d,[2],[1],info)
    call require(info == 0 .and. all(b%shape == [2000000000,2000000000]), &
                 'contract very large sparse shape')
    call require(size(b%data) == 1 .and. b%data(1) == 30.0_dp, &
                 'large sparse contraction nonzero')
    call require(all(b%coords(:,1) == [2000000000,2000000000]), &
                 'large sparse contraction coordinates')

    print *, 'sparse N-D COO reductions, broadcasting and tensor contraction: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Condition to check.
        character(len=*), intent(in) :: message !! Test explanation.
        if (.not. ok) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_ops
