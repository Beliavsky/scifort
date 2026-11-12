! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_index
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, valid_coo_nd, &
        coo_nd_from_triplets, coo_nd_get, coo_nd_to_dense_flat
    use scifort_sparse_coo_nd_index, only : coo_nd_selector, coo_nd_all, &
        coo_nd_at, coo_nd_slice, coo_nd_indices, coo_nd_select, coo_nd_gather
    implicit none
    type(coo_nd_array) :: a, b, empty, large
    type(coo_nd_selector) :: selectors(3)
    real(dp), allocatable :: flat(:)
    integer :: info

    a = coo_nd_from_triplets([3,4,5], reshape([ &
        1,1,1, 3,2,4, 3,2,4, 2,4,5, 2,3,1, 1,4,3, 3,4,5], [3,7]), &
        [2.0_dp,7.0_dp,-2.0_dp,0.0_dp,9.0_dp,-3.0_dp,11.0_dp],info)
    call require(info == sparse_success, 'input COO')
    b = coo_nd_select(a,[coo_nd_at(3),coo_nd_slice(2,4,2),coo_nd_slice(5,2,-1)],info)
    call require(info == 0 .and. valid_coo_nd(b), 'integer and reversed slice')
    call require(all(b%shape == [2,4]) .and. size(b%data) == 3, 'sliced shape and stored duplicates')
    call require(coo_nd_get(b,[1,2]) == 5.0_dp, 'duplicates summed on lookup')
    call require(coo_nd_get(b,[2,1]) == 11.0_dp, 'negative stride coordinates')
    call require(all(b%data == [7.0_dp,-2.0_dp,11.0_dp]), 'source data order retained')

    b = coo_nd_select(a,[coo_nd_at(3),coo_nd_slice(2,4,2),coo_nd_slice(5,2,-1)], &
                      info,keepdims=.true.)
    call require(info == 0 .and. all(b%shape == [1,2,4]), 'keepdims scalar axis')
    call require(coo_nd_get(b,[1,1,2]) == 5.0_dp, 'keepdims positions')
    b = coo_nd_select(a,[coo_nd_at(3),coo_nd_at(2),coo_nd_at(4)],info)
    call require(info == 0 .and. all(b%shape == [1]), 'all scalar selectors use rank-1 singleton')
    call require(size(b%data) == 2 .and. coo_nd_get(b,[1]) == 5.0_dp, 'all scalar duplicates')

    b = coo_nd_select(a,[coo_nd_indices([3,1,3]),coo_nd_indices([2,2,4]), &
                         coo_nd_at(4)],info)
    call require(info == 0 .and. all(b%shape == [3,3]), 'orthogonal indexing shape')
    call require(size(b%data) == 8, 'repeated index vectors replicate entries')
    call require(coo_nd_get(b,[1,1]) == 5.0_dp .and. coo_nd_get(b,[3,2]) == 5.0_dp, &
                 'orthogonal duplicate positions')
    call require(coo_nd_get(b,[2,1]) == 0.0_dp, 'unselected coordinates remain zero')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_slice(stride=-1),coo_nd_at(5)],info)
    call require(info == 0 .and. all(b%shape == [3,4]), 'implicit reversed full slice')
    call require(size(b%data) == 2, 'stored explicit zero preserved')
    call require(any(b%data == 0.0_dp), 'explicit zero survived slice')
    call require(coo_nd_get(b,[3,1]) == 11.0_dp, 'reversed full slice indices')

    b = coo_nd_gather(a,reshape([ &
        3,2,4, 2,4,5, 1,1,1, 3,2,4, 1,2,2], [3,5]),info)
    call require(info == 0 .and. all(b%shape == [5]), 'pairwise gather shape')
    call require(size(b%data) == 6, 'gather retains explicit zeros and duplicates')
    call require(coo_nd_get(b,[1]) == 5.0_dp .and. coo_nd_get(b,[2]) == 0.0_dp, &
                 'gather duplicate sum and explicit zero')
    call require(coo_nd_get(b,[3]) == 2.0_dp .and. coo_nd_get(b,[4]) == 5.0_dp, &
                 'gather repeated queries')
    flat = coo_nd_to_dense_flat(b,info)
    call require(info == 0 .and. all(flat == [5.0_dp,0.0_dp,2.0_dp,5.0_dp,0.0_dp]), &
                 'gather dense output')

    b = coo_nd_select(a,[coo_nd_all(),coo_nd_slice(4,1),coo_nd_all()],info)
    call require(info == 0 .and. all(b%shape == [3,0,5]), 'empty positive-stride slice')
    call require(size(b%data) == 0, 'empty section has no entries')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_slice(1,4,-1),coo_nd_all()],info)
    call require(info == 0 .and. b%shape(2) == 0, 'empty reverse slice')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_indices([integer::]),coo_nd_all()],info)
    call require(info == 0 .and. all(b%shape == [3,0,5]), 'empty index vector')
    b = coo_nd_gather(a,reshape([integer::],[3,0]),info)
    call require(info == 0 .and. all(b%shape == [0]), 'empty pointwise gather')

    b = coo_nd_select(a,[coo_nd_at(0),coo_nd_all(),coo_nd_all()],info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'invalid scalar selector')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_indices([1,5]),coo_nd_all()],info)
    call require(info == sparse_invalid_input, 'invalid index array')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_slice(stride=0),coo_nd_all()],info)
    call require(info == sparse_invalid_input, 'zero stride rejected')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_slice(0,3),coo_nd_all()],info)
    call require(info == sparse_invalid_input, 'out of bounds slice rejected')
    b = coo_nd_select(a,[coo_nd_all(),coo_nd_all()],info)
    call require(info == sparse_invalid_input, 'selector rank mismatch')
    b = coo_nd_select(a,[coo_nd_indices([3,1,3]),coo_nd_indices([2,2,4]), &
                         coo_nd_at(4)],info,max_nnz=7_int64)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(b), 'selection output cap')
    b = coo_nd_gather(a,reshape([3,2,4],[3,1]),info,max_nnz=1_int64)
    call require(info == sparse_invalid_input, 'gather output cap')
    b = coo_nd_gather(a,reshape([4,1,1],[3,1]),info)
    call require(info == sparse_invalid_input, 'gather out-of-bounds')

    empty = coo_nd_from_triplets([0,4,2],reshape([integer::],[3,0]),[real(dp)::],info)
    selectors = [coo_nd_all(),coo_nd_slice(),coo_nd_all()]
    b = coo_nd_select(empty,selectors,info)
    call require(info == 0 .and. all(b%shape == [0,4,2]) .and. size(b%data) == 0, &
                 'zero-length source with full-axis selectors')
    b = coo_nd_select(empty,[coo_nd_slice(stride=-1),coo_nd_all(),coo_nd_all()],info)
    call require(info == 0 .and. all(b%shape == [0,4,2]), 'reversed empty axis')
    b = coo_nd_select(empty,[coo_nd_indices([integer::]),coo_nd_all(),coo_nd_all()],info)
    call require(info == 0 .and. size(b%data) == 0, 'empty index vector on empty source')

    large = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([2000000000,2000000000],[2,1]),[15.0_dp],info)
    b = coo_nd_select(large,[coo_nd_slice(1999999999,2000000000), &
                             coo_nd_at(2000000000)],info)
    call require(info == 0 .and. all(b%shape == [2]), 'large sparse slice without dense allocation')
    call require(coo_nd_get(b,[2]) == 15.0_dp, 'large sparse slice coordinate')
    b = coo_nd_select(large,[coo_nd_all(),coo_nd_all()],info)
    call require(info == 0 .and. all(b%shape == [2000000000,2000000000]), &
                 'large sparse all-axes view')
    call require(all(b%coords(:,1) == [2000000000,2000000000]), &
                 'large sparse coordinate retained')

    print *, 'sparse N-D COO indexing and slicing: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Condition to check.
        character(len=*), intent(in) :: message !! Test message.
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_index
