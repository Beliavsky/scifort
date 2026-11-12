! SPDX-License-Identifier: MIT
program test_sparse_assembly
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_to_dense, csr_eye, &
        valid_csr, sparse_success, sparse_invalid_input, csr_hstack, csr_vstack, &
        csr_block_diag, csr_kron, csr_get, csr_select, csr_slice
    implicit none
    type(csr_matrix) :: a, b, c, r, invalid, empty
    type(csr_matrix) :: blocks(2), none(0)
    real(dp), allocatable :: dense(:,:), expect(:,:)
    real(dp) :: x(2,3), y(2,2), val
    integer :: info

    x = reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,3])
    y = reshape([5.0_dp,0.0_dp,6.0_dp,7.0_dp],[2,2])
    a = csr_from_dense(x)
    b = csr_from_dense(y)
    blocks = [a,b]

    c = csr_hstack(blocks,info)
    call require(info == 0 .and. valid_csr(c), 'hstack valid')
    allocate(expect(2,5))
    expect(:,1:3) = x
    expect(:,4:5) = y
    call equal(c,expect,'hstack values')
    c = csr_hstack([a],info)
    call equal(c,x,'hstack singleton')
    c = csr_hstack(none,info)
    call require(info == sparse_invalid_input .and. .not. valid_csr(c), 'hstack empty invalid')
    c = csr_hstack([a,csr_from_dense(reshape([1.0_dp],[1,1]))],info)
    call require(info == sparse_invalid_input, 'hstack wrong shape')
    c = csr_hstack([a,invalid],info)
    call require(info == sparse_invalid_input, 'hstack invalid matrix')

    b = csr_from_dense(reshape([8.0_dp, 9.0_dp, 10.0_dp],[1,3]))
    c = csr_vstack([a,b],info)
    call require(info == 0, 'vstack valid')
    deallocate(expect)
    allocate(expect(3,3))
    expect(1:2,:) = x
    expect(3,:) = [8.0_dp,9.0_dp,10.0_dp]
    call equal(c,expect,'vstack values')
    c = csr_vstack([a],info)
    call equal(c,x,'vstack singleton')
    c = csr_vstack(none,info)
    call require(info == sparse_invalid_input, 'vstack empty invalid')
    c = csr_vstack([a,csr_from_dense(y)],info)
    call require(info == sparse_invalid_input, 'vstack wrong shape')

    b = csr_from_dense(y)
    c = csr_block_diag([a,b],info)
    call require(info == 0, 'block diag valid')
    deallocate(expect)
    allocate(expect(4,5))
    expect = 0.0_dp
    expect(:2,:3) = x
    expect(3:,4:) = y
    call equal(c,expect,'block diag values')
    c = csr_block_diag(none,info)
    call require(info == 0 .and. valid_csr(c), 'block diag empty valid')
    call require(c%n_rows == 0 .and. c%n_cols == 0, 'block diag empty shape')
    c = csr_block_diag([a,invalid],info)
    call require(info == sparse_invalid_input, 'block diag invalid matrix')

    c = csr_kron(a,b,info)
    call require(info == sparse_success .and. valid_csr(c), 'kron valid')
    deallocate(expect)
    allocate(expect(4,6))
    expect = 0.0_dp
    expect(:2,:2) = x(1,1)*y
    expect(3:4,:2) = x(2,1)*y
    expect(:2,3:4) = x(1,2)*y
    expect(3:4,3:4) = x(2,2)*y
    expect(:2,5:6) = x(1,3)*y
    expect(3:4,5:6) = x(2,3)*y
    call equal(c,expect,'kron values')
    c = csr_kron(a,invalid,info)
    call require(info == sparse_invalid_input, 'kron invalid input')
    b = csr_from_dense(reshape([ieee_value(0.0_dp,ieee_positive_inf)],[1,1]),info)
    call require(info == sparse_invalid_input, 'nonfinite constructor rejected')
    a = csr_from_coo(1,1,[1],[1],[huge(1.0_dp)])
    b = csr_from_coo(1,1,[1],[1],[2.0_dp])
    c = csr_kron(a,b,info)
    call require(info == sparse_invalid_input .and. .not. valid_csr(c), 'kron overflow rejected')

    a = csr_from_coo(3,4,[1,1,2,3,3],[1,3,2,2,4], &
        [1.0_dp,0.0_dp,2.0_dp,3.0_dp,4.0_dp])
    val = csr_get(a,1,3,info)
    call require(info == 0 .and. val == 0.0_dp, 'get stored zero')
    val = csr_get(a,3,4,info)
    call require(info == 0 .and. val == 4.0_dp, 'get nonzero')
    val = csr_get(a,2,3,info)
    call require(info == 0 .and. val == 0.0_dp, 'get missing zero')
    val = csr_get(a,4,3,info)
    call require(info == sparse_invalid_input, 'get invalid index')

    c = csr_select(a,[3,1,3],[4,2,2,1],info)
    call require(info == 0 .and. valid_csr(c), 'select reorder')
    deallocate(expect)
    allocate(expect(3,4))
    expect(1,:) = [4.0_dp,3.0_dp,3.0_dp,0.0_dp]
    expect(2,:) = [0.0_dp,0.0_dp,0.0_dp,1.0_dp]
    expect(3,:) = expect(1,:)
    call equal(c,expect,'select duplicate rows and columns')
    c = csr_select(a,[1],[3,3],info)
    call require(info == 0 .and. size(c%data) == 2, 'select preserves duplicated stored zeros')
    call require(all(c%data == 0.0_dp), 'select stored zeros')
    c = csr_select(a,[integer::],[1],info)
    call require(info == 0 .and. c%n_rows == 0 .and. c%n_cols == 1, 'select zero rows')
    c = csr_select(a,[1],[integer::],info)
    call require(info == 0 .and. c%n_rows == 1 .and. c%n_cols == 0, 'select zero columns')
    c = csr_select(a,[1,0],[1],info)
    call require(info == sparse_invalid_input, 'select invalid row')
    c = csr_select(a,[1],[5],info)
    call require(info == sparse_invalid_input, 'select invalid col')

    c = csr_slice(a,1,3,2,4,info)
    call require(info == 0, 'slice full positive')
    dense = csr_to_dense(a)
    call equal(c,dense(:,2:4),'slice submatrix')
    c = csr_slice(a,3,1,4,1,info,row_step=-1,col_step=-2)
    call require(info == 0, 'slice reverse')
    call equal(c,dense(3:1:-1,4:1:-2),'slice reversed strides')
    c = csr_slice(a,1,3,1,4,info,row_step=2,col_step=3)
    call equal(c,dense(1:3:2,1:4:3),'slice positive strides')
    c = csr_slice(a,4,3,1,4,info)
    call require(info == 0 .and. c%n_rows == 0, 'empty row slice')
    c = csr_slice(a,1,3,5,4,info)
    call require(info == 0 .and. c%n_cols == 0, 'empty col slice')
    c = csr_slice(a,1,3,1,4,info,row_step=0)
    call require(info == sparse_invalid_input, 'zero stride invalid')
    c = csr_slice(a,0,3,1,4,info)
    call require(info == sparse_invalid_input, 'invalid row bound')

    a = csr_from_coo(2,1000000000,[1,2],[999999999,1000000000],[7.0_dp,8.0_dp])
    c = csr_select(a,[2,1],[1000000000,999999999,1000000000],info)
    call require(info == 0 .and. c%n_cols == 3 .and. size(c%data) == 3, 'very wide sparse select')
    call require(csr_get(c,1,1) == 8.0_dp .and. csr_get(c,2,2) == 7.0_dp, &
        'very wide sparse select values')
    a = csr_from_coo(3,4,[1,1,2,3,3],[1,3,2,2,4], &
        [1.0_dp,0.0_dp,2.0_dp,3.0_dp,4.0_dp])

    empty = csr_from_coo(0,4,[integer::],[integer::],[real(dp)::])
    c = csr_slice(empty,1,0,1,4,info)
    call require(info == 0 .and. c%n_rows == 0 .and. c%n_cols == 4, 'empty source slice')
    c = csr_block_diag([empty,a],info)
    call require(info == 0 .and. c%n_rows == 3 .and. c%n_cols == 8, 'block diag zero rows')
    c = csr_kron(empty,a,info)
    call require(info == 0 .and. c%n_rows == 0 .and. c%n_cols == 16, 'kron zero rows')

    a = csr_eye(30000)
    c = csr_hstack([a,a],info)
    call require(info == 0 .and. c%n_rows == 30000 .and. c%n_cols == 60000, &
        'large sparse hstack shape')
    call require(size(c%data) == 60000 .and. valid_csr(c), 'large sparse hstack storage')
    c = csr_vstack([a,a],info)
    call require(info == 0 .and. c%n_rows == 60000 .and. c%n_cols == 30000, &
        'large sparse vstack shape')
    call require(size(c%data) == 60000 .and. valid_csr(c), 'large sparse vstack storage')
    c = csr_block_diag([a,a],info)
    call require(info == 0 .and. c%n_rows == 60000 .and. c%n_cols == 60000, &
        'large sparse direct sum shape')
    call require(size(c%data) == 60000 .and. valid_csr(c), 'large sparse direct sum storage')

    a = csr_from_coo(1,2,[1,1],[1,2],[0.0_dp,2.0_dp])
    b = csr_from_coo(1,2,[1,1],[1,2],[3.0_dp,0.0_dp])
    c = csr_kron(a,b,info)
    call require(info == 0 .and. size(c%data) == 4, 'kron preserves stored zeros')
    call require(all(c%data == [0.0_dp,0.0_dp,6.0_dp,0.0_dp]), &
        'kron structural zero values')

    print '(a)', 'sparse assembly/indexing native tests passed'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine equal(sparse, target, message)
        type(csr_matrix), intent(in) :: sparse
        real(dp), intent(in) :: target(:,:)
        character(*), intent(in) :: message
        real(dp), allocatable :: result(:,:)
        call require(valid_csr(sparse), message//' valid')
        call require(sparse%n_rows == size(target,1), message//' row shape')
        call require(sparse%n_cols == size(target,2), message//' column shape')
        result = csr_to_dense(sparse)
        call require(all(abs(result-target) < 1.0e-12_dp), message//' entries')
    end subroutine equal
end program test_sparse_assembly
