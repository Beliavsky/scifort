! SPDX-License-Identifier: MIT
! General sparse block assembly and transactional scatter assignment tests.
program test_sparse_block_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_from_dense, csr_to_dense, &
        csr_eye, csr_bmat, csr_block, csr_set, csr_assign, csr_assign_coo, &
        csr_get, valid_csr, sparse_success, sparse_invalid_input
    implicit none
    type(csr_matrix) :: a, b, c, d, broken, target, sparse_source, blocks(2,3)
    real(dp) :: dense(3,4), expected(3,4)
    integer :: info, i

    a = csr_from_coo(2,2,[1,1,2],[1,2,2],[2.0_dp,0.0_dp,3.0_dp])
    b = csr_from_dense(reshape([4.0_dp,5.0_dp],[2,1]))
    c = csr_from_dense(reshape([6.0_dp,7.0_dp,8.0_dp],[1,3]))
    d = csr_from_dense(reshape([9.0_dp,10.0_dp,11.0_dp],[1,3]))
    blocks = csr_matrix()
    blocks(1,1) = a
    blocks(1,3) = b
    blocks(2,2) = c
    ! Need widths(2)=3, but first row lacks column 2.  Need block row 2
    ! height 1 and other occupied column 2 only; row 1 height 2.
    target = csr_bmat(blocks,info)
    call require(info == 0 .and. valid_csr(target), 'bmat absent inferred dimensions')
    call require(target%n_rows == 3 .and. target%n_cols == 6, 'bmat shape')
    call require(size(target%data) == 8, 'bmat preserves zero records')
    call require(csr_get(target,1,2) == 0.0_dp, 'bmat structural zero value')
    call require(any(target%indices(target%indptr(1):target%indptr(2)-1) == 2), 'bmat structural zero stored')
    deallocate_check: block
        real(dp) :: ref(3,6)
        ref = 0.0_dp
        ref(:2,:2) = csr_to_dense(a)
        ref(:2,6:6) = csr_to_dense(b)
        ref(3,3:5) = [6.0_dp,7.0_dp,8.0_dp]
        call equal(target,ref,'bmat values')
        target = csr_block(blocks,info)
        call equal(target,ref,'block alias')
    end block deallocate_check

    blocks(2,3) = d
    target = csr_bmat(blocks,info)
    call require(info == sparse_invalid_input, 'inconsistent block column widths')
    blocks(2,3) = csr_matrix()
    blocks(2,2) = csr_matrix()
    target = csr_bmat(blocks,info)
    call require(info == sparse_invalid_input, 'empty block row requires size')
    target = csr_bmat(blocks,info,row_sizes=[2,4],col_sizes=[2,3,1])
    call require(info == 0 .and. target%n_rows == 6 .and. target%n_cols == 6, 'sized all-empty row')
    call require(size(target%data) == 5, 'all-empty row no nnz')
    target = csr_bmat(blocks,info,row_sizes=[3,4],col_sizes=[2,3,1])
    call require(info == sparse_invalid_input, 'explicit height mismatch')
    target = csr_bmat(blocks,info,row_sizes=[2,4],col_sizes=[2,-1,1])
    call require(info == sparse_invalid_input, 'negative block size invalid')
    target = csr_bmat(blocks,info,row_sizes=[2],col_sizes=[2,3,1])
    call require(info == sparse_invalid_input, 'block size count mismatch')
    broken = csr_matrix()
    broken%n_rows = 2
    blocks(1,1) = broken
    target = csr_bmat(blocks,info,row_sizes=[2,1],col_sizes=[2,3,1])
    call require(info == sparse_invalid_input, 'malformed block is not absent')

    ! Entire block grid absent: explicit dimensions are required.
    blocks = csr_matrix()
    target = csr_bmat(blocks,info)
    call require(info == sparse_invalid_input, 'fully empty grid ambiguous')
    target = csr_bmat(blocks,info,row_sizes=[0,3],col_sizes=[4,0,2])
    call require(info == 0 .and. target%n_rows == 3 .and. target%n_cols == 6, 'all absent with given sizes')
    call require(size(target%data) == 0 .and. valid_csr(target), 'all absent canonical')
    target = csr_bmat(reshape([csr_matrix()],[1,1]),info,row_sizes=[1],col_sizes=[2])
    call require(info == 0 .and. target%n_rows == 1 .and. target%n_cols == 2, 'single absent sized')
    target = csr_bmat(reshape([csr_matrix()],[0,0]),info)
    call require(info == 0 .and. valid_csr(target), 'zero-size grid')
    target = csr_bmat(reshape([csr_matrix()],[0,1]),info)
    call require(info == sparse_invalid_input, 'degenerate nonsquare block grid')

    a = csr_from_coo(3,4,[1,1,2,3],[1,4,2,3],[1.0_dp,0.0_dp,2.0_dp,3.0_dp])
    call csr_assign_coo(a,[3,1,2,3,3],[2,1,2,2,4],[4.0_dp,7.0_dp,0.0_dp,8.0_dp,9.0_dp],info)
    call require(info == 0 .and. valid_csr(a), 'coordinate replacement and insert')
    dense = 0.0_dp
    dense(1,1) = 7.0_dp
    dense(2,2) = 0.0_dp
    dense(3,2) = 8.0_dp
    dense(3,3) = 3.0_dp
    dense(3,4) = 9.0_dp
    call equal(a,dense,'coordinate assignment values')
    call require(size(a%data) == 6, 'replaced zeros kept in structure')
    call csr_assign_coo(a,[1],[0],[42.0_dp],info)
    call require(info == sparse_invalid_input, 'bad col rejected')
    call equal(a,dense,'bad column leaves target unchanged')
    call csr_assign_coo(a,[1],[1],[ieee_value(0.0_dp,ieee_quiet_nan)],info)
    call require(info == sparse_invalid_input, 'NaN update rejected')
    call equal(a,dense,'NaN leaves target unchanged')
    call csr_assign_coo(a,[1,2],[1],[1.0_dp],info)
    call require(info == sparse_invalid_input, 'mismatched arrays rejected')
    call equal(a,dense,'mismatched arrays leave target unchanged')
    call csr_assign_coo(a,[integer::],[integer::],[real(dp)::],info)
    call require(info == 0, 'empty updates valid')

    call csr_set(a,1,3,0.0_dp,info)
    call require(info == 0 .and. size(a%data) == 7, 'set zero inserts stored zero')
    call csr_set(a,3,4,6.0_dp,info)
    call require(info == 0 .and. csr_get(a,3,4) == 6.0_dp, 'set existing element')
    call csr_set(a,4,1,0.0_dp,info)
    call require(info == sparse_invalid_input, 'set invalid index')

    a = csr_from_dense(reshape([(real(i,dp),i=1,12)],[3,4]))
    expected = csr_to_dense(a)
    call csr_assign(a,[3,1],[4,2],reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]),info)
    call require(info == 0, 'dense Cartesian assignment')
    expected(3,4) = 1.0_dp
    expected(1,4) = 2.0_dp
    expected(3,2) = 3.0_dp
    expected(1,2) = 4.0_dp
    call equal(a,expected,'dense assigned entries')
    call csr_assign(a,[2,2],[1,1],reshape([5.0_dp,6.0_dp,7.0_dp,8.0_dp],[2,2]),info)
    expected(2,1) = 8.0_dp
    call equal(a,expected,'repeated indices last-write-wins Fortran column major')
    call csr_assign(a,[1,3],[1,3],0.0_dp,info)
    expected(1,1) = 0.0_dp
    expected(3,1) = 0.0_dp
    expected(1,3) = 0.0_dp
    expected(3,3) = 0.0_dp
    call equal(a,expected,'scalar broadcast')
    sparse_source = csr_from_coo(2,2,[1],[2],[0.0_dp])
    call csr_assign(a,[2,3],[2,4],sparse_source,info)
    expected(2,2) = 0.0_dp
    expected(3,2) = 0.0_dp
    expected(2,4) = 0.0_dp
    expected(3,4) = 0.0_dp
    call equal(a,expected,'sparse source overwrites implicit zeros')
    call csr_assign(a,[1],[1],reshape([5.0_dp,6.0_dp],[2,1]),info)
    call require(info == sparse_invalid_input, 'dense shape mismatch')
    call equal(a,expected,'dense shape mismatch leaves unchanged')
    call csr_assign(a,[1],[1],sparse_source,info)
    call require(info == sparse_invalid_input, 'sparse shape mismatch')
    call equal(a,expected,'sparse shape mismatch leaves unchanged')
    call csr_assign(a,[integer::],[1],0.0_dp,info)
    call require(info == 0, 'empty Cartesian assignment accepted')
    call csr_assign(a,[5],[integer::],0.0_dp,info)
    call require(info == sparse_invalid_input, 'invalid row checked even for empty Cartesian product')
    call csr_assign(a,[integer::],[0],0.0_dp,info)
    call require(info == sparse_invalid_input, 'invalid column checked even for empty Cartesian product')
    call equal(a,expected,'empty selection validation leaves target unchanged')

    ! Huge logical column count: updates should use O(nnz+nupdate) space.
    a = csr_from_coo(2,1000000000,[1],[999999999],[11.0_dp])
    call csr_assign_coo(a,[2,1,2],[1000000000,1,1000000000], &
        [17.0_dp,19.0_dp,23.0_dp],info)
    call require(info == 0 .and. size(a%data) == 3, 'very wide sparse assign')
    call require(csr_get(a,2,1000000000) == 23.0_dp, 'very wide last-write-wins')

    a = csr_eye(30000)
    call csr_assign_coo(a,[30000,1],[1,30000],[5.0_dp,6.0_dp],info)
    call require(info == 0 .and. size(a%data) == 30002, 'large sparse assign')
    call require(valid_csr(a) .and. csr_get(a,30000,1) == 5.0_dp, 'large assignment is canonical')

    print '(a)', 'sparse bmat/scatter assignment native tests passed'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine equal(mat, reference, message)
        type(csr_matrix), intent(in) :: mat
        real(dp), intent(in) :: reference(:,:)
        character(*), intent(in) :: message
        real(dp), allocatable :: actual(:,:)
        call require(valid_csr(mat), message//' CSR validity')
        call require(mat%n_rows == size(reference,1) .and. mat%n_cols == size(reference,2), message//' shape')
        actual = csr_to_dense(mat)
        call require(all(actual == reference), message//' values')
    end subroutine equal
end program test_sparse_block_assignment
