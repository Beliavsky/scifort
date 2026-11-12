! SPDX-License-Identifier: MIT
program test_sparse_formats_norm
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, dia_matrix, bsr_matrix, &
        csr_from_dense, csr_from_coo, csr_to_dense, csr_to_csc, &
        csr_to_dia, dia_to_csr, dia_to_dense, dia_from_data, dia_matvec, valid_dia, &
        csr_to_bsr, bsr_to_csr, bsr_from_dense, bsr_to_dense, bsr_matvec, valid_bsr, &
        sparse_success, sparse_invalid_input
    use scifort_sparse_linalg, only : norm, norm_axis, sparse_linalg_success, sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a, r, empty, mat
    type(csc_matrix) :: c
    type(dia_matrix) :: d, d2
    type(bsr_matrix) :: b
    real(dp), allocatable :: dense(:,:), vec(:), cols(:), rows(:)
    real(dp) :: q(4,6), x(6), infp, infn, val, diagdata(3,6)
    integer :: i, j, info

    q = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, &
        0.0_dp, 3.0_dp, 4.0_dp, 0.0_dp, &
        5.0_dp, 0.0_dp, 0.0_dp, 6.0_dp, &
        0.0_dp, 0.0_dp, 7.0_dp, 0.0_dp, &
        0.0_dp, 8.0_dp, 0.0_dp, 0.0_dp, &
        9.0_dp, 0.0_dp, 0.0_dp, 10.0_dp], [4,6])
    x = [(real(i,dp),i=1,6)]
    a = csr_from_dense(q)
    c = csr_to_csc(a)
    d = csr_to_dia(a,info)
    call require(info == sparse_success .and. valid_dia(d), 'DIA conversion')
    call require(size(d%offsets) > 0, 'DIA diagonals present')
    dense = dia_to_dense(d,info)
    call require(info == 0 .and. maxval(abs(dense-q)) < 1.0e-14_dp, 'DIA to dense')
    r = dia_to_csr(d,info)
    call require(info == 0 .and. maxval(abs(csr_to_dense(r)-q)) < 1.0e-14_dp, 'DIA to CSR')
    vec = dia_matvec(d,x,info)
    call require(info == 0 .and. maxval(abs(vec-matmul(q,x))) < 1.0e-14_dp, 'DIA matvec')
    b = csr_to_bsr(a,2,3,info)
    call require(info == sparse_success .and. valid_bsr(b), 'BSR conversion')
    dense = bsr_to_dense(b,info)
    call require(info == 0 .and. maxval(abs(dense-q)) < 1.0e-14_dp, 'BSR to dense')
    r = bsr_to_csr(b,info)
    call require(info == 0 .and. maxval(abs(csr_to_dense(r)-q)) < 1.0e-14_dp, 'BSR to CSR')
    vec = bsr_matvec(b,x,info)
    call require(info == 0 .and. maxval(abs(vec-matmul(q,x))) < 1.0e-14_dp, 'BSR matvec')
    b = bsr_from_dense(q,2,3,info)
    call require(info == 0 .and. valid_bsr(b), 'BSR dense constructor')

    val = norm(a,info=info)
    call require(info == 0 .and. abs(val-sqrt(sum(q*q))) < 1.0e-12_dp, 'CSR Fro norm')
    val = norm(c,ord='fro',info=info)
    call require(info == 0 .and. abs(val-sqrt(sum(q*q))) < 1.0e-12_dp, 'CSC Fro norm')
    val = norm(a,1,info)
    call require(info == 0 .and. abs(val-maxval(sum(abs(q),dim=1))) < 1.0e-12_dp, 'matrix 1 norm')
    val = norm(c,-1,info)
    call require(info == 0 .and. abs(val-minval(sum(abs(q),dim=1))) < 1.0e-12_dp, 'matrix -1 norm')
    infp = ieee_value(0.0_dp,ieee_positive_inf)
    infn = ieee_value(0.0_dp,ieee_negative_inf)
    val = norm(a,infp,info)
    call require(info == 0 .and. abs(val-maxval(sum(abs(q),dim=2))) < 1.0e-12_dp, 'matrix inf norm')
    val = norm(c,infn,info)
    call require(info == 0 .and. abs(val-minval(sum(abs(q),dim=2))) < 1.0e-12_dp, 'matrix -inf norm')
    val = norm(a,2,info)
    call require(info == 0 .and. val >= 0.0_dp .and. val <= sqrt(sum(q*q))+1.0e-12_dp, 'spectral 2 norm')
    cols = norm_axis(a,1,info=info)
    call require(info == 0 .and. maxval(abs(cols-sqrt(sum(q*q,dim=1)))) < 1.0e-12_dp, 'column 2 norms')
    rows = norm_axis(c,2,ord=1.0_dp,info=info)
    call require(info == 0 .and. maxval(abs(rows-sum(abs(q),dim=2))) < 1.0e-12_dp, 'row 1 norms')
    cols = norm_axis(a,1,ord=0,info=info)
    call require(info == 0 .and. all(cols == real(count(q /= 0.0_dp,dim=1),dp)), 'column zero norms')
    rows = norm_axis(a,2,ord=infn,info=info)
    call require(info == 0 .and. all(rows == 0.0_dp), 'row -infinity norms')
    cols = norm_axis(a,1,ord=infp,info=info)
    call require(info == 0 .and. maxval(abs(cols-maxval(abs(q),dim=1))) == 0.0_dp, 'column infinity norms')
    cols = norm_axis(a,1,ord=-1,info=info)
    call require(info == 0 .and. all(cols == 0.0_dp), 'negative vector order with zeros')
    rows = norm_axis(a,3,info=info)
    call require(info == sparse_linalg_invalid_input .and. size(rows) == 0, 'axis validation')
    val = norm(a,ord=ieee_value(0.0_dp,ieee_quiet_nan),info=info)
    call require(info == sparse_linalg_invalid_input, 'NaN norm order rejected')
    val = norm(a,ord='nuclear',info=info)
    call require(info == sparse_linalg_invalid_input, 'unsupported matrix norm')
    val = norm(a,ord=0,info=info)
    call require(info == sparse_linalg_invalid_input, 'matrix zero norm disallowed')

    ! Duplicate coordinates coalesce; explicit CSR zeros become DIA padding.
    a = csr_from_coo(4,6,[1,1,4,4],[2,2,3,3],[3.0_dp,-3.0_dp,2.0_dp,1.0_dp])
    d = csr_to_dia(a)
    call require(maxval(abs(dia_to_dense(d)-csr_to_dense(a))) == 0.0_dp, 'coalesced DIA')
    b = csr_to_bsr(a,2,2)
    call require(valid_bsr(b), 'stored zero BSR')
    call require(maxval(abs(bsr_to_dense(b)-csr_to_dense(a))) == 0.0_dp, 'stored zero BSR values')

    ! Explicit SciPy-aligned diagonal data and out-of-domain padding.
    diagdata = 0.0_dp
    diagdata(1,:) = 5.0_dp ! offset +1, used columns 2..6
    diagdata(2,:) = 3.0_dp ! offset -2, used columns 1..2
    diagdata(3,:) = 1.0_dp ! offset 0, used columns 1..4
    d2 = dia_from_data(4,6,[1,-2,0],diagdata,info)
    call require(info == 0 .and. valid_dia(d2), 'sorted DIA constructor')
    call require(all(d2%offsets == [-2,0,1]), 'DIA sorted offsets')
    dense = dia_to_dense(d2)
    call require(dense(1,2) == 5.0_dp .and. dense(3,1) == 3.0_dp .and. &
        dense(2,2) == 1.0_dp .and. dense(1,1) == 1.0_dp, 'DIA offset alignment')
    d2 = dia_from_data(4,6,[0,0,2],diagdata,info)
    call require(info == sparse_invalid_input, 'duplicate DIA offsets rejected')
    b = csr_to_bsr(a,3,2,info)
    call require(info == sparse_invalid_input .and. .not. valid_bsr(b), 'BSR divisibility check')
    b = csr_to_bsr(a,0,2,info)
    call require(info == sparse_invalid_input, 'BSR block size check')

    empty = csr_from_dense(reshape([real(dp)::],[0,0]))
    d = csr_to_dia(empty,info)
    call require(info == 0 .and. valid_dia(d), 'empty DIA')
    r = dia_to_csr(d,info)
    call require(info == 0 .and. r%n_cols == 0 .and. r%n_rows == 0, 'empty DIA to CSR')
    b = csr_to_bsr(empty,2,3,info)
    call require(info == 0 .and. valid_bsr(b), 'empty BSR')
    r = bsr_to_csr(b,info)
    call require(info == 0 .and. r%n_cols == 0 .and. r%n_rows == 0, 'empty BSR to CSR')
    val = norm(empty,info=info)
    call require(info == 0 .and. val == 0.0_dp, 'empty Fro norm')
    val = norm(empty,1,info)
    call require(info == sparse_linalg_invalid_input, 'empty 1 norm rejected')

    ! Large sparse strip matrix: no dense storage in conversions or matvec.
    a = csr_from_coo(18000,18000,[1,18000],[1,18000],[2.0_dp,3.0_dp])
    d = csr_to_dia(a,info)
    call require(info == 0 .and. size(d%offsets) == 1, 'large narrow-band DIA')
    b = csr_to_bsr(a,2,2,info)
    call require(info == 0 .and. size(b%indices) == 2, 'large sparse BSR')
    val = norm(a,1,info)
    call require(info == 0 .and. val == 3.0_dp, 'large sparse norm')

    print *, 'sparse norms, DIA and BSR: PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_formats_norm
