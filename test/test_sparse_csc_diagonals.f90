! SPDX-License-Identifier: MIT
program test_sparse_csc_diagonals
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, sparse_success, sparse_invalid_input, &
        csr_from_coo, csr_to_dense, csc_from_coo, csc_from_dense, csc_to_dense, &
        csc_to_csr, csr_to_csc, valid_csc, valid_csr, csc_matvec, csc_transpose, &
        csr_diagonal, csr_diags, csr_setdiag, csr_tril, csr_triu
    implicit none
    type(csr_matrix) :: a, b, d, old
    type(csc_matrix) :: c, t, invalid
    real(dp), allocatable :: x(:), y(:,:), v(:)
    integer :: status, i, n

    a = csr_from_coo(3,4,[3,1,2,1,3,1,2],[4,2,1,2,4,4,3], &
        [2.0_dp,3.0_dp,-1.0_dp,-2.0_dp,5.0_dp,0.0_dp,1.0_dp],status)
    call require(status == sparse_success,'base COO')
    c = csr_to_csc(a,status)
    call require(status == sparse_success .and. valid_csc(c),'CSR to CSC')
    call require(all(c%indptr == [1,2,3,4,6]),'CSC column pointers')
    call require(all(c%indices == [2,1,2,1,3]),'CSC sorted row indices')
    b = csc_to_csr(c,status)
    call require(status == sparse_success .and. valid_csr(b),'CSC to CSR')
    call require(all(b%indptr == a%indptr) .and. all(b%indices == a%indices), 'round trip structure')
    call require(all(b%data == a%data),'round trip entries and zeros')
    c = csc_from_coo(3,4,[3,1,2,1,3,1,2],[4,2,1,2,4,4,3], &
        [2.0_dp,3.0_dp,-1.0_dp,-2.0_dp,5.0_dp,0.0_dp,1.0_dp],status)
    call require(status == sparse_success .and. valid_csc(c),'CSC from unsorted COO')
    call require(all(csc_to_dense(c) == csr_to_dense(a)),'CSC dense values')
    c = csc_from_dense(csr_to_dense(a),status)
    call require(status == sparse_success .and. size(c%data) == 4,'CSC from dense omits zero')
    x = csc_matvec(c,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],status)
    call require(status == sparse_success .and. all(x == [2.0_dp,2.0_dp,28.0_dp]),'CSC matvec')
    t = csc_transpose(c,status)
    call require(status == sparse_success .and. valid_csc(t),'CSC transpose')
    call require(all(csc_to_dense(t) == transpose(csr_to_dense(a))),'CSC transpose values')
    t = csc_transpose(t,status)
    call require(all(csc_to_dense(t) == csr_to_dense(a)),'CSC double transpose')
    invalid = c
    invalid%indices(2) = 0
    call require(.not. valid_csc(invalid),'invalid CSC index')
    b = csc_to_csr(invalid,status)
    call require(status == sparse_invalid_input,'reject invalid CSC')
    x = csc_matvec(c,[1.0_dp],status)
    call require(status == sparse_invalid_input .and. size(x) == 0,'invalid CSC vector')

    x = csr_diagonal(a,status=status)
    call require(status == sparse_success .and. all(x == [0.0_dp,0.0_dp,0.0_dp]),'main diagonal')
    x = csr_diagonal(a,offset=1,status=status)
    call require(status == sparse_success .and. all(x == [1.0_dp,1.0_dp,7.0_dp]),'upper diagonal')
    x = csr_diagonal(a,offset=3,status=status)
    call require(status == sparse_success .and. size(x) == 1 .and. all(x == 0.0_dp), 'long offset')
    x = csr_diagonal(a,offset=-50,status=status)
    call require(status == sparse_success .and. size(x) == 0,'outside diagonal')

    ! Three packed diagonals: -1 has 3, 0 has 4, +2 has 3 positions.
    d = csr_diags(4,5,[-1,0,2],[1,4,8,11], &
        [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp, &
         8.0_dp,9.0_dp,0.0_dp],status)
    call require(status == sparse_success .and. valid_csr(d),'packed diagonals')
    call require(size(d%data) == 10,'packed diagonals preserve explicit zeros')
    x = csr_diagonal(d,offset=-1,status=status)
    call require(all(x == [1.0_dp,2.0_dp,3.0_dp]),'lower diagonal values')
    x = csr_diagonal(d,status=status)
    call require(all(x == [4.0_dp,5.0_dp,6.0_dp,7.0_dp]),'main diagonal values')
    x = csr_diagonal(d,offset=2,status=status)
    call require(all(x == [8.0_dp,9.0_dp,0.0_dp]),'upper diagonal values')
    b = csr_diags(4,5,[0,0],[1,5,9],[real(dp):: 1,2,3,4,5,6,7,8],status)
    call require(status == sparse_invalid_input,'reject duplicate offsets')
    b = csr_diags(4,5,[1],[1,2],[1.0_dp],status)
    call require(status == sparse_invalid_input,'reject inconsistent packed lengths')
    b = csr_diags(4,5,[0],[2,6],[real(dp):: 1,2,3,4,5],status)
    call require(status == sparse_invalid_input,'reject bad diagonal pointer start')

    old = d
    call csr_setdiag(d,[9.0_dp,0.0_dp,11.0_dp,12.0_dp],status=status)
    call require(status == sparse_success .and. valid_csr(d),'setdiag success')
    x = csr_diagonal(d)
    call require(all(x == [9.0_dp,0.0_dp,11.0_dp,12.0_dp]),'setdiag overwrites values')
    call require(any(d%data == 0.0_dp),'setdiag explicitly stores zero')
    call csr_setdiag(d,[2.0_dp],status=status)
    call require(status == sparse_invalid_input,'setdiag rejects incomplete diagonal')
    x = csr_diagonal(d)
    call require(all(x == [9.0_dp,0.0_dp,11.0_dp,12.0_dp]),'invalid setdiag preserves input')
    call csr_setdiag(d,[ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,0.0_dp,0.0_dp],status=status)
    call require(status == sparse_invalid_input,'setdiag rejects NaN')
    b = csr_triu(old,offset=2,status=status)
    call require(status == sparse_success .and. valid_csr(b),'upper triangle')
    call require(size(b%data) == 3 .and. all(csr_to_dense(b) == reshape(&
        [0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
         8.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,9.0_dp,0.0_dp,0.0_dp, &
         0.0_dp,0.0_dp,0.0_dp,0.0_dp],[4,5])), 'upper triangular values')
    b = csr_tril(old,offset=-1,status=status)
    call require(status == sparse_success .and. valid_csr(b),'lower triangle')
    x = csr_diagonal(b,-1)
    call require(all(x == [1.0_dp,2.0_dp,3.0_dp]) .and. size(b%data) == 3, 'lower triangular values')

    b = csr_from_coo(0,5,[integer::],[integer::],[real(dp)::],status)
    c = csr_to_csc(b,status)
    call require(status == sparse_success .and. valid_csc(c),'empty CSC rows')
    call require(all(shape(csc_to_dense(c)) == [0,5]),'empty CSC shape')
    t = csc_transpose(c,status)
    call require(status == sparse_success .and. valid_csc(t),'empty CSC transpose')
    call require(all(shape(csc_to_dense(t)) == [5,0]),'empty CSC transposed shape')
    b = csr_diags(0,5,[0],[1,1],[real(dp)::],status)
    call require(status == sparse_success .and. valid_csr(b),'zero-row diagonal construction')
    b = csr_from_coo(5,0,[integer::],[integer::],[real(dp)::],status)
    c = csr_to_csc(b,status)
    call require(status == sparse_success .and. valid_csc(c),'empty CSC columns')

    ! Sparse linear-size test for a long rectangular matrix.
    n = 25000
    b = csr_from_coo(n,n+5,[(i,i=1,n)],[(i+2,i=1,n)],[(2.0_dp,i=1,n)],status)
    c = csr_to_csc(b,status)
    call require(status == sparse_success .and. valid_csc(c),'large CSR-to-CSC')
    b = csc_to_csr(c,status)
    call require(status == sparse_success .and. size(b%data) == n,'large CSC-to-CSR')
    x = csr_diagonal(b,2,status)
    call require(size(x) == n .and. all(x == 2.0_dp),'large diagonal extraction')
    print '(a)', 'test_sparse_csc_diagonals: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! True if test passed.
        character(len=*), intent(in) :: message !! Diagnostic on failure.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_csc_diagonals
