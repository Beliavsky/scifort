! SPDX-License-Identifier: MIT
program test_sparse_dok_lil
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_to_dense, csc_to_dense, &
        valid_csr, valid_csc, sparse_success, sparse_invalid_input
    use scifort_sparse_coo, only : coo_matrix, coo_from_triplets
    use scifort_sparse_dok_lil
    implicit none
    type(dok_matrix) :: dok, dok2, invalid_dok
    type(lil_matrix) :: lil, lil2, invalid_lil
    type(coo_matrix) :: coo
    type(csr_matrix) :: csr
    type(csc_matrix) :: csc
    real(dp), allocatable :: dense(:,:), y(:), x(:)
    real(dp) :: v, old, expected(37,43), nan, inf
    integer :: ierr, t, i, j, n

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    inf = ieee_value(0.0_dp, ieee_positive_inf)
    invalid_dok = dok_empty(-1,4,ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_dok(invalid_dok), 'invalid DOK shape')
    invalid_lil = lil_empty(2,-2,ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_lil(invalid_lil), 'invalid LIL shape')
    call dok_set(invalid_dok,1,1,3.0_dp,ierr)
    call require(ierr == sparse_invalid_input, 'invalid DOK mutation')
    call lil_set(invalid_lil,1,1,3.0_dp,ierr)
    call require(ierr == sparse_invalid_input, 'invalid LIL mutation')

    dok = dok_empty(0,5,ierr)
    call require(ierr == sparse_success .and. valid_dok(dok), 'empty DOK 0x5')
    lil = lil_empty(0,5,ierr)
    call require(ierr == sparse_success .and. valid_lil(lil), 'empty LIL 0x5')
    dense = dok_to_dense(dok,ierr)
    call require(ierr == 0 .and. all(shape(dense) == [0,5]), 'empty DOK dense shape')
    dense = lil_to_dense(lil,ierr)
    call require(ierr == 0 .and. all(shape(dense) == [0,5]), 'empty LIL dense shape')
    csr = lil_to_csr(lil,ierr)
    call require(ierr == 0 .and. valid_csr(csr), 'empty LIL CSR')
    csr = dok_to_csr(dok,ierr)
    call require(ierr == 0 .and. valid_csr(csr), 'empty DOK CSR')
    dok = dok_empty(3,0,ierr)
    lil = lil_empty(3,0,ierr)
    call require(valid_dok(dok) .and. valid_lil(lil), 'empty column dimension')
    y = dok_matvec(dok,[real(dp)::],ierr)
    call require(ierr == 0 .and. all(y == 0.0_dp) .and. size(y) == 3, 'DOK empty matvec')
    y = lil_matvec(lil,[real(dp)::],ierr)
    call require(ierr == 0 .and. all(y == 0.0_dp) .and. size(y) == 3, 'LIL empty matvec')

    ! Five canonical coordinates; two are cancellation zeros, one is explicit zero.
    coo = coo_from_triplets(3,4,[2,1,1,1,3,3,2,2], &
        [4,2,2,3,1,1,4,3], &
        [2.0_dp,5.0_dp,-5.0_dp,0.0_dp,7.0_dp,-7.0_dp,1.0_dp,4.0_dp],ierr)
    call require(ierr == 0, 'test COO input')
    dok = dok_from_coo(coo,ierr)
    call require(ierr == 0 .and. valid_dok(dok), 'COO -> DOK')
    lil = lil_from_coo(coo,ierr)
    call require(ierr == 0 .and. valid_lil(lil), 'COO -> LIL')
    call require(dok_nnz(dok) == 5 .and. lil_nnz(lil) == 5, 'duplicate canonical zeros preserved')
    dense = dok_to_dense(dok,ierr)
    call require(ierr == 0 .and. all(shape(dense) == [3,4]), 'DOK dense shape')
    call require(all(dense == lil_to_dense(lil)), 'COO DOK/LIL dense comparison')
    call require(dense(2,3) == 4.0_dp .and. dense(2,4) == 3.0_dp, 'duplicate COO values')
    call require(count(dense /= 0.0_dp) == 2, 'COO cancellation values')
    call require(dok_get(dok,3,1,ierr) == 0.0_dp .and. ierr == 0, 'DOK get explicit zero')
    call require(lil_get(lil,3,1,ierr) == 0.0_dp .and. ierr == 0, 'LIL get explicit zero')
    csr = dok_to_csr(dok,ierr)
    call require(ierr == 0 .and. valid_csr(csr), 'DOK -> CSR')
    call require(size(csr%data) == 5, 'DOK -> CSR includes stored zeros')
    call require(all(csr_to_dense(csr) == dense), 'DOK -> CSR numeric content')
    csr = lil_to_csr(lil,ierr)
    call require(ierr == 0 .and. valid_csr(csr), 'LIL -> CSR')
    call require(size(csr%data) == 5, 'LIL -> CSR includes stored zeros')
    csc = lil_to_csc(lil,ierr)
    call require(ierr == 0 .and. valid_csc(csc), 'LIL -> CSC')
    call require(all(csc_to_dense(csc) == dense), 'LIL -> CSC dense content')
    csc = dok_to_csc(dok,ierr)
    call require(ierr == 0 .and. valid_csc(csc), 'DOK -> CSC')
    call require(all(csc_to_dense(csc) == dense), 'DOK -> CSC dense content')
    dok2 = dok_from_csc(csc,ierr)
    call require(ierr == 0 .and. dok_nnz(dok2) == 5, 'CSC -> DOK preserves zeros')
    lil2 = lil_from_csc(csc,ierr)
    call require(ierr == 0 .and. lil_nnz(lil2) == 5, 'CSC -> LIL preserves zeros')
    dok2 = dok_from_csr(csr,ierr)
    call require(ierr == 0 .and. dok_nnz(dok2) == 5, 'CSR -> DOK preserves zeros')
    lil2 = lil_from_csr(csr,ierr)
    call require(ierr == 0 .and. lil_nnz(lil2) == 5, 'CSR -> LIL preserves zeros')
    dok2 = dok_from_lil(lil,ierr)
    call require(ierr == 0 .and. dok_nnz(dok2) == 5, 'LIL -> DOK preserves zeros')
    lil2 = lil_from_dok(dok,ierr)
    call require(ierr == 0 .and. lil_nnz(lil2) == 5, 'DOK -> LIL preserves zeros')
    call require(all(dok_to_dense(dok2) == dense), 'LIL -> DOK values')
    call require(all(lil_to_dense(lil2) == dense), 'DOK -> LIL values')

    dok2 = dok_transpose(dok,ierr)
    call require(ierr == 0 .and. dok_nnz(dok2) == 5, 'DOK transpose zeros')
    call require(all(shape(dok_to_dense(dok2)) == [4,3]), 'DOK transpose shape')
    call require(all(dok_to_dense(dok2) == transpose(dense)), 'DOK transpose values')
    lil2 = lil_transpose(lil,ierr)
    call require(ierr == 0 .and. lil_nnz(lil2) == 5, 'LIL transpose zeros')
    call require(all(lil_to_dense(lil2) == transpose(dense)), 'LIL transpose values')

    call dok_set(dok,1,2,0.0_dp,ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 4, 'DOK zero assignment removes explicit zero')
    call lil_set(lil,1,2,0.0_dp,ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 4, 'LIL zero assignment removes explicit zero')
    call dok_eliminate_zeros(dok,ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 2, 'DOK zero elimination')
    call lil_eliminate_zeros(lil,ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 2, 'LIL zero elimination')
    call dok_add(dok,2,3,-4.0_dp,ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 1, 'DOK additive cancellation')
    call lil_add(lil,2,3,-4.0_dp,ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 1, 'LIL additive cancellation')
    call dok_delete(dok,2,4,ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 0, 'DOK delete last')
    call lil_delete(lil,2,4,ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 0, 'LIL delete last')
    call dok_delete(dok,2,4,ierr)
    call lil_delete(lil,2,4,ierr)
    call require(ierr == 0, 'delete absent entry is harmless')

    ! Deterministic mutation stress: force DOK rehash and deletion tombstones,
    ! and exercise LIL ordered insertions/removals at the start and middle of rows.
    dok = dok_empty(37,43,ierr)
    lil = lil_empty(37,43,ierr)
    expected = 0.0_dp
    do t = 1, 2600
        i = modulo(17*t+13,37)+1
        j = modulo(37*t+11,43)+1
        old = expected(i,j)
        select case (mod(t,7))
        case (0,1,2)
            v = real(modulo(13*t,17)-8,dp)
            expected(i,j) = v
            call dok_set(dok,i,j,v,ierr)
            call require(ierr == 0, 'DOK randomized set')
            call lil_set(lil,i,j,v,ierr)
            call require(ierr == 0, 'LIL randomized set')
        case (3,4)
            v = real(modulo(7*t,5)-2,dp)
            expected(i,j) = old+v
            call dok_add(dok,i,j,v,ierr)
            call require(ierr == 0, 'DOK randomized add')
            call lil_add(lil,i,j,v,ierr)
            call require(ierr == 0, 'LIL randomized add')
        case default
            expected(i,j) = 0.0_dp
            call dok_delete(dok,i,j,ierr)
            call require(ierr == 0, 'DOK randomized delete')
            call lil_delete(lil,i,j,ierr)
            call require(ierr == 0, 'LIL randomized delete')
        end select
        call require(dok_get(dok,i,j,ierr) == expected(i,j) .and. ierr == 0, &
            'DOK randomized get')
        call require(lil_get(lil,i,j,ierr) == expected(i,j) .and. ierr == 0, &
            'LIL randomized get')
        if (mod(t,130) /= 0) cycle
        call require(valid_dok(dok) .and. valid_lil(lil), 'randomized structural invariants')
        n = count(expected /= 0.0_dp)
        call require(dok_nnz(dok) == n .and. lil_nnz(lil) == n, 'randomized stored count')
        call require(all(dok_to_dense(dok) == expected), 'randomized DOK dense equality')
        call require(all(lil_to_dense(lil) == expected), 'randomized LIL dense equality')
        csr = dok_to_csr(dok,ierr)
        call require(ierr == 0 .and. valid_csr(csr), 'randomized DOK CSR')
        call require(all(csr_to_dense(csr) == expected), 'randomized DOK CSR values')
        csr = lil_to_csr(lil,ierr)
        call require(ierr == 0 .and. valid_csr(csr), 'randomized LIL CSR')
        call require(all(csr_to_dense(csr) == expected), 'randomized LIL CSR values')
        x = [ (real(modulo(j*7,13)-5,dp), j=1,43) ]
        y = dok_matvec(dok,x,ierr)
        call require(ierr == 0 .and. all(abs(y-matmul(expected,x)) < 1.0e-10_dp), &
            'randomized DOK matvec')
        y = lil_matvec(lil,x,ierr)
        call require(ierr == 0 .and. all(abs(y-matmul(expected,x)) < 1.0e-10_dp), &
            'randomized LIL matvec')
    end do
    call dok_clear(dok,ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 0 .and. valid_dok(dok), 'DOK clear')
    call lil_clear(lil,ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 0 .and. valid_lil(lil), 'LIL clear')

    ! Sparse, large-index shapes: no dense matrix materialized.
    dok = dok_empty(120000,130000,ierr)
    lil = lil_empty(120000,130000,ierr)
    call require(ierr == 0, 'large sparse shapes')
    call dok_set(dok,120000,130000,2.0_dp,ierr)
    call lil_set(lil,120000,130000,2.0_dp,ierr)
    call require(ierr == 0, 'large shape insertion')
    call require(dok_get(dok,120000,130000,ierr) == 2.0_dp, 'large DOK get')
    call require(lil_get(lil,120000,130000,ierr) == 2.0_dp, 'large LIL get')
    csr = lil_to_csr(lil,ierr)
    call require(ierr == 0 .and. valid_csr(csr) .and. size(csr%data) == 1, &
        'large LIL to CSR')
    csr = dok_to_csr(dok,ierr)
    call require(ierr == 0 .and. valid_csr(csr) .and. size(csr%data) == 1, &
        'large DOK to CSR')

    ! Invalid coordinates, NaN/Inf and overflow must not modify the input.
    call dok_set(dok,0,2,3.0_dp,ierr)
    call require(ierr == sparse_invalid_input .and. dok_nnz(dok) == 1, 'bad DOK row')
    call lil_set(lil,1,0,3.0_dp,ierr)
    call require(ierr == sparse_invalid_input .and. lil_nnz(lil) == 1, 'bad LIL column')
    call dok_set(dok,1,1,nan,ierr)
    call require(ierr == sparse_invalid_input .and. dok_nnz(dok) == 1, 'DOK NaN rejected')
    call lil_set(lil,1,1,inf,ierr)
    call require(ierr == sparse_invalid_input .and. lil_nnz(lil) == 1, 'LIL Inf rejected')
    v = dok_get(dok,120001,1,ierr)
    call require(ierr == sparse_invalid_input .and. v == 0.0_dp, 'DOK invalid get')
    v = lil_get(lil,1,130001,ierr)
    call require(ierr == sparse_invalid_input .and. v == 0.0_dp, 'LIL invalid get')
    call dok_set(dok,1,1,huge(1.0_dp),ierr)
    call lil_set(lil,1,1,huge(1.0_dp),ierr)
    call require(ierr == 0, 'large finite DOK/LIL values accepted')
    call dok_add(dok,1,1,huge(1.0_dp),ierr)
    call require(ierr == sparse_invalid_input, 'DOK sum overflow rejected')
    call require(dok_get(dok,1,1) == huge(1.0_dp), 'DOK overflow input unchanged')
    call lil_add(lil,1,1,huge(1.0_dp),ierr)
    call require(ierr == sparse_invalid_input, 'LIL sum overflow rejected')
    call require(lil_get(lil,1,1) == huge(1.0_dp), 'LIL overflow input unchanged')
    y = dok_matvec(dok,[real(dp)::],ierr)
    call require(ierr == sparse_invalid_input .and. size(y) == 0, 'DOK matvec bad length')
    y = lil_matvec(lil,[real(dp)::],ierr)
    call require(ierr == sparse_invalid_input .and. size(y) == 0, 'LIL matvec bad length')
    dok2 = dok_from_dense(reshape([huge(1.0_dp)],[1,1]),ierr)
    y = dok_matvec(dok2,[2.0_dp],ierr)
    call require(ierr == sparse_invalid_input .and. size(y) == 0, 'DOK matvec overflow')
    lil2 = lil_from_dense(reshape([huge(1.0_dp)],[1,1]),ierr)
    y = lil_matvec(lil2,[2.0_dp],ierr)
    call require(ierr == sparse_invalid_input .and. size(y) == 0, 'LIL matvec overflow')

    ! Assignment copies are independent; constructors skip dense zeros.
    dok = dok_from_dense(reshape([0.0_dp,2.0_dp,-3.0_dp,0.0_dp],[2,2]),ierr)
    call require(ierr == 0 .and. dok_nnz(dok) == 2, 'DOK from dense')
    lil = lil_from_dense(reshape([0.0_dp,2.0_dp,-3.0_dp,0.0_dp],[2,2]),ierr)
    call require(ierr == 0 .and. lil_nnz(lil) == 2, 'LIL from dense')
    dok2 = dok
    lil2 = lil
    call dok_delete(dok2,2,1,ierr)
    call lil_delete(lil2,2,1,ierr)
    call require(dok_nnz(dok2) == 1 .and. dok_nnz(dok) == 2, 'DOK copy independent')
    call require(lil_nnz(lil2) == 1 .and. lil_nnz(lil) == 2, 'LIL copy independent')
    invalid_dok = dok_from_dense(reshape([nan],[1,1]),ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_dok(invalid_dok), &
        'DOK invalid dense NaN')
    invalid_lil = lil_from_dense(reshape([inf],[1,1]),ierr)
    call require(ierr == sparse_invalid_input .and. .not. valid_lil(invalid_lil), &
        'LIL invalid dense Inf')

    print *, 'sparse DOK/LIL: PASS'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok !! Condition required for the test to pass.
        character(len=*), intent(in) :: message !! Failure context.
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_dok_lil
