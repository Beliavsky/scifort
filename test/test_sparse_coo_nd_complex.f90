! SPDX-License-Identifier: MIT
program test_sparse_coo_nd_complex
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, valid_coo_nd
    use scifort_sparse_coo_nd_complex
    use scifort_sparse_coo_nd_complex_ops
    implicit none
    type(coo_nd_complex_array) :: a, b, c, d, v, saved
    type(coo_nd_array) :: r
    complex(dp), allocatable :: flat(:)
    complex(dp) :: input(2,3,2), output(2,3,2), total, nanval
    integer :: info

    a = coo_nd_complex_from_triplets([2,3,2], reshape([ &
        1,1,1, 1,1,1, 2,3,2, 1,2,1 ],[3,4]), &
        [cmplx(1.0_dp,2.0_dp,dp),cmplx(-1.0_dp,3.0_dp,dp), &
         cmplx(4.0_dp,-2.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp)],info)
    call require(info == sparse_success .and. valid_coo_nd_complex(a),'construct with duplicates')
    call require(size(a%data) == 4,'duplicate/zero storage count')
    call require(close(coo_nd_complex_get(a,[1,1,1],info),cmplx(0.0_dp,5.0_dp,dp)), &
        'duplicate gather')
    call require(info == sparse_success,'duplicate get status')
    saved = a
    call coo_nd_complex_sum_duplicates(a,info)
    call require(info == sparse_success .and. coo_nd_complex_is_canonical(a),'canonicalize complex COO')
    call require(size(a%data) == 3,'coalesce duplicate values')
    call require(coo_nd_complex_get(a,[1,2,1]) == cmplx(0.0_dp,0.0_dp,dp), &
        'preserve stored zeros during coalesce')
    call coo_nd_complex_eliminate_zeros(a,info)
    call require(info == sparse_success .and. size(a%data) == 2, 'remove explicit zeros')
    call require(close(coo_nd_complex_get(a,[2,3,2]),cmplx(4.0_dp,-2.0_dp,dp)), &
        'second nonzero')
    b = coo_nd_complex_conjugate(saved,info)
    call require(info == sparse_success .and. size(b%data) == 4,'conjugate preserves entries')
    call require(close(coo_nd_complex_get(b,[1,1,1]),cmplx(0.0_dp,-5.0_dp,dp)), &
        'conjugation of duplicate sum')
    b = coo_nd_complex_transpose(saved,[3,1,2],info)
    call require(info == sparse_success .and. all(b%shape == [2,2,3]),'axis permutation')
    call require(close(coo_nd_complex_get(b,[2,2,3]),cmplx(4.0_dp,-2.0_dp,dp)), &
        'axis permutation values')
    b = coo_nd_complex_reshape(saved,[3,2,2],order='C',status=info)
    call require(info == sparse_success .and. all(b%shape == [3,2,2]),'C-order reshape')
    b = coo_nd_complex_reshape(saved,[3,2,2],order='F',status=info)
    call require(info == sparse_success .and. all(b%shape == [3,2,2]),'F-order reshape')
    b = coo_nd_complex_reshape(saved,[3,2,3],status=info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd_complex(b), &
        'bad reshape rejected')

    input = cmplx(0.0_dp,0.0_dp,dp)
    input(1,1,1) = cmplx(0.0_dp,5.0_dp,dp)
    input(2,3,2) = cmplx(4.0_dp,-2.0_dp,dp)
    b = coo_nd_complex_from_dense(input,info)
    call require(info == sparse_success .and. size(b%data) == 2,'assumed-rank dense constructor')
    call coo_nd_complex_to_dense(b,output,info)
    call require(info == sparse_success .and. all(output == input),'assumed-rank dense roundtrip')
    flat = coo_nd_complex_to_dense_flat(saved,info,'C')
    call require(info == sparse_success .and. size(flat) == 12,'flat dense conversion')
    b = coo_nd_complex_from_dense_flat([2,3,2],flat,info,'C')
    call require(info == sparse_success,'from C-ordered dense flat')
    call require(close(coo_nd_complex_get(b,[1,1,1]),input(1,1,1)),'C-order flat values')
    b = coo_nd_complex_from_dense_flat([2,3,2],flat,info,'F')
    call require(info == sparse_success,'from F-ordered flat')
    flat = coo_nd_complex_to_dense_flat(saved,info,max_elements=3_int64)
    call require(info == sparse_invalid_input,'dense allocation cap')

    c = coo_nd_complex_from_triplets([2,2],reshape([1,2, 2,1],[2,2]), &
        [cmplx(2.0_dp,1.0_dp,dp),cmplx(-1.0_dp,3.0_dp,dp)],info)
    b = coo_nd_complex_hermitian(c,info)
    call require(info == sparse_success .and. all(b%shape == [2,2]), 'Hermitian shape')
    call require(close(coo_nd_complex_get(b,[2,1]),cmplx(2.0_dp,-1.0_dp,dp)), &
        'Hermitian transpose conjugates')
    b = coo_nd_complex_hermitian(saved,info)
    call require(info == sparse_invalid_input,'Hermitian restricted to rank two')

    r = coo_nd_from_triplets([2,2],reshape([1,2, 1,2],[2,2]),[1.0_dp,-2.0_dp],info)
    b = coo_nd_complex_from_real(r,info)
    call require(info == sparse_success .and. valid_coo_nd_complex(b),'real COO interoperability')
    r = coo_nd_complex_to_real(b,info)
    call require(info == sparse_success .and. valid_coo_nd(r),'real conversion')
    r = coo_nd_complex_to_real(c,info)
    call require(info == sparse_invalid_input .and. .not. valid_coo_nd(r), &
        'nonzero imaginary components rejected during real conversion')

    total = coo_nd_complex_sum_all(saved,info)
    call require(info == sparse_success .and. close(total,cmplx(4.0_dp,3.0_dp,dp)), &
        'sum all entries')
    b = coo_nd_complex_sum_axis(saved,1,status=info)
    call require(info == sparse_success .and. all(b%shape == [3,2]),'axis reduction')
    call require(close(coo_nd_complex_get(b,[1,1]),cmplx(0.0_dp,5.0_dp,dp)), &
        'axis reduction complex aggregation')
    b = coo_nd_complex_sum_axes(saved,[1,3],keepdims=.true.,status=info)
    call require(info == sparse_success .and. all(b%shape == [1,3,1]),'keepdims reduction')

    v = coo_nd_complex_from_triplets([3],reshape([2,3],[1,2]), &
        [cmplx(2.0_dp,1.0_dp,dp),cmplx(-1.0_dp,2.0_dp,dp)],info)
    b = coo_nd_complex_broadcast_to(v,[2,3],info)
    call require(info == sparse_success .and. size(b%data) == 4,'broadcast complex entries')
    b = coo_nd_complex_add(v,v,info)
    call require(info == sparse_success .and. &
        close(coo_nd_complex_get(b,[2]),cmplx(4.0_dp,2.0_dp,dp)),'complex addition')
    b = coo_nd_complex_multiply(v,v,info)
    call require(info == sparse_success .and. &
        close(coo_nd_complex_get(b,[2]),cmplx(3.0_dp,4.0_dp,dp)),'complex square')
    b = coo_nd_complex_multiply(v,v,info,max_nnz=1_int64)
    call require(info == sparse_invalid_input,'multiplication entry limit')
    b = coo_nd_complex_tensordot(v,v,[1],[1],info)
    call require(info == sparse_success .and. all(b%shape == [1]),'full complex contraction')
    call require(close(coo_nd_complex_get(b,[1]),cmplx(0.0_dp,0.0_dp,dp)), &
        'tensordot is not conjugated')
    total = coo_nd_complex_vdot(v,v,info)
    call require(info == sparse_success .and. close(total,cmplx(10.0_dp,0.0_dp,dp)), &
        'vdot is conjugated')
    b = coo_nd_complex_tensordot(v,v,[integer::],[integer::],info)
    call require(info == sparse_success .and. all(b%shape == [3,3]) .and. size(b%data) == 4, &
        'outer complex product')
    b = coo_nd_complex_tensordot(v,v,[integer::],[integer::],info,max_nnz=3_int64)
    call require(info == sparse_invalid_input,'contraction cap')
    total = coo_nd_complex_vdot(v,c,info)
    call require(info == sparse_invalid_input,'vdot shape mismatch')

    c = coo_nd_complex_from_triplets([2000000000,2000000000], &
        reshape([2000000000,1],[2,1]),[cmplx(2.0_dp,1.0_dp,dp)],info)
    d = coo_nd_complex_from_triplets([2000000000,2000000000], &
        reshape([1,2000000000],[2,1]),[cmplx(1.0_dp,2.0_dp,dp)],info)
    b = coo_nd_complex_tensordot(c,d,[2],[1],info)
    call require(info == sparse_success .and. all(b%shape == [2000000000,2000000000]), &
        'large sparse complex contraction shape')
    call require(size(b%data) == 1 .and. close(b%data(1),cmplx(0.0_dp,5.0_dp,dp)), &
        'large sparse complex contraction value')

    nanval = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    b = coo_nd_complex_from_triplets([1],reshape([1],[1,1]),[nanval],info)
    call require(info == sparse_invalid_input,'nonfinite real component rejected')
    nanval = cmplx(0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),dp)
    b = coo_nd_complex_from_triplets([1],reshape([1],[1,1]),[nanval],info)
    call require(info == sparse_invalid_input,'nonfinite imaginary component rejected')
    c = coo_nd_complex_from_triplets([1],reshape([1,1],[1,2]), &
        [cmplx(huge(1.0_dp),0.0_dp,dp),cmplx(huge(1.0_dp),0.0_dp,dp)],info)
    saved = c
    call coo_nd_complex_sum_duplicates(c,info)
    call require(info == sparse_invalid_input .and. all(c%data == saved%data), &
        'overflow leaves original COO unchanged')
    total = coo_nd_complex_sum_all(c,info)
    call require(info == sparse_invalid_input,'complex sum overflow detected')

    print '(a)', 'PASS test_sparse_coo_nd_complex'
contains
    pure logical function close(x,y)
        complex(dp), intent(in) :: x !! Compared value.
        complex(dp), intent(in) :: y !! Reference value.
        close = abs(x-y) < 1.0e-12_dp
    end function close

    subroutine require(ok, description)
        logical, intent(in) :: ok !! Checked condition.
        character(len=*), intent(in) :: description !! Descriptive failure context.
        if (.not. ok) then
            print '(a)', 'FAIL test_sparse_coo_nd_complex: '//description
            error stop 1
        end if
    end subroutine require
end program test_sparse_coo_nd_complex
