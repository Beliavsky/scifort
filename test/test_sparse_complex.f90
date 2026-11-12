! SPDX-License-Identifier: MIT
program test_sparse_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success, sparse_invalid_input, csr_from_dense
    use scifort_sparse_coo_nd_complex, only : coo_nd_complex_array, valid_coo_nd_complex
    use scifort_sparse_complex
    use scifort_sparse_linalg_complex
    implicit none
    type(csr_complex_matrix) :: a, b, c
    type(csc_complex_matrix) :: ac, cc
    type(coo_nd_complex_array) :: coo
    complex(dp) :: dense(3,3), rhs(3), truth(3), z, bad(1,1), square(3,3)
    complex(dp), allocatable :: result(:), array(:,:), column(:,:), y(:)
    integer :: info, iter
    real(dp) :: residual

    a = csr_complex_from_coo(3,3,[2,1,1,2,1,3],[2,1,1,3,3,1], &
        [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp),cmplx(1.0_dp,-3.0_dp,dp), &
         cmplx(0.0_dp,-1.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.5_dp,1.0_dp,dp)],info)
    call require(info == sparse_success .and. valid_csr_complex(a),'CSR construction')
    call require(size(a%data) == 5,'duplicate merge and explicit zeros')
    call require(a%indptr(1) == 1 .and. a%indptr(4) == 6,'one-based pointers')
    dense = cmplx(0.0_dp,0.0_dp,dp)
    dense(1,1) = cmplx(3.0_dp,-2.0_dp,dp)
    dense(1,3) = cmplx(0.0_dp,0.0_dp,dp)
    dense(2,3) = cmplx(0.0_dp,-1.0_dp,dp)
    dense(3,1) = cmplx(0.5_dp,1.0_dp,dp)
    array = csr_complex_to_dense(a,info)
    call require(info == sparse_success .and. maxval(abs(array-dense)) < 1.0e-13_dp,'CSR to dense')
    ac = csr_complex_to_csc(a,info)
    call require(info == sparse_success .and. valid_csc_complex(ac),'CSR to CSC')
    call require(size(ac%data) == 5,'CSC retains stored zeros')
    array = csc_complex_to_dense(ac,info)
    call require(info == sparse_success .and. maxval(abs(array-dense)) < 1.0e-13_dp,'CSC to dense')
    b = csc_complex_to_csr(ac,info)
    call require(info == sparse_success .and. all(b%data == a%data),'CSC to CSR')
    call require(all(b%indices == a%indices) .and. all(b%indptr == a%indptr), 'CSR roundtrip structure')
    coo = csr_complex_to_coo_nd(a,info)
    call require(info == sparse_success .and. valid_coo_nd_complex(coo),'CSR to rank-2 COO')
    b = csr_complex_from_coo_nd(coo,info)
    call require(info == sparse_success .and. all(b%data == a%data), 'COO conversion')
    b = csr_complex_transpose(a,info)
    array = csr_complex_to_dense(b)
    call require(info == sparse_success .and. maxval(abs(array-transpose(dense))) < 1.0e-13_dp, 'transpose')
    b = csr_complex_hermitian(a,info)
    array = csr_complex_to_dense(b)
    call require(info == sparse_success .and. &
        maxval(abs(array-transpose(conjg(dense)))) < 1.0e-13_dp, 'Hermitian transpose')
    cc = csc_complex_hermitian(ac,info)
    array = csc_complex_to_dense(cc)
    call require(info == sparse_success .and. &
        maxval(abs(array-transpose(conjg(dense)))) < 1.0e-13_dp, 'CSC Hermitian transpose')
    cc = csc_complex_transpose(ac,info)
    array = csc_complex_to_dense(cc)
    call require(info == sparse_success .and. maxval(abs(array-transpose(dense))) < 1.0e-13_dp,'CSC transpose')
    b = csr_complex_conjugate(a,info)
    array = csr_complex_to_dense(b)
    call require(info == sparse_success .and. maxval(abs(array-conjg(dense))) < 1.0e-13_dp,'conjugation')
    result = csr_complex_matvec(a,[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(1.0_dp,-2.0_dp,dp),cmplx(3.0_dp,1.0_dp,dp)],info)
    call require(info == sparse_success .and. size(result) == 3,'CSR matvec status')
    call require(maxval(abs(result-matmul(dense,[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(1.0_dp,-2.0_dp,dp),cmplx(3.0_dp,1.0_dp,dp)]))) < 1.0e-13_dp,'CSR matvec')
    y = csc_complex_matvec(ac,[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(1.0_dp,-2.0_dp,dp),cmplx(3.0_dp,1.0_dp,dp)],info)
    call require(info == sparse_success .and. maxval(abs(y-result)) < 1.0e-13_dp,'CSC matvec')
    column = csr_complex_matmat(a,reshape([cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.0_dp,1.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp), &
        cmplx(2.0_dp,0.0_dp,dp),cmplx(-1.0_dp,1.0_dp,dp),cmplx(1.0_dp,2.0_dp,dp)],[3,2]),info)
    call require(info == sparse_success .and. size(column,2) == 2,'CSR matrix-matrix product')
    result = csr_complex_matvec(a,[cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.0_dp,1.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp)])
    call require(maxval(abs(column(:,1)-result)) < 1.0e-13_dp,'CSR matmat first column')
    array = csr_complex_to_dense(a,info,max_elements=3_int64)
    call require(info == sparse_invalid_input .and. size(array) == 0,'dense allocation cap')
    a%indices(1) = 0
    call require(.not. valid_csr_complex(a),'bad indices rejected')
    a%indices(1) = 1
    bad(1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    b = csr_complex_from_dense(bad,info)
    call require(info == sparse_invalid_input .and. .not. valid_csr_complex(b),'reject NaN input')

    ! Three-by-three Hermitian strictly diagonally dominant positive definite.
    square = cmplx(0.0_dp,0.0_dp,dp)
    square(1,1) = cmplx(5.0_dp,0.0_dp,dp)
    square(2,2) = cmplx(6.0_dp,0.0_dp,dp)
    square(3,3) = cmplx(7.0_dp,0.0_dp,dp)
    square(1,2) = cmplx(1.0_dp,2.0_dp,dp)
    square(2,1) = conjg(square(1,2))
    square(2,3) = cmplx(0.5_dp,-1.0_dp,dp)
    square(3,2) = conjg(square(2,3))
    a = csr_complex_from_dense(square,info)
    ac = csr_complex_to_csc(a)
    truth = [cmplx(1.0_dp,0.5_dp,dp),cmplx(-2.0_dp,3.0_dp,dp),cmplx(0.7_dp,-1.0_dp,dp)]
    rhs = matmul(square,truth)
    call cg_complex(a,rhs,y,info,rtol=1.0e-12_dp,jacobi=.true.,iterations=iter,residual_norm=residual)
    call require(info == 0 .and. iter > 0 .and. residual < 1.0e-10_dp,'CG convergence')
    call require(maxval(abs(y-truth)) < 1.0e-10_dp, 'CG accuracy')
    call cg_complex(ac,rhs,y,info,rtol=1.0e-12_dp,jacobi=.false.)
    call require(info == 0 .and. maxval(abs(y-truth)) < 1.0e-10_dp, 'CG CSC wrapper')
    call cg_complex(a,rhs,y,info,rtol=1.0e-12_dp,x0=truth,iterations=iter)
    call require(info == 0 .and. iter == 0,'CG exact x0')

    ! Nonsymmetric, genuinely complex and diagonally dominant.
    square(1,1) = cmplx(4.0_dp,1.0_dp,dp)
    square(2,2) = cmplx(5.0_dp,-0.5_dp,dp)
    square(3,3) = cmplx(6.0_dp,1.0_dp,dp)
    square(1,3) = cmplx(-0.3_dp,0.5_dp,dp)
    square(3,1) = cmplx(0.9_dp,0.3_dp,dp)
    a = csr_complex_from_dense(square)
    ac = csr_complex_to_csc(a)
    rhs = matmul(square,truth)
    call bicgstab_complex(a,rhs,y,info,rtol=1.0e-12_dp,jacobi=.true., &
        iterations=iter,residual_norm=residual)
    call require(info == 0 .and. iter > 0 .and. residual < 1.0e-10_dp,'BiCGSTAB convergence')
    call require(maxval(abs(y-truth)) < 1.0e-10_dp,'BiCGSTAB accuracy')
    call bicgstab_complex(ac,rhs,y,info,rtol=1.0e-12_dp,jacobi=.false.)
    call require(info == 0 .and. maxval(abs(y-truth)) < 1.0e-10_dp,'BiCGSTAB CSC')
    call gmres_complex(a,rhs,y,info,rtol=1.0e-12_dp,jacobi=.true., &
        restart=2,maxiter=30,iterations=iter,residual_norm=residual)
    call require(info == 0 .and. iter > 0 .and. residual < 1.0e-10_dp,'restarted GMRES convergence')
    call require(maxval(abs(y-truth)) < 1.0e-10_dp,'GMRES accuracy')
    call gmres_complex(ac,rhs,y,info,rtol=1.0e-12_dp,restart=3,iterations=iter)
    call require(info == 0 .and. maxval(abs(y-truth)) < 1.0e-10_dp,'GMRES CSC')
    call gmres_complex(a,rhs,y,info,rtol=1.0e-12_dp,x0=truth,iterations=iter)
    call require(info == 0 .and. iter == 0,'GMRES exact initial guess')
    call gmres_complex(a,rhs,y,info,maxiter=0)
    call require(info > 0,'GMRES iteration limit')
    call gmres_complex(a,rhs,y,info,restart=0)
    call require(info == sparse_complex_linalg_invalid_input,'GMRES invalid restart')
    call bicgstab_complex(a,rhs,y,info,rtol=-1.0_dp)
    call require(info == sparse_complex_linalg_invalid_input,'invalid tolerance')
    call cg_complex(a,rhs,y,info,x0=[cmplx(1.0_dp,0.0_dp,dp)])
    call require(info == sparse_complex_linalg_invalid_input,'invalid initial guess')
    square(2,2) = cmplx(0.0_dp,0.0_dp,dp)
    a = csr_complex_from_dense(square)
    call gmres_complex(a,rhs,y,info,jacobi=.true.)
    call require(info == sparse_complex_linalg_bad_diagonal,'zero Jacobi diagonal')

    b = csr_complex_from_dense(reshape([cmplx(0.0_dp,0.0_dp,dp)],[1,1]),info)
    call require(info == sparse_success .and. valid_csr_complex(b),'zero matrix')
    b = csr_complex_from_coo(0,5,[integer::],[integer::],[complex(dp)::],info)
    call require(info == sparse_success .and. valid_csr_complex(b),'empty rectangle')
    print *, 'test_sparse_complex: PASS'

contains

    subroutine require(ok, message)
        logical, intent(in) :: ok !! True if assertion holds.
        character(*), intent(in) :: message !! Diagnostic on failure.
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require

end program test_sparse_complex
