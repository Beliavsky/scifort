! SPDX-License-Identifier: MIT
! Compile with src/kinds.f90, sparse dependencies and sparse_linalg_complex_funm.f90.
program complex_funm_doc_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
    use scifort_sparse_linalg_complex_funm, only : funm_multiply_krylov_complex
    implicit none
    type(csr_complex_matrix) :: a
    complex(dp) :: dense(2,2), b(2), expected(2), square(2,2)
    complex(dp), allocatable :: y(:)
    integer :: info

    dense(1,:) = [cmplx(1.0_dp,0.3_dp,dp),cmplx(0.0_dp,1.0_dp,dp)]
    dense(2,:) = [cmplx(-0.5_dp,0.2_dp,dp),cmplx(0.1_dp,-0.2_dp,dp)]
    a = csr_complex_from_dense(dense)
    b = [cmplx(2.0_dp,1.0_dp,dp),cmplx(-1.0_dp,0.5_dp,dp)]
    call funm_multiply_krylov_complex(a,b,quadratic,y,info)
    if (info /= 0) error stop 'complex funm Krylov failed'
    square = matmul(dense,dense)
    expected = matmul(square,b)
    if (maxval(abs(y-expected)) > 1.0e-12_dp) error stop 'incorrect quadratic action'
    print '(a)', 'complex_funm_doc_example: PASS'
contains
    subroutine quadratic(h,value,status)
        complex(dp), intent(in) :: h(:,:) !! Small dense projected matrix.
        complex(dp), intent(out) :: value(:,:) !! h squared.
        integer, intent(out) :: status !! Success flag.
        value = matmul(h,h)
        status = 0
    end subroutine quadratic
end program complex_funm_doc_example
