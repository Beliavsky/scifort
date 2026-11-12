! SPDX-License-Identifier: MIT
program complex_eigs_shift_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_coo
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex
    implicit none
    type(csr_complex_matrix) :: a
    complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: residual_norms(:)
    integer :: info
    a = csr_complex_from_coo(6,6,[1,2,3,4,5,6,2], &
        [1,2,3,4,5,6,5], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.1_dp,dp), &
         cmplx(3.0_dp,-0.2_dp,dp),cmplx(4.0_dp,0.3_dp,dp), &
         cmplx(5.0_dp,-0.1_dp,dp),cmplx(6.0_dp,0.2_dp,dp), &
         cmplx(0.2_dp,0.1_dp,dp)])
    call eigs_complex(a,2,eigenvalues,eigenvectors,info, &
        sigma=cmplx(3.3_dp,-0.1_dp,dp),which='LM',ncv=6, &
        tol=1.0e-10_dp,residual_norms=residual_norms)
    if (info /= 0) error stop 'shifted eigs did not converge'
    print *, eigenvalues
    print *, residual_norms
end program complex_eigs_shift_example
