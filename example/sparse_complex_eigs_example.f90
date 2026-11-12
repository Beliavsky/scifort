! SPDX-License-Identifier: MIT
program sparse_complex_eigs_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_coo
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex
    implicit none
    type(csr_complex_matrix) :: a
    complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: residual_norms(:)
    integer :: info

    a = csr_complex_from_coo(5,5,[1,2,3,4,5,1],[1,2,3,4,5,3], &
        [cmplx(1.0_dp,0.0_dp,dp), cmplx(2.0_dp,0.5_dp,dp), &
         cmplx(3.0_dp,-0.2_dp,dp), cmplx(4.0_dp,1.0_dp,dp), &
         cmplx(5.0_dp,-1.0_dp,dp), cmplx(0.3_dp,0.8_dp,dp)])
    call eigs_complex(a,2,eigenvalues,eigenvectors,info,which='LM', &
        ncv=5,tol=1.0e-10_dp,residual_norms=residual_norms)
    if (info /= 0) error stop 'complex eigs failed to converge'
    print '(a)', 'Selected eigenvalues:'
    print *, eigenvalues
    print '(a)', 'True residual norms:'
    print *, residual_norms
end program sparse_complex_eigs_example
