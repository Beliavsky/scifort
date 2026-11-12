! SPDX-License-Identifier: MIT
program linalg_pinvh_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_pinvh, only : pinvh
    implicit none

    complex(dp) :: a(3,3), ap(3,3)
    real(dp) :: real_matrix(3,3), real_inverse(3,3), residual
    integer :: status, rank

    a = cmplx(0.0_dp,0.0_dp,dp)
    a(1,1) = cmplx(-4.0_dp,0.0_dp,dp)
    a(2,2) = cmplx(2.0_dp,0.0_dp,dp)
    a(2,1) = cmplx(1.0_dp,-0.5_dp,dp)
    a(1,2) = conjg(a(2,1))
    ! The third row and column vanish, so A is singular of rank two.
    call pinvh(a,ap,status,rank)
    if (status /= linalg_status_success) error stop 'complex pinvh failure'
    residual = maxval(abs(matmul(a,matmul(ap,a))-a))
    print '(a,i0)', 'Hermitian pseudoinverse rank: ',rank
    print '(a,es12.4)', 'A*A_pinv*A reconstruction error: ',residual

    real_matrix = 0.0_dp
    real_matrix(1,1) = -4.0_dp
    real_matrix(2,2) = 2.0_dp
    real_matrix(3,3) = 1.0e-4_dp
    call pinvh(real_matrix,real_inverse,status,rank,rtol=1.0e-3_dp)
    if (status /= linalg_status_success) error stop 'real pinvh failure'
    print '(a,i0)', 'Truncated real pseudoinverse rank: ',rank
    print '(a,es12.4)', 'Reciprocal of retained negative eigenvalue: ',real_inverse(1,1)
end program linalg_pinvh_example
