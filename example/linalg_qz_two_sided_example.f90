! SPDX-License-Identifier: MIT
! Generalized complex Schur factors using inverse-free QZ iterations.
program linalg_qz_two_sided_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_qz, only : qz
    implicit none
    real(dp) :: a(3,3), b(3,3), err_a, err_b
    complex(dp) :: aa(3,3), bb(3,3), q(3,3), z(3,3), alpha(3), beta(3)
    integer :: status, i

    a = 0.0_dp
    b = 0.0_dp
    do i=1,3
        b(i,i) = 1.0_dp + 0.2_dp*real(i,dp)
    end do
    ! Two nonreal generalized eigenvalues plus a real third eigenvalue.
    a(1,2) = -2.0_dp
    a(2,1) = 2.0_dp
    a(3,3) = 3.0_dp
    a(1,3) = 0.1_dp

    call qz(a,b,aa,bb,q,z,status,alpha,beta,method='iterative')
    if (status /= linalg_status_success) error stop 'iterative generalized QZ'
    err_a = maxval(abs(matmul(q,matmul(aa,transpose(conjg(z))))-cmplx(a,kind=dp)))
    err_b = maxval(abs(matmul(q,matmul(bb,transpose(conjg(z))))-cmplx(b,kind=dp)))
    print '(a,es13.5)', 'A reconstruction error: ', err_a
    print '(a,es13.5)', 'B reconstruction error: ', err_b
    do i=1,3
        if (abs(beta(i)) > 0.0_dp) then
            print '(a,i0,a,2es16.7)', 'lambda(',i,') = ', &
                real(alpha(i)/beta(i),dp),aimag(alpha(i)/beta(i))
        else
            print '(a,i0,a)', 'lambda(',i,') = infinity'
        end if
    end do
    if (max(err_a,err_b) > 1.0e-10_dp) error stop 'QZ reconstruction'
end program linalg_qz_two_sided_example
