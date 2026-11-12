! SPDX-License-Identifier: MIT
! Generalized Schur factors for a regular pencil with individually singular A,B.
program linalg_qz_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_qz, only : qz
    implicit none

    real(dp) :: a(3,3), b(3,3), error_a, error_b
    complex(dp) :: aa(3,3), bb(3,3), q(3,3), z(3,3), alpha(3), beta(3)
    integer :: status, j

    a = 0.0_dp
    b = 0.0_dp
    a(1,1) = 2.0_dp
    a(3,3) = 3.0_dp
    a(1,2) = 0.2_dp
    b(2,2) = 1.0_dp
    b(3,3) = 2.0_dp
    b(2,3) = -0.3_dp

    call qz(a,b,aa,bb,q,z,status,alpha,beta)
    if (status /= linalg_status_success) error stop 'generalized Schur failure'
    error_a = maxval(abs(matmul(q,matmul(aa,transpose(conjg(z))))-&
        cmplx(a,kind=dp)))/max(1.0_dp,maxval(abs(a)))
    error_b = maxval(abs(matmul(q,matmul(bb,transpose(conjg(z))))-&
        cmplx(b,kind=dp)))/max(1.0_dp,maxval(abs(b)))
    print '(a,es12.4)', 'A reconstruction: ', error_a
    print '(a,es12.4)', 'B reconstruction: ', error_b
    print '(a)', 'Generalized eigenvalues (alpha/beta):'
    do j=1,3
        if (abs(beta(j)) < 1.0e-12_dp) then
            print '(a,i0,a)', '  ',j, ': infinity'
        else
            print '(a,i0,a,2es16.7)', '  ',j, ': ', &
                real(alpha(j)/beta(j),dp), aimag(alpha(j)/beta(j))
        end if
    end do
    if (max(error_a,error_b) > 1.0e-10_dp) error stop 'QZ reconstruction'
end program linalg_qz_example
