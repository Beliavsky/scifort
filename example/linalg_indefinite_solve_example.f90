! SPDX-License-Identifier: MIT
program linalg_indefinite_solve_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_ldl, only : ldl
    use scifort_linalg_indefinite_solve, only : ldl_solve, solve_indefinite
    implicit none

    real(dp) :: a(3,3), lu(3,3), d(3,3), b(3,2), x(3,2), truth(3,2)
    complex(dp) :: ah(2,2), bh(2), xh(2), expected(2)
    integer :: perm(3), status

    ! An indefinite matrix with a 2x2 pivot and multiple right-hand sides.
    a(1,:) = [0.0_dp,2.0_dp,0.0_dp]
    a(2,:) = [2.0_dp,0.0_dp,1.0_dp]
    a(3,:) = [0.0_dp,1.0_dp,3.0_dp]
    truth(:,1) = [1.0_dp,2.0_dp,3.0_dp]
    truth(:,2) = [-1.0_dp,3.0_dp,2.0_dp]
    b = matmul(a,truth)

    call ldl(a,lu,d,perm,status)
    if (status /= linalg_status_success) error stop 'LDL factorization failed'
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_success) error stop 'LDL solve failed'
    print '(a,es13.5)', 'Reused LDL factors, max solve error: ', maxval(abs(x-truth))

    ! A convenience interface factors the matrix and solves in one call.
    call solve_indefinite(a,b,x,status,lower=.false.)
    if (status /= linalg_status_success) error stop 'direct solve failed'
    print '(a,es13.5)', 'Upper-triangle direct solve, max error: ', maxval(abs(x-truth))

    ! Hermitian (not complex symmetric): use the conjugate-transpose system.
    ah = cmplx(0.0_dp,0.0_dp,dp)
    ah(1,2) = cmplx(2.0_dp,1.0_dp,dp)
    ah(2,1) = conjg(ah(1,2))
    expected = [cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp)]
    bh = matmul(ah,expected)
    call solve_indefinite(ah,bh,xh,status)
    if (status /= linalg_status_success) error stop 'complex Hermitian solve failed'
    print '(a,es13.5)', 'Complex Hermitian solution error: ', maxval(abs(xh-expected))
end program linalg_indefinite_solve_example
