! SPDX-License-Identifier: MIT
program linalg_balance_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_balance, only : matrix_balance
    implicit none

    real(dp) :: a(4,4), b(4,4), t(4,4), factors(4), error
    complex(dp) :: ac(4,4), bc(4,4), tc(4,4)
    integer :: i, j, status, perm(4)

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp
    a(4,4) = 4.0_dp
    a(2,3) = 2.0_dp**20
    a(3,2) = 2.0_dp**(-20)
    a(3,4) = 6.0_dp
    a(4,3) = 0.3_dp
    a(1,4) = 1.0_dp
    call matrix_balance(a,b,status,t,scale_factors=factors,permutation=perm)
    if (status /= linalg_status_success) error stop 'real balance failed'
    error = maxval(abs(matmul(a,t)-matmul(t,b)))
    print '(A,ES12.4)', 'real A*T - T*B residual: ',error
    print '(A,4(1X,I0))', 'permutation:',perm
    print '(A,4(1X,ES11.3))', 'scale factors:',factors

    ac = cmplx(a,0.11_dp*a,dp)
    call matrix_balance(ac,bc,status,tc)
    if (status /= linalg_status_success) error stop 'complex balance failed'
    error = maxval(abs(matmul(ac,tc)-matmul(tc,bc)))
    print '(A,ES12.4)', 'complex A*T - T*B residual: ',error
end program linalg_balance_example
