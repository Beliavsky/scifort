! SPDX-License-Identifier: MIT
program linalg_riccati_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_riccati, only : solve_continuous_are, solve_discrete_are
    implicit none
    real(dp) :: a(2,2), b(2,1), q(2,2), r(1,1), x(2,2), residual
    integer :: status, steps

    a = 0.0_dp
    a(1,1) = 0.4_dp
    a(1,2) = 0.3_dp
    a(2,1) = -0.1_dp
    a(2,2) = -0.8_dp
    b(:,1) = [1.0_dp, 0.5_dp]
    q = 0.0_dp
    q(1,1) = 1.0_dp
    q(2,2) = 1.0_dp
    r = 1.0_dp

    call solve_continuous_are(a,b,q,r,x,status,iterations=steps,residual=residual)
    if (status /= linalg_status_success) error stop 'CARE failed'
    print '(a,es12.4,a,i0)', 'CARE residual = ', residual, ', sign iterations = ', steps
    print '(a,2f13.8)', 'CARE X row 1: ', x(1,:)

    a = 0.25_dp * a
    call solve_discrete_are(a,b,q,r,x,status,iterations=steps,residual=residual)
    if (status /= linalg_status_success) error stop 'DARE failed'
    print '(a,es12.4,a,i0)', 'DARE residual = ', residual, ', doubling iterations = ', steps
    print '(a,2f13.8)', 'DARE X row 1: ', x(1,:)
end program linalg_riccati_example
