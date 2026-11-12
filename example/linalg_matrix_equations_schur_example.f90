! SPDX-License-Identifier: MIT
program linalg_matrix_equations_schur_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_matrix_equations, only : solve_sylvester, &
        solve_continuous_lyapunov, solve_discrete_lyapunov
    implicit none
    integer, parameter :: n = 24
    real(dp) :: a(n,n), b(n,n), q(n,n), x(n,n), target(n,n)
    real(dp) :: difference
    integer :: i, j, status

    a = 0.0_dp
    b = 0.0_dp
    do i = 1, n
        a(i,i) = -1.0_dp-0.01_dp*real(i,dp)
        b(i,i) = 2.0_dp+0.015_dp*real(i,dp)
        if (i < n) a(i,i+1) = 0.1_dp
    end do
    do j = 1, n
        do i = 1, n
            target(i,j) = cos(0.12_dp*real(i+j,dp))
        end do
    end do
    q = matmul(a,target)+matmul(target,b)
    call solve_sylvester(a,b,q,x,status,method='schur')
    if (status /= linalg_status_success) error stop 'Sylvester solve failed'
    difference = maxval(abs(x-target))
    print '(a,es12.4)', 'Sylvester max error: ', difference

    q = matmul(a,target)+matmul(target,transpose(a))
    call solve_continuous_lyapunov(a,q,x,status) ! 'auto' uses Schur for n>20
    if (status /= linalg_status_success) error stop 'continuous Lyapunov failed'
    print '(a,es12.4)', 'Continuous Lyapunov max error: ', maxval(abs(x-target))

    a = 0.20_dp*a
    q = target-matmul(matmul(a,target),transpose(a))
    call solve_discrete_lyapunov(a,q,x,status,method='schur')
    if (status /= linalg_status_success) error stop 'discrete Lyapunov failed'
    print '(a,es12.4)', 'Discrete Lyapunov max error: ', maxval(abs(x-target))
end program linalg_matrix_equations_schur_example
