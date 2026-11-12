! SPDX-License-Identifier: MIT
program linalg_matrix_equations_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_matrix_equations, only : solve_continuous_lyapunov
    implicit none
    real(dp) :: a(2, 2)
    real(dp) :: q(2, 2)
    real(dp) :: x(2, 2)
    integer :: status

    a = reshape([-1.0_dp, 0.0_dp, 0.25_dp, -2.0_dp], [2, 2])
    q = 0.0_dp
    q(1, 1) = -1.0_dp
    q(2, 2) = -1.0_dp
    call solve_continuous_lyapunov(a, q, x, status)
    if (status /= linalg_status_success) error stop 'Lyapunov solve failed'
    print '(a,es12.4)', 'max residual: ', maxval(abs(matmul(a, x) + &
        matmul(x, transpose(a)) - q))
    print '(a)', 'X ='
    print '(2f13.7)', x(1, :), x(2, :)
end program linalg_matrix_equations_example
