! SPDX-License-Identifier: MIT
program linalg_matrix_functions_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_matrix_functions, only : sqrtm, signm
    implicit none

    real(dp) :: a(2, 2), s(2, 2), residual
    complex(dp) :: root(2, 2)
    integer :: status

    a = 0.0_dp
    a(1, 1) = 4.0_dp
    a(2, 2) = -9.0_dp
    call sqrtm(a, root, status, residual)
    if (status /= linalg_status_success) error stop 'sqrtm failed'
    print '(a,es12.4)', 'sqrtm residual = ', residual
    print '(a,2f12.6)', 'sqrtm A(2,2) = ', root(2, 2)

    a = 0.0_dp
    a(1, 1) = -2.0_dp
    a(1, 2) = 1.0_dp
    a(2, 2) = 3.0_dp
    call signm(a, s, status, residual=residual)
    if (status /= linalg_status_success) error stop 'signm failed'
    print '(a,es12.4)', 'signm residual = ', residual
    print '(a,f12.6)', 'signm A(1,2) = ', s(1, 2)
end program linalg_matrix_functions_example
