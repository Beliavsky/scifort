! SPDX-License-Identifier: MIT
program linalg_logarithm_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_logarithm, only : logm, fractional_matrix_power
    implicit none

    real(dp) :: a(2,2), error
    complex(dp) :: logarithm(2,2), cube_root(2,2)
    integer :: status, roots, terms

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(1,2) = 3.0_dp
    a(2,1) = 1.0_dp
    a(2,2) = 4.0_dp

    call logm(a, logarithm, status, residual=error, square_roots=roots, terms=terms)
    if (status /= linalg_status_success) error stop 'logm failed'
    print '(a,es12.4)', 'relative exp(log(A))-A residual: ', error
    print '(a,i0,a,i0)', 'number of square roots: ', roots, ', atanh series terms: ', terms
    print '(a,2f13.7)', 'log(A) first row: ', real(logarithm(1,:), dp)

    call fractional_matrix_power(a, 1.0_dp/3.0_dp, cube_root, status)
    if (status /= linalg_status_success) error stop 'fractional_matrix_power failed'
    print '(a,2f13.7)', 'A**(1/3) first row: ', real(cube_root(1,:), dp)
    print '(a,es12.4)', 'cube reconstruction: ', &
        maxval(abs(matmul(cube_root, matmul(cube_root, cube_root))-cmplx(a, kind=dp)))
end program linalg_logarithm_example
