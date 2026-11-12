! SPDX-License-Identifier: MIT
program test_linalg_matrix_equations
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_singular, &
        linalg_status_invalid_shape, linalg_status_no_convergence
    use scifort_linalg_matrix_equations, only : solve_sylvester, &
        solve_continuous_lyapunov, solve_discrete_lyapunov
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: ar(2, 2), br(3, 3), xr(2, 3), qr(2, 3), gotr(2, 3)
    real(dp) :: a3(3, 3), q3(3, 3), x3(3, 3), true3(3, 3)
    real(dp) :: singular(1, 1), one(1, 1), result(1, 1)
    complex(dp) :: ac(2, 2), bc(3, 3), xc(2, 3), qc(2, 3), gotc(2, 3)
    complex(dp) :: a2c(2, 2), q2c(2, 2), x2c(2, 2), true2c(2, 2)
    real(dp) :: swap_a(2, 2), zero_b(1, 1), rhs_swap(2, 1), solution_swap(2, 1)
    real(dp) :: too_large(21, 21), large_result(21, 21)
    integer :: status, i

    ar = reshape([3.0_dp, 0.5_dp, 0.1_dp, 4.0_dp], [2,2])
    br = 0.0_dp
    do i = 1, 3
        br(i, i) = real(i+4, dp)
    end do
    xr = reshape([(real(i, dp) / 7.0_dp, i=1,6)], [2,3])
    qr = matmul(ar, xr) + matmul(xr, br)
    call solve_sylvester(ar, br, qr, gotr, status)
    if (status /= linalg_status_success .or. maxval(abs(gotr - xr)) > 1.0e-12_dp) &
        error stop 'real Sylvester'

    ac = cmplx(ar, 0.15_dp, dp)
    bc = cmplx(br, -0.1_dp, dp)
    xc = cmplx(xr, 0.3_dp * xr, dp)
    qc = matmul(ac, xc) + matmul(xc, bc)
    call solve_sylvester(ac, bc, qc, gotc, status)
    if (status /= linalg_status_success .or. maxval(abs(gotc - xc)) > 1.0e-12_dp) &
        error stop 'complex Sylvester'

    a3 = 0.0_dp
    a3(1,1) = -1.0_dp
    a3(2,2) = -2.0_dp
    a3(3,3) = -3.0_dp
    a3(2,1) = 0.35_dp
    true3 = reshape([(real(i, dp) / 11.0_dp, i=1,9)], [3,3])
    q3 = matmul(a3, true3) + matmul(true3, transpose(a3))
    call solve_continuous_lyapunov(a3, q3, x3, status)
    if (status /= linalg_status_success .or. maxval(abs(x3 - true3)) > 1.0e-12_dp) &
        error stop 'real continuous Lyapunov'

    a3 = 0.15_dp * a3
    q3 = true3 - matmul(matmul(a3, true3), transpose(a3))
    call solve_discrete_lyapunov(a3, q3, x3, status)
    if (status /= linalg_status_success .or. maxval(abs(x3 - true3)) > 1.0e-12_dp) &
        error stop 'real discrete Lyapunov'

    a2c = reshape([cmplx(-2.0_dp,0.5_dp,dp), cmplx(0.2_dp,0.1_dp,dp), &
        cmplx(0.0_dp,0.2_dp,dp), cmplx(-1.0_dp,-0.3_dp,dp)], [2,2])
    true2c = reshape([cmplx(0.4_dp,0.0_dp,dp), cmplx(0.1_dp,-0.3_dp,dp), &
        cmplx(0.1_dp,0.3_dp,dp), cmplx(0.8_dp,0.0_dp,dp)], [2,2])
    q2c = matmul(a2c, true2c) + matmul(true2c, transpose(conjg(a2c)))
    call solve_continuous_lyapunov(a2c, q2c, x2c, status)
    if (status /= linalg_status_success .or. maxval(abs(x2c - true2c)) > 1.0e-12_dp) &
        error stop 'complex continuous Lyapunov'

    a2c = 0.2_dp * a2c
    q2c = true2c - matmul(matmul(a2c, true2c), transpose(conjg(a2c)))
    call solve_discrete_lyapunov(a2c, q2c, x2c, status)
    if (status /= linalg_status_success .or. maxval(abs(x2c - true2c)) > 1.0e-12_dp) &
        error stop 'complex discrete Lyapunov'

    ! Force a row interchange in the pivoted linear solve.
    swap_a = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    zero_b = 0.0_dp
    rhs_swap(:,1) = [2.0_dp, 1.0_dp]
    call solve_sylvester(swap_a, zero_b, rhs_swap, solution_swap, status)
    if (status /= linalg_status_success .or. &
        maxval(abs(solution_swap(:,1) - [1.0_dp,2.0_dp])) > 1.0e-12_dp) &
        error stop 'row pivoting in Sylvester equation'

    too_large = 0.0_dp
    call solve_discrete_lyapunov(too_large, too_large, large_result, status, method='direct')
    if (status /= linalg_status_invalid_shape) error stop 'oversized equation not rejected'

    one(1,1) = 1.0_dp
    singular(1,1) = -1.0_dp
    call solve_sylvester(one, singular, one, result, status)
    if (status /= linalg_status_singular) error stop 'singular Sylvester not detected'
    call solve_discrete_lyapunov(one, one, result, status)
    if (status /= linalg_status_singular) error stop 'singular discrete not detected'
    singular = 0.0_dp
    call solve_continuous_lyapunov(singular, one, result, status)
    if (status /= linalg_status_singular) error stop 'singular continuous not detected'
    call solve_sylvester(ar, br, qr, x3, status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid shape not detected'
    ar(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_sylvester(ar, br, qr, gotr, status)
    if (status /= linalg_status_no_convergence) error stop 'NaN input not detected'
    print '(a)', 'matrix equation native tests passed'
end program test_linalg_matrix_equations
