! SPDX-License-Identifier: MIT
program test_linalg_expm_condition
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_expm_condition, only : expm_cond
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(3,3), x, scalar(1,1), invalid(2,3)
    complex(dp) :: c(3,3), z
    integer :: status, count, j

    a = 0.0_dp
    call expm_cond(a, x, status, iterations=count)
    if (status /= linalg_status_success .or. x /= 0.0_dp .or. count /= 0) &
        error stop 'zero matrix'
    a = 0.0_dp
    do j = 1, 3
        a(j,j) = -2.5_dp
    end do
    call expm_cond(a, x, status, iterations=count)
    if (status /= linalg_status_success .or. abs(x-2.5_dp) > 1.0e-10_dp) &
        error stop 'scalar multiple of identity'
    a = 0.0_dp
    do j = 1, 3
        a(j,j) = -300.0_dp
    end do
    call expm_cond(a, x, status)
    if (status /= linalg_status_success .or. abs(x-300.0_dp) > 1.0e-10_dp) &
        error stop 'large negative scalar shift'
    do j = 1, 3
        a(j,j) = 600.0_dp
    end do
    call expm_cond(a, x, status)
    if (status /= linalg_status_success .or. abs(x-600.0_dp) > 1.0e-10_dp) &
        error stop 'large positive scalar shift'
    ! Diagonal matrix with a very widely separated spectrum.
    a = 0.0_dp
    a(1,1) = -1000.0_dp
    a(2,2) = 1000.0_dp
    a(3,3) = -2.0_dp
    call expm_cond(a, x, status)
    if (status /= linalg_status_success .or. x <= 1000.0_dp) &
        error stop 'diagonal with extreme eigenvalues'
    scalar(1,1) = 1.5_dp
    call expm_cond(scalar, x, status)
    if (status /= linalg_status_success .or. abs(x-1.5_dp) > 1.0e-12_dp) &
        error stop 'one by one real'
    c = cmplx(0.0_dp, 0.0_dp, dp)
    z = cmplx(0.8_dp, -1.2_dp, dp)
    do j = 1, 3
        c(j,j) = z
    end do
    call expm_cond(c, x, status)
    if (status /= linalg_status_success .or. abs(x-abs(z)) > 1.0e-10_dp) &
        error stop 'complex scalar multiple of identity'
    a = 0.0_dp
    a(1,2) = 0.7_dp
    a(2,3) = -0.8_dp
    a(3,1) = -0.2_dp
    call expm_cond(a, x, status, tol=1.0e-8_dp)
    if (status /= linalg_status_success .or. x <= 0.0_dp) &
        error stop 'nonnormal real'
    call expm_cond(cmplx(a, kind=dp), x, status)
    if (status /= linalg_status_success .or. x <= 0.0_dp) &
        error stop 'nonnormal complex'
    invalid = 0.0_dp
    call expm_cond(invalid, x, status)
    if (status /= linalg_status_invalid_shape .or. x /= 0.0_dp) &
        error stop 'invalid shape'
    call expm_cond(a, x, status, tol=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'invalid tolerance'
    call expm_cond(a, x, status, max_iterations=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid iterations'
    a(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call expm_cond(a, x, status)
    if (status /= linalg_status_invalid_shape) error stop 'nonfinite input'
    a(1,1) = 0.0_dp
    call expm_cond(a, x, status, tol=1.0e-16_dp, max_iterations=1)
    if (status /= linalg_status_no_convergence .or. x /= 0.0_dp) &
        error stop 'iteration exhaustion'
    print '(a)', 'PASS test_linalg_expm_condition'
end program test_linalg_expm_condition
