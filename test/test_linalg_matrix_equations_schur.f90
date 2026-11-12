! SPDX-License-Identifier: MIT
program test_linalg_matrix_equations_schur
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_matrix_equations, only : solve_sylvester, &
        solve_continuous_lyapunov, solve_discrete_lyapunov, &
        matrix_equation_schur_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: m = 35, n = 29, k = 34
    real(dp) :: a(m,m), b(n,n), q(m,n), x(m,n), truth(m,n)
    real(dp) :: undersized(m-1,n), beyond(41,41), rejected(41,41)
    real(dp) :: empty(0,0)
    complex(dp) :: ac(k,k), qc(k,k), xc(k,k), truthc(k,k)
    complex(dp) :: ac_scalar(1,1), q_scalar(1,1), x_scalar(1,1)
    integer :: status, i, j

    if (matrix_equation_schur_max_order < 40) error stop 'Schur order limit too low'
    a = 0.0_dp
    b = 0.0_dp
    do i = 1, m
        a(i,i) = 1.0_dp + 0.017_dp*real(i,dp)
        if (i < m) a(i,i+1) = 0.12_dp
    end do
    do j = 1, n
        b(j,j) = 2.0_dp + 0.01_dp*real(j,dp)
        if (j < n) b(j,j+1) = -0.21_dp
    end do
    do j = 1, n
        do i = 1, m
            truth(i,j) = sin(0.13_dp*real(i,dp)+0.07_dp*real(j,dp))
        end do
    end do
    q = matmul(a, truth) + matmul(truth,b)
    ! Both the default auto-selection and explicit Schur method handle >400 unknowns.
    call solve_sylvester(a,b,q,x,status)
    if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 2.0e-11_dp) &
        error stop 'auto Schur Sylvester'
    call solve_sylvester(a,b,q,x,status,method='schur')
    if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 2.0e-11_dp) &
        error stop 'explicit Schur Sylvester'
    call solve_sylvester(a,b,q,x,status,method='direct')
    if (status /= linalg_status_invalid_shape) error stop 'direct cap not enforced'
    call solve_sylvester(a,b,q,undersized,status,method='schur')
    if (status /= linalg_status_invalid_shape) error stop 'bad Sylvester output shape'
    call solve_sylvester(a,b,q,x,status,method='unknown')
    if (status /= linalg_status_invalid_shape) error stop 'bad method'

    ac = cmplx(0.0_dp, 0.0_dp,dp)
    do i = 1, k
        ac(i,i) = cmplx(-1.0_dp-0.015_dp*real(i,dp), 0.006_dp*real(i,dp),dp)
        if (i < k) ac(i,i+1) = cmplx(0.15_dp,-0.02_dp,dp)
    end do
    do j = 1, k
        do i = 1, k
            truthc(i,j) = cmplx(cos(0.1_dp*real(i,dp)), sin(0.07_dp*real(j,dp)),dp)
        end do
    end do
    qc = matmul(ac,truthc)+matmul(truthc,transpose(conjg(ac)))
    call solve_continuous_lyapunov(ac,qc,xc,status)
    if (status /= linalg_status_success .or. maxval(abs(xc-truthc)) > 2.0e-11_dp) &
        error stop 'complex continuous Schur Lyapunov'
    ac = ac*0.12_dp
    qc = truthc-matmul(matmul(ac,truthc),transpose(conjg(ac)))
    call solve_discrete_lyapunov(ac,qc,xc,status,method='schur')
    if (status /= linalg_status_success .or. maxval(abs(xc-truthc)) > 2.0e-11_dp) &
        error stop 'complex discrete Schur Lyapunov'
    call solve_discrete_lyapunov(ac,qc,xc,status)
    if (status /= linalg_status_success .or. maxval(abs(xc-truthc)) > 2.0e-11_dp) &
        error stop 'auto discrete Schur Lyapunov'

    ! True resonance must remain a failure rather than produce a division by zero.
    ac_scalar(1,1) = cmplx(0.0_dp,1.0_dp,dp)
    q_scalar(1,1) = cmplx(1.0_dp,0.0_dp,dp)
    call solve_continuous_lyapunov(ac_scalar,q_scalar,x_scalar,status,method='schur')
    if (status /= linalg_status_singular) error stop 'continuous resonance'
    ac_scalar(1,1) = cmplx(1.0_dp,0.0_dp,dp)
    call solve_discrete_lyapunov(ac_scalar,q_scalar,x_scalar,status,method='schur')
    if (status /= linalg_status_singular) error stop 'discrete resonance'
    call solve_sylvester(ac_scalar,-ac_scalar,q_scalar,x_scalar,status,method='schur')
    if (status /= linalg_status_singular) error stop 'Sylvester resonance'

    beyond = 0.0_dp
    call solve_discrete_lyapunov(beyond,beyond,rejected,status,method='schur')
    if (status /= linalg_status_invalid_shape) error stop 'Schur order cap not enforced'
    call solve_discrete_lyapunov(empty,empty,empty,status)
    if (status /= linalg_status_invalid_shape) error stop 'empty matrix not rejected'
    a(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call solve_sylvester(a,b,q,x,status,method='schur')
    if (status /= linalg_status_no_convergence) error stop 'nonfinite data not rejected'

    print '(a)', 'Schur matrix-equation native tests passed'
end program test_linalg_matrix_equations_schur
