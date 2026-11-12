! SPDX-License-Identifier: MIT
! Optional coordinate-scaled and adaptive matrix-free Jacobian products.
program test_optimize_root_krylov_adaptive_steps
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: ill_x(2) = [1.0e13_dp, 0.1_dp]
    real(dp), parameter :: exact_x(4) = [1.0_dp, -0.5_dp, 2.0_dp, 0.75_dp]
    character(len=8), parameter :: methods(6) = [character(len=8) :: &
        'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres', 'tfqmr']
    type(vector_root_result) :: res
    real(dp) :: probe_dir(2), jv(2), jv_exact(2)
    logical :: probe_ok
    integer :: failures, i

    failures = 0
    probe_dir = [1.0_dp, 1.0_dp]
    jv_exact = [1.0e-12_dp + exp(0.1_dp), 0.5e-12_dp + cos(0.1_dp)]

    ! SciPy-compatible global scaling is unsuitable when one component of
    ! x is 1e13 and another is 0.1: it probes outside the residual domain.
    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='global', rdiff=1.0e-4_dp)
    call check(.not. probe_ok .and. res%nfev == 2, &
        'global directional step leaves the small-variable domain')

    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='componentwise', rdiff=1.0e-4_dp)
    call check(probe_ok .and. res%nfev == 3, &
        'componentwise forward steps cost one sample per nonzero component')
    call check(maxval(abs(jv - jv_exact)) < 1.0e-4_dp, &
        'componentwise forward directional derivative')

    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='adaptive', rdiff=1.0e-4_dp)
    call check(probe_ok .and. res%nfev == 3, &
        'adaptive mode detects ill-scaled active components')
    call check(maxval(abs(jv - jv_exact)) < 1.0e-4_dp, &
        'adaptive forward directional derivative')

    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='adaptive', jacobian_diff_scheme='3-point', &
        rdiff=1.0e-4_dp)
    call check(probe_ok .and. res%nfev == 5, &
        'adaptive centered steps cost two samples per active component')
    call check(maxval(abs(jv - jv_exact)) < 1.0e-8_dp, &
        'centered component-scaled directional derivative accuracy')

    ! An inactive huge component must not determine a nonzero coordinate step.
    probe_dir = [0.0_dp, 1.0_dp]
    jv_exact = [exp(0.1_dp), cos(0.1_dp)]
    res = newton_krylov(ill_scaled, ill_x, inner_solver=probe, &
        jacobian_step_mode='adaptive', jacobian_diff_scheme='3-point', &
        rdiff=1.0e-4_dp)
    call check(probe_ok .and. res%nfev == 3, &
        'newton_krylov forwards adaptive scaling for a sparse direction')
    call check(maxval(abs(jv - jv_exact)) < 1.0e-8_dp, &
        'inactive components do not contaminate Jacobian action')

    ! With huge residual offsets, the global norm normalization can make
    ! directional steps too small to change the function even though x moves.
    probe_dir = [1.0_dp, 1.0_dp]
    jv_exact = [1.0_dp, 2.0_dp]
    res = root(large_offset, [0.5_dp, -0.25_dp], method='krylov', &
        inner_solver=probe, jacobian_step_mode='global', rdiff=1.0e-2_dp)
    call check(res%nfev == 2 .and. maxval(abs(jv - jv_exact)) > 0.5_dp, &
        'global difference loses accuracy against a very large residual offset')
    res = root(large_offset, [0.5_dp, -0.25_dp], method='krylov', &
        inner_solver=probe, jacobian_step_mode='adaptive', rdiff=1.0e-2_dp)
    call check(probe_ok .and. res%nfev == 3, &
        'adaptive mode detects a globally undersized residual-normalized step')
    call check(maxval(abs(jv - jv_exact)) < 3.0e-2_dp, &
        'coordinate-scaled differences recover derivatives with large offsets')

    ! All-zero vectors cost no residual evaluations and are exactly zero.
    probe_dir = 0.0_dp
    jv_exact = 0.0_dp
    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='componentwise')
    call check(probe_ok .and. res%nfev == 1 .and. all(jv == 0.0_dp), &
        'coordinatewise zero direction is exact and free')

    ! On well-scaled systems adaptive mode retains the legacy global stencil.
    probe_dir = [1.0_dp, 1.0_dp]
    jv_exact = [2.0_dp, 3.0_dp]
    res = root(well_scaled, [0.5_dp, -0.25_dp], method='krylov', &
        inner_solver=probe, jacobian_step_mode='adaptive')
    call check(probe_ok .and. res%nfev == 2, &
        'adaptive mode retains one-evaluation global stencil when scaled')
    call check(maxval(abs(jv - jv_exact)) < 1.0e-7_dp, &
        'well-scaled global directional derivative')
    res = root(well_scaled, [0.5_dp, -0.25_dp], method='krylov', &
        inner_solver=probe, jacobian_step_mode='adaptive', &
        jacobian_diff_scheme='3-point')
    call check(probe_ok .and. res%nfev == 3, &
        'adaptive mode retains two-evaluation centered stencil when scaled')

    ! Extremely small coordinatewise steps are raised above the rounding floor.
    res = root(well_scaled, [0.5_dp, -0.25_dp], method='krylov', &
        inner_solver=probe, jacobian_step_mode='componentwise', rdiff=1.0e-30_dp)
    call check(probe_ok .and. res%nfev == 3, &
        'componentwise mode avoids unrepresentable tiny coordinate steps')

    ! The same option must reach all six built-in Krylov inner solvers.
    do i = 1, size(methods)
        res = root(linear_four, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            method='krylov', inner_method=trim(methods(i)), &
            jacobian_step_mode='componentwise', jacobian_diff_scheme='3-point', &
            line_search='none', inner_maxiter=30, fatol=1.0e-8_dp)
        call check(res%success .and. maxval(abs(res%fun)) < 1.0e-8_dp, &
            trim(methods(i))//' works with componentwise centered products')
        call check(maxval(abs(res%x - exact_x)) < 1.0e-7_dp, &
            trim(methods(i))//' solution matches linear-system reference')
    end do

    probe_dir = [1.0_dp, 1.0_dp]
    res = root(one_sided_domain, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='componentwise', rdiff=1.0e-4_dp)
    call check(.not. probe_ok .and. res%nfev == 3, &
        'a nonfinite second component sample is counted and rejected')
    res = root(one_sided_domain, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='componentwise', jacobian_diff_scheme='3-point', &
        rdiff=1.0e-4_dp)
    call check(.not. probe_ok .and. res%nfev == 4, &
        'a nonfinite centered component sample is counted and rejected')

    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='componentwise', rdiff=huge(1.0_dp))
    call check(.not. probe_ok .and. res%nfev == 1, &
        'an excessive componentwise step fails without overflow or extra calls')
    res = root(ill_scaled, ill_x, method='krylov', inner_solver=probe, &
        jacobian_step_mode='global', rdiff=huge(1.0_dp))
    call check(.not. res%success .and. res%nfev == 1 .and. &
        index(res%message, 'exceeds representable range') > 0, &
        'an excessive global step fails without overflowing')

    res = root(ill_scaled, ill_x, method='broyden1', jacobian_step_mode='adaptive')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') > 0, &
        'non-Krylov methods reject jacobian_step_mode')
    res = root(ill_scaled, ill_x, method='krylov', jacobian_step_mode='unknown')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'must be global') > 0, &
        'unrecognized step mode is rejected before evaluating the residual')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_adaptive_steps: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_adaptive_steps: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! true when test passes
        character(len=*), intent(in) :: label !! description of the test

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine ill_scaled(x, f)
        real(dp), intent(in) :: x(:) !! point with widely separated scales
        real(dp), intent(out) :: f(:) !! smooth nonlinear residual in its domain

        if (abs(x(2)) > 30.0_dp) then
            f = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            f(1) = x(1) / 1.0e12_dp + exp(x(2)) - 10.0_dp
            f(2) = 0.5_dp * x(1) / 1.0e12_dp + sin(x(2)) - 4.0_dp
        end if
    end subroutine ill_scaled

    subroutine one_sided_domain(x, f)
        real(dp), intent(in) :: x(:) !! variables on a one-sided domain
        real(dp), intent(out) :: f(:) !! finite residual only for x(2) <= 0.1

        if (x(2) > 0.1_dp) then
            f = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            f = [x(1) / 1.0e13_dp + 1.0_dp, sin(x(2)) + 1.0_dp]
        end if
    end subroutine one_sided_domain

    subroutine large_offset(x, f)
        real(dp), intent(in) :: x(:) !! moderate variables with large residual offsets
        real(dp), intent(out) :: f(:) !! residual with small variable dependence

        f = [1.0e12_dp + x(1), 2.0e12_dp + 2.0_dp * x(2)]
    end subroutine large_offset

    subroutine well_scaled(x, f)
        real(dp), intent(in) :: x(:) !! well-scaled independent variables
        real(dp), intent(out) :: f(:) !! linear residual

        f = [2.0_dp * x(1), 3.0_dp * x(2)]
    end subroutine well_scaled

    subroutine linear_four(x, f)
        real(dp), intent(in) :: x(:) !! four linear-system variables
        real(dp), intent(out) :: f(:) !! symmetric positive-definite residual

        f = [4.0_dp * (x(1) - exact_x(1)), &
             5.0_dp * (x(2) - exact_x(2)), &
             6.0_dp * (x(3) - exact_x(3)), &
             7.0_dp * (x(4) - exact_x(4))]
    end subroutine linear_four

    subroutine probe(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian vector product
        procedure(vector_root_linear_operator) :: preconditioner !! unused inverse
        real(dp), intent(in) :: rhs(:) !! Newton residual supplied to inner solver
        real(dp), intent(in) :: rtol !! requested relative inner tolerance
        integer, intent(in) :: maxiter !! inner iteration limit
        real(dp), intent(out) :: solution(:) !! no solve performed in this test
        logical, intent(out) :: ok !! deliberately fail after probing the product

        call matvec(probe_dir, jv, probe_ok)
        solution = 0.0_dp
        ok = .false.
    end subroutine probe

end program test_optimize_root_krylov_adaptive_steps
