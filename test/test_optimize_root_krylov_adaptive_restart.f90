! SPDX-License-Identifier: MIT
! Test adaptive LGMRES restart growth, stagnation rescue, and API validation.
program test_optimize_root_krylov_adaptive_restart
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), target(n), x0(n)
    type(vector_root_result) :: fixed, capped, adaptive, default_cap, &
        explicit_cap, via_newton, nonlinear, bad
    integer :: i
    integer :: failures
    integer :: problem
    integer :: nfev_calls
    integer :: jvp_calls
    integer :: precond_calls
    integer :: fixed_jvp
    integer :: adaptive_jvp

    failures = 0
    problem = 1
    x0 = 0.0_dp
    target = 1.0_dp

    ! GMRES(1) cannot progress on a skew-symmetric 2x2 operator: A*v
    ! is orthogonal to v. Adaptive growth rescues the stalled cycle by
    ! using GMRES(2), while a ceiling of one retains the old failure.
    call clear_counts()
    fixed = root(residual, x0(:2), method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_maxiter=8, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-12_dp, maxiter=5)
    call check(.not. fixed%success, 'fixed width cannot resolve skew case')
    call check(maxval(abs(fixed%fun)) > 0.1_dp, 'fixed width has large residual')
    call check(fixed%nfev == nfev_calls, 'fixed residual evaluation accounting')
    fixed_jvp = jvp_calls

    call clear_counts()
    capped = root(residual, x0(:2), method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., inner_restart_max=1, &
        inner_maxiter=8, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-12_dp, maxiter=5)
    call check(.not. capped%success, 'cap of one prevents restart growth')
    call check(all(capped%x == fixed%x), 'cap of one preserves fixed-width solution')
    call check(jvp_calls == fixed_jvp, 'cap of one preserves Jacobian work')

    call clear_counts()
    adaptive = root(residual, x0(:2), method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., inner_restart_max=2, &
        inner_maxiter=8, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-12_dp, maxiter=5)
    call check(adaptive%success, 'adaptive restart resolves skew case')
    call check(maxval(abs(adaptive%x - [1.0_dp, 2.0_dp])) < 1.0e-12_dp, &
        'adaptive restart finds expected skew root')
    call check(jvp_calls >= 3 .and. jvp_calls <= 8, &
        'adaptive restart respects finite Arnoldi budget')
    call check(adaptive%nfev == nfev_calls, 'adaptive skew residual accounting')

    call clear_counts()
    via_newton = newton_krylov(residual, x0(:2), inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., inner_restart_max=2, &
        inner_maxiter=8, outer_k=0, jacobian_product=jvp, &
        line_search='none', f_tol=1.0e-12_dp, maxiter=5)
    call check(via_newton%success, 'newton_krylov forwards adaptive restart')
    call check(maxval(abs(via_newton%fun)) < 1.0e-12_dp, &
        'newton_krylov skew residual')

    ! A poorly conditioned diagonal problem exercises more than two
    ! restart widths. It should benefit from growth without depending on
    ! particular iteration counts or tiny floating-point differences.
    problem = 2
    target = 1.0_dp
    x0 = 0.0_dp
    do i = 1, n
        a(i) = 2.0_dp**(0.6_dp * real(i - 1, dp))
    end do
    call clear_counts()
    fixed = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_maxiter=36, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-8_dp, maxiter=3)
    call check(maxval(abs(fixed%fun)) > 1.0e-3_dp, &
        'fixed one-column restart stagnates on scaled system')

    call clear_counts()
    adaptive = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., inner_restart_max=4, &
        inner_maxiter=36, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-8_dp, maxiter=3)
    call check(maxval(abs(adaptive%fun)) < 1.0e-5_dp, &
        'adaptive growth reduces scaled residual')
    call check(adaptive%nfev == nfev_calls, &
        'scaled problem residual evaluation accounting')
    adaptive_jvp = jvp_calls

    call clear_counts()
    default_cap = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., &
        inner_maxiter=36, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-8_dp, maxiter=3)
    call check(all(default_cap%x == adaptive%x), &
        'default growth ceiling is four times initial width')
    call check(jvp_calls == adaptive_jvp, 'default ceiling preserves callback count')

    call clear_counts()
    explicit_cap = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=1, inner_restart_adaptive=.true., inner_restart_max=4, &
        inner_maxiter=36, outer_k=0, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-8_dp, maxiter=3)
    call check(all(explicit_cap%x == adaptive%x), 'explicit restart cap deterministic')

    ! Analytic Jv, preconditioning, current-Jacobian image pruning, and a
    ! genuinely changing Jacobian must compose with adaptive restarts.
    problem = 3
    do i = 1, n
        target(i) = 0.8_dp + 0.1_dp * sin(real(i, dp))
    end do
    x0 = target + 0.2_dp
    call clear_counts()
    nonlinear = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=2, inner_restart_adaptive=.true., inner_restart_max=6, &
        inner_maxiter=20, outer_k=4, outer_image_prune_auto=.true., &
        inner_m=preconditioner, jacobian_product=jvp, &
        inner_orthogonalization='selective', line_search='none', &
        fatol=1.0e-9_dp, maxiter=15)
    call check(nonlinear%success, 'nonlinear preconditioned solve converges')
    call check(maxval(abs(nonlinear%x - target)) < 1.0e-8_dp, &
        'nonlinear preconditioned solution accuracy')
    call check(precond_calls > 0 .and. jvp_calls > 0, &
        'preconditioner and analytic Jv dispatched')
    call check(nonlinear%nfev == nfev_calls, 'nonlinear residual accounting')

    ! Without an analytic Jv, the wider cycles must also track all extra
    ! finite-difference residual evaluations correctly.
    call clear_counts()
    nonlinear = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=2, inner_restart_adaptive=.true., inner_restart_max=8, &
        inner_maxiter=24, outer_k=3, inner_orthogonalization='selective', &
        line_search='none', fatol=1.0e-8_dp, maxiter=15)
    call check(nonlinear%success, 'adaptive finite-difference solve converges')
    call check(maxval(abs(nonlinear%x - target)) < 1.0e-7_dp, &
        'adaptive finite-difference solution accuracy')
    call check(nonlinear%nfev == nfev_calls .and. nfev_calls > 2, &
        'adaptive finite-difference callback accounting')

    ! Invalid options must fail before invoking any residual callback.
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart_adaptive=.true.)
    call rejected(bad, 'adaptive restart without restart width')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_restart=2, inner_maxiter=10, inner_restart_adaptive=.true.)
    call rejected(bad, 'adaptive option on gmres')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=10, inner_maxiter=10, inner_restart_adaptive=.true.)
    call rejected(bad, 'adaptive option on unrestarted solver')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=2, inner_maxiter=10, inner_restart_max=4)
    call rejected(bad, 'cap without adaptive option')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=2, inner_maxiter=10, inner_restart_adaptive=.true., &
        inner_restart_max=1)
    call rejected(bad, 'cap below initial width')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=2, inner_maxiter=10, inner_restart_adaptive=.true., &
        inner_restart_max=11)
    call rejected(bad, 'cap above total budget')
    call clear_counts()
    bad = root(residual, x0, method='hybr', inner_restart_adaptive=.true.)
    call rejected(bad, 'adaptive option on non-Krylov root')
    call clear_counts()
    bad = root(residual, x0, method='hybr', inner_restart_max=4)
    call rejected(bad, 'restart cap on non-Krylov root')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_adaptive_restart: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_adaptive_restart: PASS'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition !! true when the assertion passes
        character(len=*), intent(in) :: description !! failure diagnostic
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine

    subroutine rejected(r, description)
        type(vector_root_result), intent(in) :: r !! expected invalid result
        character(len=*), intent(in) :: description !! rejection diagnostic
        call check(.not. r%success, trim(description)//' reports failure')
        call check(nfev_calls == 0 .and. r%nfev == 0, &
            trim(description)//' rejected before residual evaluation')
    end subroutine

    subroutine clear_counts()
        nfev_calls = 0
        jvp_calls = 0
        precond_calls = 0
    end subroutine

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: f(:) !! residual vector
        nfev_calls = nfev_calls + 1
        select case (problem)
        case (1)
            f(1) = -(x(2) - 2.0_dp)
            f(2) = x(1) - 1.0_dp
        case (2)
            f = a * (x - target)
        case (3)
            f = a * (x - target) + 0.01_dp * (x*x - target*target)
        end select
    end subroutine

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! point at which the Jacobian is evaluated
        real(dp), intent(in) :: vec(:) !! product direction
        real(dp), intent(out) :: product(:) !! analytic Jacobian-vector product
        logical, intent(out) :: ok !! true for valid product
        jvp_calls = jvp_calls + 1
        select case (problem)
        case (1)
            product(1) = -vec(2)
            product(2) = vec(1)
        case (2)
            product = a * vec
        case (3)
            product = a * vec + 0.02_dp * x * vec
        end select
        ok = .true.
    end subroutine

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! approximate inverse-Jacobian action
        logical, intent(out) :: ok !! true for valid product
        precond_calls = precond_calls + 1
        product = vec / (a + 0.02_dp * x0)
        ok = .true.
    end subroutine
end program test_optimize_root_krylov_adaptive_restart
