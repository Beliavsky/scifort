! SPDX-License-Identifier: MIT
! Evolving nonsymmetric rank-two preconditioners for genuinely nonlinear roots.
program test_optimize_root_krylov_rank_two_nonlinear_base
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 24
    integer, parameter :: refresh_intervals(3) = [2, 3, 4]
    real(dp) :: a(n), u(n,2), v(n,2), target(n), x0(n), target_value(n)
    real(dp) :: last_point(n), initial_diagonal(n), max_diagonal_change
    integer :: i, failures, base_calls, forward_calls, transpose_calls
    integer :: changed_points, fail_at, invalid_at, biased_at
    integer :: k, expected_calls
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        a(i) = 1.0_dp + 0.035_dp*real(i,dp)
        u(i,1) = 0.11_dp*sin(0.63_dp*real(i,dp)) + 0.05_dp
        u(i,2) = 0.13_dp*cos(0.48_dp*real(i,dp)) - 0.03_dp
        v(i,1) = 0.17_dp*cos(0.37_dp*real(i,dp)) + 0.06_dp
        v(i,2) = 0.19_dp*sin(0.31_dp*real(i,dp)) - 0.04_dp
        target(i) = 0.045_dp*cos(0.27_dp*real(i,dp))
        x0(i) = 0.6_dp*sin(0.23_dp*real(i,dp)) + 0.18_dp
    end do
    call nonlinear_action(target,target_value)

    ! The diagonal callback is refreshed alongside J(x).  The nonlinear
    ! cubic term changes D(x), and both rank-two factors change with x.
    call reset()
    fit = root(residual, x0, method='krylov', fatol=1.0e-11_dp, &
        maxiter=24, inner_method='gmres', inner_maxiter=20, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=1, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%success .and. maxval(abs(fit%fun)) < 1.0e-10_dp, &
        'nonlinear root converges with a varying known diagonal')
    call check(maxval(abs(fit%x-target)) < 1.0e-9_dp, &
        'nonlinear solution matches independently specified target')
    call check(base_calls >= 3 .and. changed_points >= 2 .and. &
        max_diagonal_change > 0.1_dp, &
        'multiple changing accepted iterates refresh the base diagonal')
    call check(forward_calls >= 4*base_calls .and. transpose_calls >= 3*base_calls, &
        'each dynamic refresh performs forward and transpose structure checks')

    ! A zero refresh interval deliberately holds the first diagonal and
    ! correction fixed even while nonlinear Newton steps change x.
    call reset()
    fit = root(residual, x0, method='krylov', fatol=1.0e-11_dp, &
        maxiter=24, inner_method='gmres', inner_maxiter=20, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=0, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%success .and. maxval(abs(fit%fun)) < 1.0e-10_dp, &
        'cached initial inverse still solves this nonlinear problem')
    call check(base_calls == 1 .and. changed_points == 0, &
        'zero update interval evaluates the diagonal callback once')

    ! Refresh intervals longer than one must follow accepted Newton steps,
    ! without refreshing on the final convergence check.  Unlike the
    ! zero-interval case, each scheduled refresh re-evaluates D(x).
    do k = 1, size(refresh_intervals)
        call reset()
        fit = root(residual, x0, method='krylov', fatol=1.0e-11_dp, &
            maxiter=24, inner_method='gmres', inner_maxiter=20, &
            inner_m_diag_auto=.true., &
            inner_m_diag_update_interval=refresh_intervals(k), &
            inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp, jacobian_transpose_product=jtv)
        call check(fit%success .and. maxval(abs(fit%fun)) < 1.0e-10_dp, &
            'nonlinear solve converges with an intermediate refresh interval')
        expected_calls = 1 + max(0,fit%nit-2)/refresh_intervals(k)
        call check(base_calls == expected_calls, &
            'base callbacks follow the selected Newton iteration schedule')
        call check(changed_points == max(0,base_calls-1), &
            'each scheduled callback uses a newly accepted nonlinear iterate')
    end do

    ! A delayed failure must be observed on the next scheduled refresh,
    ! rather than discarded or silently replaced by the cached base.
    call reset()
    fail_at = 2
    fit = root(residual, x0, method='krylov', fatol=1.0e-11_dp, &
        maxiter=24, inner_method='gmres', inner_maxiter=20, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=3, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(.not. fit%success .and. base_calls == 2 .and. &
        index(fit%message,'base callback failed') > 0, &
        'a delayed base failure terminates at the next scheduled refresh')

    ! The newton_krylov wrapper must forward the callback and refresh policy.
    call reset()
    fit = newton_krylov(residual, x0, f_tol=1.0e-11_dp, &
        maxiter=24, inner_method='gmres', inner_maxiter=20, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=1, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%success .and. maxval(abs(fit%fun)) < 1.0e-10_dp, &
        'newton_krylov forwards changing diagonal callback')
    call check(base_calls >= 3 .and. changed_points >= 2, &
        'newton_krylov refreshes D(x) on distinct accepted iterates')

    ! When a refreshed base diagonal is inaccurate, held-out Jacobian
    ! products reject the low-rank fit and use the exact diagonal instead.
    call reset()
    biased_at = 2
    fit = root(residual, x0, method='krylov', fatol=1.0e-11_dp, &
        maxiter=24, inner_method='gmres', inner_maxiter=20, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=1, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%success .and. maxval(abs(fit%fun)) < 1.0e-10_dp, &
        'an incorrect diagonal at a later iterate does not prevent convergence')
    call check(base_calls >= 3 .and. forward_calls >= n+4*base_calls, &
        'a later inaccurate base causes exact-diagonal fallback')

    ! Failed or invalid callback output is an explicit solver failure;
    ! it must not silently reuse a stale diagonal at a later point.
    call reset()
    fail_at = 2
    fit = root(residual, x0, method='krylov', maxiter=10, &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. .not. fit%success .and. base_calls == 2, &
        'failed base callback at second refresh terminates the solve')
    call check(index(fit%message,'base callback failed') > 0, &
        'failed diagonal callback reports its specific cause')

    call reset()
    invalid_at = 1
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. .not. fit%success .and. base_calls == 1, &
        'nonfinite callback diagonal is rejected')
    call check(index(fit%message,'invalid diagonal') > 0, &
        'nonfinite callback diagonal has a useful failure message')

    call reset()
    invalid_at = 2
    fit = root(residual, x0, method='krylov', maxiter=10, &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. .not. fit%success .and. base_calls == 2, &
        'invalid diagonal at a later refresh terminates the solve')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=a, &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0 .and. base_calls == 0, &
        'a fixed diagonal and a diagonal callback are mutually exclusive')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
        jacobian_product=jvp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'nonsymmetric diagonal callback requires a transpose product')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'a non-Krylov solver rejects the new callback option')

    if (failures > 0) error stop 1
    print '(a)', 'nonlinear evolving rank-two base tests passed'

contains
    subroutine reset()
        base_calls = 0
        forward_calls = 0
        transpose_calls = 0
        changed_points = 0
        fail_at = 0
        invalid_at = 0
        biased_at = 0
        last_point = x0
        initial_diagonal = 0.0_dp
        max_diagonal_change = 0.0_dp
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: description !! failure description
        if (.not. condition) then
            print '(a)', 'FAILED: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine nonlinear_action(x, values)
        real(dp), intent(in) :: x(:) !! current nonlinear vector
        real(dp), intent(out) :: values(:) !! unshifted function value
        real(dp) :: s1, s2
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        values = a*(x + 0.3_dp*x**3) + u(:,1)*sin(s1) + &
            u(:,2)*(s2+0.15_dp*s2**3)
    end subroutine nonlinear_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(out) :: fvec(:) !! residual at x
        call nonlinear_action(x,fvec)
        fvec = fvec - target_value
    end subroutine residual

    subroutine base_diagonal(x, diagonal, ok)
        real(dp), intent(in) :: x(:) !! current accepted nonlinear iterate
        real(dp), intent(out) :: diagonal(:) !! uncoupled diagonal D(x)
        logical, intent(out) :: ok !! whether D(x) could be computed
        base_calls = base_calls + 1
        if (maxval(abs(x-last_point)) > 1.0e-8_dp) changed_points = changed_points+1
        last_point = x
        diagonal = a*(1.0_dp+0.9_dp*x*x)
        if (base_calls == 1) initial_diagonal = diagonal
        max_diagonal_change = max(max_diagonal_change, &
            maxval(abs(diagonal-initial_diagonal)))
        if (base_calls == biased_at) diagonal = 1.15_dp*diagonal
        ok = base_calls /= fail_at
        if (base_calls == invalid_at) diagonal(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine base_diagonal

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear evaluation point
        real(dp), intent(in) :: direction(:) !! forward Jacobian direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! valid forward product
        real(dp) :: s1, s2
        forward_calls = forward_calls + 1
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+0.9_dp*x*x)*direction + &
            u(:,1)*cos(s1)*dot_product(v(:,1),direction) + &
            u(:,2)*(1.0_dp+0.45_dp*s2*s2)*dot_product(v(:,2),direction)
        ok = .true.
    end subroutine jvp

    subroutine jtv(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear evaluation point
        real(dp), intent(in) :: direction(:) !! transpose Jacobian direction
        real(dp), intent(out) :: product(:) !! transpose Jacobian-vector product
        logical, intent(out) :: ok !! valid transpose product
        real(dp) :: s1, s2
        transpose_calls = transpose_calls + 1
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+0.9_dp*x*x)*direction + &
            v(:,1)*cos(s1)*dot_product(u(:,1),direction) + &
            v(:,2)*(1.0_dp+0.45_dp*s2*s2)*dot_product(u(:,2),direction)
        ok = .true.
    end subroutine jtv
end program test_optimize_root_krylov_rank_two_nonlinear_base
