! SPDX-License-Identifier: MIT
! Diagonal inverse-Jacobian preconditioning through both public entry points.
program test_optimize_root_krylov_diagonal_preconditioner
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), invdiag(n), target(n), x0(n), bad_diag(n)
    real(dp) :: signed_a(n), signed_invdiag(n)
    type(vector_root_result) :: by_array, by_callback, via_newton, baseline, invalid
    character(len=8), parameter :: methods(6) = [character(len=8) :: &
        'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres', 'tfqmr']
    integer :: i, j, failures, f_calls, jv_calls, precond_calls
    logical :: signed_problem, changing_problem

    failures = 0
    signed_problem = .false.
    changing_problem = .false.
    do i = 1, n
        a(i) = 10.0_dp ** (-8.0_dp * real(i - 1, dp) / real(n - 1, dp))
        target(i) = 0.8_dp + 0.1_dp * sin(real(i, dp))
    end do
    invdiag = 1.0_dp / a
    x0 = 0.0_dp

    ! Check all six built-in solvers. A positive diagonal J and M=J**(-1)
    ! remove the eight-order-of-magnitude scale spread in the inner solve.
    do j = 1, size(methods)
        call reset_counters()
        by_array = root(residual, x0, method='krylov', &
            inner_method=trim(methods(j)), inner_maxiter=24, outer_k=0, &
            inner_m_diag=invdiag, jacobian_product=jvp, &
            inner_tfqmr_right_preconditioned=(trim(methods(j)) == 'tfqmr'), &
            line_search='none', fatol=1.0e-12_dp, maxiter=8)
        call check(by_array%success, trim(methods(j))//' diagonal preconditioned convergence')
        call check(maxval(abs(by_array%x - target)) < 1.0e-8_dp, &
            trim(methods(j))//' diagonal preconditioned accuracy')
        call check(by_array%nfev == f_calls, trim(methods(j))//' residual accounting')
        call check(jv_calls > 0, trim(methods(j))//' matrix-free operator used')
    end do

    ! Identical explicit procedure and diagonal-array actions must agree on
    ! the exact same GMRES/LGMRES operations without hidden residual calls.
    do j = 1, 2
        call reset_counters()
        by_array = root(residual, x0, method='krylov', &
            inner_method=trim(methods(j)), inner_maxiter=24, inner_restart=3, &
            inner_m_diag=invdiag, jacobian_product=jvp, &
            line_search='none', fatol=1.0e-12_dp, maxiter=8)
        call check(by_array%success, 'array-preconditioned restarted solver')
        call reset_counters()
        by_callback = root(residual, x0, method='krylov', &
            inner_method=trim(methods(j)), inner_maxiter=24, inner_restart=3, &
            inner_m=manual_preconditioner, jacobian_product=jvp, &
            line_search='none', fatol=1.0e-12_dp, maxiter=8)
        call check(by_callback%success, 'callback-preconditioned restarted solver')
        call check(all(by_callback%x == by_array%x), &
            'callback and array give identical restarted solutions')
        call check(by_callback%nfev == by_array%nfev .and. &
            by_callback%nit == by_array%nit, 'callback and array have identical work counts')
        call check(precond_calls > 0, 'callback preconditioner exercised')
    end do

    call reset_counters()
    via_newton = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=24, inner_restart=3, outer_k=2, &
        inner_m_diag=invdiag, jacobian_product=jvp, line_search='none', &
        f_tol=1.0e-12_dp, maxiter=8)
    call check(via_newton%success, 'newton_krylov forwards diagonal preconditioner')
    call check(maxval(abs(via_newton%x - target)) < 1.0e-8_dp, &
        'newton_krylov diagonal preconditioner accuracy')
    call check(via_newton%nfev == f_calls, 'newton_krylov callback accounting')

    ! The diagonal is only an approximation when J changes with x.  Its
    ! application stays static while the analytic Jacobian product updates.
    changing_problem = .true.
    call reset_counters()
    by_array = root(residual, target + 0.3_dp, method='krylov', &
        inner_method='lgmres', inner_maxiter=24, inner_restart=4, outer_k=3, &
        inner_m_diag=invdiag, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-11_dp, maxiter=10)
    call check(by_array%success, 'static diagonal also handles changing Jacobians')
    call check(maxval(abs(by_array%x - target)) < 1.0e-8_dp, &
        'changing-Jacobian preconditioner accuracy')
    call check(by_array%nfev == f_calls .and. jv_calls > 1, &
        'changing-Jacobian callback accounting')
    changing_problem = .false.

    ! A signed preconditioner is valid for nonsymmetric GMRES, but not for
    ! MINRES, where the diagonal must be positive definite.
    signed_problem = .true.
    signed_a = a
    signed_a(2:n:2) = -signed_a(2:n:2)
    signed_invdiag = 1.0_dp / signed_a
    call reset_counters()
    by_array = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_m_diag=signed_invdiag, jacobian_product=jvp, outer_k=0, &
        line_search='none', fatol=1.0e-12_dp, maxiter=8)
    call check(by_array%success, 'signed diagonal accepted by gmres')
    call check(maxval(abs(by_array%x - target)) < 1.0e-8_dp, &
        'signed diagonal gmres accuracy')
    signed_problem = .false.

    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=invdiag(:n-1))
    call rejected(invalid, 'wrong diagonal size')
    bad_diag = invdiag
    bad_diag(n) = 0.0_dp
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=bad_diag)
    call rejected(invalid, 'zero diagonal entry')
    bad_diag(n) = ieee_value(0.0_dp, ieee_quiet_nan)
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=bad_diag)
    call rejected(invalid, 'nan diagonal entry')
    bad_diag(n) = ieee_value(0.0_dp, ieee_positive_inf)
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=bad_diag)
    call rejected(invalid, 'infinite diagonal entry')
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=invdiag, &
        inner_m=manual_preconditioner)
    call rejected(invalid, 'two preconditioners conflict')
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=invdiag, &
        inner_m_update=update_preconditioner)
    call rejected(invalid, 'static diagonal and update hook conflict')
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_method='minres', &
        inner_m_diag=signed_invdiag)
    call rejected(invalid, 'minres requires positive preconditioner')
    call reset_counters()
    invalid = root(residual, x0, method='hybr', inner_m_diag=invdiag)
    call rejected(invalid, 'diagonal option only available for krylov')

    ! Overflow in a user-directed preconditioner action must be reported as
    ! a callback failure instead of raising an IEEE overflow exception.
    call reset_counters()
    invalid = root(residual, x0, method='krylov', inner_m_diag=invdiag, &
        inner_solver=overflow_probe, jacobian_product=jvp)
    call check(.not. invalid%success .and. invalid%nfev == f_calls, &
        'unrepresentable diagonal action safely fails')

    ! The existing no-preconditioner calling convention is unaffected.
    call reset_counters()
    baseline = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=3, inner_restart=2, outer_k=0, &
        jacobian_product=jvp, line_search='none', fatol=1.0e-12_dp, maxiter=3)
    call check(baseline%nfev == f_calls, 'unpreconditioned callback accounting')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_diagonal_preconditioner: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_diagonal_preconditioner: PASS'
contains
    subroutine check(condition, name)
        logical, intent(in) :: condition !! assertion value
        character(len=*), intent(in) :: name !! assertion diagnostic
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL '//name
        end if
    end subroutine check

    subroutine rejected(result, name)
        type(vector_root_result), intent(in) :: result !! expected rejected result
        character(len=*), intent(in) :: name !! invalid configuration
        call check(.not. result%success .and. result%nfev == 0 .and. f_calls == 0, &
            trim(name)//' rejected before evaluating the residual')
    end subroutine rejected

    subroutine reset_counters()
        f_calls = 0
        jv_calls = 0
        precond_calls = 0
    end subroutine reset_counters

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! residual
        f_calls = f_calls + 1
        if (signed_problem) then
            fvec = signed_a(1:size(x)) * (x - target(1:size(x)))
        else if (changing_problem) then
            fvec = a(1:size(x)) * (x - target(1:size(x))) * (1.0_dp + 0.1_dp*x)
        else
            fvec = a(1:size(x)) * (x - target(1:size(x)))
        end if
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! exact Jacobian action
        logical, intent(out) :: ok !! callback status
        jv_calls = jv_calls + 1
        if (signed_problem) then
            product = signed_a(1:size(x)) * vec
        else if (changing_problem) then
            product = a(1:size(x)) * (1.0_dp + 0.2_dp*x - 0.1_dp*target(1:size(x))) * vec
        else
            product = a(1:size(x)) * vec
        end if
        ok = .true.
    end subroutine jvp

    subroutine manual_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! approximate inverse-Jacobian action
        logical, intent(out) :: ok !! preconditioner status
        precond_calls = precond_calls + 1
        product = invdiag(1:size(vec)) * vec
        ok = .true.
    end subroutine manual_preconditioner

    subroutine overflow_probe(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! matrix-free operator
        procedure(vector_root_linear_operator) :: preconditioner !! inverse Jacobian
        real(dp), intent(in) :: rhs(:) !! supplied linear right-hand side
        real(dp), intent(in) :: rtol !! target linear tolerance
        integer, intent(in) :: maxiter !! iteration limit
        real(dp), intent(out) :: solution(:) !! unneeded solution
        logical, intent(out) :: ok !! false on failed probe
        real(dp) :: huge_direction(size(rhs)), product(size(rhs))
        logical :: action_ok

        huge_direction = huge(1.0_dp)
        call preconditioner(huge_direction, product, action_ok)
        call check(.not. action_ok, 'diagonal action detects multiplication overflow')
        solution = 0.0_dp
        ok = .false.
    end subroutine overflow_probe

    subroutine update_preconditioner(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: fvec(:) !! accepted residual
    end subroutine update_preconditioner
end program test_optimize_root_krylov_diagonal_preconditioner
