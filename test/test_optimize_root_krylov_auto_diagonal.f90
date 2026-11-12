! SPDX-License-Identifier: MIT
! Automatic Jacobian-diagonal preconditioning for matrix-free nonlinear solves.
program test_optimize_root_krylov_auto_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 10
    real(dp) :: scale(n), target(n), x0(n), invdiag(n)
    real(dp) :: observed(4)
    real(dp), parameter :: special_diag(4) = [0.0_dp, 1.0e-12_dp, 1.0_dp, -1.0e-4_dp]
    type(vector_root_result) :: automatic, fixed, updated, failed, finite_difference
    integer :: j, f_calls, jv_calls, solver_calls, failures
    logical :: special_mode, fail_jvp, signed_mode

    failures = 0
    special_mode = .false.
    signed_mode = .false.
    fail_jvp = .false.
    do j = 1, n
        scale(j) = 10.0_dp**(-4.0_dp * real(j - 1, dp) / real(n - 1, dp))
        target(j) = 0.5_dp + 0.1_dp * sin(real(j, dp))
    end do
    x0 = target + 0.25_dp
    invdiag = 1.0_dp / (scale * (1.0_dp + 0.1_dp * (x0 - target)))

    call reset()
    automatic = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=32, inner_restart=4, outer_k=2, &
        inner_m_diag_auto=.true., jacobian_product=jvp, &
        line_search='none', fatol=1.0e-12_dp, maxiter=14)
    call check(automatic%success, 'automatically preconditioned coupled solve converges')
    call check(maxval(abs(automatic%x - target)) < 1.0e-7_dp, 'coupled solution accuracy')
    call check(automatic%nfev == f_calls, 'analytic preconditioner has no extra residual calls')
    call check(jv_calls >= n, 'analytic probes used to estimate the diagonal')

    ! Refreshing only at the initial iterate should reproduce an explicitly
    ! supplied inverse diagonal on the same nonlinear problem.
    call reset()
    fixed = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=32, inner_restart=4, outer_k=2, inner_m_diag=invdiag, &
        jacobian_product=jvp, line_search='none', fatol=1.0e-12_dp, maxiter=14)
    call check(fixed%success, 'static diagonal comparison converges')
    call reset()
    automatic = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=32, inner_restart=4, outer_k=2, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=0, &
        jacobian_product=jvp, line_search='none', fatol=1.0e-12_dp, maxiter=14)
    call check(automatic%success, 'initial-only automatic diagonal converges')
    call check(maxval(abs(automatic%x - fixed%x)) < 1.0e-12_dp, &
        'initial-only auto diagonal agrees with explicit array')
    call check(automatic%nfev == fixed%nfev, 'analytic diagonal estimation does not evaluate residuals')

    call reset()
    updated = newton_krylov(residual, x0, inner_method='gmres', &
        inner_maxiter=32, inner_restart=4, outer_k=0, &
        inner_m_diag_auto=.true., inner_m_diag_update_interval=2, &
        jacobian_product=jvp, line_search='none', f_tol=1.0e-12_dp, maxiter=14)
    call check(updated%success, 'newton_krylov forwards automatic diagonal options')
    call check(maxval(abs(updated%x - target)) < 1.0e-7_dp, 'updated diagonal accuracy')
    call check(updated%nfev == f_calls, 'newton_krylov analytic diagonal residual accounting')

    ! In exact-iteration mode only the diagonal-estimation Jv calls are
    ! made: the custom inner solver uses the preconditioner directly.
    call reset()
    updated = newton_krylov(residual, x0, iter=3, inner_method='gmres', &
        inner_solver=quasi_newton, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=1, jacobian_product=jvp, line_search='none')
    call check(solver_calls == 3 .and. jv_calls == 3*n, 'refresh once per nonlinear solve step')
    call check(updated%nfev == f_calls .and. f_calls == 4, 'every-iteration update accounting')
    call reset()
    updated = newton_krylov(residual, x0, iter=3, inner_method='gmres', &
        inner_solver=quasi_newton, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=2, jacobian_product=jvp, line_search='none')
    call check(solver_calls == 3 .and. jv_calls == 2*n, 'refresh interval two')
    call reset()
    updated = newton_krylov(residual, x0, iter=3, inner_method='gmres', &
        inner_solver=quasi_newton, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=0, jacobian_product=jvp, line_search='none')
    call check(solver_calls == 3 .and. jv_calls == n, 'initial-only refresh interval')

    ! The built-in MINRES preconditioner must be positive definite, even
    ! when the analytic Jacobian has signed diagonal entries.
    special_mode = .true.
    signed_mode = .true.
    call reset()
    updated = root(residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='minres', inner_m_diag_auto=.true., &
        inner_maxiter=12, jacobian_product=jvp, outer_k=0, &
        line_search='none', fatol=1.0e-12_dp, maxiter=7)
    call check(updated%success, 'automatic MINRES diagonal is positive definite')
    call check(maxval(abs(updated%x - 1.0_dp)) < 1.0e-9_dp, 'signed symmetric MINRES solution')
    special_mode = .false.
    signed_mode = .false.

    ! A custom inner solver probes the actual generated diagonal.  The zero
    ! entry must use identity scaling; the tiny nonzero one is capped at 1e4.
    special_mode = .true.
    call reset()
    failed = root(residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_floor=1.0e-4_dp, &
        jacobian_product=jvp, inner_solver=probe_preconditioner)
    call check(.not. failed%success .and. solver_calls == 1, 'custom solver sees generated diagonal')
    call check(maxval(abs(observed - [1.0_dp, 1.0e4_dp, 1.0_dp, -1.0e4_dp])) < 1.0e-6_dp, &
        'zero and tiny diagonals have safe scaling')
    call check(failed%nfev == 1 .and. failed%nfev == f_calls, 'analytic probe does not evaluate residual')
    special_mode = .false.

    ! Without an analytic Jv callback, forward estimation uses one residual
    ! evaluation per coordinate, and centered estimation uses two.
    call reset()
    finite_difference = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., inner_solver=probe_preconditioner)
    call check(.not. finite_difference%success .and. solver_calls == 1, 'forward difference probe executed')
    call check(finite_difference%nfev == 1 + n .and. finite_difference%nfev == f_calls, &
        'forward finite-difference evaluation accounting')
    call reset()
    finite_difference = root(residual, x0, method='krylov', &
        inner_m_diag_auto=.true., jacobian_diff_scheme='3-point', &
        inner_solver=probe_preconditioner)
    call check(.not. finite_difference%success .and. solver_calls == 1, 'centered difference probe executed')
    call check(finite_difference%nfev == 1 + 2*n .and. finite_difference%nfev == f_calls, &
        'centered finite-difference evaluation accounting')

    ! A failed analytic derivative evaluation stops immediately and does not
    ! lead to an inner solve using a partially initialized diagonal.
    fail_jvp = .true.
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        jacobian_product=jvp, inner_solver=probe_preconditioner)
    call check(.not. failed%success .and. solver_calls == 0 .and. failed%nfev == 1, &
        'bad Jacobian product aborts estimation safely')
    fail_jvp = .false.

    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag=invdiag)
    call rejected(failed, 'manual array and auto conflict')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m=identity_preconditioner)
    call rejected(failed, 'manual callback and auto conflict')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_update=ignored_update)
    call rejected(failed, 'user update and auto conflict')
    call reset()
    failed = root(residual, x0, method='hybr', inner_m_diag_auto=.true.)
    call rejected(failed, 'auto preconditioner only supported for krylov')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_update_interval=1)
    call rejected(failed, 'interval needs auto')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=-1)
    call rejected(failed, 'negative update interval rejected')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_floor=0.0_dp)
    call rejected(failed, 'zero floor rejected')
    call reset()
    failed = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_floor=1.1_dp)
    call rejected(failed, 'floor greater than one rejected')

    if (failures /= 0) then
        print '(a,i0)', 'FAIL auto diagonal regression: ', failures
        error stop 1
    end if
    print '(a)', 'PASS automatic Krylov diagonal preconditioner'
contains
    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print '(a)', 'FAIL '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine rejected(res, description)
        type(vector_root_result), intent(in) :: res
        character(len=*), intent(in) :: description
        call check(.not. res%success .and. res%nfev == 0 .and. f_calls == 0, description)
    end subroutine rejected

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        real(dp) :: err(size(x))
        integer :: i, next
        f_calls = f_calls + 1
        if (special_mode .and. signed_mode) then
            fvec = [-2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp] * (x - 1.0_dp)
        else if (special_mode) then
            fvec = special_diag * (x - 1.0_dp)
            do i = 1, size(x)
                next = mod(i, size(x)) + 1
                fvec(i) = fvec(i) + 0.01_dp * (x(next) - 1.0_dp)
            end do
        else
            err = x - target(1:size(x))
            do i = 1, size(x)
                next = mod(i, size(x)) + 1
                fvec(i) = scale(i) * (err(i) + 0.05_dp*err(i)**2 + 0.04_dp*err(next))
            end do
        end if
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        integer :: i, next
        jv_calls = jv_calls + 1
        ok = .not. fail_jvp
        product = 0.0_dp
        if (.not. ok) return
        if (special_mode .and. signed_mode) then
            product = [-2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp] * vec
        else if (special_mode) then
            do i = 1, size(x)
                next = mod(i, size(x)) + 1
                product(i) = special_diag(i)*vec(i) + 0.01_dp*vec(next)
            end do
        else
            do i = 1, size(x)
                next = mod(i, size(x)) + 1
                product(i) = scale(i)*((1.0_dp+0.1_dp*(x(i)-target(i)))*vec(i) + &
                    0.04_dp*vec(next))
            end do
        end if
    end subroutine jvp

    subroutine probe_preconditioner(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec
        procedure(vector_root_linear_operator) :: preconditioner
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(in) :: rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        real(dp) :: probe(size(rhs)), result(size(rhs))
        logical :: precond_ok
        solver_calls = solver_calls + 1
        probe = 1.0_dp
        call preconditioner(probe, result, precond_ok)
        if (size(rhs) == 4) observed = result
        call check(precond_ok, 'generated preconditioner action succeeds')
        solution = 0.0_dp
        ok = .false.
    end subroutine probe_preconditioner

    subroutine quasi_newton(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec
        procedure(vector_root_linear_operator) :: preconditioner
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(in) :: rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        solver_calls = solver_calls + 1
        call preconditioner(rhs, solution, ok)
    end subroutine quasi_newton

    subroutine identity_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        product = vec
        ok = .true.
    end subroutine identity_preconditioner

    subroutine ignored_update(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: fvec(:)
    end subroutine ignored_update
end program test_optimize_root_krylov_auto_diagonal
