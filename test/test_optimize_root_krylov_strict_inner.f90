! SPDX-License-Identifier: BSD-3-Clause
! Focused checks for opt-in true Krylov residual certification.
program test_optimize_root_krylov_strict_inner
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    character(len=8), parameter :: methods(6) = &
        [character(len=8) :: 'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres', 'tfqmr']
    real(dp), parameter :: diag(4) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    real(dp), parameter :: target(4) = [2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp]
    type(vector_root_result) :: res
    integer :: failures
    integer :: j
    integer :: updates
    integer :: oscillate_calls
    real(dp) :: second_inner_rtol

    failures = 0
    updates = 0
    oscillate_calls = 0
    second_inner_rtol = 0.0_dp

    ! The strict option works with all six built-in solvers when the
    ! requested linear accuracy is achieved; the analytic product does not
    ! add residual evaluations.
    do j = 1, size(methods)
        res = root(diagonal_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            method='krylov', inner_method=trim(methods(j)), &
            inner_check_residual=.true., jacobian_product=diagonal_jvp, &
            inner_maxiter=60, fatol=1.0e-10_dp, line_search='none')
        if (.not. res%success) print *, 'DEBUG ',trim(methods(j)),': ',trim(res%message),res%nfev
        call check(res%success .and. maxval(abs(res%x - target)) < 1.0e-9_dp, &
            trim(methods(j))//' certified built-in step')
        call check(res%nfev >= 2 .and. res%nfev <= 4, &
            trim(methods(j))//' analytic verification adds no residual evaluations')
    end do

    ! Custom solvers may claim success without meeting their requested
    ! tolerance. Verification must catch this before either the outer update
    ! or the next nonlinear residual evaluation.
    res = root(identity_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_solver=half_solver, jacobian_product=identity_jvp, &
        inner_check_residual=.true., inner_m_update=record_update, &
        line_search='none')
    call check(.not. res%success .and. res%nfev == 1 .and. &
        index(res%message, 'exceeds tolerance') /= 0, &
        'reject inaccurate custom solve before step acceptance')
    call check(updates == 0 .and. all(res%x == 0.0_dp), &
        'rejected step does not update preconditioner or nonlinear iterate')

    ! The default preserves SciPy's accepted incomplete-iterate behavior.
    res = root(identity_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_solver=half_solver, jacobian_product=identity_jvp, &
        maxiter=1, line_search='none')
    call check(res%nfev == 2 .and. all(res%x == 0.5_dp * [1.0_dp, -2.0_dp]), &
        'default retains incomplete custom inner iterate')

    ! A correct custom solver with an analytic J*v incurs no additional fun
    ! evaluations, while finite-difference verification counts every probe.
    res = newton_krylov(identity_residual, [0.0_dp, 0.0_dp], &
        inner_solver=identity_solver, jacobian_product=identity_jvp, &
        inner_check_residual=.true., line_search='none', f_tol=1.0e-12_dp)
    call check(res%success .and. res%nfev == 2, &
        'newton_krylov analytic verified solve and evaluation accounting')
    res = root(identity_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_solver=identity_solver, inner_check_residual=.true., &
        fatol=1.0e-12_dp, line_search='none')
    call check(res%success .and. res%nfev == 3, &
        'two-point verification adds one counted residual evaluation')
    res = root(identity_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_solver=identity_solver, inner_check_residual=.true., &
        jacobian_diff_scheme='3-point', fatol=1.0e-12_dp, line_search='none')
    call check(res%success .and. res%nfev == 4, &
        'three-point verification adds two counted residual evaluations')

    ! Check both huge and tiny finite residuals with IEEE exceptions enabled:
    ! the verification must scale before subtracting, and the forcing-term
    ! update must not square an extreme norm before taking its ratio.
    res = root(scalar_residual, [1.0e180_dp], method='krylov', &
        inner_solver=identity_solver, jacobian_product=identity_jvp, &
        inner_check_residual=.true., line_search='none', fatol=0.0_dp)
    call check(res%success .and. res%x(1) == 0.0_dp, &
        'large finite norm does not overflow forcing-term update')
    res = root(scalar_residual, [1.0e-180_dp], method='krylov', &
        inner_solver=identity_solver, jacobian_product=identity_jvp, &
        inner_check_residual=.true., line_search='none', fatol=0.0_dp)
    call check(res%success .and. res%x(1) == 0.0_dp, &
        'tiny finite norm does not underflow forcing-term ratio')

    ! The overflow-safe ratio must not silently alter normal SciPy forcing
    ! terms: when the residual norm stays constant, eta_a = gamma = 0.9.
    res = root(scalar_residual, [1.0_dp], method='krylov', &
        inner_solver=oscillating_solver, maxiter=2, &
        line_search='none', fatol=0.0_dp)
    call check(oscillate_calls == 2 .and. &
        abs(second_inner_rtol - 0.9_dp) < 1.0e-12_dp, &
        'nondecreasing ordinary residual preserves SciPy forcing term')

    res = root(scalar_residual, [1.0e308_dp], method='krylov', &
        inner_solver=identity_solver, jacobian_product=negative_jvp, &
        inner_check_residual=.true., line_search='none')
    call check(.not. res%success .and. res%nfev == 1 .and. &
        index(res%message, 'exceeds tolerance') /= 0, &
        'large opposite-signed products rejected without overflow')

    res = root(identity_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_solver=identity_solver, jacobian_product=failing_jvp, &
        inner_check_residual=.true., line_search='none')
    call check(.not. res%success .and. res%nfev == 1 .and. &
        index(res%message, 'verification failed') /= 0, &
        'failed analytic verification propagates')

    res = root(identity_residual, [0.0_dp, 0.0_dp], method='broyden1', &
        inner_check_residual=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') /= 0, &
        'non-Krylov method rejects strict inner checking')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_strict_inner: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_strict_inner: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected result
        character(len=*), intent(in) :: label !! test name

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine diagonal_residual(x, f)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: f(:) !! residual

        f = diag * (x - target)
    end subroutine diagonal_residual

    subroutine diagonal_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! success

        product = diag * vec
        ok = .true.
    end subroutine diagonal_jvp

    subroutine identity_residual(x, f)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: f(:) !! residual

        f = x - [1.0_dp, -2.0_dp]
    end subroutine identity_residual

    subroutine scalar_residual(x, f)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: f(:) !! residual

        f = x
    end subroutine scalar_residual

    subroutine identity_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! identity product
        logical, intent(out) :: ok !! success

        product = vec
        ok = .true.
    end subroutine identity_jvp

    subroutine negative_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! negative identity product
        logical, intent(out) :: ok !! success

        product = -vec
        ok = .true.
    end subroutine negative_jvp

    subroutine failing_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! undefined when ok is false
        logical, intent(out) :: ok !! failed evaluation

        ok = .false.
    end subroutine failing_jvp

    subroutine identity_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! operator
        procedure(vector_root_linear_operator) :: preconditioner !! preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! tolerance
        integer, intent(in) :: maxiter !! iteration limit
        real(dp), intent(out) :: solution(:) !! exact identity solution
        logical, intent(out) :: ok !! success

        solution = rhs
        ok = .true.
    end subroutine identity_solver

    subroutine half_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! operator
        procedure(vector_root_linear_operator) :: preconditioner !! preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! tolerance
        integer, intent(in) :: maxiter !! iteration limit
        real(dp), intent(out) :: solution(:) !! inaccurate step
        logical, intent(out) :: ok !! solver reports success

        solution = 0.5_dp * rhs
        ok = .true.
    end subroutine half_solver

    subroutine oscillating_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! operator
        procedure(vector_root_linear_operator) :: preconditioner !! preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! inexact-Newton forcing tolerance
        integer, intent(in) :: maxiter !! inner iteration limit
        real(dp), intent(out) :: solution(:) !! step that reverses residual sign
        logical, intent(out) :: ok !! success

        oscillate_calls = oscillate_calls + 1
        if (oscillate_calls == 2) second_inner_rtol = rtol
        solution = 2.0_dp * rhs
        ok = .true.
    end subroutine oscillating_solver

    subroutine record_update(x, fvec)
        real(dp), intent(in) :: x(:) !! accepted iterate
        real(dp), intent(in) :: fvec(:) !! accepted residual

        updates = updates + 1
    end subroutine record_update

end program test_optimize_root_krylov_strict_inner
