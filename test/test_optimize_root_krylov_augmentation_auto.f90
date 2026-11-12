! SPDX-License-Identifier: MIT
! Automatic Jacobian-image pruning through root and newton_krylov.
program test_optimize_root_krylov_augmentation_auto
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: x0(4) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
    real(dp), parameter :: target(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    type(vector_root_result) :: baseline
    type(vector_root_result) :: disabled
    type(vector_root_result) :: filtered
    type(vector_root_result) :: res
    integer :: fails
    integer :: jvp_calls

    fails = 0
    jvp_calls = 0
    baseline = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none')
    disabled = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        outer_image_prune_auto=.false.)
    call check(baseline%success .and. baseline%nit > 3, 'baseline nonlinear solve')
    call check((disabled%success .eqv. baseline%success) .and. &
        disabled%nfev == baseline%nfev .and. disabled%nit == baseline%nit .and. &
        all(disabled%x == baseline%x), 'false option preserves default exactly')

    jvp_calls = 0
    filtered = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_orthogonalization='selective', outer_image_prune_auto=.true.)
    call check(filtered%success .and. maxval(abs(filtered%x - target)) < 1.0e-10_dp, &
        'auto-filtered nonlinear root converges')
    call check(jvp_calls > 2*filtered%nit, 'auto mode computes extra analytic images')

    res = newton_krylov(residual, x0, inner_method='lgmres', inner_maxiter=4, &
        outer_k=3, jacobian_product=jvp, f_tol=1.0e-10_dp, maxiter=25, &
        line_search='none', outer_image_prune_auto=.true.)
    call check(res%success .and. maxval(abs(res%x - target)) < 1.0e-10_dp, &
        'newton_krylov forwards automatic pruning')

    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, fatol=1.0e-7_dp, maxiter=25, &
        line_search='none', outer_image_prune_auto=.true.)
    call check(res%success .and. maxval(abs(res%x - target)) < 1.0e-6_dp, &
        'finite-difference Jacobian products support auto pruning')
    call check(res%nfev > filtered%nfev, 'extra finite-difference images counted')

    ! A one-step solve gets one additional numerical Jacobian product for
    ! the candidate image; two-point and three-point stencils differ in cost.
    baseline = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none')
    res = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', outer_image_prune_auto=.true.)
    call check(res%nfev == baseline%nfev + 1, 'auto forward-difference accounting')
    baseline = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', jacobian_diff_scheme='3-point')
    res = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', jacobian_diff_scheme='3-point', outer_image_prune_auto=.true.)
    call check(res%nfev == baseline%nfev + 2, 'auto centered-difference accounting')

    ! A failing callback in the extra candidate-image evaluation must
    ! fail the root solve rather than retain an undefined direction.
    jvp_calls = 0
    res = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, &
        jacobian_product=failing_jvp, outer_image_prune_auto=.true., &
        line_search='none')
    call check(.not. res%success .and. res%nfev == 1 .and. jvp_calls == 2, &
        'automatic pruning propagates image callback failure')

    ! Invalid combinations are rejected before any residual evaluations.
    res = root(residual, x0, method='hybr', outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') > 0, 'non-Krylov rejected')
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in lgmres') > 0, 'GMRES rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_k=0, outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0, 'outer_k=0 rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_image_prune_tol=1.0e-8_dp, outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'cannot be combined') > 0, 'fixed image cutoff rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_prune_tol=1.0e-8_dp, outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'cannot be combined') > 0, 'solution-space cutoff rejected')

    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_solver=unused_inner_solver, outer_image_prune_auto=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in lgmres') > 0, &
        'custom inner solver rejected')

    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_augmentation_auto: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_augmentation_auto: PASS'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! expectation
        character(len=*), intent(in) :: label !! description
        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual vector
        fvec = x*x - target*target
    end subroutine residual

    subroutine residual_one(x, fvec)
        real(dp), intent(in) :: x(:) !! one-dimensional nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual value
        fvec = x*x - 1.0_dp
    end subroutine residual_one

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: vec(:) !! perturbation direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! status
        jvp_calls = jvp_calls + 1
        product = 2.0_dp*x*vec
        ok = .true.
    end subroutine jvp

    subroutine failing_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! undefined after failure
        logical, intent(out) :: ok !! callback success

        jvp_calls = jvp_calls + 1
        if (jvp_calls > 1) then
            ok = .false.
            return
        end if
        product = 2.0_dp*x*vec
        ok = .true.
    end subroutine failing_jvp

    subroutine unused_inner_solver(matvec, preconditioner, rhs, rtol, &
            maxiter, sol, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! tolerance
        integer, intent(in) :: maxiter !! iteration cap
        real(dp), intent(out) :: sol(:) !! linear solution
        logical, intent(out) :: ok !! status

        sol = rhs
        ok = .true.
    end subroutine unused_inner_solver
end program test_optimize_root_krylov_augmentation_auto
