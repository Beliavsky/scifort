! SPDX-License-Identifier: MIT
! Regression for reusing current-iterate LGMRES Jacobian images during Arnoldi.
program test_optimize_root_krylov_image_reuse
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: x0(4) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
    real(dp), parameter :: target(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    type(vector_root_result) :: reference, cached, res, disabled
    integer :: fails, reference_calls, cached_calls, jvp_count
    integer :: fd_savings, cd_savings

    fails = 0
    jvp_count = 0

    ! Exact Jacobian products: nfev is unchanged while calls to J*v fall.
    reference_calls = 0
    cached_calls = 0
    reference = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_orthogonalization='selective', outer_image_prune_auto=.true.)
    reference_calls = jvp_count
    jvp_count = 0
    cached = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_orthogonalization='selective', outer_image_prune_auto=.true., &
        outer_image_reuse=.true.)
    cached_calls = jvp_count
    call check(reference%success .and. cached%success, 'automatic pruning converges with image reuse')
    call check(maxval(abs(cached%x - target)) < 1.0e-10_dp, 'correct nonlinear root')
    call check(maxval(abs(cached%x - reference%x)) < 1.0e-10_dp, &
        'caching matches recomputed-Jacobian root')
    call check(cached%nit == reference%nit .and. cached%nfev == reference%nfev, &
        'analytic residual evaluations and nonlinear iterations unchanged')
    call check(reference_calls > cached_calls, 'automatic pruning reduces analytic Jv calls')

    ! Explicitly false means exactly the legacy path, even without pruning.
    disabled = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        outer_image_prune_auto=.true., outer_image_reuse=.false.)
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        outer_image_prune_auto=.true.)
    call check(disabled%nit == res%nit .and. disabled%nfev == res%nfev .and. &
        all(disabled%x == res%x), 'explicitly disabled caching is bitwise compatible')

    ! Fixed image pruning and changing nonlinear Jacobians also share images.
    jvp_count = 0
    reference = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        f_tol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_m=preconditioner, outer_image_prune_tol=1.0e-8_dp, &
        inner_orthogonalization='twice')
    reference_calls = jvp_count
    jvp_count = 0
    cached = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        f_tol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_m=preconditioner, outer_image_prune_tol=1.0e-8_dp, &
        inner_orthogonalization='twice', outer_image_reuse=.true.)
    call check(reference%success .and. cached%success, 'fixed-image/preconditioned convergence')
    call check(maxval(abs(reference%x - cached%x)) < 1.0e-9_dp, &
        'fixed pruning with preconditioning preserves solution')
    call check(reference_calls > jvp_count, 'fixed image pruning reduces Jv calls')

    ! Central differences must save two evaluations for each reused image;
    ! forward differences save one. No images persist across Newton steps.
    reference = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=3, &
        line_search='none', outer_image_prune_auto=.true.)
    cached = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=3, &
        line_search='none', outer_image_prune_auto=.true., outer_image_reuse=.true.)
    fd_savings = reference%nfev - cached%nfev
    call check(fd_savings > 0, 'forward differences save residual evaluations')
    call check(reference%nit == cached%nit .and. &
        maxval(abs(reference%x - cached%x)) < 1.0e-7_dp, &
        'forward differences retain nonlinear behavior')

    reference = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=3, &
        line_search='none', jacobian_diff_scheme='3-point', &
        outer_image_prune_auto=.true.)
    cached = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=3, &
        line_search='none', jacobian_diff_scheme='3-point', &
        outer_image_prune_auto=.true., outer_image_reuse=.true.)
    cd_savings = reference%nfev - cached%nfev
    call check(cd_savings > 0 .and. cd_savings == 2 * fd_savings, &
        'centered differences save two evaluations per old image')
    call check(reference%nit == cached%nit .and. &
        maxval(abs(reference%x - cached%x)) < 1.0e-7_dp, &
        'centered differences retain nonlinear behavior')

    ! Preflight validation must precede the first residual evaluation.
    res = root(residual, x0, method='hybr', outer_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') > 0, 'non-Krylov rejected')
    res = root(residual, x0, method='krylov', outer_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires Jacobian-image pruning') > 0, &
        'caching without image pruning rejected')
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_image_prune_auto=.true., outer_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0, &
        'unsupported solver rejected before evaluation')
    res = root(residual, x0, method='krylov', &
        outer_image_prune_tol=0.0_dp, outer_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0, &
        'zero pruning tolerance does not enable caching')
    res = root(residual, x0, method='krylov', outer_image_reuse=.false.)
    call check(res%success, 'disabled caching allowed without pruning')

    if (fails > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_image_reuse: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_image_reuse: PASS'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! expectation
        character(len=*), intent(in) :: label !! test name
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
        real(dp), intent(in) :: x(:) !! one-dimensional iterate
        real(dp), intent(out) :: fvec(:) !! residual value
        fvec = x*x - 1.0_dp
    end subroutine residual_one

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! iterate at which J is applied
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! callback succeeded
        jvp_count = jvp_count + 1
        product = 2.0_dp*x*vec
        ok = .true.
    end subroutine jvp

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector
        real(dp), intent(out) :: product(:) !! simple inverse scaling
        logical, intent(out) :: ok !! callback succeeded
        product = 0.25_dp*vec
        ok = .true.
    end subroutine preconditioner
end program test_optimize_root_krylov_image_reuse
