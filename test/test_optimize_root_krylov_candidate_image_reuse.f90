! SPDX-License-Identifier: MIT
! Verify opt-in LGMRES candidate Jacobian-image reconstruction from Arnoldi products.
program test_optimize_root_krylov_candidate_image_reuse
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 8
    real(dp), parameter :: x0(n) = [3.0_dp, 3.2_dp, 3.4_dp, 3.6_dp, &
        3.8_dp, 4.0_dp, 4.2_dp, 4.4_dp]
    real(dp), parameter :: target(n) = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp, &
        1.4_dp, 1.5_dp, 1.6_dp, 1.7_dp]
    type(vector_root_result) :: legacy, cached, res, both, explicit_false
    integer :: fails, jvp_calls, legacy_calls, candidate_calls, old_calls, both_calls
    integer :: residual_calls

    fails = 0
    jvp_calls = 0
    residual_calls = 0

    legacy = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_auto=.true., &
        jacobian_product=jvp, line_search='none', maxiter=25, &
        fatol=1.0e-11_dp, inner_orthogonalization='selective')
    legacy_calls = jvp_calls
    call check(legacy%nfev == residual_calls, 'legacy nfev matches residual callbacks')
    call check(legacy%success, 'legacy converges')

    jvp_calls = 0
    residual_calls = 0
    cached = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_auto=.true., &
        jacobian_product=jvp, line_search='none', maxiter=25, &
        fatol=1.0e-11_dp, inner_orthogonalization='selective', &
        outer_candidate_image_reuse=.true.)
    candidate_calls = jvp_calls
    call check(cached%success .and. cached%nfev == residual_calls, &
        'candidate reuse converges and counts residual callbacks')
    call check(maxval(abs(cached%x - legacy%x)) < 1.0e-10_dp, &
        'candidate image reconstruction preserves nonlinear solution')
    call check(cached%nit == legacy%nit .and. cached%nfev == legacy%nfev, &
        'candidate reconstruction preserves iteration and residual counts')
    call check(candidate_calls < legacy_calls, 'reconstructing candidate saves Jv calls')
    call check(legacy_calls - candidate_calls >= 2, &
        'candidate savings occur on multiple nonlinear iterations')

    jvp_calls = 0
    explicit_false = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_auto=.true., &
        jacobian_product=jvp, line_search='none', maxiter=25, &
        fatol=1.0e-11_dp, inner_orthogonalization='selective', &
        outer_candidate_image_reuse=.false.)
    call check(jvp_calls == legacy_calls .and. &
        all(explicit_false%x == legacy%x) .and. &
        explicit_false%nfev == legacy%nfev, 'disabled candidate mode exactly retains default')

    ! The new optimization composes with cached retained-vector images.
    jvp_calls = 0
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_auto=.true., &
        jacobian_product=jvp, line_search='none', maxiter=25, &
        fatol=1.0e-11_dp, inner_orthogonalization='twice', &
        outer_image_reuse=.true.)
    old_calls = jvp_calls
    jvp_calls = 0
    both = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_auto=.true., &
        jacobian_product=jvp, line_search='none', maxiter=25, &
        fatol=1.0e-11_dp, inner_orthogonalization='twice', &
        outer_image_reuse=.true., outer_candidate_image_reuse=.true.)
    both_calls = jvp_calls
    call check(res%success .and. both%success, 'combined reuse converges')
    call check(maxval(abs(res%x - both%x)) < 1.0e-10_dp, &
        'combined reuse matches previously cached solution')
    call check(both_calls < old_calls .and. both_calls < candidate_calls, &
        'combined reuse saves additional Jacobian actions')

    ! The new option also works through newton_krylov and fixed image pruning.
    jvp_calls = 0
    legacy = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_tol=1.0e-8_dp, &
        jacobian_product=jvp, line_search='none', maxiter=25, f_tol=1.0e-11_dp)
    legacy_calls = jvp_calls
    jvp_calls = 0
    cached = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=6, outer_k=3, outer_image_prune_tol=1.0e-8_dp, &
        jacobian_product=jvp, line_search='none', maxiter=25, f_tol=1.0e-11_dp, &
        outer_candidate_image_reuse=.true.)
    call check(legacy%success .and. cached%success, 'newton_krylov fixed pruning')
    call check(jvp_calls < legacy_calls, 'newton_krylov reduces analytic Jv calls')
    call check(maxval(abs(legacy%x - cached%x)) < 1.0e-10_dp, &
        'fixed pruning agrees with directly evaluated candidate')

    residual_calls = 0
    res = root(residual, x0, method='hybr', outer_candidate_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. residual_calls == 0, &
        'non-Krylov usage rejected before residual evaluation')
    res = root(residual, x0, method='krylov', jacobian_product=jvp, &
        outer_candidate_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires Jacobian-image pruning') > 0, &
        'candidate image reuse requires pruning')
    res = root(residual, x0, method='krylov', outer_image_prune_auto=.true., &
        outer_candidate_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires analytic jacobian_product') > 0, &
        'candidate image reuse requires analytic Jacobian')
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_image_prune_auto=.true., jacobian_product=jvp, &
        outer_candidate_image_reuse=.true.)
    call check(.not. res%success .and. res%nfev == 0, 'only built-in LGMRES allowed')
    res = root(residual, x0, method='krylov', outer_image_prune_auto=.true., &
        outer_candidate_image_reuse=.false.)
    call check(res%success, 'disabled candidate mode needs no analytic callback')

    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_candidate_image_reuse: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_candidate_image_reuse: PASS'
    print '(a,i0,a,i0,a,i0)', 'automatic Jv calls: default=', legacy_calls, &
        ', candidate-only=', candidate_calls, ', combined=', both_calls

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
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual vector
        real(dp) :: delta(n)

        residual_calls = residual_calls + 1
        delta = x - target
        fvec = (x*x - target*target) + 0.002_dp * sum(delta) * delta
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: vec(:) !! Krylov direction
        real(dp), intent(out) :: product(:) !! analytic Jacobian action
        logical, intent(out) :: ok !! success
        real(dp) :: delta(n)

        jvp_calls = jvp_calls + 1
        delta = x - target
        product = 2.0_dp * x * vec + 0.002_dp * &
            (sum(delta) * vec + sum(vec) * delta)
        ok = .true.
    end subroutine jvp
end program test_optimize_root_krylov_candidate_image_reuse
