! SPDX-License-Identifier: MIT
! Regression for memory-bounded LGMRES restarts and persistent augmentation.
program test_optimize_root_krylov_lgmres_restart
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none
    integer, parameter :: n = 36
    real(dp) :: x0(n), target(n)
    type(vector_root_result) :: baseline, equal_restart, restarted, via_newton, &
        pruned, finite_difference, bad
    integer :: i, failures, nfev_calls, jvp_calls, precond_calls
    integer :: original_nfev, original_jvp, original_precond
    integer :: jvp_fail_after

    failures = 0
    do i = 1, n
        target(i) = 0.8_dp + 0.25_dp * sin(0.23_dp * real(i, dp))
    end do
    x0 = target + 0.55_dp

    call clear_counts()
    baseline = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=25, outer_k=6, inner_orthogonalization='selective', &
        jacobian_product=jvp, inner_m=preconditioner, &
        line_search='none', fatol=1.0e-10_dp, maxiter=30)
    call verify(baseline, 'original LGMRES')
    original_nfev = nfev_calls
    original_jvp = jvp_calls
    original_precond = precond_calls

    call clear_counts()
    equal_restart = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=25, inner_restart=25, outer_k=6, &
        inner_orthogonalization='selective', jacobian_product=jvp, &
        inner_m=preconditioner, line_search='none', fatol=1.0e-10_dp, maxiter=30)
    call verify(equal_restart, 'full-width LGMRES')
    call check(all(equal_restart%x == baseline%x), 'full-width path unchanged')
    call check(nfev_calls == original_nfev, 'full-width nfev unchanged')
    call check(jvp_calls == original_jvp, 'full-width Jv count unchanged')
    call check(precond_calls == original_precond, 'full-width preconditioner count unchanged')

    call clear_counts()
    restarted = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=25, inner_restart=5, outer_k=6, &
        inner_orthogonalization='selective', jacobian_product=jvp, &
        inner_m=preconditioner, inner_check_residual=.true., &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(restarted, 'restarted LGMRES')
    call check(maxval(abs(restarted%x - baseline%x)) < 1.0e-8_dp, &
        'restarted solution agrees with original')

    call clear_counts()
    via_newton = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_restart=4, inner_maxiter=25, outer_k=6, &
        jacobian_product=jvp, inner_m=preconditioner, &
        inner_orthogonalization='twice', line_search='none', &
        f_tol=1.0e-10_dp, maxiter=35)
    call verify(via_newton, 'newton_krylov forwarding')

    call clear_counts()
    pruned = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=25, inner_restart=4, outer_k=6, &
        outer_image_prune_auto=.true., outer_image_reuse=.true., &
        outer_candidate_image_reuse=.true., outer_candidate_cache_max_mib=1, &
        outer_image_basis_reuse=.true., jacobian_product=jvp, &
        inner_m=preconditioner, inner_orthogonalization='selective', &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(pruned, 'restarted image-pruned LGMRES')

    call clear_counts()
    finite_difference = root(residual, x0, method='krylov', &
        inner_method='lgmres', inner_maxiter=25, inner_restart=5, &
        outer_k=6, line_search='none', fatol=1.0e-9_dp, maxiter=35)
    call check(finite_difference%success, 'restarted finite-difference convergence')
    call check(maxval(abs(finite_difference%x - target)) < 1.0e-7_dp, &
        'finite-difference root accuracy')
    call check(finite_difference%nfev == nfev_calls, &
        'finite-difference residual evaluation accounting')

    ! The outer iteration limit is one, so there can be at most seven
    ! ordinary Arnoldi steps plus one residual Jv check per restart cycle.
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=7, inner_restart=3, outer_k=2, &
        jacobian_product=jvp, line_search='none', &
        fatol=1.0e-40_dp, maxiter=1)
    call check(jvp_calls >= 1 .and. jvp_calls <= 14, &
        'restart respects ordinary Arnoldi budget')
    call check(bad%nfev == nfev_calls, 'budget test counts residual evaluations')

    ! Failure of the operator inside a restarted cycle must be propagated.
    call clear_counts()
    jvp_fail_after = 2
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=20, inner_restart=3, outer_k=2, &
        jacobian_product=jvp, line_search='none', maxiter=10)
    call check(.not. bad%success, 'restarted Jv callback failure is reported')
    call check(bad%nfev == nfev_calls, 'callback failure counts residual evaluations')

    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=0)
    call rejected(bad, 'zero restart')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_restart=-1)
    call rejected(bad, 'negative restart')
    call clear_counts()
    bad = root(residual, x0, method='krylov', inner_method='bicgstab', &
        inner_restart=3)
    call rejected(bad, 'unsupported method')
    call clear_counts()
    bad = root(residual, x0, method='hybr', inner_restart=3)
    call rejected(bad, 'non-Krylov method')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_lgmres_restart: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_lgmres_restart: PASS'
contains
    subroutine clear_counts()
        nfev_calls = 0
        jvp_calls = 0
        precond_calls = 0
        jvp_fail_after = huge(0)
    end subroutine

    subroutine check(condition, label)
        logical, intent(in) :: condition !! true to pass
        character(len=*), intent(in) :: label !! diagnostic message
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine

    subroutine rejected(r, label)
        type(vector_root_result), intent(in) :: r !! invalid result
        character(len=*), intent(in) :: label !! diagnostic
        call check(.not. r%success, trim(label)//' reports failure')
        call check(nfev_calls == 0 .and. r%nfev == 0, &
            trim(label)//' rejected before residual evaluation')
    end subroutine

    subroutine verify(r, label)
        type(vector_root_result), intent(in) :: r !! solver result
        character(len=*), intent(in) :: label !! case name
        call check(r%success, trim(label)//' converges: '//trim(r%message))
        call check(all(ieee_is_finite(r%x)), trim(label)//' has finite solution')
        call check(maxval(abs(r%x - target)) < 1.0e-8_dp, &
            trim(label)//' solves expected system')
        call check(r%nfev == nfev_calls, trim(label)//' counts residual evaluations')
        call check(jvp_calls > 0 .and. precond_calls > 0, &
            trim(label)//' uses Jv and preconditioner')
    end subroutine

    subroutine apply_linear(v, product)
        real(dp), intent(in) :: v(:) !! direction
        real(dp), intent(out) :: product(:) !! linear image
        integer :: j
        product = 0.0_dp
        do j = 1, n
            product(j) = (2.1_dp + 0.03_dp * real(j, dp)) * v(j)
            if (j > 1) product(j) = product(j) - 0.16_dp * v(j - 1)
            if (j < n) product(j) = product(j) + 0.11_dp * v(j + 1)
        end do
        product = product + 0.03_dp * sum(v) / real(n, dp)
    end subroutine

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(out) :: f(:) !! nonlinear residual
        nfev_calls = nfev_calls + 1
        call apply_linear(x - target, f)
        f = f + 0.02_dp * (x*x - target*target)
    end subroutine

    subroutine jvp(x, v, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: v(:) !! direction
        real(dp), intent(out) :: product(:) !! J(x)*v
        logical, intent(out) :: ok !! true for successful application
        jvp_calls = jvp_calls + 1
        if (jvp_calls >= jvp_fail_after) then
            product = 0.0_dp
            ok = .false.
            return
        end if
        call apply_linear(v, product)
        product = product + 0.04_dp * x * v
        ok = all(ieee_is_finite(product))
    end subroutine

    subroutine preconditioner(v, product, ok)
        real(dp), intent(in) :: v(:) !! right hand side
        real(dp), intent(out) :: product(:) !! approximate inverse
        logical, intent(out) :: ok !! success
        integer :: j
        precond_calls = precond_calls + 1
        do j = 1, n
            product(j) = v(j) / (2.1_dp + 0.03_dp * real(j, dp) + 0.04_dp * x0(j))
        end do
        ok = all(ieee_is_finite(product))
    end subroutine
end program test_optimize_root_krylov_lgmres_restart
