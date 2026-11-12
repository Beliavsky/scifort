! SPDX-License-Identifier: MIT
! Public-entry regression for current-Jacobian retained-image basis reuse.
program test_optimize_root_krylov_image_basis_reuse
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none
    integer, parameter :: n = 8
    real(dp) :: target(n)
    real(dp) :: x0(n)
    type(vector_root_result) :: old_result, cached_result, wrapper_result, bad_result
    integer :: fails
    integer :: i
    integer :: k
    integer :: calls
    integer :: jvp_calls
    logical :: automatic
    logical :: retained_reuse

    fails = 0
    do i = 1, n
        target(i) = 0.7_dp + 0.2_dp * real(i, dp)
        x0(i) = target(i) + 0.4_dp + 0.05_dp * real(i, dp)
    end do
    do k = 1, 4
        automatic = mod(k, 2) == 0
        retained_reuse = k > 2
        jvp_calls = 0
        old_result = root(residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=8, outer_k=2, &
            outer_image_prune_tol=merge(0.0_dp, 1.0e-8_dp, automatic), &
            outer_image_prune_auto=automatic, outer_image_reuse=retained_reuse, &
            outer_candidate_image_reuse=.true., jacobian_product=jvp, &
            fatol=1.0e-12_dp, line_search='none', maxiter=40)
        cached_result = root(residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=8, outer_k=2, &
            outer_image_prune_tol=merge(0.0_dp, 1.0e-8_dp, automatic), &
            outer_image_prune_auto=automatic, outer_image_reuse=retained_reuse, &
            outer_candidate_image_reuse=.true., outer_image_basis_reuse=.true., &
            jacobian_product=jvp, fatol=1.0e-12_dp, &
            line_search='none', maxiter=40)
        call check(old_result%success .and. cached_result%success, &
            'analytic solvers converge')
        call check(cached_result%nit >= 4, 'multiple steps exercise image insertion')
        call check(old_result%nfev == cached_result%nfev .and. &
            old_result%nit == cached_result%nit, 'analytic work counts agree')
        call check(maxval(abs(cached_result%x - old_result%x)) < 1.0e-9_dp, &
            'analytic roots agree')
        call check(maxval(abs(cached_result%fun)) < 1.0e-11_dp, &
            'cached root has small true residual')
    end do
    old_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_image_prune_tol=1.0e-8_dp, &
        fatol=1.0e-8_dp, line_search='none', maxiter=40)
    cached_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_image_prune_tol=1.0e-8_dp, &
        outer_image_basis_reuse=.true., &
        fatol=1.0e-8_dp, line_search='none', maxiter=40)
    call check(old_result%success .and. cached_result%success, &
        'finite-difference solvers converge')
    call check(old_result%nfev == cached_result%nfev .and. &
        old_result%nit == cached_result%nit, 'finite-difference counts agree')
    call check(maxval(abs(cached_result%x - old_result%x)) < 1.0e-7_dp, &
        'finite-difference roots agree')
    wrapper_result = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_image_prune_tol=1.0e-8_dp, &
        outer_image_basis_reuse=.true., jacobian_product=jvp, &
        f_tol=1.0e-12_dp, line_search='none', maxiter=40)
    call check(wrapper_result%success, 'newton_krylov forwards new option')
    call check(maxval(abs(wrapper_result%fun)) < 1.0e-11_dp, &
        'newton_krylov cached root is accurate')

    ! Configuration errors must be reported without evaluating residuals.
    calls = 0
    bad_result = root(residual, x0, method='hybr', outer_image_basis_reuse=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'non-Krylov rejected')
    bad_result = root(residual, x0, method='krylov', outer_image_basis_reuse=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'missing image pruning rejected')
    bad_result = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_image_prune_tol=1.0e-8_dp, outer_image_basis_reuse=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'GMRES rejected')
    bad_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_k=0, outer_image_prune_tol=1.0e-8_dp, &
        outer_image_basis_reuse=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'zero retention rejected')
    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_image_basis_reuse: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_image_basis_reuse: PASS'
contains
    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! variable vector
        real(dp), intent(out) :: f(:) !! nonlinear residual
        calls = calls + 1
        f = x*x - target*target + 0.05_dp * sum(x - target)
    end subroutine residual

    subroutine jvp(x, vec, result_vec, ok)
        real(dp), intent(in) :: x(:) !! current linearization point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: result_vec(:) !! exact J*v
        logical, intent(out) :: ok !! success indicator
        jvp_calls = jvp_calls + 1
        result_vec = 2.0_dp*x*vec + 0.05_dp * sum(vec)
        ok = .true.
    end subroutine jvp

    subroutine check(condition, description)
        logical, intent(in) :: condition !! required condition
        character(len=*), intent(in) :: description !! diagnostic label
        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_root_krylov_image_basis_reuse
