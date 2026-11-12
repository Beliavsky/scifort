! SPDX-License-Identifier: MIT
! Bound raw Arnoldi image storage without approximating candidate Jacobian products.
program test_optimize_root_krylov_candidate_cache_budget
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 8000
    real(dp), allocatable :: target(:), starting_point(:)
    type(vector_root_result) :: unbounded, limited, sufficient, res
    integer :: failed, jvp_calls, unlimited_calls, limited_calls
    integer :: residual_calls, k
    logical :: retained_cache, nonlinear_mode

    allocate(target(n), starting_point(n))
    target = 1.0_dp
    starting_point = 2.0_dp
    failed = 0
    nonlinear_mode = .false.

    do k = 1, 2
        retained_cache = k == 2
        jvp_calls = 0
        residual_calls = 0
        unbounded = root(residual, starting_point, method='krylov', &
            inner_method='lgmres', inner_maxiter=32, outer_k=3, &
            outer_image_prune_auto=.true., outer_image_reuse=retained_cache, &
            outer_candidate_image_reuse=.true., jacobian_product=jvp, &
            fatol=1.0e-12_dp, maxiter=10, line_search='none')
        unlimited_calls = jvp_calls
        call check(unbounded%success .and. unbounded%nfev == residual_calls, &
            'unbounded cache convergence and accounting')

        jvp_calls = 0
        residual_calls = 0
        sufficient = root(residual, starting_point, method='krylov', &
            inner_method='lgmres', inner_maxiter=32, outer_k=3, &
            outer_image_prune_auto=.true., outer_image_reuse=retained_cache, &
            outer_candidate_image_reuse=.true., outer_candidate_cache_max_mib=8, &
            jacobian_product=jvp, fatol=1.0e-12_dp, maxiter=10, line_search='none')
        call check(all(unbounded%x == sufficient%x) .and. &
            jvp_calls == unlimited_calls, 'sufficient budget preserves exact cached path')
        call check(unbounded%nit == sufficient%nit .and. &
            unbounded%nfev == sufficient%nfev, 'unbounded work-count compatibility')

        jvp_calls = 0
        residual_calls = 0
        limited = root(residual, starting_point, method='krylov', &
            inner_method='lgmres', inner_maxiter=32, outer_k=3, &
            outer_image_prune_auto=.true., outer_image_reuse=retained_cache, &
            outer_candidate_image_reuse=.true., outer_candidate_cache_max_mib=1, &
            jacobian_product=jvp, fatol=1.0e-12_dp, maxiter=10, line_search='none')
        limited_calls = jvp_calls
        call check(limited%success .and. limited%nfev == residual_calls, &
            'bounded-cache convergence and nfev accounting')
        call check(maxval(abs(limited%x - unbounded%x)) < 1.0e-12_dp, &
            'exact recomputed image agrees with cached candidate image')
        call check(limited%nit == unbounded%nit .and. &
            limited%nfev == unbounded%nfev, 'bounded-cache nonlinear work counts')
        call check(limited_calls > unlimited_calls, 'small budget makes direct analytic Jv call')
    end do

    jvp_calls = 0
    res = newton_krylov(residual, starting_point, inner_method='lgmres', &
        inner_maxiter=32, outer_k=3, outer_image_prune_auto=.true., &
        outer_candidate_image_reuse=.true., outer_candidate_cache_max_mib=1, &
        jacobian_product=jvp, f_tol=1.0e-12_dp, maxiter=10, line_search='none')
    call check(res%success .and. maxval(abs(res%x - target)) < 1.0e-12_dp, &
        'newton_krylov bounded-cache forwarding')

    ! With image reuse disabled, the first inner solve fits in 1 MiB,
    ! but a later accepted augmentation direction raises the projected
    ! cache capacity above the limit. Exercise switching after allocation.
    nonlinear_mode = .true.
    do k = 1, n
        target(k) = 1.0_dp + 0.1_dp * sin(0.01_dp * real(k, dp))
        starting_point(k) = target(k) + 0.3_dp + 0.05_dp * cos(0.02_dp * real(k, dp))
    end do
    jvp_calls = 0
    unbounded = root(residual, starting_point, method='krylov', &
        inner_method='lgmres', inner_maxiter=16, outer_k=4, &
        outer_image_prune_auto=.true., outer_candidate_image_reuse=.true., &
        jacobian_product=jvp, fatol=1.0e-9_dp, maxiter=30, line_search='none')
    unlimited_calls = jvp_calls
    jvp_calls = 0
    limited = root(residual, starting_point, method='krylov', &
        inner_method='lgmres', inner_maxiter=16, outer_k=4, &
        outer_image_prune_auto=.true., outer_candidate_image_reuse=.true., &
        outer_candidate_cache_max_mib=1, jacobian_product=jvp, &
        fatol=1.0e-9_dp, maxiter=30, line_search='none')
    call check(unbounded%success .and. limited%success, &
        'cache cap switch with a changing Jacobian')
    call check(maxval(abs(unbounded%x - limited%x)) < 1.0e-8_dp .and. &
        limited%nfev == unbounded%nfev, 'late switch preserves nonlinear root')
    call check(jvp_calls > unlimited_calls, &
        'later growing history triggers the budget fallback')
    nonlinear_mode = .false.

    residual_calls = 0
    res = root(residual, starting_point, method='hybr', &
        outer_candidate_cache_max_mib=1)
    call check(.not. res%success .and. residual_calls == 0 .and. res%nfev == 0, &
        'non-Krylov use rejected before residual evaluation')
    res = root(residual, starting_point, method='krylov', &
        outer_candidate_cache_max_mib=-1)
    call check(.not. res%success .and. residual_calls == 0 .and. res%nfev == 0, &
        'negative memory cap rejected')
    res = root(residual, starting_point, method='krylov', &
        outer_candidate_cache_max_mib=1)
    call check(.not. res%success .and. residual_calls == 0 .and. res%nfev == 0, &
        'memory cap without candidate reuse rejected')

    if (failed > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_candidate_cache_budget: FAIL ', failed
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_candidate_cache_budget: PASS'
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: fvec(:) !! residual at current state
        residual_calls = residual_calls + 1
        if (nonlinear_mode) then
            fvec = x*x - target*target
        else
            fvec = x - target
        end if
    end subroutine residual

    subroutine jvp(x, v, product, ok)
        real(dp), intent(in) :: x(:) !! current state (Jacobian is constant)
        real(dp), intent(in) :: v(:) !! Krylov direction
        real(dp), intent(out) :: product(:) !! exact Jacobian image
        logical, intent(out) :: ok !! callback validity
        jvp_calls = jvp_calls + 1
        if (nonlinear_mode) then
            product = 2.0_dp * x * v
        else
            product = v
        end if
        ok = .true.
    end subroutine jvp

    subroutine check(condition, label)
        logical, intent(in) :: condition !! required test condition
        character(len=*), intent(in) :: label !! readable assertion label
        if (.not. condition) then
            failed = failed + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check
end program test_optimize_root_krylov_candidate_cache_budget
