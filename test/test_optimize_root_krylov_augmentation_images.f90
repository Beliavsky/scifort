! SPDX-License-Identifier: MIT
! Public root tests for opt-in current-Jacobian LGMRES augmentation filtering.
program test_optimize_root_krylov_augmentation_images
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: x0(4) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
    real(dp), parameter :: expected(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    type(vector_root_result) :: baseline
    type(vector_root_result) :: zero_tol
    type(vector_root_result) :: filtered
    type(vector_root_result) :: res
    integer :: fails
    integer :: jvp_calls

    fails = 0
    jvp_calls = 0
    baseline = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none')
    zero_tol = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        outer_image_prune_tol=0.0_dp)
    call check(baseline%success .and. baseline%nit > 3, 'baseline takes multiple steps')
    call check((zero_tol%success .eqv. baseline%success) .and. &
        zero_tol%nit == baseline%nit .and. zero_tol%nfev == baseline%nfev .and. &
        all(zero_tol%x == baseline%x), 'zero tolerance exactly preserves default')

    jvp_calls = 0
    filtered = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_orthogonalization='selective', outer_image_prune_tol=1.0e-8_dp)
    call check(filtered%success .and. maxval(abs(filtered%x - expected)) < 1.0e-10_dp, &
        'image-filtered LGMRES converges')
    call check(jvp_calls > 2*filtered%nit, 'image filtering invokes additional products')

    res = newton_krylov(residual, x0, inner_method='lgmres', inner_maxiter=4, &
        outer_k=3, jacobian_product=jvp, f_tol=1.0e-10_dp, maxiter=25, &
        line_search='none', outer_image_prune_tol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%x - expected)) < 1.0e-10_dp, &
        'newton_krylov forwards image pruning')

    res = newton_krylov(residual, x0, inner_method='lgmres', inner_maxiter=4, &
        outer_k=3, jacobian_product=jvp, inner_m=diagonal_preconditioner, &
        inner_orthogonalization='twice', f_tol=1.0e-10_dp, maxiter=25, &
        line_search='none', outer_image_prune_tol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%x - expected)) < 1.0e-10_dp, &
        'preconditioned and reorthogonalized image pruning converges')

    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, fatol=1.0e-7_dp, maxiter=25, &
        line_search='none', outer_image_prune_tol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%x - expected)) < 1.0e-6_dp, &
        'finite-difference products support image filtering')
    call check(res%nfev > filtered%nfev, &
        'finite-difference image checks are counted in nfev')

    ! Exactly one extra J*v per first step, using two-point differences;
    ! a centered three-point stencil requires exactly two extra samples.
    baseline = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none')
    filtered = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', outer_image_prune_tol=1.0e-8_dp)
    call check(filtered%nfev == baseline%nfev + 1, &
        'one extra forward-difference residual sample for candidate image')
    baseline = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', jacobian_diff_scheme='3-point')
    filtered = root(residual_one, [2.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=1, outer_k=1, maxiter=1, &
        line_search='none', jacobian_diff_scheme='3-point', &
        outer_image_prune_tol=1.0e-8_dp)
    call check(filtered%nfev == baseline%nfev + 2, &
        'two extra centered-difference samples for candidate image')

    ! Callback failure specifically on the *additional* image check must
    ! not silently retain a direction or accept the nonlinear solution.
    jvp_calls = 0
    res = root(residual_one, [2.0_dp], method='krylov', inner_method='lgmres', &
        inner_maxiter=1, outer_k=1, jacobian_product=failing_jvp, &
        line_search='none', outer_image_prune_tol=1.0e-8_dp)
    call check(.not. res%success .and. res%nfev == 1 .and. jvp_calls == 2, &
        'image callback failure propagates')

    ! Bad options are rejected before evaluating fun.
    res = root(residual, x0, method='krylov', outer_image_prune_tol=-1.0_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'finite and in [0, 1)') /= 0, 'negative tolerance rejected')
    res = root(residual, x0, method='krylov', outer_image_prune_tol=1.0_dp)
    call check(.not. res%success .and. res%nfev == 0, 'unit tolerance rejected')
    res = root(residual, x0, method='krylov', &
        outer_image_prune_tol=ieee_value(1.0_dp, ieee_quiet_nan))
    call check(.not. res%success .and. res%nfev == 0, 'NaN tolerance rejected')
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_image_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in lgmres') /= 0, 'GMRES rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_k=0, outer_image_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0, 'outer_k=0 rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_solver=unused_inner_solver, outer_image_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0, 'custom inner solver rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_image_prune_tol=1.0e-6_dp, outer_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'cannot be combined') /= 0, 'two filters rejected')
    res = root(residual, x0, method='broyden1', outer_image_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') /= 0, 'non-Krylov rejected')

    if (fails > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_augmentation_images: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_augmentation_images: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected result
        character(len=*), intent(in) :: label !! test description

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual values

        fvec = x*x - expected*expected
    end subroutine residual

    subroutine residual_one(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual values

        fvec = x*x - 1.0_dp
    end subroutine residual_one

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! success

        jvp_calls = jvp_calls + 1
        product = 2.0_dp*x*vec
        ok = .true.
    end subroutine jvp

    subroutine failing_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! success

        jvp_calls = jvp_calls + 1
        ok = jvp_calls /= 2
        product = 2.0_dp*x*vec
    end subroutine failing_jvp

    subroutine diagonal_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! input vector
        real(dp), intent(out) :: product(:) !! inverse approximation
        logical, intent(out) :: ok !! success

        product = vec / 8.0_dp
        ok = .true.
    end subroutine diagonal_preconditioner

    subroutine unused_inner_solver(matvec, preconditioner, rhs, rtol, &
            maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! current Jacobian
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! inner right-hand side
        real(dp), intent(in) :: rtol !! target relative tolerance
        integer, intent(in) :: maxiter !! inner iteration limit
        real(dp), intent(out) :: solution(:) !! approximate step
        logical, intent(out) :: ok !! whether solve succeeded

        solution = rhs
        ok = .true.
    end subroutine unused_inner_solver

end program test_optimize_root_krylov_augmentation_images
