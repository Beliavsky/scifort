! SPDX-License-Identifier: BSD-3-Clause
! Public LGMRES regression for opt-in redundant outer-vector pruning.
program test_optimize_root_krylov_augmentation
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

    fails = 0
    baseline = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none')
    zero_tol = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        outer_prune_tol=0.0_dp)
    call check(baseline%success .and. baseline%nit > 3, &
        'default LGMRES takes multiple nonlinear steps')
    call check(zero_tol%success .eqv. baseline%success, &
        'zero tolerance preserves default success flag')
    call check(zero_tol%nit == baseline%nit .and. zero_tol%nfev == baseline%nfev .and. &
        all(zero_tol%x == baseline%x), 'zero tolerance retains bitwise default behavior')

    filtered = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=4, outer_k=3, jacobian_product=jvp, &
        fatol=1.0e-10_dp, maxiter=25, line_search='none', &
        inner_orthogonalization='selective', outer_prune_tol=1.0e-8_dp)
    call check(filtered%success .and. maxval(abs(filtered%x - expected)) < 1.0e-10_dp, &
        'filtered, reorthogonalized LGMRES converges')

    res = newton_krylov(residual, x0, inner_method='lgmres', inner_maxiter=4, &
        outer_k=3, jacobian_product=jvp, f_tol=1.0e-10_dp, maxiter=25, &
        line_search='none', outer_prune_tol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%x - expected)) < 1.0e-10_dp, &
        'newton_krylov forwards outer_prune_tol')

    res = newton_krylov(residual, x0, inner_method='lgmres', inner_maxiter=4, &
        outer_k=3, jacobian_product=jvp, inner_m=diagonal_preconditioner, &
        inner_orthogonalization='twice', f_tol=1.0e-10_dp, maxiter=25, &
        line_search='none', outer_prune_tol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%x - expected)) < 1.0e-10_dp, &
        'filtered LGMRES works with preconditioning and two-pass Arnoldi')

    ! Invalid options must be reported before calling the residual function.
    res = root(residual, x0, method='krylov', outer_prune_tol=-1.0_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'finite and in [0, 1)') /= 0, 'negative tolerance rejected')
    res = root(residual, x0, method='krylov', outer_prune_tol=1.0_dp)
    call check(.not. res%success .and. res%nfev == 0, 'unit tolerance rejected')
    res = root(residual, x0, method='krylov', &
        outer_prune_tol=ieee_value(1.0_dp, ieee_quiet_nan))
    call check(.not. res%success .and. res%nfev == 0, 'nonfinite tolerance rejected')
    res = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in lgmres') /= 0, &
        'unsupported GMRES inner method rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_k=0, outer_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0, &
        'no-augmentation case rejected')
    res = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_solver=unused_inner_solver, outer_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in lgmres') /= 0, &
        'user-provided inner solver rejected')
    res = root(residual, x0, method='broyden1', outer_prune_tol=1.0e-6_dp)
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') /= 0, &
        'non-Krylov method rejected')

    if (fails > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_augmentation: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_augmentation: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! success criterion
        character(len=*), intent(in) :: label !! description

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual

        fvec = x*x - expected*expected
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! analytic Jacobian-vector product
        logical, intent(out) :: ok !! callback result

        product = 2.0_dp*x*vec
        ok = .true.
    end subroutine jvp

    subroutine diagonal_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! linear-system vector
        real(dp), intent(out) :: product(:) !! scaled approximate inverse
        logical, intent(out) :: ok !! valid result

        product = vec / 8.0_dp
        ok = .true.
    end subroutine diagonal_preconditioner

    subroutine unused_inner_solver(matvec, preconditioner, rhs, rtol, &
            maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! linear-system right-hand side
        real(dp), intent(in) :: rtol !! inner relative stopping tolerance
        integer, intent(in) :: maxiter !! inner iteration limit
        real(dp), intent(out) :: solution(:) !! computed step
        logical, intent(out) :: ok !! successful step

        solution = rhs
        ok = .true.
    end subroutine unused_inner_solver

end program test_optimize_root_krylov_augmentation
