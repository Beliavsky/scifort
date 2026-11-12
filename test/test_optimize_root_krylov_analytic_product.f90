! SPDX-License-Identifier: MIT
! Analytic Krylov Jacobian-vector actions and cancellation/exception behavior.
program test_optimize_root_krylov_analytic_product
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: diag(4) = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: target(4) = [0.25_dp, -1.0_dp, 0.8_dp, 1.5_dp]
    real(dp), parameter :: curved_target(2) = [1.3_dp, 1.2_dp]
    character(len=8), parameter :: methods(6) = [character(len=8) :: &
        'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres', 'tfqmr']
    type(vector_root_result) :: r
    integer :: failures, jvp_calls, j
    real(dp) :: measured(2), measured_default(2)
    logical :: probe_ok

    failures = 0
    jvp_calls = 0
    measured = 0.0_dp
    measured_default = 0.0_dp

    ! Even a nominally representable global perturbation can disappear in
    ! subtraction against a residual offset of O(1e12).
    r = root(offset_residual, [0.1_dp, -0.2_dp], method='krylov', &
        inner_solver=probe_jv, rdiff=1.0e-8_dp)
    measured_default = measured
    call check(probe_ok .and. r%nfev == 2, &
        'baseline finite-difference callback uses a residual evaluation')
    call check(maxval(abs(measured_default - [3.0_dp, -2.0_dp])) > 0.5_dp, &
        'offset causes significant cancellation in finite differences')

    jvp_calls = 0
    r = root(offset_residual, [0.1_dp, -0.2_dp], method='krylov', &
        inner_solver=probe_jv, jacobian_product=offset_jvp)
    call check(probe_ok .and. r%nfev == 1, &
        'analytic callback has no additional residual evaluations')
    call check(jvp_calls == 1 .and. all(measured == [3.0_dp, -2.0_dp]), &
        'analytic J*v avoids cancellation and skips zero directions')

    ! All built-in methods must receive the same analytic matrix-free product.
    do j = 1, size(methods)
        jvp_calls = 0
        r = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            method='krylov', inner_method=trim(methods(j)), &
            jacobian_product=linear_jvp, inner_maxiter=30, &
            line_search='none', fatol=1.0e-12_dp)
        call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-10_dp, &
            trim(methods(j))//' analytic product converged')
        call check(jvp_calls > 0 .and. r%nfev <= 3, &
            trim(methods(j))//' products are not counted as fun evaluations')
    end do

    ! Verify the callback receives the evolving nonlinear iterate, not x0.
    jvp_calls = 0
    r = newton_krylov(curved_residual, [1.0_dp, 1.0_dp], &
        inner_method='gmres', jacobian_product=curved_jvp, &
        line_search='none', f_tol=1.0e-12_dp)
    if (.not. r%success) print *, 'DEBUG curved:',r%message, r%nfev,r%nit,r%x,r%fun,jvp_calls
    call check(r%success .and. maxval(abs(r%x - curved_target)) < 1.0e-10_dp, &
        'newton_krylov refreshes analytic Jacobian at each nonlinear iterate')
    call check(jvp_calls > 1 .and. r%nfev < 15, &
        'nonlinear analytic callback count separate from objective calls')

    ! Composes with true-residual checking and the new right-TFQMR path.
    jvp_calls = 0
    r = newton_krylov(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        inner_method='tfqmr', jacobian_product=linear_jvp, &
        inner_m=diagonal_inverse, inner_verify_residual=.true., &
        inner_tfqmr_right_preconditioned=.true., inner_maxiter=30, &
        line_search='none', f_tol=1.0e-12_dp)
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-10_dp, &
        'analytic product with preconditioned verified TFQMR')

    ! The low-level custom solver also receives the analytic operator.
    jvp_calls = 0
    r = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_solver=custom_inner, &
        jacobian_product=linear_jvp, fatol=1.0e-12_dp, line_search='none')
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-10_dp, &
        'custom inner solver consumes exact product')
    call check(jvp_calls > 0 .and. r%nfev <= 3, &
        'custom inner solver analytic calls leave nfev unchanged')

    ! Failed and nonfinite analytic products must fail the inner solve rather
    ! than accept uninitialized or nonfinite outputs.
    jvp_calls = 0
    r = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', jacobian_product=failed_jvp)
    call check(.not. r%success .and. r%nfev == 1 .and. jvp_calls > 0 .and. &
        index(r%message, 'inversion failed') > 0, &
        'false ok from analytic callback propagates')

    jvp_calls = 0
    r = newton_krylov(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        jacobian_product=nan_jvp)
    call check(.not. r%success .and. r%nfev == 1 .and. jvp_calls > 0 .and. &
        index(r%message, 'inversion failed') > 0, &
        'nonfinite analytic product rejected')

    jvp_calls = 0
    r = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='broyden1', jacobian_product=linear_jvp)
    call check(.not. r%success .and. r%nfev == 0 .and. &
        index(r%message, 'requires method=krylov') > 0, &
        'non-Krylov methods reject analytic product')
    r = newton_krylov(curved_residual, [1.0_dp, 1.0_dp], &
        jacobian_product=curved_jvp, rdiff=1.0e-5_dp)
    call check(.not. r%success .and. r%nfev == 0 .and. &
        index(r%message, 'cannot be combined') > 0, &
        'analytic callback cannot silently override rdiff')
    r = root(curved_residual, [1.0_dp, 1.0_dp], method='krylov', &
        jacobian_product=curved_jvp, jacobian_diff_scheme='3-point')
    call check(.not. r%success .and. r%nfev == 0 .and. &
        index(r%message, 'cannot be combined') > 0, &
        'analytic callback cannot silently override difference stencil')
    r = root(curved_residual, [1.0_dp, 1.0_dp], method='krylov', &
        jacobian_product=curved_jvp, jacobian_step_mode='adaptive')
    call check(.not. r%success .and. r%nfev == 0 .and. &
        index(r%message, 'cannot be combined') > 0, &
        'analytic callback cannot silently override difference step mode')
    call check(jvp_calls == 0, 'invalid options do not invoke analytic callback')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_analytic_product: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_analytic_product: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! test condition
        character(len=*), intent(in) :: description !! test label

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

    subroutine offset_residual(x, f)
        real(dp), intent(in) :: x(:) !! independent variables
        real(dp), intent(out) :: f(:) !! residual with a large constant offset

        f(1) = 1.0e12_dp + x(1)
        f(2) = 1.0e12_dp - 2.0_dp * x(2)
    end subroutine offset_residual

    subroutine offset_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! exact Jacobian product
        logical, intent(out) :: ok !! valid result

        jvp_calls = jvp_calls + 1
        product = [vec(1), -2.0_dp * vec(2)]
        ok = .true.
    end subroutine offset_jvp

    subroutine probe_jv(matvec, precond, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! sampled Jacobian action
        procedure(vector_root_linear_operator) :: precond !! unused inverse Jacobian
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! relative tolerance
        integer, intent(in) :: maxiter !! iteration cap
        real(dp), intent(out) :: solution(:) !! intentionally unused
        logical, intent(out) :: ok !! deliberately reject after measurement

        real(dp) :: zero(2)
        logical :: zero_ok

        call matvec([0.0_dp, 0.0_dp], zero, zero_ok)
        call check(zero_ok .and. all(zero == 0.0_dp), 'zero direction is exact')
        call matvec([3.0_dp, 1.0_dp], measured, probe_ok)
        solution = 0.0_dp
        ok = .false.
    end subroutine probe_jv

    subroutine linear_residual(x, f)
        real(dp), intent(in) :: x(:) !! independent variables
        real(dp), intent(out) :: f(:) !! linear residual

        f = diag * (x - target)
    end subroutine linear_residual

    subroutine linear_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! diagonal Jacobian action
        logical, intent(out) :: ok !! success flag

        jvp_calls = jvp_calls + 1
        product = diag * vec
        ok = .true.
    end subroutine linear_jvp

    subroutine curved_residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: f(:) !! nonlinear residual

        f = x*x - curved_target*curved_target
    end subroutine curved_residual

    subroutine curved_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: vec(:) !! product direction
        real(dp), intent(out) :: product(:) !! Jacobian at current x times vec
        logical, intent(out) :: ok !! validity

        jvp_calls = jvp_calls + 1
        product = 2.0_dp * x * vec
        ok = .true.
    end subroutine curved_jvp

    subroutine diagonal_inverse(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! preconditioner input
        real(dp), intent(out) :: product(:) !! exact inverse action
        logical, intent(out) :: ok !! validity

        product = vec / diag
        ok = .true.
    end subroutine diagonal_inverse

    subroutine custom_inner(matvec, precond, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! analytic operator
        procedure(vector_root_linear_operator) :: precond !! inverse action
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! tolerance
        integer, intent(in) :: maxiter !! iteration limit
        real(dp), intent(out) :: solution(:) !! linear-system solution
        logical, intent(out) :: ok !! success

        real(dp) :: sample(size(rhs))
        logical :: valid

        call matvec([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], sample, valid)
        ok = valid
        if (.not. ok) return
        call check(maxval(abs(sample - [diag(1), 0.0_dp, 0.0_dp, 0.0_dp])) < &
            1.0e-12_dp, 'custom solver sees analytic J*v')
        solution = rhs / diag
    end subroutine custom_inner

    subroutine failed_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! undefined on failure
        logical, intent(out) :: ok !! failure status

        jvp_calls = jvp_calls + 1
        ok = .false.
    end subroutine failed_jvp

    subroutine nan_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! nonfinite action
        logical, intent(out) :: ok !! misleading user success flag

        jvp_calls = jvp_calls + 1
        product = ieee_value(0.0_dp, ieee_quiet_nan)
        ok = .true.
    end subroutine nan_jvp

end program test_optimize_root_krylov_analytic_product
