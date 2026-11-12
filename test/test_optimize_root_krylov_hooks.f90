! SPDX-License-Identifier: MIT
! Matrix-free Krylov preconditioner, adaptive update, and inner-solver tests.
program test_optimize_root_krylov_hooks
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : newton_krylov, root, vector_root_result
    implicit none

    real(dp), parameter :: diagonal(4) = [1.0_dp, 10.0_dp, 100.0_dp, 1000.0_dp]
    real(dp), parameter :: rhs(4) = [2.0_dp, -3.0_dp, 4.0_dp, -5.0_dp]
    integer :: failures
    integer :: apply_calls
    integer :: update_calls
    integer :: custom_calls
    integer :: jv_calls
    integer :: identity_calls
    real(dp) :: current_x(2)
    type(vector_root_result) :: r

    failures = 0
    apply_calls = 0
    update_calls = 0
    custom_calls = 0
    jv_calls = 0
    identity_calls = 0

    r = root(scaled_linear, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='gmres', inner_maxiter=1, &
        line_search='none', rdiff=1.0e-6_dp, fatol=1.0e-7_dp, inner_m=diag_inverse)
    call check(r%success, 'preconditioned one-step GMRES converged')
    call check(r%nit == 2 .and. r%nfev == 4, &
        'GMRES matches SciPy 1.17 iteration and evaluation counts')
    call check(maxval(abs(r%x - rhs / diagonal)) < 1.0e-7_dp, &
        'preconditioned GMRES solution')
    call check(apply_calls > 0, 'GMRES invokes right preconditioner')

    apply_calls = 0
    r = newton_krylov(scaled_linear, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        inner_method='lgmres', inner_maxiter=1, line_search='none', &
        rdiff=1.0e-6_dp, f_tol=1.0e-7_dp, inner_m=diag_inverse)
    call check(r%success, 'preconditioned LGMRES converged')
    call check(r%nit == 2 .and. r%nfev == 3, &
        'LGMRES matches SciPy 1.17 iteration and evaluation counts')
    call check(maxval(abs(r%x - rhs / diagonal)) < 1.0e-7_dp, &
        'preconditioned LGMRES solution')
    call check(apply_calls > 0, 'LGMRES invokes right preconditioner')

    current_x = [1.0_dp, 1.0_dp]
    apply_calls = 0
    update_calls = 0
    r = newton_krylov(quadratic, [1.0_dp, 1.0_dp], rdiff=1.0e-6_dp, &
        inner_maxiter=2, inner_m=adaptive_inverse, inner_m_update=update_inverse, &
        line_search='none', f_tol=1.0e-8_dp, outer_k=1)
    call check(r%success, 'adaptive nonlinear preconditioner converged')
    call check(r%nit == 7 .and. r%nfev == 16, &
        'adaptive LGMRES matches SciPy 1.17 work counts')
    call check(maxval(abs(r%x - [2.0_dp, 3.0_dp])) < 1.0e-7_dp, &
        'adaptive preconditioner solution')
    call check(apply_calls > 0, 'adaptive inverse invoked')
    call check(update_calls == r%nit - 1 .and. update_calls > 1, &
        'preconditioner update called once per accepted step')
    call check(maxval(abs(current_x - r%x)) < 1.0e-12_dp, &
        'preconditioner receives latest accepted iterate')

    custom_calls = 0
    jv_calls = 0
    apply_calls = 0
    r = root(scaled_linear, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='user-defined', inner_solver=custom_linear, &
        inner_m=diag_inverse, inner_maxiter=4, rdiff=1.0e-6_dp, &
        line_search='none', fatol=1.0e-7_dp)
    call check(r%success, 'custom inner linear solver converged')
    call check(maxval(abs(r%x - rhs / diagonal)) < 1.0e-7_dp, &
        'custom inner linear solution')
    call check(custom_calls > 0 .and. jv_calls > 0 .and. apply_calls > 0, &
        'custom solver received both Jacobian product and preconditioner')
    call check(r%nfev == 1 + jv_calls + custom_calls, &
        'custom Jacobian evaluations counted exactly')

    custom_calls = 0
    r = root(scaled_linear, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='tfqmr', inner_solver=custom_linear, &
        inner_verify_residual=.true.)
    call check(.not. r%success .and. &
        index(r%message, 'requires the built-in tfqmr') > 0, &
        'verification option rejects an overriding custom solver')
    call check(custom_calls == 0, 'invalid option fails before custom callback')

    identity_calls = 0
    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], &
        inner_solver=custom_identity, inner_maxiter=2, line_search='none', &
        f_tol=1.0e-7_dp)
    call check(r%success .and. identity_calls > 0, &
        'custom inner solver gets identity preconditioner by default')

    r = root(scaled_linear, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='broyden1', inner_m=diag_inverse)
    call check(.not. r%success .and. index(r%message, 'require method=krylov') > 0, &
        'reject preconditioner for a non-Krylov root method')

    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], inner_m=failing_operator)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'failure from preconditioner propagates')

    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], inner_solver=failing_solver)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'failure from user-defined inner solver propagates')

    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], inner_solver=nan_solver)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'nonfinite custom inner solver output rejected')

    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], inner_method='unknown')
    call check(.not. r%success .and. index(r%message, 'unsupported Krylov') > 0, &
        'built-in solver name validated without a custom callback')

    r = newton_krylov(identity_system, [0.0_dp, 0.0_dp], inner_m=nan_operator)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'nonfinite preconditioner output rejected')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root_krylov_hooks: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_root_krylov_hooks: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! condition expected to be true
        character(len=*), intent(in) :: description !! assertion label

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine scaled_linear(x, fvec)
        real(dp), intent(in) :: x(:) !! current trial iterate
        real(dp), intent(out) :: fvec(:) !! diagonally scaled residual

        fvec = diagonal * x - rhs
    end subroutine scaled_linear

    subroutine diag_inverse(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! inverse diagonal action
        logical, intent(out) :: ok !! action succeeded

        apply_calls = apply_calls + 1
        product = vec / diagonal
        ok = .true.
    end subroutine diag_inverse

    subroutine quadratic(x, fvec)
        real(dp), intent(in) :: x(:) !! current trial iterate
        real(dp), intent(out) :: fvec(:) !! nonlinear residual

        fvec = x**2 - [4.0_dp, 9.0_dp]
    end subroutine quadratic

    subroutine adaptive_inverse(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! approximate inverse Jacobian action
        logical, intent(out) :: ok !! action succeeded

        apply_calls = apply_calls + 1
        product = vec / (2.0_dp * current_x)
        ok = .true.
    end subroutine adaptive_inverse

    subroutine update_inverse(x, fvec)
        real(dp), intent(in) :: x(:) !! just-accepted iterate
        real(dp), intent(in) :: fvec(:) !! residual at x

        update_calls = update_calls + 1
        current_x = x
        if (maxval(abs(fvec - (x**2 - [4.0_dp, 9.0_dp]))) > 1.0e-12_dp) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: adaptive inverse update residual mismatch'
        end if
    end subroutine update_inverse

    subroutine custom_linear(matvec, preconditioner, b, tol, maxiter, sol, ok)
        use scifort_optimize, only : vector_root_linear_operator
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector operation
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(in) :: tol !! relative tolerance request
        integer, intent(in) :: maxiter !! allowed inner iterations
        real(dp), intent(out) :: sol(:) !! approximate inner solution
        logical, intent(out) :: ok !! status of inner solve

        real(dp) :: product(size(b)), correction(size(b))
        logical :: mat_ok

        custom_calls = custom_calls + 1
        call preconditioner(b, sol, ok)
        if (.not. ok) return
        call matvec(sol, product, mat_ok)
        if (.not. mat_ok) then
            ok = .false.
            return
        end if
        jv_calls = jv_calls + 1
        call preconditioner(b - product, correction, ok)
        if (ok) sol = sol + correction
        if (tol < 0.0_dp .or. tol > 1.0_dp .or. maxiter /= 4) ok = .false.
    end subroutine custom_linear

    subroutine identity_system(x, fvec)
        real(dp), intent(in) :: x(:) !! current trial iterate
        real(dp), intent(out) :: fvec(:) !! residual

        fvec = x - [1.0_dp, 2.0_dp]
    end subroutine identity_system

    subroutine custom_identity(matvec, preconditioner, b, tol, maxiter, sol, ok)
        use scifort_optimize, only : vector_root_linear_operator
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector operation
        procedure(vector_root_linear_operator) :: preconditioner !! identity action
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(in) :: tol !! requested relative tolerance
        integer, intent(in) :: maxiter !! inner iteration budget
        real(dp), intent(out) :: sol(:) !! solution
        logical, intent(out) :: ok !! solver status

        real(dp) :: probe(size(b))
        logical :: mat_ok

        call preconditioner(b, sol, ok)
        if (.not. ok) return
        call matvec(sol, probe, mat_ok)
        ok = mat_ok .and. maxval(abs(probe - b)) < 1.0e-5_dp
        if (ok .and. maxiter == 2 .and. tol >= 0.0_dp) identity_calls = identity_calls + 1
    end subroutine custom_identity

    subroutine failing_operator(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! input vector
        real(dp), intent(out) :: product(:) !! transformed vector
        logical, intent(out) :: ok !! failure flag

        product = 0.0_dp
        ok = .false.
    end subroutine failing_operator

    subroutine nan_operator(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! input vector
        real(dp), intent(out) :: product(:) !! transformed vector
        logical, intent(out) :: ok !! success flag

        product = ieee_value(0.0_dp, ieee_quiet_nan)
        ok = .true.
    end subroutine nan_operator

    subroutine nan_solver(matvec, preconditioner, b, tol, maxiter, sol, ok)
        use scifort_optimize, only : vector_root_linear_operator
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector operation
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(in) :: tol !! relative tolerance
        integer, intent(in) :: maxiter !! allowed iterations
        real(dp), intent(out) :: sol(:) !! nonfinite output vector
        logical, intent(out) :: ok !! misleading success status

        sol = ieee_value(0.0_dp, ieee_quiet_nan)
        ok = .true.
    end subroutine nan_solver

    subroutine failing_solver(matvec, preconditioner, b, tol, maxiter, sol, ok)
        use scifort_optimize, only : vector_root_linear_operator
        procedure(vector_root_linear_operator) :: matvec !! Jacobian-vector operation
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(in) :: tol !! relative tolerance
        integer, intent(in) :: maxiter !! allowed iterations
        real(dp), intent(out) :: sol(:) !! output vector
        logical, intent(out) :: ok !! failure status

        ok = .false.
    end subroutine failing_solver

end program test_optimize_root_krylov_hooks
