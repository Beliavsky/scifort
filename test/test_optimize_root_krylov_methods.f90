! SPDX-License-Identifier: MIT
! Public root and newton_krylov dispatch for the four additional inner methods.
program test_optimize_root_krylov_methods
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    real(dp), parameter :: target(4) = [1.0_dp, -2.0_dp, 3.0_dp, -1.0_dp]
    real(dp) :: a(4, 4), rhs(4), first(4)
    type(vector_root_result) :: r
    integer :: failures, j, preconditioner_calls
    integer, parameter :: scipy_nfev(3) = [9, 10, 9]
    integer, parameter :: scipy_curved_nfev(3) = [34, 30, 37]
    integer, parameter :: scipy_curved_nit(3) = [9, 8, 9]
    character(len=8), parameter :: methods(3) = [character(len=8) :: &
        'bicgstab', 'cgs', 'tfqmr']

    failures = 0
    preconditioner_calls = 0
    a = 0.0_dp
    a(1, :) = [4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a(2, :) = [0.0_dp, 3.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 0.0_dp, 2.0_dp, 1.0_dp]
    a(4, :) = [1.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
    rhs = matmul(a, target)
    first = 0.0_dp

    do j = 1, size(methods)
        r = root(linear_residual, first, method='krylov', &
            inner_method=trim(methods(j)), inner_maxiter=30, &
            line_search='none', fatol=1.0e-7_dp, rdiff=1.0e-6_dp)
        call check(r%success, trim(methods(j))//' public root converged')
        call check(maxval(abs(r%x - target)) < 1.0e-6_dp, &
            trim(methods(j))//' nonlinear root solution')
        call check(r%nit == 2 .and. r%nfev == scipy_nfev(j), &
            trim(methods(j))//' SciPy 1.17 iteration and evaluation counts')

        r = newton_krylov(linear_residual, first, inner_method=trim(methods(j)), &
            inner_maxiter=30, line_search='none', f_tol=1.0e-7_dp, rdiff=1.0e-6_dp)
        call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-6_dp, &
            trim(methods(j))//' direct newton_krylov solution')
    end do

    ! Opt-in TFQMR verifies J*s = f on the original unpreconditioned system.
    ! Its extra matrix-free product is intentionally not included in the
    ! SciPy-compatible reference counts checked above.
    r = root(linear_residual, first, method='krylov', &
        inner_method='tfqmr', inner_maxiter=30, inner_verify_residual=.true., &
        line_search='none', fatol=1.0e-7_dp, rdiff=1.0e-6_dp)
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-6_dp, &
        'public root with verified TFQMR converges')
    call check(r%nfev > scipy_nfev(3), &
        'verified TFQMR counts additional Jacobian-vector products')

    r = newton_krylov(linear_residual, first, inner_method='tfqmr', &
        inner_verify_residual=.true., inner_maxiter=30, &
        line_search='none', f_tol=1.0e-7_dp, rdiff=1.0e-6_dp)
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-6_dp, &
        'direct newton_krylov forwards verified TFQMR')

    r = root(linear_residual, first, method='krylov', &
        inner_method='tfqmr', inner_maxiter=1, inner_verify_residual=.true., &
        line_search='none', fatol=1.0e-7_dp, rdiff=1.0e-6_dp)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'verified TFQMR reports failure when inner budget is inadequate')

    r = root(linear_residual, first, method='krylov', &
        inner_method='bicgstab', inner_verify_residual=.true.)
    call check(.not. r%success .and. &
        index(r%message, 'requires the built-in tfqmr') > 0, &
        'verified TFQMR option rejects other Krylov methods')

    r = root(linear_residual, first, method='broyden1', &
        inner_verify_residual=.true.)
    call check(.not. r%success .and. &
        index(r%message, 'requires method=krylov') > 0, &
        'verified TFQMR option rejects non-Krylov methods')

    a = 0.0_dp
    a(1, :) = [4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a(2, :) = [1.0_dp, -3.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
    a(4, :) = [0.0_dp, 0.0_dp, 1.0_dp, -4.0_dp]
    rhs = matmul(a, target)
    r = root(linear_residual, first, method='krylov', inner_method='minres', &
        inner_maxiter=30, line_search='none', fatol=1.0e-7_dp, rdiff=1.0e-6_dp)
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-6_dp, &
        'MINRES public root symmetric indefinite system')
    call check(r%nit == 2 .and. r%nfev == 6, &
        'MINRES SciPy 1.17 iteration and evaluation counts')

    preconditioner_calls = 0
    r = newton_krylov(linear_residual, first, inner_method='minres', &
        inner_maxiter=30, line_search='none', f_tol=1.0e-7_dp, &
        inner_m=half_inverse, rdiff=1.0e-6_dp)
    call check(r%success .and. maxval(abs(r%x - target)) < 1.0e-6_dp, &
        'MINRES public root symmetric preconditioner')
    call check(preconditioner_calls > 0, 'public MINRES preconditioner invoked')

    r = newton_krylov(linear_residual, first, inner_method='minres', &
        inner_maxiter=30, inner_m=negative_inverse, line_search='none', &
        f_tol=1.0e-7_dp)
    call check(.not. r%success .and. index(r%message, 'inversion failed') > 0, &
        'MINRES invalid preconditioner rejected by root wrapper')

    ! A genuinely nonlinear system checks repeated Jacobian reconstruction.
    do j = 1, size(methods)
        r = root(quadratic_residual, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
            method='krylov', inner_method=trim(methods(j)), &
            inner_maxiter=30, line_search='none', fatol=1.0e-8_dp, &
            rdiff=1.0e-6_dp)
        call check(r%success, trim(methods(j))//' curved root converges')
        call check(maxval(abs(r%x - [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])) < 1.0e-7_dp, &
            trim(methods(j))//' curved root solution')
        call check(r%nit == scipy_curved_nit(j) .and. &
            r%nfev == scipy_curved_nfev(j), &
            trim(methods(j))//' curved root SciPy 1.17 work counts')
    end do
    r = root(quadratic_residual, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
        method='krylov', inner_method='minres', inner_maxiter=30, &
        line_search='none', fatol=1.0e-8_dp, rdiff=1.0e-6_dp)
    call check(r%success .and. maxval(abs(r%x - &
        [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])) < 1.0e-7_dp, &
        'MINRES curved symmetric-Jacobian root')
    call check(r%nit == 9 .and. r%nfev == 27, &
        'MINRES curved root SciPy 1.17 work counts')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root_krylov_methods: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_root_krylov_methods: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! assertion
        character(len=*), intent(in) :: description !! failure description

        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//description
        end if
    end subroutine check

    subroutine linear_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual

        fvec = matmul(a, x) - rhs
    end subroutine linear_residual

    subroutine quadratic_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear trial point
        real(dp), intent(out) :: fvec(:) !! vector of quadratic residuals

        fvec = x**2 - [4.0_dp, 9.0_dp, 16.0_dp, 25.0_dp]
    end subroutine quadratic_residual

    subroutine half_inverse(v, w, ok)
        real(dp), intent(in) :: v(:) !! input
        real(dp), intent(out) :: w(:) !! scaled inverse result
        logical, intent(out) :: ok !! preconditioner status

        preconditioner_calls = preconditioner_calls + 1
        w = 0.5_dp * v
        ok = .true.
    end subroutine half_inverse

    subroutine negative_inverse(v, w, ok)
        real(dp), intent(in) :: v(:) !! input
        real(dp), intent(out) :: w(:) !! invalid inverse result
        logical, intent(out) :: ok !! preconditioner status

        w = -v
        ok = .true.
    end subroutine negative_inverse

end program test_optimize_root_krylov_methods
